#!/usr/bin/env python3
"""Integration oracle for achieved film sizing, dose convergence and shifts."""
import argparse,hashlib,json,math,os,re,subprocess,sys,time
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/'decks/subcmos2014'))
from validate import read

def cut(mesh,x,material,species):
 coords,tris,nodes=mesh;points={}
 for mat,vs in tris:
  if mat!=material:continue
  for ia,ib in zip(vs,vs[1:]+vs[:1]):
   a,b=coords[ia],coords[ib]
   if a[0]==b[0]:
    if a[0]==x:
     points[a[1]]=nodes[ia,mat].get(species,0);points[b[1]]=nodes[ib,mat].get(species,0)
   elif min(a[0],b[0])<=x<=max(a[0],b[0]):
    f=(x-a[0])/(b[0]-a[0]);y=a[1]+f*(b[1]-a[1]);points[y]=(1-f)*nodes[ia,mat].get(species,0)+f*nodes[ib,mat].get(species,0)
 y=np.array(sorted(points));c=np.array([points[z] for z in y]);return float(np.trapezoid(c,y)*1e-4),float(max(c))

def metrics(mesh,mat,species):
 coords,tris,nodes=mesh;area=inventory=maxedge=0
 for material,vs in tris:
  a,b,c=[coords[v] for v in vs];signed=((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))/2
  assert signed>0,'nonpositive element'
  if material!=mat:continue
  area+=signed;inventory+=signed*sum(nodes[v,mat].get(species,0) for v in vs)/3
  for ia,ib in zip(vs,vs[1:]+vs[:1]):maxedge=max(maxedge,math.dist(coords[ia],coords[ib]))
 return dict(area_um2=area,inventory_per_cm=inventory*1e-8,max_edge_um=maxedge)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args();binary=args.binary.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
 env={**os.environ,'ASAN_OPTIONS':'detect_leaks=0:halt_on_error=1','UBSAN_OPTIONS':'halt_on_error=1'}
 for key,name in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[key]=str(binary.parent.parent/'data'/name)
 results=[]
 for name,h,step in [('flat20',.02,None),('flat10',.01,None),('flat5',.005,None),('step_left',.01,.975),('step_center',.01,1.),('step_right',.01,1.025)]:
  d=out/name;d.mkdir();deck='''option quiet
line x loc=0 spac=0.2 tag=left
line x loc=2 spac=0.2 tag=right
line y loc=0 spac=0.1 tag=top
line y loc=2 spac=0.2 tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
boundary exposed xlo=left xhi=right ylo=top yhi=top
boundary backside xlo=left xhi=right ylo=bottom yhi=bottom
init boron conc=1e15 orient=100
structure outfile=initial.str
'''
  if step is not None:deck+=f'deposit oxide thick=0.15 divisions=3 space=0.05\netch oxide left p1.x={step}\n'
  deck+=f'deposit poly thick=0.5 divisions=10 space={h}\nstructure outfile=before.str\nimplant phosphorus dose=2e16 energy=50 pearson\nstructure outfile=after.str\nquit\n'
  (d/'input.in').write_text(deck);start=time.monotonic()
  with (d/'run.log').open('w') as log:r=subprocess.run([str(binary)],input=deck,text=True,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600)
  assert r.returncode==0,(name,r.returncode)
  log=(d/'run.log').read_text();assert not re.search(r'ERROR:|panic:|AddressSanitizer:|runtime error:',log),name
  before=read(d/'before.str');after=read(d/'after.str');initial=read(d/'initial.str');pm=metrics(before,4,3);sm=metrics(before,3,5);im=metrics(initial,3,5)
  assert pm['max_edge_um']<=math.sqrt(2)*h*(1+1e-5),(name,pm)
  assert abs(sm['area_um2']/im['area_um2']-1)<1e-6
  assert abs(sm['inventory_per_cm']/im['inventory_per_cm']-1)<1e-6
  cuts={str(x):dict(zip(('sheet_cm2','peak_cm3'),cut(after,x,4,3))) for x in (.25,.5,1.5,1.75)}
  # Flat cuts away from the rounded step must recover the incident dose.
  assert all(abs(c['sheet_cm2']/2e16-1)<.015 for c in cuts.values()),(name,cuts)
  record=dict(name=name,h_um=h,step_x_um=step,elapsed_seconds=time.monotonic()-start,poly=pm,substrate=sm,cuts=cuts,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest());results.append(record);(out/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(record,flush=True)
 # Differences between actual h and h/2 are checked, rather than node counts.
 for x in ('0.25','0.5','1.5','1.75'):
  assert abs(results[1]['cuts'][x]['sheet_cm2']/results[2]['cuts'][x]['sheet_cm2']-1)<.01
 print('PASS: achieved spacing, positive elements, substrate inventory, profile convergence and shifted steps')
if __name__=='__main__':main()
