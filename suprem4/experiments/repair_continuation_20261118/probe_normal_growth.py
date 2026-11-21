#!/usr/bin/env python3
"""Measure achieved through-oxide spacing during planar thermal growth."""
import argparse,hashlib,json,os,subprocess,time,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'decks/subcmos2014'))
from validate import read
from audit_replay import inspect
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
p.add_argument('--grid-um',type=float,required=True);p.add_argument('--normal-layers',action='store_true')
a=p.parse_args();out=a.output.resolve();out.mkdir(exist_ok=False);binary=a.binary.resolve()
if a.grid_um<=0:p.error('positive oxide spacing required')
deck='''option quiet
line x loc=0 spac=0.1 tag=left
line x loc=1 spac=0.1 tag=right
line y loc=0 spac=0.01 tag=top
line y loc=0.1 spac=0.02
line y loc=1 spac=0.1 tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
boundary exposed xlo=left xhi=right ylo=top yhi=top
init boron conc=1e15 orient=100
method fermi vertical grid.oxide=GRID
diffuse time=45 temp=900 dryo2
structure outfile=grown.str
quit
'''.replace('GRID',str(a.grid_um))
(out/'input.in').write_text(deck)
env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
env.pop('SUPREM4_NORMAL_LAYERS',None)
if a.normal_layers:env['SUPREM4_NORMAL_LAYERS']='1'
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(binary.parent.parent/'data'/n)
digest=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
r={'binary':str(binary),'binary_sha256':digest(binary),'input_sha256':digest(out/'input.in'),'requested_normal_spacing_nm':a.grid_um*1000,'runtime_options':{k:v for k,v in env.items() if k.startswith('SUPREM4_')},'accepted_full_cmos':False,
   'data_sha256':{n:digest(binary.parent.parent/'data'/n) for n in ('modelrc','suprem.uk','sup4gs.imp')}}
t=time.monotonic()
with (out/'run.log').open('w') as log:
 try:r['exit']=subprocess.run([str(binary)],input=deck,text=True,stdout=log,stderr=subprocess.STDOUT,cwd=out,env=env,timeout=180).returncode
 except subprocess.TimeoutExpired:r['exit']=124
r['elapsed_seconds']=time.monotonic()-t
r['spacing_reset_reported']='resetting' in (out/'run.log').read_text()
r['achieved_requested_edge_spacing']=False
r['audit_errors']=[]
if (out/'grown.str').exists():
 coords,tris,nodes=read(out/'grown.str');cuts=[]
 for x in [.025,.275,.525,.775,.975]:
  ys=[]
  for mat,vs in tris:
   if mat!=1:continue
   pts=[coords[v] for v in vs]
   for b,c in zip(pts,pts[1:]+pts[:1]):
    if b[0]!=c[0] and min(b[0],c[0])<=x<=max(b[0],c[0]):ys.append(b[1]+(c[1]-b[1])*(x-b[0])/(c[0]-b[0]))
  ys=sorted(set(ys));cuts.append({'x_um':x,'oxide_thickness_nm':1000*(ys[-1]-ys[0]),'cut_max_gap_nm':1000*max(v-u for u,v in zip(ys,ys[1:])),'intersections':len(ys)})
 r['cuts']=cuts;r['points']=len(coords)
 r['maximum_oxide_edge_normal_span_nm']=1000*max(abs(coords[u][1]-coords[v][1]) for mat,vs in tris if mat==1 for u,v in zip(vs,vs[1:]+vs[:1]))
 r['achieved_requested_edge_spacing']=not r['spacing_reset_reported'] and r['maximum_oxide_edge_normal_span_nm']<=r['requested_normal_spacing_nm']*(1+1e-5)
 try:r['independent_audit']=inspect(out/'grown.str')
 except (OSError,ValueError,KeyError) as error:r['audit_errors'].append(str(error))
r['passed']=r['exit']==0 and r['achieved_requested_edge_spacing'] and not r['audit_errors']
(out/'result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
raise SystemExit(not r['passed'])
