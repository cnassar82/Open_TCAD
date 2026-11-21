#!/usr/bin/env python3
"""Analytic oxide/poly/oxide stack at half, normal, and double mesh density."""
import argparse,hashlib,json,os,re,subprocess,sys
from pathlib import Path
from partition_audit import audit_partition
from run_oxidation_perturbations import audit
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'decks/subcmos2014'))
from validate import read,gate_dielectric,intervals
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
env={**os.environ,'SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide','SUPREM4_PROCESS_SURFACE':'1','SUPREM4_MESH_PATCH_REMOVAL':'1'}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(binary.parent.parent/'data'/n)
def areas(path):
 mats={int(v[1]):int(v[2]) for line in path.read_text().splitlines() if (v:=line.split()) and v[0]=='r'}
 result={}
 for r,area in audit(path)['area_by_region_um2'].items():result[mats[r]]=result.get(mats[r],0)+area
 return result
results=[]
for divisions in (5,10,20):
 d=out/str(divisions);d.mkdir();r={'divisions':divisions,'errors':[]}
 deck=f'''option quiet
line x loc=0 spac={2/divisions} tag=left
line x loc=10 spac={2/divisions} tag=right
line y loc=0 spac={.2/divisions} tag=top
line y loc=2 spac={1/divisions} tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
boundary exposed xlo=left xhi=right ylo=top yhi=top
boundary backside xlo=left xhi=right ylo=bottom yhi=bottom
initialize boron conc=1e15 orient=100
deposit oxide thick=0.02 divisions={divisions}
deposit poly thick=0.5 divisions={divisions}
etch poly left p1.x=4
etch poly right p1.x=6
deposit oxide thick=0.5 divisions={divisions}
structure outfile=before.str
etch oxide dry thick=1.5
structure outfile=after.str
quit
'''
 (d/'input.in').write_text(deck)
 try:
  with (d/'run.log').open('w') as log:
   r['returncode']=subprocess.run([str(binary)],input=deck,text=True,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=120).returncode
  if r['returncode'] or re.search(r'errors detected|ERROR:|runtime error:|AddressSanitizer:|suprem4 panic:',(d/'run.log').read_text(),re.I):raise ValueError('simulator error')
  r['partition']=audit_partition(d/'after.str');before=areas(d/'before.str');after=areas(d/'after.str');r['areas_um2']=after
  for mat in (3,4):
   if abs(after[mat]-before[mat])>1e-12*max(1,before[mat]):raise ValueError(f'non-target material {mat} changed')
  # Use the represented rectangle, including deposition's float projection,
  # as the exact area reference; the nominal value is 0.04 square microns.
  source=read(d/'before.str');coords,tris,_=source
  poly_x=[coords[v][0] for mat,vs in tris if mat==4 for v in vs]
  expected=(max(poly_x)-min(poly_x))*(intervals(source,5,3)[0][0]-intervals(source,5,4)[0][1])
  r['expected_buried_area_um2']=expected
  # Etch snap-rounds its temporary cutter: a mask edge may be routed through
  # an existing vertex that the Boolean's ON predicate cannot distinguish from
  # it (per axis 2*FLT_EPSILON*max|coord|) rather than creating a sliver
  # beside it. The area bound is that band times the cut height, per side.
  thickness=intervals(source,5,3)[0][0]-intervals(source,5,4)[0][1]
  band=2*2.0**-23*max(abs(x) for x in poly_x)
  tolerance=max(1e-12,2*band*thickness)
  r['buried_area_tolerance_um2']=tolerance
  if abs(after.get(1,0)-expected)>tolerance:raise ValueError(f'buried oxide area {after.get(1,0)} != {expected} um2 (tolerance {tolerance})')
  mesh=read(d/'after.str');r['gate_oxide']=gate_dielectric(mesh)
  for x in (3,7):
   if intervals(mesh,x,1):raise ValueError('unmasked oxide was not cleared')
 except (OSError,ValueError,KeyError,subprocess.TimeoutExpired) as e:r['errors'].append(str(e))
 r['passed']=not r['errors'];results.append(r);print('PASS' if r['passed'] else 'FAIL',divisions,r['errors'],flush=True)
report={'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'results':results,'passed':all(r['passed'] for r in results)}
(out/'results.json').write_text(json.dumps(report,indent=2)+'\n');raise SystemExit(not report['passed'])
