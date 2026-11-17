#!/usr/bin/env python3
"""Exercise etch polygon growth at continue/done and spacing-buffer ownership."""
import argparse, hashlib, json, os, re, subprocess
from pathlib import Path
from partition_audit import audit_partition
from run_oxidation_perturbations import audit
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
data=binary.parent.parent/'data';env=os.environ.copy()
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
# These masks lie outside the mesh. Their geometry must not alter the substrate.
# Each has explicit per-vertex spacing, including through each capacity change.
def polygon(count):
 points=[]
 for i in range(count):
  t=4*i/count
  if t<1:x,y=20+10*t,-3
  elif t<2:x,y=30,-3+(t-1)
  elif t<3:x,y=30-10*(t-2),-2
  else:x,y=20,-2-(t-3)
  points.append((x,y))
 return ''.join(f'etch oxide {"start" if i==0 else "done" if i==count-1 else "continue"} x={x:.17g} y={y:.17g} spac=0.2\n' for i,(x,y) in enumerate(points))
deck='''option quiet
line x loc=0 spac=0.5 tag=left
line x loc=12 spac=0.5 tag=right
line y loc=0 spac=0.1 tag=top
line y loc=2 spac=0.2 tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
boundary exposed xlo=left xhi=right ylo=top yhi=top
boundary backside xlo=left xhi=right ylo=bottom yhi=bottom
initialize boron conc=1e15 orient=100
deposit oxide thick=0.5 divisions=10 space=0.04
'''+polygon(501)+polygon(1003)+'''etch oxide dry thick=0.49
structure outfile=after.str
quit
'''
(out/'input.in').write_text(deck);r={'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'polygon_sizes':[501,1003],'errors':[]}
with (out/'input.in').open('rb') as stdin,(out/'run.log').open('wb') as log:
 try:r['returncode']=subprocess.run([str(binary)],cwd=out,env=env,stdin=stdin,stdout=log,stderr=subprocess.STDOUT,timeout=90).returncode
 except subprocess.TimeoutExpired:r['returncode']=None;r['errors'].append('timeout')
text=(out/'run.log').read_text(errors='replace')
if re.search(r'errors detected|ERROR:|runtime error:|AddressSanitizer:|suprem4 panic:',text,re.I):r['errors'].append('simulator error')
try:
 r['partition']=audit_partition(out/'after.str');r['mesh']=audit(out/'after.str')
 regions={int(v[1]):int(v[2]) for line in (out/'after.str').read_text().splitlines() if (v:=line.split()) and v[0]=='r'}
 r['oxide_area_um2']=sum(area for region,area in r['mesh']['area_by_region_um2'].items() if regions[region]==1)
 if abs(r['oxide_area_um2']-.12)>2e-5:r['errors'].append('unexpected oxide area after external masks and dry etch')
except (OSError,ValueError,KeyError) as e:r['errors'].append(str(e))
r['passed']=r['returncode']==0 and not r['errors'];(out/'result.json').write_text(json.dumps(r,indent=2)+'\n')
print('PASS' if r['passed'] else 'FAIL',r['errors']);raise SystemExit(not r['passed'])
