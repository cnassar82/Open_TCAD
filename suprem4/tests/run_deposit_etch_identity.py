#!/usr/bin/env python3
"""Unmet construction gate: a planar deposit/equal dry etch must cancel.

This intentionally fails engines which manufacture a residual oxide film.
No tolerance deletes a finite material region to satisfy the test.
"""
import argparse, hashlib, json, os, subprocess
from pathlib import Path
from run_oxidation_perturbations import audit
from partition_audit import audit_partition

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve()
out.mkdir(parents=True,exist_ok=False);data=binary.parent.parent/'data'
env={**os.environ,'SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide'}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
results=[]
for label,y in [('zero',0.),('shift_down',.1),('shift_up',-.427)]:
 d=out/label;d.mkdir()
 deck=f'''option quiet
line x loc=0 spac=0.5 tag=left
line x loc=12 spac=0.5 tag=right
line y loc={y:.17g} spac=0.1 tag=top
line y loc={y+2:.17g} spac=0.2 tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
boundary exposed xlo=left xhi=right ylo=top yhi=top
boundary backside xlo=left xhi=right ylo=bottom yhi=bottom
initialize boron conc=1e15 orient=100
structure outfile=initial.str
deposit oxide thick=0.5 divisions=10 space=0.04
structure outfile=deposited.str
etch oxide dry thick=0.5
structure outfile=etched.str
quit
'''
 (d/'input.in').write_text(deck)
 r={'case':label,'errors':[],'deck_sha256':hashlib.sha256(deck.encode()).hexdigest()}
 with (d/'input.in').open('rb') as stdin,(d/'simulator.log').open('wb') as log:
  try:r['returncode']=subprocess.run([str(binary)],cwd=d,env=env,stdin=stdin,stdout=log,stderr=subprocess.STDOUT,timeout=30).returncode
  except subprocess.TimeoutExpired:r['returncode']=None;r['errors'].append('timeout')
 try:
  r['partition']=audit_partition(d/'etched.str');r['mesh']=audit(d/'etched.str')
  materials={int(v[1]):int(v[2]) for line in (d/'etched.str').read_text().splitlines() if (v:=line.split()) and v[0]=='r'}
  r['residual_oxide_area_um2']=sum(area for reg,area in r['mesh']['area_by_region_um2'].items() if materials[reg]==1)
  if r['residual_oxide_area_um2']!=0:r['errors'].append('equal-depth planar construction leaves a finite oxide film')
 except (OSError,ValueError,KeyError) as error:r['errors'].append(str(error))
 r['passed']=r['returncode']==0 and not r['errors'];results.append(r)
 print(('PASS' if r['passed'] else 'FAIL'),label,r.get('residual_oxide_area_um2'))
report={'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'results':results,'passed':all(r['passed'] for r in results)}
(out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
raise SystemExit(not report['passed'])
