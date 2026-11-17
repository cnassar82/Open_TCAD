#!/usr/bin/env python3
"""Construction geometry: cancellation, retained films, restart and invalidation.

All etches use the stated depth; no area-based material deletion is allowed.
The sloped control must retain physical oxide even for equal nominal depths.
"""
import argparse, hashlib, json, os, re, subprocess, struct
from pathlib import Path
from run_oxidation_perturbations import audit
from partition_audit import audit_partition
from inventory_audit import inventories

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
data=binary.parent.parent/'data'
env={**os.environ,'SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide','SUPREM4_PROCESS_SURFACE':'1'}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
float32=lambda x:struct.unpack('f',struct.pack('f',x))[0]

def run(d,name,deck):
 (d/(name+'.in')).write_text(deck)
 with (d/(name+'.in')).open('rb') as stdin,(d/(name+'.log')).open('wb') as log:
  result=subprocess.run([str(binary)],cwd=d,env=env,stdin=stdin,stdout=log,stderr=subprocess.STDOUT,timeout=30)
 text=(d/(name+'.log')).read_text(errors='replace')
 if result.returncode or re.search(r'errors detected|ERROR:|Input mesh error|AddressSanitizer:|runtime error:|suprem4 panic:',text,re.I):
  raise ValueError(f'{name} simulator error: exit={result.returncode}')
 return text

def initial(y,width):
 return f'''option quiet
line x loc=0 spac=0.5 tag=left
line x loc={width} spac=0.5 tag=right
line y loc={y:.17g} spac=0.1 tag=top
line y loc={y+2:.17g} spac=0.2 tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
boundary exposed xlo=left xhi=right ylo=top yhi=top
boundary backside xlo=left xhi=right ylo=bottom yhi=bottom
initialize boron conc=1e15 orient=100
structure outfile=initial.str
'''

cases=[]
for y in (0.,.1,-.427):
 for width in (12,100):
  for depth in (.5,.49,.49999,.51):cases.append((y,width,depth,True,False))
for y in (0.,.1,-.427):cases.append((y,12,.5,False,False))
cases.append((-.427,12,.5,True,True))
results=[]
for index,(y,width,depth,restart,stale) in enumerate(cases):
 d=out/f'case_{index:02}';d.mkdir();r={'y_um':y,'width_um':width,'etch_um':depth,'restart':restart,'stale_control':stale,'errors':[]}
 try:
  deposit='deposit oxide thick=0.5 divisions=10 space=0.04\nstructure outfile=deposited.str\n'
  etch=f'etch oxide dry thick={depth}\nstructure outfile=etched.str\nquit\n'
  if restart:
   run(d,'deposit',initial(y,width)+deposit+'quit\n')
   if stale:
    # Change the actual surface without altering its construction descriptor.
    lines=[]
    for line in (d/'deposited.str').read_text().splitlines():
     s=line.split()
     if s and s[0]=='c':s[3]=str(float(s[3])+.02);line=' '.join(s)
     lines.append(line)
    (d/'changed.str').write_text('\n'.join(lines)+'\n');source='changed.str'
   else:source='deposited.str'
   log=run(d,'etch',f'option quiet\nstructure infile={source}\n'+etch)
  else:log=run(d,'direct',initial(y,width)+deposit+etch)
  used='dry etch uses construction' in log;r['used_construction']=used
  if used==stale:raise ValueError('construction validity guard did not select the expected path')
  r['partition']=audit_partition(d/'etched.str');r['mesh']=audit(d/'etched.str')
  regions={int(s[1]):int(s[2]) for line in (d/'etched.str').read_text().splitlines() if (s:=line.split()) and s[0]=='r'}
  area=sum(x for k,x in r['mesh']['area_by_region_um2'].items() if regions[k]==1);r['oxide_area_um2']=area
  if not stale:
   expected=width*max(0.,float32(.5)-float32(depth));r['expected_area_um2']=expected
   if expected==0 and area!=0:raise ValueError('a canceled planar layer still exists')
   if expected>0 and (area<=0 or abs(area-expected)>2e-5):raise ValueError('real retained oxide was lost or distorted')
  before=inventories(d/('changed.str' if stale else 'initial.str'),species=(5,))
  after=inventories(d/'etched.str',species=(5,))
  error=abs(after[3,5]-before[3,5])/before[3,5];r['silicon_boron_relative_error']=error
  if error>1e-12:raise ValueError('etch changed substrate chemical inventory')
 except (OSError,ValueError,KeyError,subprocess.TimeoutExpired) as e:r['errors'].append(str(e))
 r['passed']=not r['errors'];results.append(r);print('PASS' if r['passed'] else 'FAIL',index,r['errors'],flush=True)
# A true sloped conformal film has greater vertical thickness than its normal
# deposition thickness. Equal-depth directional etching must preserve it.
d=out/'sloped_control';d.mkdir();r={'case':'sloped_control','errors':[]}
try:
 run(d,'initial',initial(0.,12)+'quit\n');lines=[]
 for line in (d/'initial.str').read_text().splitlines():
  s=line.split()
  if s and s[0]=='c':s[3]=str(float(s[3])+.1*float(s[2]));line=' '.join(s)
  lines.append(line)
 (d/'sloped.str').write_text('\n'.join(lines)+'\n')
 run(d,'sloped','option quiet\nstructure infile=sloped.str\ndeposit oxide thick=0.5 divisions=10 space=0.04\netch oxide dry thick=0.5\nstructure outfile=etched.str\nquit\n')
 r['partition']=audit_partition(d/'etched.str')
 import sys
 sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'decks/subcmos2014'))
 from validate import read,intervals
 spans=intervals(read(d/'etched.str'),6,1);thickness=sum(b-a for a,b in spans)
 expected=.5*((1+.1**2)**.5-1);r.update(thickness_um=thickness,expected_thickness_um=expected)
 if thickness<=0 or abs(thickness-expected)>2e-6:raise ValueError('sloped physical residual was erased or distorted')
except (OSError,ValueError,KeyError,subprocess.TimeoutExpired) as e:r['errors'].append(str(e))
r['passed']=not r['errors'];results.append(r);print('PASS' if r['passed'] else 'FAIL','sloped',r['errors'])
report={'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'runtime_options':{k:env[k] for k in ('SUPREM4_PROCESS_SURFACE','SUPREM4_DEPOSITION_GEOMETRY','SUPREM4_MATERIAL_DOMAINS')},'results':results,'passed':all(r['passed'] for r in results)}
(out/'results.json').write_text(json.dumps(report,indent=2)+'\n');raise SystemExit(not report['passed'])
