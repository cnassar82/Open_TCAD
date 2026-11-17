#!/usr/bin/env python3
"""Captured contact etch: buried film, open contacts, depth/origin perturbations."""
import argparse,gzip,hashlib,json,os,re,subprocess,sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from partition_audit import audit_partition
from run_oxidation_perturbations import audit
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'decks/subcmos2014'))
from validate import read,gate_dielectric,intervals
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
fixture=Path(__file__).parent/'fixtures/rit_nmos_contact_shielding'
raw=gzip.decompress((fixture/'before.str.gz').read_bytes())
assert hashlib.sha256(raw).hexdigest()==json.loads((fixture/'fixture.json').read_text())['sha256']
env=os.environ.copy();data=binary.parent.parent/'data'
for key,name in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[key]=str(data/name)
def run(case):
 label,shift,depth=case;d=out/label;d.mkdir();lines=[]
 for line in raw.decode().splitlines():
  v=line.split()
  if v and v[0]=='c':v[2]=format(float(v[2])+shift,'.17g');line=' '.join(v)
  lines.append(line)
 (d/'input.str').write_text('\n'.join(lines)+'\n')
 deck=f'option quiet\nstructure infile=input.str\netch oxide dry thick={depth}\netch photoresist all\nstructure outfile=after.str\nquit\n'
 (d/'input.in').write_text(deck)
 result={'case':label,'shift_um':shift,'etch_depth_um':depth,'errors':[]}
 try:
  with (d/'input.in').open('rb') as stdin,(d/'run.log').open('wb') as log:
   result['returncode']=subprocess.run([str(binary)],stdin=stdin,stdout=log,stderr=subprocess.STDOUT,cwd=d,env=env,timeout=180).returncode
  if result['returncode'] or re.search(r'errors detected|ERROR:|runtime error:|AddressSanitizer:|suprem4 panic:',(d/'run.log').read_text(),re.I):raise ValueError('simulator error')
  result['mesh']=audit(d/'after.str');result['partition']=audit_partition(d/'after.str')
  mesh=read(d/'after.str');result['gate_oxide']=gate_dielectric(mesh)
  result['open_contacts']={}
  for name,x,mat in [('body',1.25,3),('source',4.25,3),('gate',6,4),('drain',7.75,3)]:
   x+=shift;top=intervals(mesh,x,mat)[0][0]
   if any(lo<top and hi>lo for lo,hi in intervals(mesh,x,1)):
    raise ValueError(f'oxide blocks {name} contact')
   result['open_contacts'][name]={'x_um':x,'conductor_top_um':top}
 except (OSError,ValueError,KeyError,subprocess.TimeoutExpired) as e:result['errors'].append(str(e))
 result['passed']=not result['errors'];print('PASS' if result['passed'] else 'FAIL',label,result['errors'],flush=True)
 return result
cases=[('baseline',0,1.5),('left_2nm',-.002,1.5),('right_2nm',.002,1.5),('depth_minus_1pct',0,1.485),('depth_plus_1pct',0,1.515)]
with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,cases))
report={'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'fixture_sha256':hashlib.sha256(raw).hexdigest(),'results':results,'passed':all(r['passed'] for r in results)}
(out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
raise SystemExit(not report['passed'])
