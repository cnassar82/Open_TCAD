#!/usr/bin/env python3
"""Native-oxide strip construction on a complete exposed boundary loop."""
import argparse,hashlib,json,os,re,subprocess
from pathlib import Path
from partition_audit import audit_partition
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);data=binary.parent.parent/'data';env=os.environ.copy()
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
def run(name,deck):
 (out/(name+'.in')).write_text(deck)
 with (out/(name+'.in')).open('rb') as f,(out/(name+'.log')).open('wb') as log:
  r=subprocess.run([str(binary)],cwd=out,env=env,stdin=f,stdout=log,stderr=subprocess.STDOUT,timeout=20)
 if r.returncode or re.search(r'errors detected|ERROR:|runtime error:|AddressSanitizer:|suprem4 panic:',(out/(name+'.log')).read_text(errors='replace'),re.I):raise ValueError(name+' simulator failure')
r={'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'errors':[]}
try:
 run('initial','''option quiet
line x loc=0 spac=0.1 tag=left
line x loc=1 spac=0.1 tag=right
line y loc=0 spac=0.1 tag=top
line y loc=1 spac=0.1 tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
boundary exposed xlo=left xhi=right ylo=top yhi=top
initialize boron conc=1e15 orient=100
structure outfile=initial.str
quit
''')
 lines=[]
 for line in (out/'initial.str').read_text().splitlines():
  s=line.split()
  if s and s[0]=='t':
   for i in range(6,9):
    if int(s[i])<0:s[i]='-1022'
   line=' '.join(s)
  lines.append(line)
 (out/'island.str').write_text('\n'.join(lines)+'\n')
 # A tiny physical interval isolates seeding from significant oxide growth.
 run('seed','option quiet\nmethod fermi vertical grid.oxide=0.04\nstructure infile=island.str\ndiffuse time=0.000000001 temp=825 dryo2\nstructure outfile=after.str\nquit\n')
 r['partition']=audit_partition(out/'after.str')
 # The silicon has no remaining exposed face; every side received oxide.
 regions={int(s[1]):int(s[2]) for line in (out/'after.str').read_text().splitlines() if (s:=line.split()) and s[0]=='r'}
 oxide=0
 for line in (out/'after.str').read_text().splitlines():
  s=line.split()
  if s and s[0]=='t':
   if regions[int(s[2])]==3 and '-1022' in s[6:9]:raise ValueError('uncoated silicon edge')
   if regions[int(s[2])]==1:oxide+=1
 if not oxide:raise ValueError('no oxide seeded')
 r['oxide_triangles']=oxide
except (OSError,ValueError,KeyError,subprocess.TimeoutExpired) as e:r['errors'].append(str(e))
r['passed']=not r['errors'];(out/'result.json').write_text(json.dumps(r,indent=2)+'\n');print('PASS' if r['passed'] else 'FAIL',r['errors']);raise SystemExit(not r['passed'])
