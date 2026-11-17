#!/usr/bin/env python3
"""A failed implant must preserve geometry and dopant values for later commands."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
a=p.parse_args(); binary=a.binary.resolve(); data=binary.parent.parent/'data'
work=Path(tempfile.mkdtemp(prefix='suprem4_implant_failure_'))
text='''option quiet
line x loc=0 spac=0.25 tag=left
line x loc=1 spac=0.25 tag=right
line y loc=0 spac=0.05 tag=top
line y loc=1 spac=0.2 tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
initialize boron conc=1e15 orient=100
structure outfile=before.str
implant boron dose=1e13 energy=60 pearson
structure outfile=after.str
quit
'''
env=os.environ.copy()
for key,name in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc')]:
 env[key]=str(data/name)
env['IMPFILE']=str(work/'deliberately_missing.imp')
r=subprocess.run([str(binary)],input=text,text=True,cwd=work,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
(work/'suprem4.log').write_text(r.stdout)
assert 'ERROR: implant aborted' in r.stdout, r.stdout
assert 'panic' not in r.stdout.lower(),r.stdout
# imp_select may add an active species; compare the geometry and original B value.
def read(name):
 coords={};values={};species=[];tri=[]
 for l in (work/name).read_text().splitlines():
  s=l.split()
  if not s:continue
  if s[0]=='c':coords[s[1]]=s[2:]
  elif s[0]=='t':tri.append(s[1:])
  elif s[0]=='s':species=list(map(int,s[2:]))
  elif s[0]=='n':values[tuple(s[1:3])]=float(s[3+species.index(5)])
 return coords,tri,values
assert read('before.str')==read('after.str'),f'failed implant changed mesh or boron: {work}'
print(f'PASS: failed implant preserves coordinates, connectivity and boron; {work}')
