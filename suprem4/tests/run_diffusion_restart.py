#!/usr/bin/env python3
"""Independent diffusion calls must retain the method's configured first step."""
import argparse,hashlib,json,os,re,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path)
a=p.parse_args();binary=a.binary.resolve();data=binary.parent.parent/'data'
work=(a.output or Path(tempfile.mkdtemp(prefix='rit_diffusion_restart_'))).resolve();work.mkdir(parents=True,exist_ok=True)
if any(work.iterdir()):p.error('populated output')
text='''option quiet
mode one.dim
line x loc=0 spac=0.01 tag=top
line x loc=1 spac=0.1 tag=bottom
region silicon xlo=top xhi=bottom
boundary exposed xlo=top xhi=top
boundary backside xlo=bottom xhi=bottom
init boron conc=1e15 orient=100
implant phosphorus dose=1e13 energy=60 pearson
method fermi init.time=0.125
diffuse time=0.01 temp=800 nitrogen
diffuse time=0.01 temp=800 nitrogen continue
diffuse time=0.01 temp=800 nitrogen
structure outfile=after.str
quit
'''
(work/'restart.in').write_text(text);env=os.environ.copy()
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
with (work/'restart.in').open('rb') as stdin,(work/'restart.log').open('wb') as log:
 try:code=subprocess.run([str(binary)],stdin=stdin,stdout=log,stderr=subprocess.STDOUT,cwd=work,env=env,timeout=30).returncode
 except subprocess.TimeoutExpired:code=None
log=(work/'restart.log').read_text(errors='replace')
steps=list(map(float,re.findall(r'Solving\s+0\s+\+\s+([0-9.eE+-]+)',log)))
report={'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
        'returncode':code,'independent_first_steps_seconds':steps,'expected_seconds':[.125,.125]}
report['passed']=code==0 and steps==[.125,.125] and not re.search(r'ERROR:|runtime error:|AddressSanitizer:|suprem4 panic:',log,re.I)
(work/'result.json').write_text(json.dumps(report,indent=2)+'\n')
print(f"{'PASS' if report['passed'] else 'FAIL'} diffusion restart: {steps}; {work}")
raise SystemExit(int(not report['passed']))
