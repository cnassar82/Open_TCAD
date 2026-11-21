#!/usr/bin/env python3
"""Probe the saved PMOS failure with the original stressed oxidation law."""
import argparse,hashlib,json,os,re,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
p.add_argument('--initial',type=Path,default=HERE/'pmos_post_viscous/pmos/panic.str')
p.add_argument('--jcheck',action='store_true');p.add_argument('--staged',action='store_true');p.add_argument('--trust',action='store_true')
p.add_argument('--timeout',type=float,default=120)
a=p.parse_args();b=a.binary.resolve();d=a.output.resolve();d.mkdir(exist_ok=False);initial=a.initial.resolve()
if a.jcheck and a.trust:p.error('this diagnostic is in the Newton path')
if a.staged and 'getenv("SUPREM4_FE_STAGED_HOMOTOPY")' not in (b.parent/'oxide/viscous.c').read_text():p.error('candidate does not support the staged option')
deck=f'option chat\noxide Vc=300 Vr=30 Vd=25 stress.dep=true\nmethod fermi viscous oxide.rel=1e-2 grid.oxide=0.04\nstructure infile={initial}\ndiffuse time=0.00166666666666667 temp=900 dryo2\nstructure outfile=oxidized.str\nquit\n';(d/'input.in').write_text(deck)
env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','SUPREM4_PROCESS_SURFACE':'1','SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide','SUPREM4_MESH_PATCH_REMOVAL':'1','SUPREM4_FE_OXIDANT_P1':'1','SUPREM4_FE_BBAR_P1':'1','SUPREM4_FE_LINEAR_TRACE':'1'}
for flag,on in [('SUPREM4_FE_JCHECK',a.jcheck),('SUPREM4_FE_STAGED_HOMOTOPY',a.staged),('SUPREM4_FE_TRUST_REGION',a.trust)]:
 env.pop(flag,None)
 if on:env[flag]='stop' if flag=='SUPREM4_FE_JCHECK' else '1'
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(b.parent.parent/'data'/n)
h=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
r={'accepted_full_cmos':False,'scope':'saved PMOS 0.1 second stressed oxidation probe','binary':str(b),'binary_sha256':h(b),'input_sha256':h(d/'input.in'),'initial_sha256':h(initial),'runtime_options':{k:v for k,v in env.items() if k.startswith('SUPREM4_')},'data_sha256':{n:h(b.parent.parent/'data'/n) for n in ('modelrc','suprem.uk','sup4gs.imp')}};t=time.monotonic()
with (d/'run.log').open('w') as log:
 try:r['exit']=subprocess.run([str(b)],input=deck,text=True,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=a.timeout).returncode
 except subprocess.TimeoutExpired:r['exit']=124
log=(d/'run.log').read_text(errors='replace');r['elapsed_seconds']=time.monotonic()-t
r['diagnostic_completed']=a.jcheck and 'stalled Jacobian diagnostic complete' in log
r['execution_completed']=r['exit']==0 and (d/'oxidized.str').is_file() and not re.search('Giving up|panic:|ERROR:',log)
(d/'execution.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
raise SystemExit(not (r['diagnostic_completed'] if a.jcheck else r['execution_completed']))
