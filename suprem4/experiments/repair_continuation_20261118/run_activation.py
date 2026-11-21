#!/usr/bin/env python3
"""Replay the captured poly implant with an unchanged activation recipe."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[3]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--initial',type=Path)
p.add_argument('--phase',choices=['first','full','remaining'],default='first')
p.add_argument('--first-minutes',type=float,default=12,
               help='bounded diagnostic duration for the first 800 C hold')
p.add_argument('--timeout',type=float,default=1800)
args=p.parse_args();out=args.output.resolve();out.mkdir(exist_ok=False)
binary=args.binary.resolve()
if args.phase=='remaining' and args.initial is None:p.error('--phase remaining requires --initial from the completed first 12-minute hold')
initial=(args.initial or ROOT/'suprem4/experiments/repair_continuation_20261118/cmos_mixed_trace/nmos/poly_implanted.str').resolve()
if args.first_minutes<=0:p.error('--first-minutes must be positive')
thermal=f'diffuse time={args.first_minutes:g} temp=800 nitrogen\n'
if args.phase in ('full','remaining'):
    thermal=(ROOT/'decks/subcmos2014/nmos.in').read_text().split('structure outfile=poly_implanted.str\n',1)[1].split('structure outfile=poly_doped.str',1)[0]
    if args.phase=='remaining':
        # A saved state does not serialize a multistep integrator's history.
        # Preserve every remaining hold/ramp, but initialize the first one
        # explicitly instead of requesting nonexistent continuation history.
        thermal=thermal.split('diffuse time=12 temp=800 nitrogen\n',1)[1]
        thermal=thermal.replace(' nitrogen continue',' nitrogen',1)
deck='option chat\nmethod fermi vertical grid.oxide=0.04\nstructure infile='+str(initial)+'\n'+thermal+'structure outfile=activated.str\nquit\n'
(out/'input.in').write_text(deck)
env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
for key,name in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[key]=str(binary.parent.parent/'data'/name)
digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
report={'binary':str(binary),'binary_sha256':digest(binary),'input_sha256':digest(out/'input.in'),
        'initial_sha256':digest(initial),'phase':args.phase,'first_minutes':args.first_minutes,'accepted_full_cmos':False,
        'runtime_options':{k:v for k,v in env.items() if k.startswith('SUPREM4_')},
        'data_sha256':{n:digest(binary.parent.parent/'data'/n) for n in ('modelrc','suprem.uk','sup4gs.imp')}}
report['checkpoint_restart']=args.phase=='remaining'
t=time.monotonic();report['timed_out']=False
with (out/'run.log').open('w') as log:
    proc=subprocess.Popen([str(binary)],stdin=subprocess.PIPE,stdout=log,stderr=subprocess.STDOUT,
                          cwd=out,env=env,text=True,start_new_session=True)
    try:proc.communicate(deck,timeout=args.timeout)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid,signal.SIGKILL);proc.communicate();report['timed_out']=True
report.update(exit=proc.returncode,elapsed_seconds=time.monotonic()-t)
report['reported_error']=bool(re.search(r'ERROR:|suprem4 panic:|runtime error:|AddressSanitizer:|Giving up',(out/'run.log').read_text(errors='replace')))
report['execution_completed']=proc.returncode==0 and not report['reported_error'] and (out/'activated.str').is_file()
(out/'execution.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ('exit','timed_out','elapsed_seconds','execution_completed')}),flush=True)
raise SystemExit(not report['execution_completed'])
