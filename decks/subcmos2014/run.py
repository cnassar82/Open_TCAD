#!/usr/bin/env python3
"""Run and audit the PDF-derived NMOS/PMOS process decks in isolated folders."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import shutil
import subprocess
import sys
import tempfile
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'suprem4/tests'))
sys.path.insert(0,str(ROOT/'suprem4/tools'))
from run_oxidation_perturbations import audit
from partition_audit import audit_partition
from inventory_audit import audit_chemistry
from str_compare import parse_str,summarize
from build_decks import CHECKPOINTS
from validate import validate


def audit_checkpoint(path, device, stage):
    return {'audit': audit(path), 'summary': summarize(parse_str(path)),
            'partition': audit_partition(path), 'chemistry': audit_chemistry(path),
            'device_checks': validate(path, device, stage)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True,help='native64 candidate produced by build_native64.py')
    p.add_argument('--devices',nargs='+',choices=['nmos','pmos'],default=['nmos','pmos'])
    p.add_argument('--through',choices=CHECKPOINTS,default='metal1')
    p.add_argument('--output',type=Path)
    p.add_argument('--timeout',type=float,default=1800)
    p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--resume',type=Path,help='previous run folder, containing device subfolders')
    p.add_argument('--after',choices=CHECKPOINTS,help='resume after this saved checkpoint')
    args=p.parse_args()
    if bool(args.resume)!=bool(args.after):p.error('--resume and --after must be supplied together')
    if args.after and CHECKPOINTS.index(args.after)>=CHECKPOINTS.index(args.through):
        p.error('--through must follow --after')
    binary=args.binary.resolve()
    if not binary.is_file():p.error(f'missing executable: {binary}')
    data=binary.parent.parent/'data'
    required_data={'suprem.uk','modelrc','sup4gs.imp'}
    missing=sorted(name for name in required_data if not (data/name).is_file())
    if missing:p.error(f'missing candidate data files: {missing} in {data}')
    build_manifest=binary.parent.parent/'build.json'
    provenance={
        'runtime_options': {key: os.environ.get(key) for key in
                            ('SUPREM4_MATERIAL_DOMAINS', 'SUPREM4_MESH_PATCH_REMOVAL',
                             'SUPREM4_MASS_LUMPING', 'SUPREM4_TIME_INTEGRATOR', 'SUPREM4_POSITIVITY_STEP',
                             'SUPREM4_DEPOSITION_GEOMETRY','SUPREM4_VERTICAL_CONSTRAINTS','SUPREM4_PROCESS_SURFACE')},
        'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
        'data_sha256':{name:hashlib.sha256((data/name).read_bytes()).hexdigest()
                       for name in sorted(required_data)},
        'build_manifest_sha256':(hashlib.sha256(build_manifest.read_bytes()).hexdigest()
                                 if build_manifest.is_file() else None),
    }
    work=(args.output or Path(tempfile.mkdtemp(prefix='suprem4_subcmos2014_'))).resolve()
    work.mkdir(parents=True,exist_ok=True)
    print(f'workdir: {work}\nbinary: {binary}',flush=True)
    expected=CHECKPOINTS[:CHECKPOINTS.index(args.through)+1]
    def run(device):
        d=work/device
        if d.exists() and any(d.iterdir()):raise ValueError(f'output folder is not empty: {d}')
        d.mkdir(exist_ok=True)
        shutil.copy2(HERE/'resist_top.dat',d/'resist_top.dat')
        text=(HERE/(device+'.in')).read_text()
        if args.through!='metal1':
            marker=f'structure outfile={args.through}.str'
            text=text.split(marker)[0]+marker+'\nquit\n'
        (d/'process.in').write_text(text)
        if args.resume:
            source=args.resume.resolve()/device
            marker=f'structure outfile={args.after}.str'
            before,remaining=text.split(marker,1)
            previous_file=source/'process.in'
            if not previous_file.exists():previous_file=source/(device+'.in')
            previous=previous_file.read_text().split(marker)[0]
            if previous!=before:raise ValueError(f'resume process prefix differs: {source}')
            for name in CHECKPOINTS[:CHECKPOINTS.index(args.after)+1]:
                saved=source/(name+'.str');audit(saved)
                shutil.copy2(saved,d/saved.name)
            settings=[line for line in before.splitlines()
                      if line.startswith(('method ','oxide ','option '))]
            text='\n'.join(settings)+f'\nstructure infile={args.after}.str\n'+remaining.lstrip()
        deck=d/(device+'.in');deck.write_text(text)
        env=os.environ.copy()
        for key,name in [('SUP4KEYFILE','suprem.uk'),('KEYFILE','suprem.uk'),
                         ('SUP4MODELRC','modelrc'),('SUP4IMPDATA','sup4gs.imp'),
                         ('IMPFILE','sup4gs.imp')]:env[key]=str(data/name)
        r={'device':device,'binary':str(binary),'through':args.through,
           'input_sha256':hashlib.sha256(text.encode()).hexdigest(),
           'provenance':provenance,
           'timed_out':False,'checkpoints':{},'errors':[]}
        if args.resume:r['resumed_from']={'folder':str(args.resume.resolve()),'after':args.after}
        t=time.monotonic()
        with deck.open('rb') as stdin,(d/'suprem4.log').open('wb') as log:
            proc=subprocess.Popen([str(binary)],stdin=stdin,stdout=log,stderr=subprocess.STDOUT,
                                  cwd=d,env=env,start_new_session=True)
            try:proc.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid,signal.SIGKILL);proc.wait();r['timed_out']=True
                r['errors'].append('runtime limit exceeded')
        r['elapsed_seconds']=time.monotonic()-t;r['returncode']=proc.returncode
        log=(d/'suprem4.log').read_text(errors='replace')
        if re.search(r'errors detected on command input|suprem4 panic:|No mesh defined|not a valid command|ERROR:|bad initial grid|Newton convergence failed|runtime error:|AddressSanitizer:|Giving up',log,re.I):
            r['errors'].append('simulator reported a command or mesh error')
        for name in expected:
            try:r['checkpoints'][name]=audit_checkpoint(d/(name+'.str'),device,name)
            except (OSError,ValueError,KeyError) as error:r['errors'].append(f'{name}: {error}')
        r['passed']=proc.returncode==0 and not r['errors']
        (d/'result.json').write_text(json.dumps(r,indent=2)+'\n')
        print(f"{'PASS' if r['passed'] else 'FAIL'} {device}: exit={proc.returncode}, "
              f"checkpoints={len(r['checkpoints'])}/{len(expected)}, {r['errors'][:2]}",flush=True)
        return r
    with ThreadPoolExecutor(max_workers=max(1,args.jobs)) as pool:results=list(pool.map(run,args.devices))
    (work/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    return int(not all(r['passed'] for r in results))

if __name__=='__main__':raise SystemExit(main())
