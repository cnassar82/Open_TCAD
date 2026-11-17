#!/usr/bin/env python3
"""Planar poly/oxide/silicon diagnostic using the unchanged poly thermal budget.

This is a same-engine dimensional/resolution comparison, not independent
physical calibration. Uniform 1e17 boron substitutes for the process well.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
from validate import read, concentration

HERE=Path(__file__).resolve().parent


def metrics(path, dimensions):
    c={};r={};nodes={};species=[];elements=[]
    for line in path.read_text().splitlines():
        p=line.split()
        if not p: continue
        if p[0]=='c':c[int(p[1])]=tuple(map(float,p[2:2+dimensions]))
        elif p[0]=='r':r[int(p[1])]=int(p[2])
        elif p[0]=='t':elements.append((r[int(p[2])],tuple(map(int,p[3:3+dimensions+1]))))
        elif p[0]=='s':species=list(map(int,p[2:]))
        elif p[0]=='n':nodes[int(p[1])+1,int(p[2])]=dict(zip(species,map(float,p[3:])))
    totals={}
    for mat,vs in elements:
        pts=[c[v] for v in vs]
        if dimensions==1: measure=abs(pts[1][0]-pts[0][0])*1e-4
        else:
            a,b,d=pts
            measure=abs((b[0]-a[0])*(d[1]-a[1])-(b[1]-a[1])*(d[0]-a[0]))*.5e-8
            measure/=1e-4  # one um width: normalize 2D inventory by width in cm
        values=[nodes[v,mat].get(3,0.) for v in vs]
        if not all(math.isfinite(v) and v>=0 for v in values):raise ValueError('invalid phosphorus')
        totals[mat]=totals.get(mat,0.)+measure*sum(values)/len(values)
    profile={}
    for depth in (.005,.01,.025,.05,.1,.25,.5):
        if dimensions==2:val=concentration(read(path),.5,depth)['phosphorus_cm3']
        else:
            val=None
            for mat,(a,b) in elements:
                if mat==3 and min(c[a][0],c[b][0])<=depth<=max(c[a][0],c[b][0]):
                    w=(depth-c[a][0])/(c[b][0]-c[a][0]);val=(1-w)*nodes[a,mat][3]+w*nodes[b,mat][3];break
            if val is None:raise ValueError('missing silicon probe')
        profile[str(depth)]=val
    return {'chemical_phosphorus_atoms_cm2_by_material':totals,
            'silicon_phosphorus_cm3_by_depth_um':profile,'elements':len(elements)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--dimensions',type=int,choices=[1,2],default=1)
    p.add_argument('--oxide-divisions',type=int,default=8)
    p.add_argument('--oxide-um',type=float,default=.004639937545)
    p.add_argument('--si-space',type=float,default=.001)
    p.add_argument('--time-step',type=float)
    p.add_argument('--output',type=Path)
    args=p.parse_args();binary=args.binary.resolve();data=binary.parent.parent/'data'
    if args.oxide_divisions<1 or args.oxide_um<=0 or args.si_space<=0 or (args.time_step is not None and args.time_step<=0):p.error('positive resolution/time controls required')
    work=(args.output or Path(tempfile.mkdtemp(prefix='rit_poly_stack_'))).resolve()
    work.mkdir(parents=True,exist_ok=True)
    if any(work.iterdir()):p.error('populated output')
    axis='x' if args.dimensions==1 else 'y'
    lines=['option quiet']
    if args.dimensions==1:lines+=['mode one.dim']
    else:lines+=['line x loc=0 spac=0.1 tag=left','line x loc=1 spac=0.1 tag=right']
    lines += [f'line {axis} loc=0 spac={args.si_space} tag=top',f'line {axis} loc=0.1 spac=0.005',
              f'line {axis} loc=0.5 spac=0.02',f'line {axis} loc=8 spac=0.25 tag=bottom']
    if args.dimensions==1:
        lines+=['region silicon xlo=top xhi=bottom','boundary exposed xlo=top xhi=top','boundary backside xlo=bottom xhi=bottom']
    else:lines+=['region silicon xlo=left xhi=right ylo=top yhi=bottom','boundary exposed xlo=left xhi=right ylo=top yhi=top']
    lines+=['init boron conc=1e17 orient=100','method fermi'+(f' formula="{args.time_step}"' if args.time_step else ''),
            f'deposit oxide thick={args.oxide_um} divisions={args.oxide_divisions}',
            'deposit poly thick=0.5 divisions=40','implant phosphorus dose=2e16 energy=50 pearson','structure outfile=implanted.str']
    thermal=(HERE/'nmos.in').read_text().split('implant phosphorus dose=2e+16 energy=50 pearson',1)[1].split('structure outfile=poly_doped.str',1)[0]
    text='\n'.join(lines)+'\n'+thermal+'structure outfile=activated.str\nquit\n'
    deck=work/'stack.in';deck.write_text(text)
    digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    env=os.environ.copy()
    for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
    report={'scope':'same-engine planar diagnostic','controls':{k:v for k,v in vars(args).items() if k not in ('binary','output')},
            'binary':str(binary),'binary_sha256':digest(binary),'deck_sha256':digest(deck),
            'data_sha256':{n:digest(data/n) for n in ('modelrc','suprem.uk','sup4gs.imp')}}
    print(f'workdir: {work}',flush=True)
    with deck.open('rb') as stdin,(work/'stack.log').open('wb') as log:
        try:r=subprocess.run([str(binary)],stdin=stdin,stdout=log,stderr=subprocess.STDOUT,cwd=work,env=env,timeout=300);report['returncode']=r.returncode
        except subprocess.TimeoutExpired:report['returncode']=None
    report['errors']=[]
    if re.search(r'ERROR:|runtime error:|AddressSanitizer:|suprem4 panic:|Newton convergence failed|Giving up', (work/'stack.log').read_text(errors='replace'),re.I):report['errors'].append('simulator error')
    for stage in ('implanted','activated'):
        try:report[stage]=metrics(work/(stage+'.str'),args.dimensions)
        except (OSError,ValueError,KeyError) as error:report['errors'].append(f'{stage}: {error}')
    report['passed']=report['returncode']==0 and not report['errors']
    (work/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS diagnostic' if report['passed'] else str(report['errors']))
    return int(not report['passed'])


if __name__=='__main__':raise SystemExit(main())
