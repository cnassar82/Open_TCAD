#!/usr/bin/env python3
"""Spacer etch that removes exactly the deposited thickness.

A poly gate on a thin gate oxide receives an oxide deposit of thickness T,
then a dry etch of T. On every flat face the etch front coincides with the
pre-deposit surface, but only to binary32 rounding. The Boolean must treat
that as exact contact: no vertex pairs within rounding of each other, no
zero-thickness residue on flat faces, the gate oxide and spacers retained,
and silicon/poly unchanged. Thickness, origin and mesh density are varied.
"""
import argparse, json, math, os, re, subprocess, sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_oxidation_perturbations import audit
from partition_audit import audit_partition
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from str_compare import parse_str

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
data=binary.parent.parent/'data'
env={**os.environ,'SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide',
     'SUPREM4_PROCESS_SURFACE':'1'}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),
            ('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
GATE_OXIDE=.02;POLY=.3

def deck(x0,y0,spacer,space):
    gl,gr=x0+1.5,x0+2.5
    return f'''option quiet
line x loc={x0:.17g} spac={space} tag=left
line x loc={x0+4:.17g} spac={space} tag=right
line y loc={y0:.17g} spac=0.02 tag=top
line y loc={y0+1:.17g} spac=0.2 tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
boundary exposed xlo=left xhi=right ylo=top yhi=top
boundary backside xlo=left xhi=right ylo=bottom yhi=bottom
initialize boron conc=1e15 orient=100
deposit oxide thick={GATE_OXIDE} divisions=2
deposit poly thick={POLY} divisions=4
etch poly start x={x0-1:.17g} y={y0-5:.17g}
etch continue x={x0-1:.17g} y={y0+1:.17g}
etch continue x={gl:.17g} y={y0+1:.17g}
etch done x={gl:.17g} y={y0-5:.17g}
etch poly start x={gr:.17g} y={y0-5:.17g}
etch continue x={gr:.17g} y={y0+1:.17g}
etch continue x={x0+5:.17g} y={y0+1:.17g}
etch done x={x0+5:.17g} y={y0-5:.17g}
structure outfile=gate.str
deposit oxide thick={spacer} divisions=6 space=0.04
etch oxide dry thick={spacer}
structure outfile=etched.str
quit
'''

def areas(s):
    by=defaultdict(float)
    for t in s['triangles']:
        (ax,ay),(bx,by_),(cx,cy)=(s['coords'][v] for v in t['verts'])
        by[s['regions'][t['region']]]+=abs((bx-ax)*(cy-ay)-(by_-ay)*(cx-ax))/2
    return by

def material_at(s,x,y):
    for t in s['triangles']:
        P=[s['coords'][v] for v in t['verts']]
        d=[(P[(i+1)%3][0]-P[i][0])*(y-P[i][1])-(P[(i+1)%3][1]-P[i][1])*(x-P[i][0]) for i in range(3)]
        if all(v>=0 for v in d) or all(v<=0 for v in d):return s['regions'][t['region']]
    return 0

results=[];failed=0
for x0 in (0.,.123456):
    for y0 in (0.,-.427):
        for spacer in (.1,.25,.5):
            for space in (.1,.05):
                name=f'x{x0}_y{y0}_t{spacer}_s{space}';d=out/name;d.mkdir()
                r={'case':name,'errors':[]}
                try:
                    (d/'input.in').write_text(deck(x0,y0,spacer,space))
                    with (d/'input.in').open('rb') as stdin,(d/'run.log').open('wb') as log:
                        rc=subprocess.run([str(binary)],cwd=d,env=env,stdin=stdin,stdout=log,
                                          stderr=subprocess.STDOUT,timeout=300).returncode
                    text=(d/'run.log').read_text(errors='replace')
                    if rc or re.search(r'suprem4 panic:|errors detected|ERROR:|runtime error:',text):
                        raise ValueError(f'simulator error exit={rc}')
                    audit(d/'etched.str');audit_partition(d/'etched.str')
                    before=parse_str(d/'gate.str');after=parse_str(d/'etched.str')
                    A0,A1=areas(before),areas(after)
                    for m,label in ((3,'silicon'),(4,'poly')):
                        if abs(A1[m]-A0[m])>1e-9*max(1,A0[m]):raise ValueError(f'{label} area changed')
                    # No two vertices within storage rounding of each other.
                    pts=sorted(after['coords'].values())
                    close=0
                    for i in range(len(pts)):
                        for j in range(i+1,len(pts)):
                            if pts[j][0]-pts[i][0]>1e-5:break
                            if math.dist(pts[i],pts[j])<1e-5:close+=1
                    r['near_vertex_pairs']=close
                    if close:raise ValueError(f'{close} vertex pairs closer than 1e-5 um')
                    gl,gr=x0+1.5,x0+2.5;top=y0-GATE_OXIDE;ptop=top-POLY
                    # Gate oxide retained under the flat field; nothing above it.
                    for x in (x0+.2,x0+3.8):
                        if material_at(after,x,y0-GATE_OXIDE/2)!=1:raise ValueError('gate oxide lost')
                        if material_at(after,x,top-.002)!=0:raise ValueError('residue on flat field')
                    # Nothing left on the poly top.
                    if material_at(after,(gl+gr)/2,ptop-.002)!=0:raise ValueError('residue on gate top')
                    # Spacers retained against both sidewalls.
                    for x in (gl-.01,gr+.01):
                        if material_at(after,x,top-POLY/2)!=1:raise ValueError('spacer missing')
                    r['oxide_area_um2']=A1[1]
                except Exception as error:
                    r['errors'].append(str(error))
                ok=not r['errors'];failed+=not ok
                print(('PASS ' if ok else 'FAIL ')+name,r['errors'],flush=True);results.append(r)
(out/'results.json').write_text(json.dumps(results,indent=1)+'\n')
sys.exit(1 if failed else 0)
