#!/usr/bin/env python3
"""Native oxide on a silicon mesa between oxide fills of every height.

The bare mesa top joins oxide at both ends. A fill above the mesa top
(rising) is joined by noding its free boundary at the film offset. A flush
or falling fill (spacer top, recessed fill) never reaches that offset; the
film must end with a positive-width end face and must not coat the old
oxide's top. Every case must cover the whole mesa top, leave the old fills
unchanged and pass the topology and partition audits.
"""
import argparse
from collections import defaultdict
import os
from pathlib import Path
import subprocess
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_oxidation_perturbations import audit
from partition_audit import audit_partition, orientation
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from str_compare import parse_str

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();data=binary.parent.parent/'data'
work=Path(tempfile.mkdtemp(prefix='suprem4_native_end_face_'));env=os.environ.copy()
env['SUPREM4_NATIVE_JOIN_CAP']='1'
for key,name in [('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc')]:env[key]=str(data/name)

# y increases into the wafer; y=0 is the mesa top, exposed to gas.
CASES={'flush':(0.,0.),'falling':(.3,.3),'rising':(-.05,-.05),
       'flush_falling':(0.,.3),'rising_falling':(-.05,.3),'shallow_fall':(.0004,.0004)}

def write_structure(path,left,right):
    coords=[(0,.5),(1,.5),(2,.5),(3,.5),(0,1),(1,1),(2,1),(3,1),
            (1,0),(2,0),(0,left),(3,right)]
    triangles=[(1,(1,2,6)),(1,(1,6,5)),(1,(9,10,3)),(1,(9,3,2)),
               (1,(2,3,7)),(1,(2,7,6)),(1,(3,4,8)),(1,(3,8,7)),
               (2,(11,9,2)),(2,(11,2,1)),(3,(10,12,4)),(3,(10,4,3))]
    for _,vs in triangles:
        assert orientation(*(coords[v-1] for v in vs))>0,vs
    regions={1:3,2:1,3:1}
    text='v SUPREM-IV.GS B.9305\nD 2 3 3\n'
    for i,(x,y) in enumerate(coords,1):text+=f'c {i} {x:.17g} {y:.17g} 0\n'
    for r in sorted(regions):text+=f'r {r} {regions[r]}\n'
    owners=defaultdict(list)
    for i,(_,vs) in enumerate(triangles,1):
        for j in range(3):owners[tuple(sorted((vs[(j+1)%3],vs[(j+2)%3])))].append(i)
    nodes=set()
    for i,(r,vs) in enumerate(triangles,1):
        neighbors=[]
        for j in range(3):
            edge=tuple(sorted((vs[(j+1)%3],vs[(j+2)%3])));adj=owners[edge]
            if len(adj)==2:neighbor=next(t for t in adj if t!=i)
            elif all(coords[v-1][1]==1 for v in edge):neighbor=-1023
            elif all(coords[v-1][0]==0 for v in edge) or all(coords[v-1][0]==3 for v in edge):neighbor=-1024
            else:neighbor=-1022
            neighbors.append(neighbor)
        text+=f't {i} {r} '+' '.join(map(str,(*vs,*neighbors)))+'\n'
        nodes.update((v-1,regions[r]) for v in vs)
    text+='s 1 5\n'
    for point,material in sorted(nodes):text+=f'n {point} {material} '+('1e15' if material==3 else '1')+'\n'
    for i in (4,5,6,7):text+=f'n {i} 0 1e5\n'
    text+='M 100 1173\n';path.write_text(text)
    return coords,triangles

def oxide_triangles(path,result):
    regions={int(s[1]):int(s[2]) for line in path.read_text().splitlines()
             if (s:=line.split()) and s[0]=='r'}
    out=[]
    for line in path.read_text().splitlines():
        s=line.split()
        if s and s[0]=='t' and regions[int(s[2])]==1:
            out.append([result['coords'][int(v)] for v in s[3:6]])
    return out

def area(t):
    (ax,ay),(bx,by),(cx,cy)=t
    return abs((bx-ax)*(cy-ay)-(by-ay)*(cx-ax))/2

def covered(triangles,point):
    return any(all(orientation(a,b,point)>=0 for a,b in zip(t,t[1:]+t[:1])) for t in triangles)

for case,(left,right) in CASES.items():
    folder=work/case;folder.mkdir()
    coords,triangles=write_structure(folder/'input.str',left,right)
    audit(folder/'input.str')
    old_oxide=oxide_triangles(folder/'input.str',parse_str(folder/'input.str'))
    old_area=sum(area(t) for t in old_oxide)
    commands=('option quiet\nmethod fermi vertical\nstructure infile=input.str\n'
              'diffuse time=1e-9 temp=800 dryo2\nstructure outfile=seed.str\nquit\n')
    r=subprocess.run([str(binary)],input=commands,text=True,cwd=folder,env=env,
                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
    (folder/'suprem4.log').write_text(r.stdout)
    assert r.returncode==0 and 'panic' not in r.stdout.lower() and 'errors detected' not in r.stdout,f'{folder}: {r.stdout[-2000:]}'
    audit(folder/'seed.str');audit_partition(folder/'seed.str')
    result=parse_str(folder/'seed.str');oxide=oxide_triangles(folder/'seed.str',result)
    # The complete mesa top, including both junction ends, is coated.
    for x in (1.0002,1.001,1.25,1.5,1.75,1.999,1.9998):
        assert covered(oxide,(x,-.0001)),(case,'uncoated mesa top',x)
    # The old fills are still oxide: sample their interiors.
    for x,y in ((.5,(left+.5)/2+.01),(2.5,(right+.5)/2+.01)):
        assert covered(oxide,(x,y)),(case,'old fill lost',x,y)
    # A flush/falling fill is not coated across its top: the film ends at the junction.
    # Old fill tops run from the domain side to the mesa corner (y=0).
    for x,top,far in ((.5,left*.5,left),(2.5,right*.5,right)):
        if far>=0:
            assert not covered(oxide,(x,top-.0005)),(case,'old oxide top coated',x)
    new_area=sum(area(t) for t in oxide)
    film=new_area-old_area
    # One mesa width (1 um) of a ~1 nm film, plus at most its end faces.
    assert .0005<film<.003,(case,'unexpected native film area',film)
    print(f'PASS native end face {case}: film area {film:.3g} um^2, {len(result["triangles"])} positive triangles',flush=True)
