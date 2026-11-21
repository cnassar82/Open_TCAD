#!/usr/bin/env python3
"""Seed oxide on bare Si, beside an oxide wedge, and across a one-edge gap."""
import argparse
from collections import defaultdict
import os
from pathlib import Path
import subprocess
import struct
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_oxidation_perturbations import audit
from partition_audit import audit_partition, orientation
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from str_compare import parse_str
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();data=binary.parent.parent/'data'
work=Path(tempfile.mkdtemp(prefix='suprem4_native_join_'));env=os.environ.copy()
for key,name in [('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc')]:env[key]=str(data/name)
coords=[(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(2,1),(3,1),
        (0,-.05),(.99,-.0005),(2.01,-.0005),(3,-.05)]
for case in ('bare','wedge','gap'):
 folder=work/case;folder.mkdir()
 triangles=[(1,(1,2,6)),(1,(1,6,5)),(1,(2,3,7)),(1,(2,7,6)),
            (1,(3,4,8)),(1,(3,8,7))]
 if case!='bare':triangles.extend([(2,(9,10,1)),(2,(10,2,1))])
 if case=='gap':triangles.extend([(3,(3,11,4)),(3,(11,12,4))])
 points=8 if case=='bare' else 10 if case=='wedge' else 12
 text='v SUPREM-IV.GS B.9305\nD 2 3 3\n'
 for i,(x,y) in enumerate(coords[:points],1):text+=f'c {i} {x:.17g} {y:.17g} 0\n'
 regions={1:3,2:1,3:1}
 for r in sorted({r for r,_ in triangles}):text+=f'r {r} {regions[r]}\n'
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
  text+=f't {i} {r} '+ ' '.join(map(str,(*vs,*neighbors)))+'\n'
  nodes.update((v-1,regions[r]) for v in vs)
 text+='s 1 5\n'
 for point,material in sorted(nodes):text+=f'n {point} {material} '+('1e15' if material==3 else '1')+'\n'
 for i in (4,5,6,7):text+=f'n {i} 0 1e5\n'
 text+='M 100 1173\n';(folder/'input.str').write_text(text);audit(folder/'input.str')
 commands='option quiet\nmethod fermi vertical\nstructure infile=input.str\ndiffuse time=1e-9 temp=800 dryo2\nstructure outfile=seed.str\nquit\n'
 r=subprocess.run([str(binary)],input=commands,text=True,cwd=folder,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
 (folder/'suprem4.log').write_text(r.stdout)
 assert r.returncode==0 and 'panic' not in r.stdout.lower() and 'errors detected' not in r.stdout,f'{folder}: {r.stdout}'
 audit(folder/'seed.str');audit_partition(folder/'seed.str');result=parse_str(folder/'seed.str')
 if case!='bare':
  for vertex in ((1,2,9,10) if case=='wedge' else (1,2,3,4,9,10,11,12)):
   expected=tuple(struct.unpack('f',struct.pack('f',v*1e-4))[0]*1e4
                  for v in coords[vertex-1])
   distance=min(max(abs(x-y) for x,y in zip(expected,actual)) for actual in result['coords'].values())
   assert distance<1e-7,(expected,distance,folder)
 if case=='bare':
  top=[p for p in result['coords'].values() if p[1]<0]
  assert top and max(p[1] for p in top)<-.00099,(folder,top)
 if case=='gap':
  # A conforming cap may span the complete gap without an interior vertex.
  # Check material coverage, rather than requiring a particular triangulation.
  regions={int(s[1]):int(s[2]) for line in (folder/'seed.str').read_text().splitlines()
           if (s:=line.split()) and s[0]=='r'}
  oxide=[]
  for line in (folder/'seed.str').read_text().splitlines():
   s=line.split()
   if s and s[0]=='t' and regions[int(s[2])]==1:
    oxide.append([result['coords'][int(v)] for v in s[3:6]])
  for x in (1.1,1.5,1.9):
   for y in (-.0001,-.0005,-.00099):
    assert any(all(orientation(a,b,(x,y))>=0 for a,b in zip(t,t[1:]+t[:1]))
               for t in oxide),(folder,'uncoated gap',x,y)
 print(f'PASS native oxide {case}: {len(result["triangles"])} positive triangles',flush=True)
print(f'workdir: {work}')
