#!/usr/bin/env python3
"""Verify that saving/reloading narrow float triangles preserves mesh and fields."""
import argparse
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_oxidation_perturbations import audit

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();data=binary.parent.parent/'data'
work=Path(tempfile.mkdtemp(prefix='suprem4_structure_roundtrip_'))
def f32(v):return struct.unpack('f',struct.pack('f',v))[0]
# Two float ULPs distinguish these sliver vertices. Six significant digits
# make all three y coordinates equal, destroying a finite triangle.
points=[(f32(8.25532e-4),f32(.23093e-4)),
        (f32(8.25563e-4),f32(.2309302e-4)),
        (f32(8.25594e-4),f32(.2309305e-4))]
# Force positive orientation independently of the chosen rounded vertices.
a0,b,c=points
if (b[0]-a0[0])*(c[1]-a0[1])-(b[1]-a0[1])*(c[0]-a0[0])<0:
 points[1],points[2]=points[2],points[1]
text='v SUPREM-IV.GS B.9305\nD 2 3 3\n'
for i,(x,y) in enumerate(points,1):text+=f'c {i} {x*1e4:.17g} {y*1e4:.17g} 0\n'
text+='r 1 3\nt 1 1 1 2 3 -1022 -1022 -1022\ns 1 5\n'
for i in range(3):text+=f'n {i} 3 {1.23456789012345e15+i*1e10:.17g}\n'
text+='M 100 1173\n'
(work/'input.str').write_text(text);audit(work/'input.str')
env=os.environ.copy()
for key,name in [('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc')]:env[key]=str(data/name)
commands='structure infile=input.str\nstructure outfile=first.str\nstructure infile=first.str\nstructure outfile=second.str\nquit\n'
r=subprocess.run([str(binary)],input=commands,text=True,cwd=work,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
(work/'suprem4.log').write_text(r.stdout)
assert r.returncode==0 and 'panic' not in r.stdout.lower(),r.stdout
for name in ['first.str','second.str']:audit(work/name)
def records(name):return [l.split() for l in (work/name).read_text().splitlines() if l.startswith(('c ','t ','n '))]
assert records('first.str')==records('second.str'),str(work)
# The material node IDs may be reordered; compare solutions by point ID.
def fields(name):return {s[1]:float(s[3]) for s in records(name) if s[0]=='n' and s[2]=='3'}
assert fields('input.str')==fields('first.str'),str(work)
print(f'PASS: narrow mesh and double-precision fields survive save/reload; {work}')
