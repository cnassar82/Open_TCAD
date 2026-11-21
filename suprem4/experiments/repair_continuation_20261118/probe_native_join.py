#!/usr/bin/env python3
"""Detect dielectric pinch-outs when a native film joins existing oxide."""
import argparse,hashlib,json,os,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE),str(ROOT/'decks/subcmos2014')]
from validate import read
from audit_replay import inspect
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--join-x',type=float,default=.5);a=p.parse_args();b=a.binary.resolve();d=a.output.resolve();d.mkdir(exist_ok=False)
if not 0<a.join_x<1:p.error('--join-x must lie inside the structure')
deck='''option quiet
line x loc=0 spac=0.1 tag=left
line x loc=0.5 spac=0.05
line x loc=1 spac=0.1 tag=right
line y loc=0 spac=0.01 tag=top
line y loc=0.1 spac=0.02
line y loc=0.5 spac=0.1 tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
boundary exposed xlo=left xhi=right ylo=top yhi=top
init boron conc=1e15 orient=100
method fermi vertical grid.oxide=0.04
deposit oxide thick=0.02 divisions=2
etch oxide left p1.x=0.5
structure outfile=before_seed.str
diffuse time=1e-9 temp=900 dryo2
structure outfile=seeded.str
deposit poly thick=0.1 divisions=10 space=0.01
structure outfile=stack.str
quit
''';deck=deck.replace('etch oxide left p1.x=0.5',f'etch oxide left p1.x={a.join_x:.17g}');(d/'input.in').write_text(deck)
env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','SUPREM4_PROCESS_SURFACE':'1','SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide','SUPREM4_MESH_PATCH_REMOVAL':'1'}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(b.parent.parent/'data'/n)
h=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();r={'accepted_full_cmos':False,'binary':str(b),'binary_sha256':h(b),'input_sha256':h(d/'input.in'),'runtime_options':{k:v for k,v in env.items() if k.startswith('SUPREM4_')},'errors':[]};t=time.monotonic()
with (d/'run.log').open('w') as log:
 try:r['exit']=subprocess.run([str(b)],input=deck,text=True,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=60).returncode
 except subprocess.TimeoutExpired:r['exit']=124
r['elapsed_seconds']=time.monotonic()-t
try:
 r['mesh_audit']=inspect(d/'stack.str');c,tris,_=read(d/'stack.str');materials={}
 for mat,vs in tris:
  for v in vs:materials.setdefault(v,set()).add(mat)
 r['poly_silicon_contacts']=[{'point':v,'xy_um':c[v],'materials':sorted(m)} for v,m in materials.items() if {3,4}<=m]
 r['continuous_barrier']=not r['poly_silicon_contacts']
except (OSError,ValueError,KeyError) as e:r['errors'].append(str(e));r['continuous_barrier']=False
r['passed']=r['exit']==0 and not r['errors'] and r['continuous_barrier'];(d/'result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
raise SystemExit(not r['passed'])
