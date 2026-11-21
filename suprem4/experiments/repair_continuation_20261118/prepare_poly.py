#!/usr/bin/env python3
"""Prepare the unchanged blanket poly implant from an audited gate-oxide cut."""
import argparse,hashlib,json,os,re,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE),str(ROOT/'decks/subcmos2014')]
from audit_replay import inspect
from validate import read,dielectric_separation
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--initial',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
b=a.binary.resolve();initial=a.initial.resolve();d=a.output.resolve();d.mkdir(exist_ok=False)
recipe=ROOT/'decks/subcmos2014/nmos.in';commands=recipe.read_text().split('structure outfile=gate_oxide.str\n',1)[1].split('structure outfile=poly_implanted.str',1)[0]
deck=f'option quiet\noxide Vc=300 Vr=30 Vd=25 stress.dep=true\nmethod fermi vertical grid.oxide=0.04\nstructure infile={initial}\n'+commands+'structure outfile=poly_implanted.str\nquit\n';(d/'input.in').write_text(deck)
env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','SUPREM4_PROCESS_SURFACE':'1','SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide','SUPREM4_MESH_PATCH_REMOVAL':'1'}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(b.parent.parent/'data'/n)
h=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();r={'accepted_full_cmos':False,'binary':str(b),'binary_sha256':h(b),'input_sha256':h(d/'input.in'),'initial_sha256':h(initial),'recipe_sha256':h(recipe),'runtime_options':{k:v for k,v in env.items() if k.startswith('SUPREM4_')},'errors':[]};t=time.monotonic()
with (d/'run.log').open('w') as log:
 try:r['exit']=subprocess.run([str(b)],input=deck,text=True,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180).returncode
 except subprocess.TimeoutExpired:r['exit']=124
r['elapsed_seconds']=time.monotonic()-t
if re.search('Giving up|panic:|ERROR:|errors detected on command input',(d/'run.log').read_text(errors='replace')):r['errors'].append('simulator error')
try:
 f=d/'poly_implanted.str';r['mesh_audit']=inspect(f);r['dielectric_separation']=dielectric_separation(read(f));r['output_sha256']=h(f)
except (OSError,ValueError,KeyError) as e:r['errors'].append(str(e))
r['passed']=r['exit']==0 and not r['errors'];(d/'result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
raise SystemExit(not r['passed'])
