from pathlib import Path
import argparse,os,subprocess,shutil,json,time,hashlib,re
root=Path('/home/chris/Open_TCAD');p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--full-locos',action='store_true');p.add_argument('--mesh-only',action='store_true');args=p.parse_args();source=root/'suprem4/experiments/implant_repair_20261117/fe_pivoted/pmos';d=args.output.resolve();d.mkdir(exist_ok=False);shutil.copy2(source/'locos_mask.str',d/'locos_mask.str');deck=(source/'pmos.in').read_text().replace('option quiet','option chat')
if not args.full_locos:deck=deck.split('diffuse time=65 temp=1100 weto2')[0]+'structure outfile=flood.str\nquit\n'
if args.full_locos:deck=deck.replace('diffuse time=65 temp=1100 weto2','structure outfile=flood.str\ndiffuse time=65 temp=1100 weto2')
if args.mesh_only:deck=deck.split('structure infile=locos_mask.str')[0]+'structure infile=locos_mask.str\nstructure outfile=quality_initial.str\nquit\n'
else:deck=deck.replace('structure infile=locos_mask.str','structure infile=locos_mask.str\nstructure outfile=quality_initial.str')
(d/'input.in').write_text(deck);b=args.binary.resolve();env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','SUPREM4_PROCESS_SURFACE':'1','SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide','SUPREM4_MESH_PATCH_REMOVAL':'1'}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp')]:env[k]=str(b.parent.parent/'data'/n)
t=time.monotonic()
timed_out=False
try:
 with (d/'run.log').open('w') as f:r=subprocess.run([str(b)],input=deck,text=True,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800)
 code=r.returncode
except subprocess.TimeoutExpired:
 timed_out=True;code=124
rec={'exit':code,'timed_out':timed_out,'runtime_options':{k:v for k,v in env.items() if k.startswith('SUPREM4_')},'source_structure_sha256':hashlib.sha256((d/'locos_mask.str').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t,'binary':str(b),'binary_sha256':hashlib.sha256(b.read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(deck.encode()).hexdigest()};(d/'execution.json').write_text(json.dumps(rec,indent=2));print(rec,flush=True)

raise SystemExit(code!=0)
