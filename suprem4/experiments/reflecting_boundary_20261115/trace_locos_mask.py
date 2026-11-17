import os,subprocess,shutil
from pathlib import Path
root=Path('/home/chris/Open_TCAD'); src=root/'decks/subcmos2014/results/astra_finite_clip_continue/pmos'
build=root/'suprem4/build/native64_yuj0z1w7'; w=root/'decks/subcmos2014/results/astra_locos_mask_trace'; w.mkdir(exist_ok=True)
shutil.copy(src/'wells.str',w/'input.str')
s=(src/'process.in').read_text().split('structure outfile=wells.str\n')[1].split('structure outfile=locos_mask.str')[0]
lines=['option quiet','method fermi vertical grid.oxide=0.04','structure infile=input.str']; index=[]
for i,l in enumerate(s.splitlines()):
 lines.append(l)
 if l.strip() and not l.startswith('#') and not ('etch continue' in l or 'etch nitride start' in l):
  lines.append(f'structure outfile=step_{i}.str');index.append((i,l))
lines+=['quit']; (w/'trace.in').write_text('\n'.join(lines)+'\n')
env=os.environ.copy();env.update(SUPREM4_MATERIAL_DOMAINS='oxide',SUPREM4_MESH_PATCH_REMOVAL='1')
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(build/'data'/n)
with (w/'trace.in').open() as f,(w/'trace.log').open('w') as log:print('exit',subprocess.run([str(build/'src/suprem4')],stdin=f,stdout=log,stderr=subprocess.STDOUT,cwd=w,env=env,timeout=180).returncode)
for i,l in index:
 f=w/f'step_{i}.str'
 if not f.exists():continue
 a=[]
 for row in f.read_text().splitlines():
  v=row.split()
  if v and v[0]=='c' and 0<abs(float(v[2]))<1e-10:a.append(v[:4])
 print(i,l,'drift',a[:3])
