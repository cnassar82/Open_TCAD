#!/usr/bin/env python3
"""Perturb captured spacer/seed meshes without moving reflecting side planes.

This is a bounded geometric robustness check, not a mesh-convergence study or
an electrical calibration. No process depth, dose, or thermal budget changes.
"""
import argparse,gzip,hashlib,json,math,os,re,struct,subprocess
from pathlib import Path
from partition_audit import audit_partition
from inventory_audit import inventories
from run_oxidation_perturbations import audit
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);data=binary.parent.parent/'data'
env={**os.environ,'SUPREM4_PROCESS_SURFACE':'1','SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide','SUPREM4_MESH_PATCH_REMOVAL':'1'}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
f32=lambda x:struct.unpack('f',struct.pack('f',x))[0]
results=[]
for case in ('nmos_construction_spacer','pmos_native_seed'):
 f=Path(__file__).resolve().parent/'fixtures'/('rit_'+case);raw=gzip.decompress((f/'before.str.gz').read_bytes());info=json.loads((f/'fixture.json').read_text());assert hashlib.sha256(raw).hexdigest()==info['sha256']
 lines=raw.decode().splitlines();assert not any(l.startswith('G ') for l in lines)
 xs=[float(s[2]) for l in lines if (s:=l.split()) and s[0]=='c'];lo,hi=min(xs),max(xs)
 for amplitude in (-.002,.002):
  d=out/(case+('_left' if amplitude<0 else '_right'));d.mkdir();r={'case':case,'interior_shift_amplitude_um':amplitude,'fixture_sha256':info['sha256'],'errors':[]}
  try:
   changed=[]
   for l in lines:
    s=l.split()
    if s and s[0]=='c':
     x=float(s[2]);shift=0 if x in (lo,hi) else amplitude*math.sin(math.pi*(x-lo)/(hi-lo))
     s[2]=format(f32((x+shift)*1e-4)*1e4,'.17g');l=' '.join(s)
    changed.append(l)
   (d/'input.str').write_text('\n'.join(changed)+'\n');r['input_partition']=audit_partition(d/'input.str')
   deck='option quiet\n'+'\n'.join(info.get('settings',[]))+'\nstructure infile=input.str\n'+info['command']+'\nstructure outfile=after.str\nquit\n';(d/'input.in').write_text(deck)
   with (d/'input.in').open('rb') as stdin,(d/'run.log').open('wb') as log:
    r['returncode']=subprocess.run([str(binary)],cwd=d,env=env,stdin=stdin,stdout=log,stderr=subprocess.STDOUT,timeout=60).returncode
   if r['returncode'] or re.search(r'errors detected|ERROR:|runtime error:|AddressSanitizer:|suprem4 panic:',(d/'run.log').read_text(errors='replace'),re.I):raise ValueError('simulator error')
   r['mesh']=audit(d/'after.str');r['partition']=audit_partition(d/'after.str')
   if case=='nmos_construction_spacer':
    before=inventories(d/'input.str');after=inventories(d/'after.str')
    r['silicon_chemical_relative_errors']={str(species):(after[3,species]-before[3,species])/before[3,species] for species in (3,5)}
    if any(abs(v)>1e-12 for v in r['silicon_chemical_relative_errors'].values()):raise ValueError('oxide spacer construction changed silicon chemical inventory')
  except (OSError,ValueError,KeyError,subprocess.TimeoutExpired) as e:r['errors'].append(str(e))
  r['passed']=not r['errors'];results.append(r);print('PASS' if r['passed'] else 'FAIL',case,amplitude,r['errors'],flush=True)
r={'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'runtime_options':{k:env[k] for k in ('SUPREM4_PROCESS_SURFACE','SUPREM4_DEPOSITION_GEOMETRY','SUPREM4_MATERIAL_DOMAINS','SUPREM4_MESH_PATCH_REMOVAL')},'cases':results,'passed':all(r['passed'] for r in results)}
(out/'results.json').write_text(json.dumps(r,indent=2)+'\n');raise SystemExit(not r['passed'])
