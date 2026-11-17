#!/usr/bin/env python3
"""Replay captured RIT process failures with an independent partition audit."""
import argparse,gzip,hashlib,json,os,re,subprocess,tempfile
from pathlib import Path
from run_oxidation_perturbations import audit
from partition_audit import audit_partition
from inventory_audit import audit_transfer
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path)
p.add_argument('--timeout',type=float,default=60,help='simulator timeout in seconds (allow extra time for sanitizer builds)')
p.add_argument('--case',choices=['nmos_locos_strip','pmos_activation_seed','pmos_ild_boundary','pad_side_boundary','pmos_kooi','pmos_spacer_clip','pmos_reflecting_normal','pmos_spacer_flat','nmos_domain_spacer','nmos_construction_spacer','pmos_native_seed','nmos_contact_shielding','vertical_constraints'],default='nmos_locos_strip')
a=p.parse_args();binary=a.binary.resolve();data=binary.parent.parent/'data'
work=(a.output or Path(tempfile.mkdtemp(prefix='rit_locos_strip_'))).resolve()
work.mkdir(parents=True,exist_ok=True)
if any(work.iterdir()):p.error('populated output')
f=Path(__file__).resolve().parent/'fixtures'/('rit_'+a.case)
raw=gzip.decompress((f/'before.str.gz').read_bytes());info=json.loads((f/'fixture.json').read_text())
assert hashlib.sha256(raw).hexdigest()==info['sha256']
(work/'input.str').write_bytes(raw)
text='option quiet\n'+'\n'.join(info.get('settings',[]))+'\nstructure infile=input.str\n'+info['command']+'\nstructure outfile=after.str\nquit\n'
(work/'strip.in').write_text(text);env=os.environ.copy()
if a.case=='vertical_constraints':
 env['SUPREM4_CONSTRAINT_TRACE']='1'
 env['SUPREM4_DIAG_LIVE_TRIANGLES']='1'
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
r={'binary':str(binary),'binary_sha256':digest(binary),'fixture_sha256':info['sha256'],
   'deck_sha256':digest(work/'strip.in'),'data_sha256':{n:digest(data/n) for n in ('modelrc','suprem.uk','sup4gs.imp')},'errors':[], 'runtime_options':{k:env.get(k) for k in ('SUPREM4_MATERIAL_DOMAINS','SUPREM4_MESH_PATCH_REMOVAL','SUPREM4_MASS_LUMPING','SUPREM4_TIME_INTEGRATOR','SUPREM4_POSITIVITY_STEP','SUPREM4_DEPOSITION_GEOMETRY','SUPREM4_VERTICAL_CONSTRAINTS','SUPREM4_PROCESS_SURFACE')}}
m=binary.parent.parent/'build.json'
if m.exists():r['build_manifest_sha256']=digest(m)
with (work/'strip.in').open('rb') as stdin,(work/'strip.log').open('wb') as log:
 try:r['returncode']=subprocess.run([str(binary)],stdin=stdin,stdout=log,stderr=subprocess.STDOUT,cwd=work,env=env,timeout=a.timeout).returncode
 except subprocess.TimeoutExpired:r['returncode']=None;r['errors'].append('timeout')
if re.search(r'errors detected|ERROR:|runtime error:|AddressSanitizer:|suprem4 panic:',(work/'strip.log').read_text(errors='replace'),re.I):r['errors'].append('simulator error')
try:r['mesh']=audit(work/'after.str');r['partition']=audit_partition(work/'after.str')
except (OSError,ValueError,KeyError) as error:r['errors'].append(str(error))
if info.get('require_gate_oxide') and (work/'after.str').exists():
 import sys
 sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'decks/subcmos2014'))
 from validate import read,gate_dielectric
 try:r['gate_oxide_material_probes']=gate_dielectric(read(work/'after.str'))
 except (ValueError,KeyError) as error:r['errors'].append(str(error))
if a.case=='vertical_constraints' and env.get('SUPREM4_VERTICAL_CONSTRAINTS')=='1':
 try:
  r['constraint_transfer']=audit_transfer(work/'constraint_before.str',work/'constraint_after.str')
  r['constraint_partition']=audit_partition(work/'constraint_after.str')
 except (OSError,ValueError,KeyError) as error:r['errors'].append(str(error))
if info.get('preserve_input_regions') and 'mesh' in r:
 before=audit(work/'input.str')['area_by_region_um2']
 after=r['mesh']['area_by_region_um2'];r['input_region_area_changes_um2']={}
 for region,value in before.items():
  change=after.get(region,0.)-value;r['input_region_area_changes_um2'][region]=change
  if abs(change)>1e-12*max(1.,abs(value)):
   r['errors'].append(f'deposition changed input region {region} area by {change} um2')
if info.get('reflecting_planes') and (work/'after.str').exists():
 coordinates={};faces=0
 for line in (work/'after.str').read_text().splitlines():
  v=line.split()
  if not v:continue
  if v[0]=='c':coordinates[int(v[1])]=tuple(map(float,v[2:4]))
  elif v[0]=='t':
   vertices=list(map(int,v[3:6]))
   for i,neighbor in enumerate(map(int,v[6:9])):
    if neighbor!=-1024:continue
    faces+=1;edge=[coordinates[vertices[(i+j)%3]] for j in (1,2)]
    if not any(all(point[axis]==value for point in edge) for axis,value in info['reflecting_planes']):
     r['errors'].append(f'reflecting face left its prescribed plane: triangle {v[1]}, face {i}, {edge}')
 r['reflecting_faces_checked']=faces
 if not faces:r['errors'].append('no reflecting faces checked')
r['passed']=r['returncode']==0 and not r['errors'];(work/'result.json').write_text(json.dumps(r,indent=2)+'\n')
print(f"{'PASS' if r['passed'] else 'FAIL'} {a.case}: {work}")
raise SystemExit(int(not r['passed']))
