"""Reproduce activation failures from saved valid states; retain diagnostics."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import os,subprocess,json,hashlib
root=Path(__file__).resolve().parents[3]
binary=root/'suprem4/build/native64_a9xbvu_4/src/suprem4'
data=binary.parent.parent/'data'
source=root/'decks/subcmos2014/results/astra_conservative_constraints_activation'
out=root/'decks/subcmos2014/results/astra_constraints_activation_diagnosis'
out.mkdir(exist_ok=False)
def run(device):
 d=out/device;d.mkdir();raw=(source/device/'step_failure.str').read_bytes();(d/'input.str').write_bytes(raw)
 (d/'input.in').write_text('option quiet\nmethod fermi vertical grid.oxide=0.04\nstructure infile=input.str\ndiffuse time=0.05 temp=825 dryo2 continue\nstructure outfile=after.str\nquit\n')
 env={**os.environ,'SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide','SUPREM4_MESH_PATCH_REMOVAL':'1','SUPREM4_VERTICAL_CONSTRAINTS':'1','SUPREM4_CONSTRAINT_TRACE':'1','SUPREM4_SOLVE_DIAG':'1','SUPREM4_MOTION_TRACE':'1'}
 for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
 with (d/'input.in').open('rb') as i,(d/'diagnostic.log').open('wb') as log:
  try:rc=subprocess.run([str(binary)],stdin=i,stdout=log,stderr=subprocess.STDOUT,cwd=d,env=env,timeout=90).returncode
  except subprocess.TimeoutExpired:rc=None
 r={'device':device,'returncode':rc,'source':str(source/device/'step_failure.str'),'source_sha256':hashlib.sha256(raw).hexdigest()}
 (d/'result.json').write_text(json.dumps(r,indent=2)+'\n');print(r,flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,['nmos','pmos']))
