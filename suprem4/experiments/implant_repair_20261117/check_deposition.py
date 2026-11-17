from pathlib import Path
import os,subprocess,json,sys,numpy as np,hashlib
from concurrent.futures import ThreadPoolExecutor
root=Path('/home/chris/Open_TCAD');out=Path(__file__).resolve().parent/'deposition_spacing';out.mkdir(exist_ok=True)
helper=root/'suprem4/experiments/presentation_audit_20261117/extract_profiles.py';ns={'__file__':str(helper)};exec(helper.read_text().split('summary={};all_profiles={}')[0],ns)
binary=root/'suprem4/build/native64_2eyknnhy/src/suprem4';data=binary.parent.parent/'data'
env={**os.environ,'SUPREM4_PROCESS_SURFACE':'1','SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide','SUPREM4_MESH_PATCH_REMOVAL':'1'}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
def run(case):
 dev,h=case;d=out/f'{dev}_h{h:g}';d.mkdir(exist_ok=False);source=root/'decks/subcmos2014/results/directional_final_fresh_20261117'/dev/'gate_oxide.str'
 deck=f'option quiet\nmethod fermi vertical grid.oxide=0.04\nstructure infile={source}\ndeposit poly thick=0.5 divisions={round(.5/h)} space={h}\nstructure outfile=before.str\nimplant phosphorus dose=2e16 energy=50 pearson\nstructure outfile=after.str\nquit\n';(d/'input.in').write_text(deck)
 with (d/'run.log').open('w') as log:r=subprocess.run([str(binary)],input=deck,text=True,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600)
 record={'device':dev,'spacing_um':h,'returncode':r.returncode,'deck_sha256':hashlib.sha256(deck.encode()).hexdigest()}
 if r.returncode==0:
  mesh=ns['read'](d/'after.str');y,v=ns['cut'](mesh,6,4);record.update(center_poly_P_sheet_cm2=float(np.trapezoid(v[:,1],y)*1e-4),center_poly_peak_P_cm3=float(max(v[:,1])),vertices=len(mesh[0]),triangles=len(mesh[1]))
 (d/'result.json').write_text(json.dumps(record,indent=2)+'\n');print(record,flush=True);return record
with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,[('nmos',.05),('nmos',.02),('nmos',.01),('pmos',.02)]))
(out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
