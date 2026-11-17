import os,json,subprocess,shutil,hashlib,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
root=Path('/home/chris/Open_TCAD');out=Path(__file__).resolve().parent/'electrical';out.mkdir(exist_ok=True)
previous=root/'suprem4/experiments/directional_etch_20261117/electrical';runtime=root/'suprem4/experiments/rit_electrical_20261116/genius_runtime';binary=runtime/'bin/genius.LINUX'
env={**os.environ,'GENIUS_DIR':str(runtime),'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
jobs=[]
for dev in ('pmos',):
 sign=1 if dev=='nmos' else -1
 base=(previous/dev/'equilibrium.in').read_text().split('EXPORT ')[0]
 jobs.append((dev,'transfer_vd3p3',base+f'SOLVE Type=DC Vscan=Drain Vstart=0 Vstep={sign*.05:g} Vstop={sign*3.3:g} out.prefix=drain_prep\nSOLVE Type=DC Vscan=Gate Vstart=0 Vstep={sign*.05:g} Vstop={sign*3.3:g} out.prefix=transfer\nEXPORT VTKFile=final.vtu\nEND\n'))
 for gate in (.75,1.5,3.3):
  name=f'output_vg{gate:g}'.replace('.','p')
  jobs.append((dev,name,base+f'SOLVE Type=DC Vscan=Gate Vstart=0 Vstep={sign*.05:g} Vstop={sign*gate:g} out.prefix=gate_prep\nSOLVE Type=DC Vscan=Drain Vstart=0 Vstep={sign*.05:g} Vstop={sign*3.3:g} out.prefix=output\nEXPORT VTKFile=final.vtu\nEND\n'))
def run(job):
 dev,name,deck=job;d=out/dev/name;d.mkdir(parents=True,exist_ok=False);shutil.copy2(previous/dev/'device.tif',d/'device.tif');(d/'device.in').write_text(deck)
 start=time.monotonic()
 with (d/'stdout.log').open('w') as log:
  try:r=subprocess.run([str(binary),'-i','device.in'],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=900);code=r.returncode
  except subprocess.TimeoutExpired:code=None
 record={'device':dev,'case':name,'returncode':code,'elapsed_seconds':time.monotonic()-start,'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'deck_sha256':hashlib.sha256(deck.encode()).hexdigest(),'input_tif_sha256':hashlib.sha256((d/'device.tif').read_bytes()).hexdigest()}
 (d/'execution.json').write_text(json.dumps(record,indent=2)+'\n');print(record,flush=True);return record
with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,jobs))
(out/'pmos_executions.json').write_text(json.dumps(results,indent=2)+'\n')
raise SystemExit(not all(r['returncode']==0 for r in results))
