import concurrent.futures,hashlib,json,os,subprocess,time
from pathlib import Path
root=Path('/home/chris/Open_TCAD')
process=root/'decks/subcmos2014/results/directional_final_fresh_20261117'
output=root/'suprem4/experiments/directional_etch_20261117/electrical'
runtime=root/'suprem4/experiments/rit_electrical_20261116/genius_runtime'
deadline=time.monotonic()+1800
print('Waiting for audited fresh process results',flush=True)
while not (process/'results.json').exists():
 if time.monotonic()>deadline:raise TimeoutError('Fresh process did not finish')
 time.sleep(2)
results=json.loads((process/'results.json').read_text())
# The process runner returns a list of per-device results.
print('Fresh process finished; checking result format',flush=True)
entries=results if isinstance(results,list) else results.get('results',results.get('devices',[]))
if isinstance(entries,dict):entries=list(entries.values())
assert len(entries)==2 and all(r.get('passed',False) for r in entries),results
subprocess.run(['python3',str(root/'decks/subcmos2014/prepare_electrical.py'),'--process-results',str(process),'--output',str(output)],check=True,cwd=root)
env={**os.environ,'GENIUS_DIR':str(runtime),'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
binary=runtime/'bin/genius.LINUX'
def run(dev):
 d=output/dev;t=time.monotonic()
 with (d/'transfer_stdout.log').open('w') as log:
  r=subprocess.run([str(binary),'-i','transfer.in'],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600)
 record={'returncode':r.returncode,'elapsed_seconds':time.monotonic()-t,'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'deck_sha256':hashlib.sha256((d/'transfer.in').read_bytes()).hexdigest()}
 (d/'transfer_execution.json').write_text(json.dumps(record,indent=2)+'\n')
 print(dev,record,flush=True)
 return r.returncode
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:codes=list(pool.map(run,['nmos','pmos']))
assert codes==[0,0],codes
