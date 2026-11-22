"""Genius DD transfer sweeps on the AdvCmos2014 (0.25 um lambda, 0.5 um poly) cuts, 7 um device active; SRIM P-in-nitride table; QF 5e11; geometric poly re-ox, 0.40 um as-etched -> 0.356 um.
Same bias protocol as the SubCmos2014 runs; dual-poly gate workfunctions."""
import os,json,subprocess,hashlib,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
root=Path('/home/chris/Open_TCAD');here=Path(__file__).resolve().parent
process=root/'decks/advcmos2014/results/polyreox_geo_20261122'
export=here/'export';runtime=root/'suprem4/experiments/rit_electrical_20261116/genius_runtime';binary=runtime/'bin/genius.LINUX'
if not export.exists():
 subprocess.run(['python3',str(root/'decks/advcmos2014/prepare_electrical.py'),'--process-results',str(process),'--output',str(export)],check=True,cwd=root/'decks/advcmos2014')
# hwloc's GL plugin (via OpenMPI) blocks connecting to an X display when
# DISPLAY is set; the solver needs no display.
env={k:v for k,v in os.environ.items() if k!='DISPLAY'}
env.update({'GENIUS_DIR':str(runtime),'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','HWLOC_COMPONENTS':'-gl'})
jobs=[]
for dev in ('nmos','pmos'):
 s=1 if dev=='nmos' else -1
 base=(export/dev/'equilibrium.in').read_text().split('EXPORT ')[0]
 for vd in (.05,3.3):
  jobs.append((dev,f'transfer_vd{vd:g}'.replace('.','p'),base+f'SOLVE Type=DC Vscan=Drain Vstart=0 Vstep={s*min(vd,.05):g} Vstop={s*vd:g} out.prefix=drain_prep\nSOLVE Type=DC Vscan=Gate Vstart=0 Vstep={s*.05:g} Vstop={s*3.3:g} out.prefix=transfer\nEXPORT VTKFile=final.vtu\nEND\n'))
 # sub-zero gate leg to see off-state / buried-channel behaviour
 jobs.append((dev,'transfer_neg_vd0p05',base+f'SOLVE Type=DC Vscan=Drain Vstart=0 Vstep={s*.01:g} Vstop={s*.05:g} out.prefix=drain_prep\nSOLVE Type=DC Vscan=Gate Vstart=0 Vstep={-s*.05:g} Vstop={-s*1.0:g} out.prefix=transfer\nEND\n'))
def run(job):
 dev,name,deck=job;d=here/dev/name;d.mkdir(parents=True,exist_ok=True)
 (d/'device.tif').write_bytes((export/dev/'device.tif').read_bytes());(d/'device.in').write_text(deck)
 t=time.monotonic()
 with (d/'stdout.log').open('w') as log:
  try:code=subprocess.run([str(binary),'-i','device.in'],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800).returncode
  except subprocess.TimeoutExpired:code=None
 rec={'device':dev,'case':name,'returncode':code,'elapsed_seconds':time.monotonic()-t,'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'deck_sha256':hashlib.sha256(deck.encode()).hexdigest()}
 (d/'execution.json').write_text(json.dumps(rec,indent=2)+'\n');print(rec,flush=True);return rec
with ThreadPoolExecutor(6) as p:res=list(p.map(run,jobs))
(here/'executions.json').write_text(json.dumps(res,indent=2)+'\n')
