"""Run bounded numerical/geometry checks and record exact source hashes."""
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
CASES=['run_geometry_weights.sh','run_reflecting_normal.sh','run_interface_measure.sh',
       'run_moving_flux.sh','run_be_time.sh','run_be_mass_term.sh','run_boundary_angle.sh',
       'run_rect_projection.sh','run_ear_clip.sh','run_quad_split.sh',
       'test_time_stage.py','test_diffusion_retry.py','run_p1_mass.sh']
def run(name):
 env={**os.environ,'ASAN_OPTIONS':'detect_leaks=0','UBSAN_OPTIONS':'halt_on_error=1'}
 if name=='run_p1_mass.sh':env['SUPREM4_MASS_LUMPING']='p1'
 command=['sh' if name.endswith('.sh') else 'python3',str(ROOT/'suprem4/tests'/name)]
 r=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env,cwd=ROOT,timeout=90)
 log=OUT/(name+'.log');log.write_bytes(r.stdout)
 return {'test':name,'returncode':r.returncode,'passed':r.returncode==0,
         'log':str(log.relative_to(ROOT)),'log_sha256':hashlib.sha256(r.stdout).hexdigest()}
if __name__=='__main__':
 with ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(run,CASES))
 (OUT/'focused_checks.json').write_text(json.dumps(results,indent=2)+'\n')
 for r in results:print('PASS' if r['passed'] else 'FAIL',r['test'])
 raise SystemExit(any(not r['passed'] for r in results))
