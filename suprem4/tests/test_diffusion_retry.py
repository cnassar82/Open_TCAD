"""Bound real solstep recursion without running a costly diffusion solve."""
from pathlib import Path
import os,subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]
class RetryTest(unittest.TestCase):
    def test_retry_success_depth_and_work_limits(self):
        with tempfile.TemporaryDirectory() as d:
            binary=Path(d)/'retry'
            subprocess.run(['gcc','-std=gnu89','-g','-O1','-fcommon','-ffunction-sections','-fdata-sections',
                '-fsanitize=address,undefined','-fno-omit-frame-pointer','-fno-pie','-no-pie',
                '-I',str(ROOT/'src/include'),str(ROOT/'tests/diffusion_retry.c'),
                str(ROOT/'src/diffuse/solve_time.c'),'-Wl,--gc-sections','-lm','-o',str(binary)],check=True)
            env={**os.environ,'ASAN_OPTIONS':'detect_leaks=0','UBSAN_OPTIONS':'halt_on_error=1'}
            for case,code,needle in [(0,0,'PASS split step'),(1,1,'depth=16'),(2,1,'trials=129'),(3,0,'PASS velocity history')]:
                run=subprocess.run([str(binary),str(case)],capture_output=True,text=True,env=env,timeout=5)
                self.assertEqual(run.returncode,code,run.stdout+run.stderr)
                self.assertIn(needle,run.stdout+run.stderr)
                self.assertNotIn('runtime error:',run.stderr)
                if case in (1,2):self.assertIn('checkpoint calls=',run.stdout)
if __name__=='__main__':unittest.main()
