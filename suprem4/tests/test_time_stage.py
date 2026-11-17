"""Failed comparison solves must never become valid timestep estimates."""
from pathlib import Path
import os
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class StageTest(unittest.TestCase):
    def test_stage_acceptance(self):
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory) / 'stages'
            source = os.environ.get('SUPREM_TIME_SOURCE', str(ROOT / 'src/diffuse/time.c'))
            subprocess.run([
                'gcc', '-std=gnu89', '-g', '-O1', '-fcommon',
                '-ffunction-sections', '-fdata-sections',
                '-fsanitize=address,undefined', '-fno-omit-frame-pointer',
                '-fno-pie', '-no-pie', '-I', str(ROOT / 'src/include'),
                str(ROOT / 'tests/time_stage.c'), source,
                '-Wl,--gc-sections', '-lm', '-o', str(binary)], check=True)
            env = {**os.environ, 'ASAN_OPTIONS': 'detect_leaks=0',
                   'UBSAN_OPTIONS': 'halt_on_error=1'}
            for scenario in range(5):
                with self.subTest(scenario=scenario):
                    run = subprocess.run([str(binary), str(scenario)],
                                         capture_output=True, text=True, env=env, timeout=5)
                    self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
                    self.assertIn('PASS', run.stdout)


if __name__ == '__main__':
    unittest.main()
