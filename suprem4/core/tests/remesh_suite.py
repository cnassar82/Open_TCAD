#!/usr/bin/env python3
"""WP3 test: remesh every .str in DIR/*.str and DIR/*/*.str with core_remesh.

PASS requires: invariants, B-rep, geometry preservation (within the
coarsening tolerance) and gas-node checks ok, no Mmg failure, and the box-
method criterion: no edge coupling (cot sum / 2, per region) below
--min-coupling (default -0.01, i.e. a boundary edge opposite > ~91 deg).

usage: remesh_suite.py CORE_REMESH DIR [DIR ...] [--out OUTDIR] [-- extra core_remesh options]
"""
import argparse
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

Q = re.compile(r'points (\d+) triangles (\d+) \| min angle ([\d.]+) max angle ([\d.]+) \| <20deg (\d+) <10deg (\d+) '
               r'>120deg (\d+) \| negative couplings (\d+) \| obtuse boundary edges (\d+) \| couplings < -0.01: (\d+) worst (\S+)')


def quality(line):
    m = Q.search(line or '')
    if not m:
        return None
    k = ['points', 'triangles', 'min', 'max', 'lt20', 'lt10', 'gt120', 'neg', 'obtuse', 'sig', 'worst']
    return {a: float(b) for a, b in zip(k, m.groups())}


def main():
    argv = sys.argv[1:]
    extra = []
    if '--' in argv:
        i = argv.index('--')
        argv, extra = argv[:i], argv[i + 1:]
    ap = argparse.ArgumentParser()
    ap.add_argument('binary')
    ap.add_argument('dirs', nargs='+')
    ap.add_argument('--out')
    ap.add_argument('--min-coupling', type=float, default=-0.01)
    a = ap.parse_args(argv)
    files = sorted({p for d in a.dirs for p in list(Path(d).glob('*.str')) + list(Path(d).glob('*/*.str'))})
    out_dir = Path(a.out) if a.out else Path(tempfile.mkdtemp())
    out_dir.mkdir(parents=True, exist_ok=True)
    fails = 0
    tot_in = tot_out = 0
    print(f'{"structure":28s} {"pts in":>7s} {"out":>7s} {"maxang in":>9s} {"out":>6s} {"neg in":>6s} {"out":>4s} '
          f'{"worst coupling in":>17s} {"out":>8s} {"sec":>5s}  checks')
    for f in files:
        name = f'{f.parent.name}/{f.name}'
        out = out_dir / f'{f.parent.name}_{f.name}'
        t0 = time.time()
        r = subprocess.run([a.binary, str(f), str(out)] + extra, capture_output=True, text=True)
        dt = time.time() - t0
        lines = r.stdout.splitlines()
        qi = quality(next((l for l in lines if l.startswith('input')), ''))
        qo = quality(next((l for l in lines if l.startswith('output')), ''))
        checks = {k: ('ok' if f'{k} ok' in r.stdout else 'FAIL') for k in ('invariants', 'brep', 'geometry', 'gas-nodes')}
        mmg_fail = re.search(r'mmg failures (\d+)', r.stdout)
        ok = r.returncode == 0 and qo and all(v == 'ok' for v in checks.values())
        ok = ok and mmg_fail and int(mmg_fail.group(1)) == 0
        ok = ok and qo['worst'] >= a.min_coupling
        fails += not ok
        if qi and qo:
            tot_in += qi['points']; tot_out += qo['points']
            print(f'{name:28s} {qi["points"]:7.0f} {qo["points"]:7.0f} {qi["max"]:9.1f} {qo["max"]:6.1f} '
                  f'{qi["neg"]:6.0f} {qo["neg"]:4.0f} {qi["worst"]:17.3g} {qo["worst"]:8.3g} {dt:5.1f}  '
                  f'{"PASS" if ok else "FAIL"} ' + ' '.join(k for k, v in checks.items() if v != 'ok'))
        else:
            print(f'{name:28s} FAIL rc={r.returncode} {r.stderr.strip()[:200]}')
        for l in lines:
            if l.startswith('NOTE') or (' ' in l and l.split()[0] in checks and l.split()[1] not in ('ok', 'FAILED')):
                print('    ' + l)
    print(f'{len(files) - fails}/{len(files)} passed; total {tot_in:.0f} -> {tot_out:.0f} points; output in {out_dir}')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
