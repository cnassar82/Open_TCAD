#!/usr/bin/env python3
"""WP6b test: Lagrangian update (core_update) of every .str given.

A file passes when all core_update checks are ok: invariants, B-rep,
geometry (B-rep within the removal tolerance), gas nodes, the SUPREM
contract (WP6a), retained (every node at a point outside the changed set
keeps its old values bit for bit) and inventory (conserved inventories
restored to 1e-10).

Two kinds of input:
  - a sequence of oxidation states (SUPREM4_CORE_DUMP=DIR writes every mesh
    the hook receives): --reference DIR/core_in_0000.str, the grid at the
    start of the oxidation, as in the hook;
  - static structures (no reference: each file is its own reference, so the
    update must change little).

usage: update_suite.py CORE_BUILD [--reference REF.str] [--events DELTA] [-j N] PATH...
  PATH: .str files or directories (scanned for *.str and */*.str)
"""
import argparse
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

CHECKS = ['invariants', 'brep', 'geometry', 'gas-nodes', 'suprem', 'retained', 'inventory']


def files_of(paths):
    out = []
    for p in map(Path, paths):
        out += sorted(p.glob('*.str')) + sorted(p.glob('*/*.str')) if p.is_dir() else [p]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('core_build')
    ap.add_argument('paths', nargs='+')
    ap.add_argument('--reference')
    ap.add_argument('-j', type=int, default=8)
    ap.add_argument('--events', help='event size (um): WP6c local repair before the update, as in the hook')
    a = ap.parse_args()
    binary = str(Path(a.core_build).resolve() / 'core_update')

    def one(f):
        cmd = [binary, str(f), '/dev/null', '--float'] + (['--reference', a.reference] if a.reference else []) + \
              (['--events', a.events, '--box'] if a.events else [])
        r = subprocess.run(cmd, capture_output=True, text=True)
        status = {}
        for line in r.stdout.splitlines():
            k = line.split(' ', 1)
            if len(k) == 2 and k[0] in CHECKS and k[1] in ('ok', 'FAILED'):
                status[k[0]] = k[1]
        upd = next((l for l in r.stdout.splitlines() if l.startswith('update ')), '')
        ok = r.returncode == 0 and all(status.get(c) == 'ok' for c in CHECKS)
        failed = [c for c in CHECKS if status.get(c) != 'ok']
        detail = '' if ok else (r.stderr.strip() or next((l for l in r.stdout.splitlines() if l.split(' ', 1)[0] in failed), ''))
        return f, ok, failed, upd, detail

    with ThreadPoolExecutor(a.j) as ex:
        results = list(ex.map(one, files_of(a.paths)))
    for f, ok, failed, upd, detail in results:
        m = re.search(r'points (\d+) -> (\d+).*new triangles (\d+) of (\d+)', upd)
        stats = f'{m[1]:>7} -> {m[2]:>7} points, {m[3]:>6}/{m[4]:<7} triangles new' if m else ''
        print(f'{"PASS" if ok else "FAIL"} {stats}  {f}' + ('' if ok else f'\n     {",".join(failed)}: {detail[:200]}'))
    n = sum(ok for _, ok, *_ in results)
    print(f'{n}/{len(results)} passed')
    return 0 if n == len(results) else 1


if __name__ == '__main__':
    sys.exit(main())
