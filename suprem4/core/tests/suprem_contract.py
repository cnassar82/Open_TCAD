#!/usr/bin/env python3
"""WP6a test: the core's SUPREM validity contract (core_check_suprem) must
agree with SUPREM itself.

For every .str file given (directories are scanned for *.str and */*.str) the
core checker classifies the file as OK / WARN (SUPREM reloads it; a strict
rule fails) / FAIL (SUPREM panics or crashes on reload, or its embedding check
rejects it), and SUPREM loads it with SUPREM4_LOAD_CHECK=1 ("structure
infile=", then its float-orientation count and mesh_embedding_check).  Each
file must agree: OK/WARN <-> SUPREM loads it and its embedding check passes,
FAIL <-> SUPREM panics, crashes or rejects the embedding.  Files the core
cannot read (ERR: e.g. non-sequential record numbers) are listed separately.

With --mutate, five corruptions of the first file are generated and checked
the same way (inverted vertex order, missing exposed gas node, gas node at an
interior point, empty region, coincident points).

usage: suprem_contract.py CORE_BUILD SUPREM_BUILD [--mutate] [-j N] PATH...
  CORE_BUILD    directory with core_check_suprem
  SUPREM_BUILD  native64 build directory (src/suprem4, data/)
"""
import argparse
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def files_of(paths):
    out = []
    for p in map(Path, paths):
        out += sorted(p.glob('*.str')) + sorted(p.glob('*/*.str')) if p.is_dir() else [p]
    return out


def core_verdict(binary, f):
    r = subprocess.run([binary, str(f)], capture_output=True, text=True)
    first = (r.stdout.splitlines() or ['ERR'])[0]
    return first.split()[0], r.stdout.strip()


def suprem_verdict(build, f):
    d = Path(build) / 'data'
    env = dict(os.environ, SUP4KEYFILE=str(d / 'suprem.uk'), KEYFILE=str(d / 'suprem.uk'),
               SUP4MODELRC=str(d / 'modelrc'), SUP4IMPDATA=str(d / 'sup4gs.imp'), IMPFILE=str(d / 'sup4gs.imp'),
               SUPREM4_LOAD_CHECK='1')
    with tempfile.TemporaryDirectory() as w:
        try:
            r = subprocess.run([str(Path(build) / 'src' / 'suprem4')], input=f'structure infile={Path(f).resolve()}\n',
                               capture_output=True, text=True, cwd=w, env=env, timeout=300)
        except subprocess.TimeoutExpired:
            return 'TIMEOUT', ''
    text = r.stdout + r.stderr
    line = next((l for l in text.splitlines() if '[load-check]' in l), '')
    if r.returncode == 0 and ' orientation 0 embedding ok' in line:
        return 'LOADS', line
    if 'embedding FAIL' in line:
        return 'REJECTS', line
    why = next((l for l in text.splitlines() if 'panic' in l or 'self test' in l), f'exit {r.returncode}')
    return 'PANIC', why


def mutate(src, out_dir):
    """Five corruptions of a valid structure; returns their paths."""
    L = Path(src).read_text().split('\n')
    tri = [i for i, l in enumerate(L) if l.startswith('t ')]
    nod = [i for i, l in enumerate(L) if l.startswith('n ')]
    crd = [i for i, l in enumerate(L) if l.startswith('c ')]
    reg = [i for i, l in enumerate(L) if l.startswith('r ')]
    made = []

    def write(name, lines):
        p = Path(out_dir) / name
        p.write_text('\n'.join(lines))
        made.append(p)

    M = L[:]; f = M[tri[100]].split(); f[3], f[4] = f[4], f[3]; M[tri[100]] = ' '.join(f) + ' '
    write('m1_inverted_order.str', M)
    gas = [i for i in nod if L[i].split()[2] == '0']
    M = L[:]; del M[gas[len(gas) // 2]]
    write('m2_missing_exposed_gas.str', M)
    exposed = set()
    for i in tri:
        f = list(map(int, L[i].split()[1:9])); v, c = f[2:5], f[5:8]
        for j in range(3):
            if c[j] in (-1022, -1023): exposed |= {v[(j + 1) % 3], v[(j + 2) % 3]}
    interior = sorted({int(L[i].split()[1]) for i in nod if int(L[i].split()[1]) + 1 not in exposed})
    p = interior[len(interior) // 2]
    ref = next(i for i in nod if int(L[i].split()[1]) == p)
    f = L[ref].split(); f[2] = '0'; M = L[:]; M.insert(ref + 1, ' '.join(f))
    write('m3_interior_gas.str', M)
    M = L[:]; M.insert(reg[-1] + 1, f'r {len(reg) + 1}   1')
    write('m4_empty_region.str', M)
    f = list(map(int, L[tri[300]].split()[1:9])); a, b = f[2], f[3]
    ca, cb = L[crd[a - 1]].split(), L[crd[b - 1]].split(); M = L[:]
    M[crd[a - 1]] = ' '.join(['c', ca[1], cb[2], cb[3]]) + '  0'
    write('m5_coincident.str', M)
    return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('core_build'); ap.add_argument('suprem_build'); ap.add_argument('paths', nargs='+')
    ap.add_argument('--mutate', action='store_true'); ap.add_argument('-j', type=int, default=8)
    a = ap.parse_args()
    a.core_build, a.suprem_build = Path(a.core_build).resolve(), Path(a.suprem_build).resolve()
    checker = str(Path(a.core_build) / 'core_check_suprem')
    files = files_of(a.paths)
    tmp = tempfile.TemporaryDirectory()
    if a.mutate:
        files += mutate(files[0], tmp.name)

    def one(f):
        c, cdet = core_verdict(checker, f)
        s, sdet = suprem_verdict(a.suprem_build, f)
        # ERR: outside the core's input format (e.g. non-sequential record numbers, which SUPREM
        # renumbers silently); reported, not counted as a disagreement
        agree = (c in ('OK', 'WARN') and s == 'LOADS') or (c == 'FAIL' and s in ('PANIC', 'REJECTS')) or c == 'ERR'
        return f, c, s, agree, cdet, sdet

    with ThreadPoolExecutor(a.j) as ex:
        results = list(ex.map(one, files))
    counts = {}
    bad = 0
    for f, c, s, agree, cdet, sdet in results:
        counts[(c, s)] = counts.get((c, s), 0) + 1
        if not agree:
            bad += 1
            print(f'DISAGREE core {c} suprem {s} {f}\n  core: {cdet[:300]}\n  suprem: {sdet[:200]}')
    for (c, s), n in sorted(counts.items()):
        print(f'core {c:4s} suprem {s:7s} {n}')
    for f, c, s, agree, cdet, sdet in results:
        if c == 'ERR':
            print(f'UNREADABLE by the core (SUPREM: {s}) {f}: {cdet[:160]}')
    print(f'{len(results) - bad}/{len(results)} agree')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
