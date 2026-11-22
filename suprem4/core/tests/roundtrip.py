#!/usr/bin/env python3
"""WP1 test: every .str read into the core Mesh and written back must be
byte-identical, and the structural invariants must hold.

usage: roundtrip.py CORE_STR_BINARY DIR [DIR ...]
Scans DIR/*.str and DIR/*/*.str (e.g. a run folder with nmos/ and pmos/)."""
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    binary, dirs = sys.argv[1], sys.argv[2:]
    files = sorted({p for d in dirs for p in list(Path(d).glob('*.str')) + list(Path(d).glob('*/*.str'))})
    failures = 0
    with tempfile.TemporaryDirectory() as tmp:
        for f in files:
            out = Path(tmp) / 'rt.str'
            r = subprocess.run([binary, str(f), '--write', str(out)], capture_output=True, text=True)
            identical = out.exists() and out.read_bytes() == f.read_bytes()
            inv = 'ok' if 'invariants ok' in r.stdout else 'FAILED'
            brep = 'ok' if 'brep ok' in r.stdout else 'FAILED'
            status = 'PASS' if identical and inv == 'ok' and brep == 'ok' else 'FAIL'
            if status == 'FAIL':
                failures += 1
            summary = r.stdout.splitlines()[0] if r.stdout else r.stderr.strip()
            brl = next((l for l in r.stdout.splitlines() if l.startswith('brep loops')), '')
            print(f'{status} {f.parent.name}/{f.name}: identical={identical} invariants={inv} brep={brep} | {brl}')
            for line in r.stdout.splitlines():
                if line.startswith(('INVARIANT', 'BREP')):
                    print('    ' + line)
    print(f'{len(files) - failures}/{len(files)} passed')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
