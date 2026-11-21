#!/usr/bin/env python3
"""Audit saved replay structures independently of simulator exit status."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'suprem4/tests'), str(ROOT / 'decks/subcmos2014')]
from inventory_audit import audit_chemistry
from partition_audit import audit_partition
from run_oxidation_perturbations import audit
from validate import read


def inspect(path):
    result = {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'topology': audit(path), 'partition': audit_partition(path),
              'chemistry': audit_chemistry(path)}
    coordinates, triangles, _ = read(path)
    quality = []
    for _, vertices in triangles:
        a, b, c = [coordinates[v] for v in vertices]
        determinant = ((b[0]-a[0])*(c[1]-a[1]) -
                       (b[1]-a[1])*(c[0]-a[0]))
        denominator = sum(math.dist(x, y)**2 for x, y in ((a,b),(b,c),(c,a)))
        quality.append(2*math.sqrt(3)*determinant/denominator)
    result['points'] = len(coordinates)
    result['minimum_quality'] = min(quality)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    parser.add_argument('--final', choices=['flood', 'locos'], default='flood')
    args = parser.parse_args()
    result = {'accepted_full_cmos': False, 'errors': [], 'structures': {}}
    for name in ('locos_mask', 'quality_initial', args.final):
        try:
            result['structures'][name] = inspect(args.folder / (name+'.str'))
        except (OSError, ValueError, KeyError) as error:
            result['errors'].append(f'{name}: {error}')
    try:
        result['execution'] = json.loads((args.folder / 'execution.json').read_text())
        if result['execution']['exit'] != 0:
            result['errors'].append('simulator did not complete successfully')
    except (OSError, ValueError, KeyError) as error:
        result['errors'].append(f'execution: {error}')
    result['replay_mesh_checks_passed'] = not result['errors']
    (args.folder / 'independent_audit.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'folder':str(args.folder), 'passed':not result['errors'],
                      'errors':result['errors']}), flush=True)
    return int(bool(result['errors']))


if __name__ == '__main__':
    raise SystemExit(main())
