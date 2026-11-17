#!/usr/bin/env python3
"""Replay unchanged poly doping/activation with per-command field inventories.

Diagnostic resumed runs only; these do not count as fresh CMOS acceptance.
Inventories integrate linear nodal chemical concentrations over triangle area
(cm^2), giving atoms/cm of out-of-plane depth. No active-species substitution.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

from validate import read, intervals, concentration

HERE = Path(__file__).resolve().parent


def inventory(path):
    coords, triangles, nodes = mesh = read(path)
    # Recover region identity separately; validate.read retains material only.
    regions = {}
    totals = {}
    for line in path.read_text().splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == 'r':
            regions[int(p[1])] = int(p[2])
        elif p[0] == 't':
            region = int(p[2]); material = regions[region]
            vs = list(map(int, p[3:6])); a, b, c = [coords[v] for v in vs]
            area = ((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])) / 2
            if area <= 0 or not math.isfinite(area):
                raise ValueError(f'invalid triangle area {area}')
            result = totals.setdefault(str(region), {'material': material, 'area_um2': 0.,
                                                     'chemical_atoms_per_cm': {'boron': 0., 'phosphorus': 0.}})
            result['area_um2'] += area
            for name, species in [('boron', 5), ('phosphorus', 3)]:
                values = [nodes[v, material].get(species, 0.) for v in vs]
                if not all(math.isfinite(x) and x >= 0 for x in values):
                    raise ValueError(f'invalid {name} field')
                result['chemical_atoms_per_cm'][name] += area * 1e-8 * sum(values) / 3
    materials = {}
    for result in totals.values():
        out = materials.setdefault(str(result['material']), {'area_um2': 0.,
                              'chemical_atoms_per_cm': {'boron': 0., 'phosphorus': 0.}})
        out['area_um2'] += result['area_um2']
        for species, value in result['chemical_atoms_per_cm'].items():
            out['chemical_atoms_per_cm'][species] += value
    probes = {}
    for name, x in [('body', 1.25), ('source', 4.25), ('gate', 6.), ('drain', 7.75)]:
        si = intervals(mesh, x, 3)
        if not si:
            raise ValueError(f'missing silicon at {name}')
        surface = si[0][0]
        poly = intervals(mesh, x, 4)
        probes[name] = {'silicon_surface_um': surface,
                        'poly_to_silicon_gap_um': surface-poly[-1][1] if poly else None,
                        'profile': {str(d): concentration(mesh, x, surface+d)
                                    for d in (.005, .01, .025, .05, .1, .25, .5)}}
    return {'regions': totals, 'materials': materials, 'probes': probes}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--resume', type=Path, required=True, help='folder containing device/gate_oxide.str')
    p.add_argument('--device', choices=['nmos', 'pmos'], required=True)
    p.add_argument('--after', choices=['locos_mask', 'gate_oxide', 'ldd', 'source_drain'], default='gate_oxide')
    p.add_argument('--through', choices=['locos', 'poly_doped', 'source_drain', 'annealed'], default='poly_doped')
    p.add_argument('--output', type=Path)
    p.add_argument('--timeout', type=float, default=900)
    args = p.parse_args()
    binary = args.binary.resolve(); data = binary.parent.parent/'data'
    work = (args.output or Path(tempfile.mkdtemp(prefix='rit_poly_trace_'))).resolve()
    work.mkdir(parents=True, exist_ok=True)
    if any(work.iterdir()):
        p.error(f'populated output: {work}')
    text = (HERE/(args.device+'.in')).read_text()
    prefix, tail = text.split(f'structure outfile={args.after}.str', 1)
    marker = f'structure outfile={args.through}.str'
    if marker not in tail:
        p.error('--through must follow --after')
    segment = tail.split(marker, 1)[0]
    source = args.resume.resolve()/args.device/(args.after+'.str')
    previous = source.parent/'process.in'
    if previous.read_text().split(f'structure outfile={args.after}.str', 1)[0] != prefix:
        p.error('resume process prefix differs')
    (work/'input.str').write_bytes(source.read_bytes())
    (work/'resist_top.dat').write_bytes((HERE/'resist_top.dat').read_bytes())
    settings = [line for line in prefix.splitlines() if line.startswith(('method ', 'oxide ', 'option '))]
    lines = settings + ['structure infile=input.str', 'structure outfile=00_input.str']
    stages = [{'name': '00_input', 'command': 'input checkpoint'}]
    for line in segment.splitlines():
        if not line.strip() or line.startswith(('#', 'structure ')):
            continue
        if re.match(r'etch (?:\S+ start|continue)\b', line):
            lines.append(line)
            continue
        name = f'{len(stages):02d}_' + line.split()[0]
        lines += [line, f'structure outfile={name}.str']
        stages.append({'name': name, 'command': line})
    lines += ['quit']
    deck = work/'trace.in'; deck.write_text('\n'.join(lines)+'\n')
    env = os.environ.copy()
    for key, name in [('KEYFILE', 'suprem.uk'), ('SUP4KEYFILE', 'suprem.uk'),
                      ('SUP4MODELRC', 'modelrc'), ('SUP4IMPDATA', 'sup4gs.imp'), ('IMPFILE', 'sup4gs.imp')]:
        env[key] = str(data/name)
    report = {'diagnostic_resumed_run': True, 'device': args.device,
              'after': args.after, 'through': args.through,
              'binary': str(binary), 'binary_sha256': digest(binary),
              'input_sha256': digest(source), 'deck_sha256': digest(deck),
              'data_sha256': {name: digest(data/name) for name in ('modelrc', 'suprem.uk', 'sup4gs.imp')},
              'stages': stages}
    manifest = binary.parent.parent/'build.json'
    if manifest.exists():
        report['build_manifest_sha256'] = digest(manifest)
    print(f'workdir: {work}', flush=True)
    with deck.open('rb') as stdin, (work/'trace.log').open('wb') as log:
        try:
            result = subprocess.run([str(binary)], stdin=stdin, stdout=log, stderr=subprocess.STDOUT,
                                    env=env, cwd=work, timeout=args.timeout)
            report['returncode'] = result.returncode
        except subprocess.TimeoutExpired:
            report['returncode'] = None
            report['timed_out'] = True
    log = (work/'trace.log').read_text(errors='replace')
    report['log_error'] = bool(re.search(r'runtime error:|AddressSanitizer:|suprem4 panic:|ERROR:|Newton convergence failed|Giving up', log, re.I))
    for stage in stages:
        try:
            stage['inventory'] = inventory(work/(stage['name']+'.str'))
            stage['status'] = 'PASS'
        except (OSError, ValueError, KeyError) as error:
            stage.update(status='FAIL', error=str(error))
    report['passed'] = report['returncode'] == 0 and not report['log_error'] and all(s['status']=='PASS' for s in stages)
    (work/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(f"{'PASS' if report['passed'] else 'FAIL'} diagnostic trace (not device acceptance)")
    return int(not report['passed'])


if __name__ == '__main__':
    raise SystemExit(main())
