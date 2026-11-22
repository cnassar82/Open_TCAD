"""Stage checks for the Advanced CMOS cuts; mesh helpers are shared with subcmos2014."""
import importlib.util
from pathlib import Path

from build_decks import GATE, CONTACTS, PROBES, FIELD

_spec = importlib.util.spec_from_file_location(
    'subcmos2014_validate', Path(__file__).resolve().parents[1] / 'subcmos2014/validate.py')
_shared = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_shared)
read, intervals, concentration = _shared.read, _shared.intervals, _shared.concentration
dielectric_separation, gate_dielectric = _shared.dielectric_separation, _shared.gate_dielectric


def validate(path, device, stage):
    mesh = read(path); result = {}
    sign = 1 if device == 'nmos' else -1          # + for n-type S/D and gate
    def surface(x): return intervals(mesh, x, 3)[0][0]
    if stage in ('recessed_oxide', 'wells', 'gate_oxide'):
        spans = intervals(mesh, sum(FIELD[1]) / 2, 1)   # middle of the inner field
        result['field_oxide_um'] = sum(b - a for a, b in spans)
        if stage == 'recessed_oxide' and not .25 < result['field_oxide_um'] < .55:
            raise ValueError(f"field oxide {result['field_oxide_um']} um, expected ~0.365")
    if stage == 'wells':
        probe = concentration(mesh, PROBES['channel'], surface(PROBES['channel']) + .5)
        if probe['net_donor_cm3'] * -sign <= 0: raise ValueError(f'wrong well polarity: {probe}')
        result['well_0p5um_below_surface'] = probe
    if stage == 'gate_oxide':
        x = PROBES['channel']; ox = intervals(mesh, x, 1)
        result['gate_oxide_um'] = ox[-1][1] - ox[-1][0] if ox else 0
        if not .006 < result['gate_oxide_um'] < .016:
            raise ValueError(f"gate oxide {result['gate_oxide_um']} um, expected ~0.010")
    if stage in ('gate', 'annealed', 'contacts', 'metal1'):
        coords, tris, _ = mesh
        poly = [coords[v] for mat, vs in tris if mat == 4 for v in vs]
        if not poly: raise ValueError('missing poly gate')
        width = max(p[0] for p in poly) - min(p[0] for p in poly)
        drawn = GATE[1] - GATE[0]
        # Poly re-oxidation (when modeled) may consume the sidewalls after the 'gate' stage.
        # geometric poly re-ox etches the poly to its post-re-ox length
        if stage == 'gate' and not .8 * drawn - 1e-4 < width < drawn + 1e-4:
            raise ValueError(f'gate width {width} um, drawn {drawn}')
        if not .5 * drawn < width < drawn + 1e-4: raise ValueError(f'gate width {width} um, drawn {drawn}')
        result['gate_length_um'] = width
        x = sum(GATE) / 2; spans = intervals(mesh, x, 4)
        result['gate_dielectric_um'] = surface(x) - spans[0][1]
        if result['gate_dielectric_um'] <= 0: raise ValueError('missing gate dielectric')
        result['gate_oxide_material_probes'] = gate_dielectric(mesh)
        if stage != 'gate':
            probe = concentration(mesh, x, sum(spans[0]) / 2, 4)
            if probe['net_donor_cm3'] * sign <= 0:
                raise ValueError(f"poly gate is not {'n' if sign > 0 else 'p'}-type: {probe}")
            result['poly_midpoint'] = probe
    if stage in ('source_drain', 'annealed', 'metal1'):
        for name, s in [('source', sign), ('drain', sign), ('body', -sign)]:
            x = PROBES[name]; probe = concentration(mesh, x, surface(x) + .03)
            if s * probe['net_donor_cm3'] <= 0: raise ValueError(f'wrong {name} polarity: {probe}')
            result[name + '_30nm_below_surface'] = probe
    if stage == 'metal1':
        _, tris, _ = mesh; parent = {}
        def root(v):
            while parent[v] != v: parent[v] = parent[parent[v]]; v = parent[v]
            return v
        for mat, vs in tris:
            if mat != 6: continue
            for v in vs: parent.setdefault(v, v)
            for v in vs[1:]: parent[root(v)] = root(vs[0])
        count = len({root(v) for v in parent})
        if count != 4: raise ValueError(f'{count} metal islands, expected body/source/gate/drain')
        result['metal_islands'] = count
        for name, (lo, hi) in CONTACTS.items():
            x = (lo + hi) / 2; metal = intervals(mesh, x, 6)
            target = intervals(mesh, x, 4 if name == 'gate' else 3)
            if not metal or not target: raise ValueError(f'missing {name} contact')
            gap = target[0][0] - metal[-1][1]
            if abs(gap) > 1e-5: raise ValueError(f'{name} contact gap {gap} um')
            result[name + '_contact_gap_um'] = gap
    return result
