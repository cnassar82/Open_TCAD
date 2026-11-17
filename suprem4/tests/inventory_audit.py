"""Independent P1 chemical inventory in SUPREM structure files.

Integrates stored nodal fields over stored triangles, grouped by material.
This is a rezoning check; it does not assert inventory is constant during a
physical process step with transport, reaction, or boundary fluxes.
"""
from collections import defaultdict
from math import fsum, isfinite
from pathlib import Path


def inventories(path, species=(3, 5)):
    coordinates, regions, nodes, triangles = {}, {}, {}, []
    streams = []
    for line in Path(path).read_text().splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == 'c':
            coordinates[int(p[1])] = tuple(map(float, p[2:4]))
        elif p[0] == 'r':
            regions[int(p[1])] = int(p[2])
        elif p[0] == 's':
            streams = list(map(int, p[2:]))
        elif p[0] == 'n':
            nodes[int(p[1])+1, int(p[2])] = dict(zip(streams, map(float, p[3:])))
        elif p[0] == 't':
            triangles.append((int(p[2]), tuple(map(int, p[3:6]))))
    terms = defaultdict(list)
    for region, vertices in triangles:
        a, b, c = [coordinates[v] for v in vertices]
        area = ((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))/2
        if not isfinite(area) or area <= 0:
            raise ValueError('invalid triangle area in inventory audit')
        material = regions[region]
        for v in vertices:
            for imp in species:
                value = nodes[v, material][imp]
                if not isfinite(value) or value < 0:
                    raise ValueError('invalid chemical concentration')
                terms[material, imp].append(area/3*value)
    if not terms:
        raise ValueError('no inventory records')
    return {key: fsum(values) for key, values in terms.items()}


def audit_transfer(before, after, tolerance=1e-12):
    old, new = inventories(before), inventories(after)
    if old.keys() != new.keys():
        raise ValueError('material or species lost during rezoning')
    errors = {f'{m}:{s}': (new[m, s]-value)/max(abs(value), 1e-300)
              for (m, s), value in old.items()}
    if any(abs(e) > tolerance for e in errors.values()):
        raise ValueError(f'rezoning changed chemical inventory: {errors}')
    return {'relative_errors_by_material_species': errors, 'tolerance': tolerance}


def audit_chemistry(path, species=(3, 5)):
    """Audit available process dopants without asserting conserved process dose.

    Physical implants, etches and boundary fluxes can change these inventories.
    Nonnegative finite chemistry and the integrated values are recorded here;
    a before/after conservation comparison belongs to a pure rezoning test.
    """
    available = []
    for line in Path(path).read_text().splitlines():
        values = line.split()
        if values and values[0] == 's':
            available = [s for s in species if s in list(map(int, values[2:]))]
            break
    if not available:
        raise ValueError('no requested process dopants present')
    return {'species_checked': available,
            'inventory_by_material_species': {
                f'{m}:{s}': value for (m, s), value in
                inventories(path, species=available).items()}}
