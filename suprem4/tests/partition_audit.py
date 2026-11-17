"""Independent geometric checks on the stored coordinates of a SUPREM mesh.

Filtered orientation predicates fall back to exact rational arithmetic. No
snapping tolerance is used: this checks the represented mesh, not its intent.
Holes require a separate reference-domain/area check.
"""
from functools import lru_cache
from fractions import Fraction
from pathlib import Path


@lru_cache(maxsize=100000)
def _exact(value):
    return Fraction(value)


def orientation(a, b, c):
    x = (b[0] - a[0]) * (c[1] - a[1])
    y = (b[1] - a[1]) * (c[0] - a[0])
    det = x - y
    # Conservative filter; the fallback preserves the exact sign of the
    # binary floating point coordinates, including almost collinear edges.
    if abs(det) > 1e-14 * (abs(x) + abs(y)):
        return (det > 0) - (det < 0)
    ax, ay = map(_exact, a)
    bx, by = map(_exact, b)
    cx, cy = map(_exact, c)
    det = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    return (det > 0) - (det < 0)


def _separated(a, b):
    return any(all(orientation(p, q, v) <= 0 for v in b)
               for p, q in zip(a, a[1:] + a[:1]))


def _bounds(points, identity):
    return (min(p[0] for p in points), max(p[0] for p in points),
            min(p[1] for p in points), max(p[1] for p in points), identity, points)


def _pairs(records, include_touch=False):
    active = []
    for record in sorted(records):
        x0, _, y0, y1 = record[:4]
        active = [r for r in active if r[1] >= x0] if include_touch else [
            r for r in active if r[1] > x0]
        for other in active:
            if include_touch:
                separated = other[3] < y0 or y1 < other[2]
            else:
                separated = other[3] <= y0 or y1 <= other[2]
            if not separated:
                yield record, other
        active.append(record)


def audit_partition(path):
    coordinates, triangles, boundary = {}, [], []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "c":
            coordinates[int(parts[1])] = tuple(map(float, parts[2:4]))
        elif parts[0] == "t":
            tid = int(parts[1])
            vertices = tuple(map(int, parts[3:6]))
            points = [coordinates[v] for v in vertices]
            if orientation(*points) <= 0:
                raise ValueError(f"nonpositive exact orientation: triangle {tid}")
            triangles.append(_bounds(points, tid))
            for i, neighbor in enumerate(map(int, parts[6:9])):
                if neighbor < 0:
                    edge = (vertices[(i + 1) % 3], vertices[(i + 2) % 3])
                    boundary.append(_bounds([coordinates[v] for v in edge], edge))
    if not triangles:
        raise ValueError("no triangles in partition")
    pairs = 0
    for first, second in _pairs(triangles):
        pairs += 1
        a, b = first[5], second[5]
        if not _separated(a, b) and not _separated(b, a):
            raise ValueError(f"overlapping triangle interiors: {first[4]}, {second[4]}")
    boundary_pairs = 0
    for first, second in _pairs(boundary, include_touch=True):
        boundary_pairs += 1
        a, b = first[5]
        c, d = second[5]
        ac, ad = orientation(a, b, c), orientation(a, b, d)
        ca, cb = orientation(c, d, a), orientation(c, d, b)
        if ac * ad < 0 and ca * cb < 0:
            raise ValueError(f"crossing boundary edges: {first[4]}, {second[4]}")
        if ac == ad == ca == cb == 0:
            axis = 0 if a[0] != b[0] else 1
            lo = max(min(a[axis], b[axis]), min(c[axis], d[axis]))
            hi = min(max(a[axis], b[axis]), max(c[axis], d[axis]))
            if hi > lo:
                raise ValueError(f"overlapping collinear boundary edges: {first[4]}, {second[4]}")
        # A boundary endpoint inside another boundary segment is an unnoded
        # contact, even if the incident material triangles do not overlap.
        for p, side, u, v in ((c, ac, a, b), (d, ad, a, b),
                               (a, ca, c, d), (b, cb, c, d)):
            if side == 0 and p != u and p != v and all(
                    min(u[k], v[k]) <= p[k] <= max(u[k], v[k]) for k in (0, 1)):
                raise ValueError(f"unnoded boundary contact: {first[4]}, {second[4]}")
    return {"triangle_pairs_checked": pairs, "boundary_pairs_checked": boundary_pairs,
            "triangles": len(triangles), "boundary_edges": len(boundary)}
