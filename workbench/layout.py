"""Project hierarchical GDSII mask geometry onto a 2D SUPREM4 process cut."""

from math import hypot, isfinite
import re

from runner import NUMBER
from suprem_cards import polygon_cards


def read_layout(path):
    import gdstk

    library = gdstk.read_gds(str(path), unit=1e-6)
    if not library.cells:
        raise ValueError("The GDS file has no cells.")
    return library


def layer_pairs(cell):
    """Collect specifications through the hierarchy without expanding every instance."""
    pairs = set()
    for item in [cell, *cell.dependencies(recursive=True)]:
        pairs.update((shape.layer, shape.datatype) for shape in item.polygons)
        for path in item.paths:
            pairs.update(zip(path.layers, path.datatypes))
    return sorted(pairs)


def merge_intervals(intervals, tolerance):
    merged = []
    for lo, hi in sorted(intervals):
        if hi - lo <= tolerance:
            continue
        if merged and lo <= merged[-1][1] + tolerance:
            merged[-1] = (merged[-1][0], max(hi, merged[-1][1]))
        else:
            merged.append((lo, hi))
    return merged


def section(polygons, start, end, *, tolerance=1e-9):
    """Return the cut length and drawn intervals in microns from the start."""
    x0, y0 = start
    length = hypot(end[0] - x0, end[1] - y0)
    if not isfinite(length) or length <= tolerance:
        raise ValueError("Choose distinct cut endpoints with a finite separation.")
    ux, uy = (end[0] - x0) / length, (end[1] - y0) / length
    intervals = []
    for polygon in polygons:
        projected = []
        for x, y in polygon.points:
            dx, dy = x - x0, y - y0
            offset = ux * dy - uy * dx
            projected.append((ux * dx + uy * dy, 0.0 if abs(offset) <= tolerance else offset))
        crosses = []
        for (x1, y1), (x2, y2) in zip(projected, projected[1:] + projected[:1]):
            if y1 == y2 == 0:
                intervals.append(tuple(sorted((x1, x2))))
            elif (y1 > 0) != (y2 > 0):
                crosses.append(x1 + (x2 - x1) * (-y1) / (y2 - y1))
        crosses.sort()
        intervals.extend(zip(crosses[::2], crosses[1::2]))
    clipped = [(max(0.0, lo), min(length, hi)) for lo, hi in intervals]
    return length, merge_intervals(clipped, tolerance)


def invert_intervals(intervals, length, tolerance):
    starts = [0.0, *(hi for _, hi in intervals)]
    ends = [*(lo for lo, _ in intervals), length]
    return [(lo, hi) for lo, hi in zip(starts, ends) if hi - lo > tolerance]


def fmt(value):
    return format(value, ".10g") if value else "0"


def mask_commands(intervals, length, origin, *, mode, spacing=0.02, material="photoresist",
                  top=-2.0, bottom=0.5, tolerance=1e-9):
    if mode == "Mesh lines":
        # Endpoints normally already have the tags from the structure wizard.
        edges = (edge for pair in intervals for edge in pair if tolerance < edge < length - tolerance)
        return "\n".join(f"line x loc={fmt(origin + x)} spac={fmt(spacing)}" for x in edges)
    return "\n".join(
        polygon_cards(material, "\n".join(
            f"{fmt(origin + x)} {fmt(y)}" for x, y in ((lo, top), (hi, top), (hi, bottom), (lo, bottom))
        ), "")
        for lo, hi in intervals
    )


X_LINE = re.compile(r"^\s*line\s+x(?:\.direction)?(?=\s|$)", re.I)
LOCATION = re.compile(r"\bloc(?:ation)?\s*=\s*(" + NUMBER + r")(?=\s|$)", re.I)


def merge_mesh_lines(deck, additions):
    """Insert edge lines into the initial mesh in the order SUPREM4 requires."""
    lines = deck.splitlines(keepends=True)
    limit = next((i for i, line in enumerate(lines)
                  if re.match(r"^\s*init(?:ialize)?(?:\s|$)", line, re.I)), len(lines))
    indexed = []
    for index, line in enumerate(lines[:limit]):
        if not X_LINE.match(line):
            continue
        location = LOCATION.search(line.partition("#")[0])
        if not location:
            raise ValueError("Use numeric X mesh locations before merging GDS edges.")
        indexed.append((index, float(location[1])))
    if len(indexed) < 2:
        raise ValueError("Set up the starting structure before adding GDS mesh lines.")
    locations = [x for _, x in indexed]
    if any(a >= b for a, b in zip(locations, locations[1:])):
        raise ValueError("The existing X mesh lines must be in increasing order.")
    insertions = {}
    for command in additions.splitlines():
        x = float(LOCATION.search(command)[1])
        if any(abs(x - old) < 1e-9 for old in locations):
            continue
        if not locations[0] < x < locations[-1]:
            raise ValueError(f"Mask edge X = {fmt(x)} µm is outside the starting X mesh.")
        index = next(i for i, location in indexed if location > x)
        insertions.setdefault(index, []).append(command + "\n")
    return "".join("".join(insertions.get(index, ())) + line for index, line in enumerate(lines))
