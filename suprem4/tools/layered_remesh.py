#!/usr/bin/env python3
"""Deck-driven layered remesh for SUPREM .str files.

Reconstructs the x/y node profile a SUPREM deck would have produced (via the
same dvpram/addint logic as src/mesh/rect.c), drops interior points into each
material polygon found in the current .str, and uses Gmsh only as a constrained
Delaunay engine that connects those points without inserting refinement Steiner
points of its own.

The existing scalar/anisotropic path in gmsh_remesh_str.py is left in place;
this module is its sibling for the deck-anchored layered case.
"""

import argparse
import math
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import gmsh_remesh_str as gmstr


EPS = 1.0e-6
MR = 1.5  # mesh ratio cap; matches rect.c

GAS_MATERIAL = 0
SiO2_MATERIAL = 1
SiNi_MATERIAL = 2
SILICON_MATERIAL = 3
POLY_MATERIAL = 4
# Materials that get deck-driven interior points. Oxide/nitride are excluded by
# default because SUPREM's etch/oxidation logic panics ("node not in any
# triangle") when adaptive remeshing has to swallow our interior points.
DEFAULT_TARGET_MATERIALS = (SILICON_MATERIAL,)
# Materials where log-interpolation of dopant nodal values is appropriate
# (multi-decade dopant profiles). Other materials use arithmetic interpolation
# so that triangle-area integrals are exactly preserved.
DEFAULT_LOG_INTERP_MATERIALS = (SILICON_MATERIAL, POLY_MATERIAL)


def parse_deck_mesh_lines(deck_path):
    """Return {"x": [...], "y": [...]} of {"loc","spac","tag"} dicts from a SUPREM deck."""
    out = {"x": [], "y": []}
    line_re = re.compile(r"\bline\s+(x|y)\b", re.IGNORECASE)
    kv_re = re.compile(r"\b(loc|spac|tag)\s*=\s*([^\s#]+)", re.IGNORECASE)
    for raw in Path(deck_path).read_text(encoding="utf-8").splitlines():
        stripped = raw.lstrip()
        if stripped.startswith("#"):
            continue
        m = line_re.search(stripped)
        if not m:
            continue
        dim = m.group(1).lower()
        vals = {k.lower(): v for k, v in kv_re.findall(stripped)}
        if "loc" not in vals:
            continue
        try:
            loc = float(vals["loc"])
        except ValueError:
            continue
        spac = None
        if "spac" in vals:
            try:
                spac = float(vals["spac"])
            except ValueError:
                pass
        out[dim].append({"loc": loc, "spac": spac, "tag": vals.get("tag")})
    out["x"].sort(key=lambda d: d["loc"])
    out["y"].sort(key=lambda d: d["loc"])
    return out


def smooth_spacings(lines):
    """Match rect.c::smooth: cap spacings by neighbour interval, then propagate the smallest."""
    n = len(lines)
    if n == 0:
        return []
    # Initial pass: fill missing spacings with MAXFLOAT, then cap by neighbour interval.
    big = float("inf")
    spa = [(item["spac"] if item["spac"] is not None else big) for item in lines]
    # rect.c keeps user spacings negative to preserve them; reproduce by tracking a "user" flag.
    user = [(item["spac"] is not None) for item in lines]
    loc = [item["loc"] for item in lines]
    for i in range(n - 1):
        lej = abs(loc[i + 1] - loc[i])
        if not user[i] and lej < spa[i]:
            spa[i] = lej
        if not user[i + 1] and lej < spa[i + 1]:
            spa[i + 1] = lej
    # Smallest-first sweep to cap neighbours by mr.
    done = [False] * n
    for _ in range(n):
        jmin = -1
        hmin = big
        for j in range(n):
            if done[j]:
                continue
            if spa[j] < hmin:
                hmin = spa[j]
                jmin = j
        if jmin < 0:
            break
        for jmax in (jmin - 1, jmin + 1):
            if 0 <= jmax < n:
                lej = abs(loc[jmax] - loc[jmin])
                hM = (hmin + (MR - 1.0) * lej) / MR
                if hM < spa[jmin]:
                    spa[jmin] = hM
        done[jmin] = True
    return spa


def dvpram(hl, hr, el):
    """1:1 port of rect.c::dvpram. Returns (nnew, ratio, first)."""
    if hl == 0 or hr == 0:
        return 0, 1.0, 1.0
    nnew = 0
    ratio = 1.0
    first = 1.0
    if hl > el:
        hl = el
    if hr > el:
        hr = el
    if abs(el - hl) <= EPS * el and abs(el - hr) <= EPS * el:
        return nnew, ratio, first
    if abs(hr - hl) <= EPS * el:
        nnew = int((el / hl) - 0.49999)
        if nnew < 0:
            nnew = 0
        first = 1.0 / (nnew + 1)
        return nnew, ratio, first
    hr2 = hl / MR + (1.0 - 1.0 / MR) * el
    hl2 = hr / MR + (1.0 - 1.0 / MR) * el
    hr = min(hr, hr2)
    hl = min(hl, hl2)
    r = (el - hl) / (el - hr)
    rn = hr / hl
    if r <= 0 or rn <= 0:
        return 0, 1.0, 1.0
    nnew = int((math.log(rn) / math.log(r)) + 0.5)
    if nnew == 0:
        return 0, 1.0, 1.0
    ratio = math.exp(math.log(rn) / nnew)
    if ratio < 1.0 / MR:
        ratio = 1.0 / MR
    if ratio > MR:
        ratio = MR
    first = (ratio - 1.0) / (math.exp((nnew + 1) * math.log(ratio)) - 1.0)
    return nnew, ratio, first


def addint(lines):
    """Replicate rect.c::addint: insert intermediate locations between user planes."""
    if not lines:
        return []
    spa = smooth_spacings(lines)
    loc = [item["loc"] for item in lines]
    rloc = [loc[0]]
    for i in range(len(lines) - 1):
        x = loc[i]
        dx = loc[i + 1] - x
        nn, r, f = dvpram(spa[i], spa[i + 1], abs(dx))
        for _ in range(nn):
            x = x + f * dx
            rloc.append(x)
            f *= r
        rloc.append(loc[i + 1])
    return rloc


def compute_dim_profile(lines):
    """Return final 1D coordinate list for one dimension's user mesh lines."""
    return addint(lines)


def material_y_top(coords, regions, triangles, material):
    """Smallest y observed for any vertex of the given material."""
    ys = []
    for _tid, region, verts in triangles:
        if regions.get(region) != material:
            continue
        ys.extend(coords[v][1] for v in verts)
    if not ys:
        return None
    return min(ys)


def split_collinear_boundary_edges(edges, coords, tol=1.0e-7):
    """Split any boundary edge (a,b) that geometrically passes through another
    boundary vertex c. Returns a new edge list with elementary segments only.

    SUPREM input meshes sometimes have an edge skipping intermediate vertices
    (e.g. one triangle spans (top, bottom) while three triangles cover the same
    area through interior chord vertices). Gmsh cannot recover the long line
    constraint without complaining, so we sub-segment it.
    """
    vertex_ids = {v for a, b in edges for v in (a, b)}
    vertex_xy = {v: coords[v] for v in vertex_ids}

    def split_one(a, b):
        ax, ay = vertex_xy[a]
        bx, by = vertex_xy[b]
        dx, dy = bx - ax, by - ay
        L2 = dx * dx + dy * dy
        if L2 <= 0.0:
            return [(a, b)]
        best = None  # (t, vertex)
        for c, (cx, cy) in vertex_xy.items():
            if c == a or c == b:
                continue
            t = ((cx - ax) * dx + (cy - ay) * dy) / L2
            if t <= tol or t >= 1.0 - tol:
                continue
            px = ax + t * dx
            py = ay + t * dy
            if (px - cx) ** 2 + (py - cy) ** 2 > tol * tol:
                continue
            if best is None or t < best[0]:
                best = (t, c)
        if best is None:
            return [(a, b)]
        return split_one(a, best[1]) + split_one(best[1], b)

    out = []
    for a, b in edges:
        out.extend(split_one(a, b))
    return out


def region_polygons(coords, regions, triangles):
    """Return {region_id: [list of loops, each a list of vertex ids]} per region.

    First loop is the outer boundary; subsequent loops are holes.
    """
    boundaries = gmstr.build_region_boundary(triangles)
    out = {}
    for region, edges in boundaries.items():
        edges = split_collinear_boundary_edges(edges, coords)
        loops = gmstr.loops_from_edges(edges)
        if not loops:
            continue
        loops.sort(key=lambda lp: abs(gmstr.signed_area(coords, lp)), reverse=True)
        out[region] = loops
    return out


def point_in_polygon_xy(point, loop_xy):
    """Ray-cast test using (x,y) tuples instead of vertex ids."""
    x, y = point
    inside = False
    n = len(loop_xy)
    for i in range(n):
        ax, ay = loop_xy[i]
        bx, by = loop_xy[(i + 1) % n]
        if (ay > y) != (by > y):
            xinters = (bx - ax) * (y - ay) / (by - ay) + ax
            if x < xinters:
                inside = not inside
    return inside


def segment_distance(point, a, b):
    """Distance from point to segment ab in 2D."""
    px, py = point
    ax, ay = a
    bx, by = b
    dx = bx - ax
    dy = by - ay
    denom = dx * dx + dy * dy
    if denom <= 0.0:
        return math.hypot(px - ax, py - ay)
    t = ((px - ax) * dx + (py - ay) * dy) / denom
    t = max(0.0, min(1.0, t))
    qx = ax + t * dx
    qy = ay + t * dy
    return math.hypot(px - qx, py - qy)


def min_polygon_edge_distance(point, loop_xy):
    return min(segment_distance(point, loop_xy[i], loop_xy[(i + 1) % len(loop_xy)])
               for i in range(len(loop_xy)))


def clip_grid_to_regions(x_coords, y_coords, y_spacings, coords, polygons,
                         regions, target_materials, boundary_xy_set,
                         sliver_factor=0.4):
    """Drop candidate (x,y) points into target-material region polygons.

    Drops points within sliver_factor * local_y_spacing of any region edge
    (including non-target-material edges) and points that coincide with an
    existing boundary vertex. Non-target materials receive no interior points
    (boundary-only triangulation).
    """
    target_region_loops = {}
    for region, loops in polygons.items():
        if regions.get(region) not in target_materials:
            continue
        outer = [coords[v] for v in loops[0]]
        holes = [[coords[v] for v in loop] for loop in loops[1:]]
        target_region_loops[region] = (outer, holes)

    # All boundary loops (every region), used as a global edge-distance filter.
    all_loops_xy = []
    for region, loops in polygons.items():
        for loop in loops:
            all_loops_xy.append([coords[v] for v in loop])

    interior_by_region = defaultdict(list)
    for x in x_coords:
        for y, sy in zip(y_coords, y_spacings):
            if (round(x, 6), round(y, 6)) in boundary_xy_set:
                continue
            tol = max(sliver_factor * sy, 1.0e-9)
            # Global edge-distance filter — also kills points sitting on a region
            # interface even when the inside/outside test mis-classifies them.
            on_boundary = False
            for loop_xy in all_loops_xy:
                if min_polygon_edge_distance((x, y), loop_xy) < tol:
                    on_boundary = True
                    break
            if on_boundary:
                continue
            for region, (outer, holes) in target_region_loops.items():
                if not point_in_polygon_xy((x, y), outer):
                    continue
                if any(point_in_polygon_xy((x, y), h) for h in holes):
                    continue
                interior_by_region[region].append((x, y))
                break
    return interior_by_region


def y_spacing_at_anchor(profile, anchor_index):
    """Local spacing at a given index of the y profile (forward difference)."""
    if anchor_index + 1 < len(profile):
        return profile[anchor_index + 1] - profile[anchor_index]
    if anchor_index > 0:
        return profile[anchor_index] - profile[anchor_index - 1]
    return 1.0


def build_silicon_y_profile(deck_y_lines, silicon_top_y, structure_bottom_y,
                            interface_buffer=0.0):
    """Deck dvpram profile anchored at silicon_top_y, downward only.

    Skips d=0 (the interface itself) and any sample within interface_buffer of
    the interface. The buffer leaves the band of silicon nearest the oxide at
    its original SUPREM-mesh density, avoiding the aspect-ratio cliff that
    triggers SUPREM's triang::chop and interp panics during downstream
    oxidation.
    """
    raw = compute_dim_profile(deck_y_lines)
    if not raw or len(raw) < 2:
        return [], []
    offsets = [y - raw[0] for y in raw]
    deck_spacings = [offsets[i + 1] - offsets[i] for i in range(len(offsets) - 1)]
    deck_spacings.append(deck_spacings[-1])
    y_out, s_out = [], []
    for d, sy in zip(offsets[1:], deck_spacings[1:]):
        if d < interface_buffer:
            continue
        y = silicon_top_y + d
        if y > structure_bottom_y + 0.5 * sy:
            break
        y_out.append(y)
        s_out.append(sy)
    return y_out, s_out


def build_above_silicon_y_profile(silicon_top_y, structure_top_y, oxide_y_spacing,
                                  interface_buffer=0.0):
    """Uniform isotropic y profile above silicon_top_y for oxide / nitride / etc.

    Mirroring the deck's near-surface 0.5 nm spacing into oxide causes SUPREM
    panics ('node not in any triangle') during downstream oxidation because the
    resulting high-aspect-ratio triangles cannot survive interface motion. A
    uniform spacing comparable to the lateral resolution gives enough integral
    accuracy without crushing the aspect ratio. Within interface_buffer of the
    silicon-oxide interface, no new oxide layered points are emitted so the
    band immediately above silicon stays at its original SUPREM density.
    """
    if oxide_y_spacing <= 0:
        return [], []
    y_out, s_out = [], []
    y = silicon_top_y - max(oxide_y_spacing, interface_buffer)
    while y > structure_top_y + 0.5 * oxide_y_spacing:
        y_out.append(y)
        s_out.append(oxide_y_spacing)
        y -= oxide_y_spacing
    y_out.reverse()
    s_out.reverse()
    return y_out, s_out


def build_x_profile_for_str(deck_x_lines, structure_x_min, structure_x_max):
    raw = compute_dim_profile(deck_x_lines)
    if not raw:
        return []
    shift = structure_x_min - raw[0]
    shifted = [x + shift for x in raw]
    return [x for x in shifted
            if structure_x_min - 1.0e-9 <= x <= structure_x_max + 1.0e-9]


def structure_bbox(coords):
    xs = [p[0] for p in coords.values()]
    ys = [p[1] for p in coords.values()]
    return (min(xs), max(xs), min(ys), max(ys))


def boundary_vertex_xy(coords, polygons, decimals=6):
    """Set of (rounded x, rounded y) tuples; default 1 nm rounding is robust to FP fuzz."""
    out = set()
    for loops in polygons.values():
        for loop in loops:
            for v in loop:
                x, y = coords[v]
                out.add((round(x, decimals), round(y, decimals)))
    return out


def classify_rectangular_polygon(loop, coords, tol=1.0e-7):
    """If `loop` has 4 bbox-corner vertices and straight axis-aligned bottom,
    left, and right sides, return (True, corners, sides). The top side may be
    a non-flat curve (e.g. silicon-oxide interface after oxidation has eaten
    into silicon). Gmsh's Transfinite Surface handles curved boundary sides
    fine — only the four corners must be at the bbox extremes.

    Returns:
      - corners = [BL, BR, TR, TL] vertex IDs at the bbox extremes
      - sides   = {"bottom": [pid,...], "right": [...], "top": [...], "left": [...]}
        each containing the vertex IDs along that side in walking order, corners
        included.
    """
    if len(loop) < 4:
        return False, None, None
    xs = [coords[v][0] for v in loop]
    ys = [coords[v][1] for v in loop]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    if x_max - x_min < tol or y_max - y_min < tol:
        return False, None, None

    # Find the "corners" as the four side extremes — i.e. the vertices at the
    # extremes of each axis-aligned side. The bbox's true min-y / max-y may sit
    # at an interior x (e.g. a stepped silicon-oxide interface dipping deepest
    # at some x between the structure edges), so we don't require bbox corners
    # to exist as polygon vertices.
    left_vs  = [v for v in loop if abs(coords[v][0] - x_min) < tol]
    right_vs = [v for v in loop if abs(coords[v][0] - x_max) < tol]
    bottom_vs = [v for v in loop if abs(coords[v][1] - y_max) < tol]
    if not (left_vs and right_vs and bottom_vs):
        return False, None, None
    # BL/TL/BR/TR refer to bbox-quadrant naming (y increases downward):
    #   BL = top-left in screen = leftmost vertex with smallest y
    #   TL = bottom-left = leftmost vertex with largest y
    #   BR = top-right = rightmost vertex with smallest y
    #   TR = bottom-right = rightmost vertex with largest y
    bl = min(left_vs, key=lambda v: coords[v][1])
    tl = max(left_vs, key=lambda v: coords[v][1])
    br = min(right_vs, key=lambda v: coords[v][1])
    tr = max(right_vs, key=lambda v: coords[v][1])

    n = len(loop)
    bl_idx = loop.index(bl)
    walk = [loop[(bl_idx + k) % n] for k in range(n)]
    if gmstr.signed_area(coords, walk) < 0:
        walk = [loop[(bl_idx - k) % n] for k in range(n)]
    try:
        ibr = walk.index(br)
        itr = walk.index(tr)
        itl = walk.index(tl)
    except ValueError:
        return False, None, None
    if not (0 < ibr < itr < itl):
        return False, None, None

    bottom = walk[0:ibr + 1]
    right  = walk[ibr:itr + 1]
    top    = walk[itr:itl + 1]
    left   = walk[itl:] + [walk[0]]

    # In our convention the bbox "bottom" side (going from BL to BR) is the
    # silicon-oxide interface — it can be a curve. The other three sides must
    # be axis-aligned (left/right vertical, top horizontal).
    for v in left:
        if abs(coords[v][0] - x_min) > tol:
            return False, None, None
    for v in right:
        if abs(coords[v][0] - x_max) > tol:
            return False, None, None
    y_back = coords[tl][1]   # silicon back y (large y in SUPREM convention)
    for v in top:
        if abs(coords[v][1] - y_back) > tol:
            return False, None, None

    return True, [bl, br, tr, tl], {"bottom": bottom, "right": right, "top": top, "left": left}


def write_layered_geo(geo_path, coords, polygons, interior_by_region,
                      regions, point_size, remesh_regions=None,
                      extra_embed_pids_by_region=None,
                      transfinite_regions=None,
                      gmsh_algo="front2d"):
    """Emit a Gmsh .geo with boundaries + embedded interior points and refinement disabled.

    interior_by_region:        {region_id: [(x, y), ...]} new point coords to embed
    extra_embed_pids_by_region: {region_id: [pid, ...]} original-coord pids to also embed

    If remesh_regions is given, only those region IDs become Plane Surfaces; other
    regions still have their boundary vertices emitted (so shared interface points
    keep consistent IDs) but Gmsh doesn't retriangulate them — they'll be merged
    back in by the caller using the original triangles.
    """
    geo_coords = dict(coords)
    next_id = (max(coords) if coords else 0) + 1
    interior_point_ids = defaultdict(list)
    for region, pts in interior_by_region.items():
        for x, y in pts:
            geo_coords[next_id] = (x, y)
            interior_point_ids[region].append(next_id)
            next_id += 1
    if extra_embed_pids_by_region:
        for region, pids in extra_embed_pids_by_region.items():
            interior_point_ids[region].extend(pids)

    line_ids = {}
    line_defs = []
    next_line = 1

    loop_defs = []
    surface_defs = []  # (surface_id, [loop_ids], region_id)
    surface_region = {}
    physical_defs = defaultdict(list)
    next_loop = 1
    next_surface = 1

    # Track line IDs traversed when walking each region's outer loop, used by
    # the transfinite-silicon path to issue Transfinite Curve directives.
    region_outer_walk = {}  # region_id -> list of (a, b, line_id) tuples in walk order

    for region in sorted(polygons):
        if remesh_regions is not None and region not in remesh_regions:
            continue
        loops = polygons[region]
        gmsh_loops = []
        for loop_index, loop in enumerate(loops):
            signed_lines = []
            walk = []
            for a, b in zip(loop, loop[1:] + loop[:1]):
                key = gmstr.edge_key(a, b)
                if key not in line_ids:
                    line_ids[key] = next_line
                    line_defs.append((next_line, key[0], key[1]))
                    next_line += 1
                lid = line_ids[key]
                signed_lines.append(lid if (a, b) == key else -lid)
                walk.append((a, b, lid))
            if loop_index == 0:
                region_outer_walk[region] = walk
            loop_id = next_loop
            next_loop += 1
            loop_defs.append((loop_id, signed_lines))
            gmsh_loops.append((loop_id, loop))

        for loop_id, loop in gmsh_loops:
            if any(other_id != loop_id
                   and abs(gmstr.signed_area(coords, other_loop)) > abs(gmstr.signed_area(coords, loop))
                   and gmstr.polygon_contains(coords, other_loop, loop)
                   for other_id, other_loop in gmsh_loops):
                continue
            holes = []
            for hole_id, hole in gmsh_loops:
                if hole_id == loop_id:
                    continue
                if gmstr.polygon_contains(coords, loop, hole):
                    holes.append(hole_id)
            surf = next_surface
            next_surface += 1
            surface_defs.append((surf, [loop_id] + holes, region))
            surface_region[surf] = region
            physical_defs[regions.get(region, -1)].append(surf)

    algo_name_to_id = {
        "meshadapt": 1, "auto": 2, "initial2d": 3, "del2d": 5, "front2d": 6,
        "delquad": 8, "quadqs": 11,
    }
    algo_id = algo_name_to_id.get(gmsh_algo, 6)
    with open(geo_path, "w", encoding="utf-8") as handle:
        handle.write('SetFactory("Built-in");\n')
        handle.write("Mesh.MeshSizeFromPoints = 0;\n")
        handle.write("Mesh.MeshSizeFromCurvature = 0;\n")
        handle.write("Mesh.MeshSizeExtendFromBoundary = 0;\n")
        handle.write(f"Mesh.MeshSizeMin = {point_size:.12g};\n")
        handle.write("Mesh.MeshSizeMax = 1.0e30;\n")
        handle.write(f"Mesh.Algorithm = {algo_id};  // {gmsh_algo}\n")
        # Quality passes: with MeshSizeFromPoints=0 these will not insert Steiner
        # points; they only smooth and edge-swap the existing triangulation,
        # raising the minimum triangle quality that SUPREM's downstream mesh
        # mutation code can tolerate.
        # NOTE: disabled (was Smoothing=5, Optimize=1, OptimizeNetgen=1) because
        # they hang on thin high-aspect-ratio oxide regions with many embedded
        # layered points. We tolerate "ugly" triangles in exchange for the
        # algorithm actually terminating; SUPREM's grid mutation will smooth
        # them downstream anyway.
        handle.write("Mesh.Smoothing = 0;\n")
        handle.write("Mesh.Optimize = 0;\n")
        handle.write("Mesh.OptimizeNetgen = 0;\n")
        for pid in sorted(geo_coords):
            x, y = geo_coords[pid]
            handle.write(f"Point({pid}) = {{{x:.12g}, {y:.12g}, 0, {point_size:.12g}}};\n")
        for lid, a, b in line_defs:
            handle.write(f"Line({lid}) = {{{a}, {b}}};\n")
        for loop_id, signed_lines in loop_defs:
            handle.write(f"Curve Loop({loop_id}) = {{{', '.join(map(str, signed_lines))}}};\n")
        for surf, loops, _region in surface_defs:
            handle.write(f"Plane Surface({surf}) = {{{', '.join(map(str, loops))}}};\n")
        for material, surfaces in sorted(physical_defs.items()):
            handle.write(f"Physical Surface({material}) = {{{', '.join(map(str, surfaces))}}};\n")
        # Resolve which regions get Transfinite (skip embedding their interior).
        transfinite_emitted = set()
        if transfinite_regions:
            for surf, _loops, region in surface_defs:
                tf = transfinite_regions.get(region)
                if tf is None:
                    continue
                # tf is {"sides": {bottom/right/top/left -> [vid,...]},
                #        "corners": [BL, BR, TR, TL]}
                sides = tf["sides"]
                corners = tf["corners"]

                def line_ids_for_side(side_vids):
                    ids = []
                    for a, b in zip(side_vids, side_vids[1:]):
                        key = gmstr.edge_key(a, b)
                        lid = line_ids.get(key)
                        if lid is not None:
                            ids.append(lid)
                    return ids

                bottom_lines = line_ids_for_side(sides["bottom"])
                top_lines    = line_ids_for_side(sides["top"])
                left_lines   = line_ids_for_side(sides["left"])
                right_lines  = line_ids_for_side(sides["right"])

                def emit_opposite_pair(rich_lines, sparse_lines, side_name):
                    """Transfinite Surface needs each opposite-side pair to have
                    matching total subdivisions. Set rich-side per-line count to
                    2 (= each segment is one subdivision); set sparse-side to
                    `total_segments + 1` so its single line matches.
                    """
                    if not rich_lines or not sparse_lines:
                        return False
                    if len(rich_lines) < len(sparse_lines):
                        rich_lines, sparse_lines = sparse_lines, rich_lines
                    if len(sparse_lines) != 1 and len(sparse_lines) != len(rich_lines):
                        return False  # incompatible distribution
                    rich_n = 2
                    if len(rich_lines) == len(sparse_lines):
                        sparse_n = 2
                    else:
                        sparse_n = len(rich_lines) + 1
                    handle.write(
                        f"Transfinite Curve {{{', '.join(str(l) for l in rich_lines)}}}"
                        f" = {rich_n} Using Progression 1;  // {side_name} rich\n"
                    )
                    handle.write(
                        f"Transfinite Curve {{{', '.join(str(l) for l in sparse_lines)}}}"
                        f" = {sparse_n} Using Progression 1;  // {side_name} sparse\n"
                    )
                    return True

                ok_x = emit_opposite_pair(top_lines, bottom_lines, "x")
                ok_y = emit_opposite_pair(left_lines, right_lines, "y")
                if not (ok_x and ok_y):
                    # Incompatible side counts; skip the Transfinite Surface
                    # declaration so Gmsh falls back to standard meshing for it.
                    continue
                handle.write(
                    f"Transfinite Surface {{{surf}}} = "
                    f"{{{', '.join(str(c) for c in corners)}}};\n"
                )
                transfinite_emitted.add(region)
        for surf, _loops, region in surface_defs:
            if region in transfinite_emitted:
                continue  # Transfinite controls all nodes; no embedded points.
            embed = interior_point_ids.get(region, [])
            if embed:
                handle.write(f"Point{{{', '.join(map(str, embed))}}} In Surface{{{surf}}};\n")
        handle.write("Coherence;\n")

    return surface_region


def build_xy_to_original_id(coords, used_vertex_ids, decimals=6):
    """Map (rounded x, rounded y) -> original vertex id, for the vertices that
    appear in the preserved triangles. Used to merge Gmsh-output silicon nodes
    back onto the same IDs at the silicon-oxide interface."""
    xy_to_id = {}
    for vid in used_vertex_ids:
        x, y = coords[vid]
        xy_to_id[(round(x, decimals), round(y, decimals))] = vid
    return xy_to_id


def merge_silicon_with_preserved(coords, regions, triangles, surface_region,
                                 gmsh_nodes, gmsh_tris, remesh_region_ids,
                                 decimals=6):
    """Combine preserved (non-remesh) triangles with Gmsh-output silicon triangles.

    Returns (coords_merged, tris_merged, gmsh_node_to_merged_id) where:
      - coords_merged: {pid -> (x, y)} for every vertex referenced by any output triangle.
      - tris_merged: [(old_region_id, (v1, v2, v3))] for the entire merged mesh.
      - gmsh_node_to_merged_id: {gmsh_node_id -> merged_pid} for downstream interpolation.

    Preserved triangles keep their original vertex IDs. Gmsh silicon nodes that
    match (within 1 nm) a vertex used by a preserved triangle are re-bound to
    that original ID, so the silicon-oxide interface is shared exactly. The
    rest get fresh IDs above max(original).
    """
    preserved_tris = []
    preserved_vertex_ids = set()
    for _tid, region, verts in triangles:
        if region in remesh_region_ids:
            continue
        preserved_tris.append((region, verts))
        preserved_vertex_ids.update(verts)

    xy_to_orig = build_xy_to_original_id(coords, preserved_vertex_ids, decimals)

    coords_merged = {vid: coords[vid] for vid in preserved_vertex_ids}
    next_id = (max(coords) if coords else 0) + 1
    gmsh_to_merged = {}
    for gid, (x, y) in gmsh_nodes.items():
        key = (round(x, decimals), round(y, decimals))
        if key in xy_to_orig:
            mapped = xy_to_orig[key]
        else:
            mapped = next_id
            coords_merged[next_id] = (x, y)
            next_id += 1
        gmsh_to_merged[gid] = mapped

    tris_merged = list(preserved_tris)
    for _eid, tags, verts in gmsh_tris:
        if len(tags) < 2 or tags[1] not in surface_region:
            continue
        region = surface_region[tags[1]]
        mapped = tuple(gmsh_to_merged[v] for v in verts)
        # CCW vs CW: keep .str convention (negative signed area means CCW
        # in screen coords; mirror the existing writer's choice).
        p0, p1, p2 = (coords_merged[v] for v in mapped)
        area2 = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])
        if area2 > 0.0:
            mapped = (mapped[0], mapped[2], mapped[1])
        tris_merged.append((region, mapped))
    return coords_merged, tris_merged, gmsh_to_merged


def write_combined_str(path, coords_merged, tris_merged, regions, node_values,
                       sol_line, implant_line, model_line):
    """Emit a SUPREM .str from the merged preserve+remesh mesh.

    node_values: {(merged_pid, material): [vals]} indexed by sol_ids order.
    """
    point_ids = {gid: i + 1 for i, gid in enumerate(sorted(coords_merged))}
    region_ids = {old: i + 1 for i, old in enumerate(sorted({r for r, _ in tris_merged}))}

    edge_count = defaultdict(int)
    for _region, verts in tris_merged:
        for a, b in ((verts[1], verts[2]), (verts[2], verts[0]), (verts[1], verts[0])):
            edge_count[gmstr.edge_key(a, b)] += 1

    used_materials = defaultdict(set)
    for region, verts in tris_merged:
        material = regions.get(region, -1)
        for v in verts:
            used_materials[v].add(material)

    xs = [p[0] for p in coords_merged.values()]
    ys = [p[1] for p in coords_merged.values()]
    xmin, xmax = min(xs), max(xs)
    ymax = max(ys)

    def boundary_code(a, b):
        if edge_count[gmstr.edge_key(a, b)] != 1:
            return -1024
        ax, ay = coords_merged[a]
        bx, by = coords_merged[b]
        tol = 1.0e-8
        if abs(ay - ymax) < tol and abs(by - ymax) < tol:
            return -1023
        if (abs(ax - xmin) < tol and abs(bx - xmin) < tol) or (
            abs(ax - xmax) < tol and abs(bx - xmax) < tol
        ):
            return -1024
        return -1022

    for old_region, verts in tris_merged:
        for a, b in ((verts[1], verts[2]), (verts[2], verts[0]), (verts[1], verts[0])):
            bc = boundary_code(a, b)
            if edge_count[gmstr.edge_key(a, b)] == 1 and bc != -1024:
                used_materials[a].add(0)
                used_materials[b].add(0)

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("v SUPREM-IV.GS B.9305\n")
        handle.write("D 2 3 3\n")
        for gid in sorted(coords_merged):
            pid = point_ids[gid]
            x, y = coords_merged[gid]
            handle.write(f"c {pid} {x:.9g} {y:.9g}  0\n")
        for old_region in sorted(region_ids):
            handle.write(f"r {region_ids[old_region]}   {regions[old_region]}\n")
        for tid, (old_region, verts) in enumerate(tris_merged, start=1):
            pv = [point_ids[v] for v in verts]
            bc = [
                boundary_code(verts[1], verts[2]),
                boundary_code(verts[2], verts[0]),
                boundary_code(verts[1], verts[0]),
            ]
            handle.write(
                f"t {tid} {region_ids[old_region]} {pv[0]} {pv[1]} {pv[2]} "
                f"{bc[0]} {bc[1]} {bc[2]} -1 -1 \n"
            )
        handle.write(f"{model_line}\n")
        if sol_line is None:
            raise RuntimeError("input .str has no solution header")
        handle.write(f"{sol_line}\n")
        for gid in sorted(coords_merged):
            pid = point_ids[gid]
            for material in sorted(used_materials[gid]):
                vals = node_values.get((gid, material))
                if vals is None:
                    continue
                handle.write(f"n {pid - 1}   {material}  ")
                handle.write(" ".join(f"{v:.6e}" for v in vals))
                handle.write("\n")
        if implant_line:
            handle.write(f"{implant_line}\n")


def build_interface_edge_map(triangles, regions):
    """For each undirected edge, return the set of materials whose triangles use it.

    Used to classify a region's boundary edges:
      * shared with another material  -> material interface (silicon-oxide, etc.)
      * unique to this region         -> external structure boundary
    """
    edge_materials = defaultdict(set)
    for _tid, region, verts in triangles:
        material = regions.get(region, -1)
        for a, b in ((verts[0], verts[1]), (verts[1], verts[2]), (verts[2], verts[0])):
            edge_materials[gmstr.edge_key(a, b)].add(material)
    return edge_materials


def polygon_signed_orientation(loop_xy):
    area = 0.0
    n = len(loop_xy)
    for i in range(n):
        ax, ay = loop_xy[i]
        bx, by = loop_xy[(i + 1) % n]
        area += ax * by - bx * ay
    return area  # > 0 CCW, < 0 CW


def generate_interface_layer_points(region, coords, polygons, regions, triangles,
                                    edge_materials, neighbour_materials,
                                    layer_distances, edge_sample_spacing,
                                    sliver_factor=0.4):
    """Generate inward-offset points along a region's boundary edges that are
    interfaces with one of `neighbour_materials`.

    layer_distances: iterable of perpendicular offsets (um) into the region;
        e.g. [0.05] for one layer, [0.05, 0.15] for two.
    edge_sample_spacing: target spacing along each edge between samples.
    """
    loops = polygons.get(region, [])
    if not loops:
        return []
    outer = [coords[v] for v in loops[0]]
    holes = [[coords[v] for v in loop] for loop in loops[1:]]
    # Inward normal sign: for a CCW outer loop with screen y-down, inward normal
    # is the right of the edge tangent; with y-up it's the left. We make it
    # robust by checking against the polygon centroid.
    cx = sum(p[0] for p in outer) / len(outer)
    cy = sum(p[1] for p in outer) / len(outer)

    # Walk every adjacent pair of vertices in every loop (outer + holes) once.
    points = []
    seen = set()
    for loop in loops:
        n = len(loop)
        for i in range(n):
            a = loop[i]
            b = loop[(i + 1) % n]
            key = gmstr.edge_key(a, b)
            if key in seen:
                continue
            seen.add(key)
            mats = edge_materials.get(key, set())
            if not (mats & neighbour_materials):
                continue
            ax, ay = coords[a]
            bx, by = coords[b]
            ex = bx - ax
            ey = by - ay
            L = math.hypot(ex, ey)
            if L < 1.0e-12:
                continue
            tx, ty = ex / L, ey / L
            # Two candidate inward normals; pick the one that points toward the centroid.
            mid_x = 0.5 * (ax + bx)
            mid_y = 0.5 * (ay + by)
            n1x, n1y = -ty, tx
            d1 = (cx - mid_x) * n1x + (cy - mid_y) * n1y
            if d1 < 0:
                n1x, n1y = -n1x, -n1y
            nsamples = max(1, int(round(L / edge_sample_spacing)))
            for k in range(nsamples):
                t = (k + 0.5) / nsamples
                bx_t = ax + t * ex
                by_t = ay + t * ey
                for d in layer_distances:
                    px = bx_t + d * n1x
                    py = by_t + d * n1y
                    if not point_in_polygon_xy((px, py), outer):
                        continue
                    if any(point_in_polygon_xy((px, py), h) for h in holes):
                        continue
                    edge_dist = min_polygon_edge_distance((px, py), outer)
                    for h in holes:
                        edge_dist = min(edge_dist, min_polygon_edge_distance((px, py), h))
                    if edge_dist < sliver_factor * d:
                        continue
                    points.append((px, py))
    return points


def diagnose_solution_drift(old_coords, old_regions, old_triangles, old_solutions,
                            new_coords, converted_tris, new_values, sol_ids):
    old_int = gmstr.material_solution_integrals(old_coords, old_regions, old_triangles, old_solutions, sol_ids)
    new_int = gmstr.material_solution_integrals(new_coords, old_regions, converted_tris, new_values, sol_ids)
    drift = {}
    for material, old_vals in old_int.items():
        if material == GAS_MATERIAL:
            continue
        for i, sol in enumerate(sol_ids):
            old_total = old_vals[i]
            new_total = new_int.get(material, [0.0] * len(sol_ids))[i]
            denom = max(abs(old_total), abs(new_total), 1.0e-30)
            drift[(material, sol)] = (new_total - old_total) / denom
    return drift


def run_layered_remesh(args):
    coords, regions, triangles, node_solutions, sol_line, implant_line, model_line, _header = gmstr.parse_str(args.input)
    deck = parse_deck_mesh_lines(args.deck)
    if not deck["x"] or not deck["y"]:
        raise RuntimeError(f"deck {args.deck} did not yield line x and y commands")

    polygons = region_polygons(coords, regions, triangles)
    if not polygons:
        raise RuntimeError("no region polygons recovered from input .str")

    xmin, xmax, ymin, ymax = structure_bbox(coords)
    silicon_top = material_y_top(coords, regions, triangles, args.silicon_material)
    if silicon_top is None:
        silicon_top = ymin

    si_y, si_sy = build_silicon_y_profile(deck["y"], silicon_top, ymax,
                                          interface_buffer=args.interface_buffer)
    ox_y, ox_sy = build_above_silicon_y_profile(silicon_top, ymin, args.oxide_y_spacing,
                                                interface_buffer=args.interface_buffer)
    x_coords = build_x_profile_for_str(deck["x"], xmin, xmax)
    if not x_coords or not si_y:
        raise RuntimeError("empty x or silicon-y profile after applying deck spacings")

    silicon_targets = {args.silicon_material}
    oxide_targets = {SiO2_MATERIAL, SiNi_MATERIAL}
    above_targets = set(args.target_materials or DEFAULT_TARGET_MATERIALS) - silicon_targets
    boundary_xy = boundary_vertex_xy(coords, polygons)

    interior_by_region = clip_grid_to_regions(
        x_coords, si_y, si_sy, coords, polygons, regions,
        silicon_targets, boundary_xy, sliver_factor=args.sliver_factor,
    )
    if above_targets and ox_y:
        upper = clip_grid_to_regions(
            x_coords, ox_y, ox_sy, coords, polygons, regions,
            above_targets, boundary_xy, sliver_factor=args.sliver_factor,
        )
        for region, pts in upper.items():
            interior_by_region.setdefault(region, []).extend(pts)

    # Resolve oxide mode (legacy --light-oxide-remesh -> preserve-interior).
    oxide_mode = args.oxide_mode
    if args.light_oxide_remesh and oxide_mode == "boundary-only":
        oxide_mode = "preserve-interior"

    def loop_is_simple(loop):
        return len(set(loop)) == len(loop)
    silicon_region_ids = {r for r in polygons if regions.get(r) in silicon_targets}
    remesh_region_ids = set(polygons)  # everything in Gmsh by default

    extra_embed_pids_by_region = defaultdict(list)
    oxide_layer_points = defaultdict(list)

    if oxide_mode == "interface-layers":
        edge_materials = build_interface_edge_map(triangles, regions)
        edge_spacing = args.oxide_interface_edge_spacing
        if edge_spacing is None and len(x_coords) >= 2:
            edge_spacing = x_coords[1] - x_coords[0]
        if edge_spacing is None or edge_spacing <= 0:
            edge_spacing = 1.0
        n_layers = max(1, args.oxide_interface_layers)
        d0 = args.oxide_interface_distance
        step = max(1.0, args.oxide_interface_step)
        layer_distances = [d0 * (step ** k) for k in range(n_layers)]
        # Neighbour materials whose interface we layer against. Silicon is the
        # primary driver (oxidation moves the Si/SiO2 boundary); poly likewise.
        neighbour_materials = {SILICON_MATERIAL, POLY_MATERIAL}
        for r in polygons:
            if regions.get(r) not in oxide_targets:
                continue
            pts = generate_interface_layer_points(
                r, coords, polygons, regions, triangles, edge_materials,
                neighbour_materials, layer_distances, edge_spacing,
                sliver_factor=args.sliver_factor,
            )
            if pts:
                oxide_layer_points[r] = pts
                interior_by_region.setdefault(r, []).extend(pts)
        print(f"oxide-mode=interface-layers: layers={layer_distances} "
              f"edge_spacing={edge_spacing:.4g}; "
              f"per-region counts=" + " ".join(f"r{r}={len(v)}"
                                               for r, v in sorted(oxide_layer_points.items())))

    elif oxide_mode == "preserve-interior":
        oxide_region_ids = set()
        skipped_pinched = []
        for r in polygons:
            if regions.get(r) not in oxide_targets:
                continue
            if all(loop_is_simple(lp) for lp in polygons[r]):
                oxide_region_ids.add(r)
            else:
                skipped_pinched.append(r)
        if skipped_pinched:
            print(f"preserve (pinched-loop fallback) regions: {skipped_pinched}")
            remesh_region_ids -= set(skipped_pinched)

        region_triangle_vertices = defaultdict(set)
        for _tid, region, verts in triangles:
            region_triangle_vertices[region].update(verts)
        for region in oxide_region_ids:
            loops = polygons.get(region, [])
            if not loops:
                continue
            loop_verts = set()
            for loop in loops:
                loop_verts.update(loop)
            outer_xy = [coords[v] for v in loops[0]]
            holes_xy = [[coords[v] for v in loop] for loop in loops[1:]]
            candidates = region_triangle_vertices[region] - loop_verts
            embed_tol = 1.0e-7
            kept = []
            for vid in candidates:
                x, y = coords[vid]
                if not point_in_polygon_xy((x, y), outer_xy):
                    continue
                if any(point_in_polygon_xy((x, y), h) for h in holes_xy):
                    continue
                edge_dist = min_polygon_edge_distance((x, y), outer_xy)
                for h in holes_xy:
                    edge_dist = min(edge_dist, min_polygon_edge_distance((x, y), h))
                if edge_dist < embed_tol:
                    continue
                kept.append(vid)
            if kept:
                extra_embed_pids_by_region[region].extend(sorted(kept))

    all_y_spacings = list(si_sy) + list(ox_sy)
    point_size = max(max(all_y_spacings), x_coords[-1] - x_coords[0]) * 5.0
    geo_path = args.out_prefix.with_suffix(".geo")
    msh_path = args.out_prefix.with_suffix(".msh")

    transfinite_regions = None
    if args.silicon_mesher == "transfinite":
        transfinite_regions = {}
        for region in silicon_region_ids:
            loops = polygons.get(region, [])
            if not loops:
                continue
            ok, corners, sides = classify_rectangular_polygon(loops[0], coords)
            if not ok:
                print(f"silicon-mesher=transfinite: region {region} is not a simple "
                      f"axis-aligned rectangle; falling back to Delaunay for it")
                continue
            transfinite_regions[region] = {
                "corners": corners,
                "sides": sides,
            }
        if transfinite_regions:
            print(f"silicon-mesher=transfinite: applying to regions {sorted(transfinite_regions)}")

    surface_region = write_layered_geo(
        geo_path, coords, polygons, interior_by_region, regions, point_size,
        remesh_regions=remesh_region_ids,
        extra_embed_pids_by_region=extra_embed_pids_by_region,
        transfinite_regions=transfinite_regions,
        gmsh_algo=args.gmsh_algo,
    )

    cmd = [
        args.gmsh, "-2", "-format", "msh2",
        "-algo", args.gmsh_algo,
        "-o", str(msh_path), str(geo_path),
    ]
    proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          timeout=args.gmsh_timeout)
    if proc.returncode != 0:
        print(proc.stdout)
        return proc.returncode

    nodes, out_tris = gmstr.parse_msh(msh_path)
    qualities = [gmstr.triangle_quality(nodes, verts) for _eid, _tags, verts in out_tris]
    print(f"geo: {geo_path}")
    print(f"msh: {msh_path}")
    print(f"deck x layers={len(x_coords)} silicon-y layers={len(si_y)} above-silicon-y layers={len(ox_y)}")
    print(f"interior points: " + " ".join(f"r{r}={len(v)}" for r, v in sorted(interior_by_region.items())))
    print(f"output nodes={len(nodes)} triangles={len(out_tris)}")
    if qualities:
        print(f"quality min={min(qualities):.6g} avg={sum(qualities)/len(qualities):.6g}")

    if args.write_str:
        log_mats = set(args.log_interp_materials or DEFAULT_LOG_INTERP_MATERIALS)
        if oxide_mode != "preserve-interior":
            # Default test6-equivalent path: all regions came from Gmsh; the
            # existing writer applies per-material conservation scaling and
            # handles gas nodes correctly.
            bulk_only_materials = None
            bulk_only_neighbours = None
            if args.conservation_scope == "bulk-only":
                bulk_only_materials = oxide_targets  # SiO2, SiNi
                bulk_only_neighbours = silicon_targets | {POLY_MATERIAL}
            gmstr.write_remeshed_str(
                args.write_str,
                nodes, out_tris, surface_region,
                regions, coords, triangles, node_solutions,
                sol_line, implant_line, model_line,
                conserve_dopants=args.conserve_dopants,
                verbose=args.verbose,
                log_materials=log_mats,
                bulk_only_materials=bulk_only_materials,
                bulk_only_neighbours=bulk_only_neighbours,
            )
            print(f"str: {args.write_str}")
            if args.verbose:
                new_coords, new_regions, new_tris, new_solutions, _sl, _il, _ml, _h = gmstr.parse_str(args.write_str)
                sol_ids = gmstr.solution_ids(sol_line)
                converted = [(t[1], t[2]) for t in new_tris]
                drift = diagnose_solution_drift(
                    coords, regions, triangles, node_solutions,
                    {pid: new_coords[pid] for pid in new_coords},
                    converted, new_solutions, sol_ids,
                )
                worst = sorted(drift.items(), key=lambda kv: -abs(kv[1]))[:8]
                print("solution drift (top 8 by |rel|):")
                for (mat, sol), rd in worst:
                    print(f"  material={mat} sol={sol} drift={rd:+.3e}")
            return 0
        coords_merged, tris_merged, gmsh_to_merged = merge_silicon_with_preserved(
            coords, regions, triangles, surface_region,
            nodes, out_tris, remesh_region_ids,
        )

        sol_ids = gmstr.solution_ids(sol_line)
        used_materials = defaultdict(set)
        for region, verts in tris_merged:
            material = regions.get(region, -1)
            for v in verts:
                used_materials[v].add(material)

        # Mirror gmsh_remesh_str.write_remeshed_str: any vertex on a single-tri
        # boundary edge that isn't a reflecting side must also have a gas
        # (material 0) node entry, otherwise SUPREM's mesh self-test fails with
        # "exposed point N has no gas node".
        edge_count = defaultdict(int)
        for _region, verts in tris_merged:
            for a, b in ((verts[1], verts[2]), (verts[2], verts[0]), (verts[1], verts[0])):
                edge_count[gmstr.edge_key(a, b)] += 1
        xs = [p[0] for p in coords_merged.values()]
        ys = [p[1] for p in coords_merged.values()]
        xmin_c, xmax_c = min(xs), max(xs)
        ymax_c = max(ys)
        bc_tol = 1.0e-8
        for _region, verts in tris_merged:
            for a, b in ((verts[1], verts[2]), (verts[2], verts[0]), (verts[1], verts[0])):
                if edge_count[gmstr.edge_key(a, b)] != 1:
                    continue
                ax, ay = coords_merged[a]
                bx, by = coords_merged[b]
                if (abs(ax - xmin_c) < bc_tol and abs(bx - xmin_c) < bc_tol) or \
                   (abs(ax - xmax_c) < bc_tol and abs(bx - xmax_c) < bc_tol):
                    continue  # reflecting sides — no gas
                used_materials[a].add(0)
                used_materials[b].add(0)

        old_lookup = gmstr.build_old_lookup(coords, regions, triangles)
        stats = defaultdict(int)
        # Vertex IDs that came from the input (preserved or interface-shared);
        # for these we copy the original solution values verbatim so oxide and
        # aluminum stay bit-identical.
        original_pids = set(coords.keys()) & set(coords_merged.keys())
        node_values = {}
        for gid, materials in used_materials.items():
            for material in materials:
                if gid in original_pids and (gid, material) in node_solutions:
                    node_values[(gid, material)] = list(node_solutions[(gid, material)])
                else:
                    node_values[(gid, material)] = gmstr.interpolate_solution(
                        coords_merged[gid], material, old_lookup, node_solutions,
                        coords, sol_ids, stats, log_materials=log_mats,
                    )

        if args.conserve_dopants:
            gmstr.conserve_dopant_integrals(
                coords, regions, triangles, node_solutions,
                coords_merged, tris_merged, node_values, sol_ids, args.verbose,
            )

        write_combined_str(
            args.write_str, coords_merged, tris_merged, regions, node_values,
            sol_line, implant_line, model_line,
        )
        print(f"str: {args.write_str}")
        print(f"merged: preserved-tris={sum(1 for r, _ in tris_merged if r not in remesh_region_ids)} "
              f"new-silicon-tris={sum(1 for r, _ in tris_merged if r in remesh_region_ids)}")
        print(f"interp stats: exact={stats['interpolated']} nearest={stats['nearest_fallback']} "
              f"preserved-original={sum(1 for (g, _m) in node_values if g in original_pids)}")

        if args.verbose:
            # Re-read for drift diagnostic.
            new_coords, new_regions, new_tris, new_solutions, _sl, _il, _ml, _h = gmstr.parse_str(args.write_str)
            converted = [(t[1], t[2]) for t in new_tris]
            drift = diagnose_solution_drift(
                coords, regions, triangles, node_solutions,
                {pid: new_coords[pid] for pid in new_coords},
                converted, new_solutions, sol_ids,
            )
            worst = sorted(drift.items(), key=lambda kv: -abs(kv[1]))[:8]
            print("solution drift (top 8 by |rel|):")
            for (mat, sol), rd in worst:
                print(f"  material={mat} sol={sol} drift={rd:+.3e}")
    return 0


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="input .str")
    parser.add_argument("--deck", type=Path, required=True, help="SUPREM deck (.in) whose line x/y commands drive the profile")
    parser.add_argument("--out-prefix", type=Path)
    parser.add_argument("--write-str", type=Path)
    parser.add_argument("--silicon-material", type=int, default=SILICON_MATERIAL)
    parser.add_argument("--target-material", type=int, action="append", dest="target_materials",
                        help="material id to receive deck-driven interior points (may repeat)")
    parser.add_argument("--oxide-y-spacing", type=float, default=0.0,
                        help="uniform y spacing (um) for material above silicon top; 0 disables (default, "
                             "since SUPREM panics during downstream etch/oxidation when oxide interior nodes "
                             "exist)")
    parser.add_argument("--interface-buffer", type=float, default=0.02,
                        help="distance (um) from the silicon-oxide interface within which no new layered "
                             "points are added on either side. Leaves the interface band at original "
                             "SUPREM-mesh density to avoid the aspect-ratio cliff that triggers triang::chop "
                             "and interp panics during downstream oxidation. Default 0.02 (20 nm).")
    parser.add_argument("--silicon-mesher", choices=("delaunay", "transfinite"), default="delaunay",
                        help="silicon meshing engine. 'delaunay' (default) uses Gmsh constrained Delaunay "
                             "with embedded layered points. 'transfinite' uses Gmsh Transfinite Surface to "
                             "produce a deterministic Cartesian-split-into-triangles grid (matches SUPREM's "
                             "own rect.c output style). Transfinite requires silicon to be an axis-aligned "
                             "rectangle; falls back to delaunay per region if not.")
    parser.add_argument("--oxide-mode", choices=("boundary-only", "interface-layers", "preserve-interior"),
                        default="boundary-only",
                        help="how to mesh oxide/nitride regions: "
                             "boundary-only=test6-equivalent, survives downstream oxidation but ~33%% oxide-B "
                             "drift at pmos.str; "
                             "interface-layers=structured layers normal to silicon-oxide interface (1-2 layers, "
                             "controlled spacing); aims to lower drift without triggering the solver segfault; "
                             "preserve-interior=embed all original oxide interior vertices, lowest drift but "
                             "SEGFAULTS DOWNSTREAM if any further oxide growth follows the remesh — only safe "
                             "for post-thermal device-extraction.")
    parser.add_argument("--oxide-interface-layers", type=int, default=1,
                        help="number of interface layers when --oxide-mode=interface-layers (default 1)")
    parser.add_argument("--oxide-interface-distance", type=float, default=0.05,
                        help="perpendicular distance (um) for the first interface layer (default 0.05)")
    parser.add_argument("--oxide-interface-step", type=float, default=2.0,
                        help="multiplier between successive interface-layer distances (default 2.0, so "
                             "layers at d, 2d, 4d, ...)")
    parser.add_argument("--oxide-interface-edge-spacing", type=float, default=None,
                        help="along-edge spacing for interface layer points; defaults to the deck x spacing")
    parser.add_argument("--light-oxide-remesh", action="store_true",
                        help="deprecated alias for --oxide-mode preserve-interior")
    parser.add_argument("--log-interp-material", type=int, action="append", dest="log_interp_materials",
                        help="material id where dopant interp uses log(C+1) (default silicon+poly); other "
                             "materials use arithmetic interpolation so triangle integrals stay exact")
    parser.add_argument("--sliver-factor", type=float, default=0.4,
                        help="drop interior points within this fraction of local y spacing of any polygon edge")
    parser.add_argument("--conserve-dopants", action=argparse.BooleanOptionalAction, default=True,
                        help="apply per-material dopant-integral scaling after interpolation (on by default)")
    parser.add_argument("--conservation-scope", choices=("material", "bulk-only"), default="material",
                        help="scaling target: material=uniform per material (default, today's behavior); "
                             "bulk-only=hold silicon-oxide interface vertex values fixed and scale only "
                             "bulk oxide vertices to match the integral, avoiding boron-enhanced "
                             "oxidation perturbation downstream.")
    parser.add_argument("--gmsh", default="gmsh")
    parser.add_argument("--gmsh-timeout", type=float, default=120.0)
    parser.add_argument("--gmsh-algo", default="front2d",
                        choices=("auto", "meshadapt", "del2d", "front2d", "delquad", "quadqs",
                                 "initial2d"),
                        help="Gmsh 2D meshing algorithm. 'del2d' (the old default) is pure "
                             "Delaunay and scales poorly with many constrained points. 'front2d' "
                             "(Frontal-Delaunay, new default) is the recommended general-purpose "
                             "gmsh algorithm: faster on dense constrained inputs and produces "
                             "better-quality triangles near interfaces. 'meshadapt' is the most "
                             "robust legacy option if 'front2d' has trouble.")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.out_prefix is None:
        args.out_prefix = args.input.with_suffix("")
    return run_layered_remesh(args)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
