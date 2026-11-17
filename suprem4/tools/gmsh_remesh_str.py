#!/usr/bin/env python3
"""Prototype SUPREM .str -> Gmsh constrained remesh export.

It builds one conforming multi-region Gmsh geometry from the material
boundaries in a .str file, meshes it, reports geometry/quality metrics, and
can write a remeshed .str with same-material interpolated node solutions.
"""

import argparse
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


LOG_INTERPOLATED_IMPURITIES = {
    2,   # As
    3,   # P
    4,   # Sb
    5,   # B
    20,  # Asa
    21,  # Pa
    22,  # Sba
    23,  # Ba
    31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48,
}


def parse_str(path):
    coords = {}
    regions = {}
    triangles = []
    node_solutions = {}
    sol_line = None
    implant_line = None
    model_line = "M 0 1.173000e+03"
    header = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "c" and len(parts) >= 4:
                coords[int(parts[1])] = (float(parts[2]), float(parts[3]))
            elif parts[0] == "r" and len(parts) >= 3:
                regions[int(parts[1])] = int(parts[2])
            elif parts[0] == "t" and len(parts) >= 6:
                triangles.append(
                    (int(parts[1]), int(parts[2]), (int(parts[3]), int(parts[4]), int(parts[5])))
                )
            elif parts[0] == "n" and len(parts) >= 4:
                # SUPREM .str files use zero-based point ids in n records,
                # while c/t records use one-based point ids.
                node_solutions[(int(parts[1]) + 1, int(parts[2]))] = [float(x) for x in parts[3:]]
            elif parts[0] == "s":
                sol_line = line.rstrip()
                header.append(line.rstrip())
            elif parts[0] == "I":
                implant_line = line.rstrip()
            elif parts[0] == "M":
                model_line = line.rstrip()
                header.append(line.rstrip())
            elif parts[0] in {"v", "D", "s", "M"}:
                header.append(line.rstrip())
    return coords, regions, triangles, node_solutions, sol_line, implant_line, model_line, header


def solution_ids(sol_line):
    if sol_line is None:
        raise RuntimeError("input .str has no solution header")
    parts = sol_line.split()
    return [int(x) for x in parts[2:]]


def edge_key(a, b):
    return (a, b) if a < b else (b, a)


def signed_area(coords, loop):
    area = 0.0
    for a, b in zip(loop, loop[1:] + loop[:1]):
        ax, ay = coords[a]
        bx, by = coords[b]
        area += ax * by - bx * ay
    return 0.5 * area


def polygon_contains(coords, outer, inner):
    x, y = coords[inner[0]]
    inside = False
    for a, b in zip(outer, outer[1:] + outer[:1]):
        ax, ay = coords[a]
        bx, by = coords[b]
        if (ay > y) != (by > y):
            xinters = (bx - ax) * (y - ay) / (by - ay) + ax
            if x < xinters:
                inside = not inside
    return inside


def build_region_boundary(triangles):
    by_region = defaultdict(list)
    for _, region, verts in triangles:
        for a, b in ((verts[0], verts[1]), (verts[1], verts[2]), (verts[2], verts[0])):
            by_region[region].append((a, b))

    boundaries = {}
    for region, edges in by_region.items():
        counts = defaultdict(int)
        oriented = {}
        for a, b in edges:
            key = edge_key(a, b)
            counts[key] += 1
            oriented[key] = (a, b)
        boundaries[region] = [oriented[key] for key, count in counts.items() if count == 1]
    return boundaries


def select_top_material_points(coords, regions, triangles, material, depth):
    span = material_y_span(coords, regions, triangles, material)
    if span is None:
        return {}
    top_y, _bottom_y = span
    limit_y = top_y + depth
    selected = defaultdict(set)
    for _tid, region, verts in triangles:
        if regions.get(region) != material:
            continue
        for point in verts:
            if coords[point][1] <= limit_y:
                selected[region].add(point)
    return selected


def measured_top_material_y_spacing(coords, regions, triangles, material, depth):
    span = material_y_span(coords, regions, triangles, material)
    if span is None:
        return None
    top_y, _bottom_y = span
    limit_y = top_y + depth
    dys = []
    for _tid, region, verts in triangles:
        if regions.get(region) != material:
            continue
        for a, b in ((verts[0], verts[1]), (verts[1], verts[2]), (verts[2], verts[0])):
            ya = coords[a][1]
            yb = coords[b][1]
            if ya > limit_y or yb > limit_y:
                continue
            dy = abs(ya - yb)
            if dy > 1.0e-9:
                dys.append(dy)
    if not dys:
        return None
    dys.sort()
    return dys[max(0, min(len(dys) - 1, len(dys) // 10))]


def containing_region(point, coords, regions, triangles, material):
    for _tid, region, verts in triangles:
        if regions.get(region) != material:
            continue
        pts = [coords[v] for v in verts]
        if point_in_tri(point, pts, 1.0e-8) is not None:
            return region
    return None


def add_top_grid_points(coords, regions, triangles, material, depth, x_size, y_top_size, y_bottom_size):
    span = material_y_span(coords, regions, triangles, material)
    if span is None:
        return {}, defaultdict(set)
    top_y, _bottom_y = span
    xs = []
    for _tid, region, verts in triangles:
        if regions.get(region) == material:
            xs.extend(coords[v][0] for v in verts)
    if not xs:
        return {}, defaultdict(set)
    xmin, xmax = min(xs), max(xs)
    x_size = max(x_size, 1.0e-6)
    y_top_size = max(y_top_size, 1.0e-6)
    y_bottom_size = max(y_bottom_size, y_top_size)

    y_levels = []
    y = top_y + y_top_size
    step = y_top_size
    limit = top_y + depth
    while y < limit - 1.0e-9 and len(y_levels) < 200:
        y_levels.append(y)
        step = min(step * 1.6, y_bottom_size)
        y += step

    new_coords = {}
    by_region = defaultdict(set)
    next_id = max(coords) + 1
    existing = list(coords.values())
    x = xmin + x_size
    while x < xmax - 1.0e-9:
        for y in y_levels:
            if any((px - x) ** 2 + (py - y) ** 2 < 1.0e-12 for px, py in existing):
                continue
            region = containing_region((x, y), coords, regions, triangles, material)
            if region is None:
                continue
            new_coords[next_id] = (x, y)
            by_region[region].add(next_id)
            next_id += 1
        x += x_size
    return new_coords, by_region


def loops_from_edges(edges):
    unused = {edge_key(a, b): (a, b) for a, b in edges}
    loops = []
    while unused:
        _, (a, b) = unused.popitem()
        loop = [a, b]
        while True:
            tail = loop[-1]
            found = None
            for key, (u, v) in unused.items():
                if u == tail:
                    found = (key, v)
                    break
                if v == tail:
                    found = (key, u)
                    break
            if found is None:
                break
            key, nxt = found
            unused.pop(key)
            if nxt == loop[0]:
                break
            loop.append(nxt)
        if len(loop) >= 3:
            loops.append(loop)
    return loops


def material_y_span(coords, regions, triangles, material):
    ys = []
    for _tid, region, verts in triangles:
        if regions.get(region) == material:
            ys.extend(coords[v][1] for v in verts)
    if not ys:
        return None
    return min(ys), max(ys)


def graded_size(y, top_y, bottom_y, top_size, bottom_size, power=1.0):
    if bottom_y <= top_y:
        return top_size
    t = (y - top_y) / (bottom_y - top_y)
    t = min(max(t, 0.0), 1.0)
    t = t ** power
    return top_size + t * (bottom_size - top_size)


def write_geo(
    coords,
    regions,
    triangles,
    geo_path,
    mesh_size,
    top_size=None,
    bottom_size=None,
    grade_top=None,
    grade_bottom=None,
    grade_power=1.0,
    silicon_profile=False,
    silicon_material=3,
    active_depth=None,
    preserve_top_material=False,
    preserve_depth=None,
    anisotropic=False,
    x_size=None,
    y_top_size=None,
    y_bottom_size=None,
    structure_top_spacing=False,
    embed_top_grid=False,
):
    geo_coords = dict(coords)
    boundaries = build_region_boundary(triangles)
    boundary_points = {p for edges in boundaries.values() for edge in edges for p in edge}
    embedded_points_by_region = defaultdict(set)
    if preserve_top_material:
        depth = active_depth if preserve_depth is None else preserve_depth
        if depth is None:
            depth = 2.0
        embedded_points_by_region = select_top_material_points(coords, regions, triangles, silicon_material, depth)
    if embed_top_grid:
        depth = active_depth if active_depth is not None else 2.0
        hx = mesh_size if x_size is None else x_size
        hy0 = y_top_size
        if hy0 is None:
            hy0 = measured_top_material_y_spacing(coords, regions, triangles, silicon_material, depth)
        if hy0 is None:
            hy0 = mesh_size
        hy1 = bottom_size if y_bottom_size is None else y_bottom_size
        new_coords, grid_by_region = add_top_grid_points(coords, regions, triangles, silicon_material, depth, hx, hy0, hy1)
        geo_coords.update(new_coords)
        for region, points in grid_by_region.items():
            embedded_points_by_region[region].update(points)
    embedded_points = {p for points in embedded_points_by_region.values() for p in points}
    used_points = sorted(boundary_points | embedded_points)
    line_ids = {}
    line_defs = []
    next_line = 1
    loop_defs = []
    surface_defs = []
    physical_defs = defaultdict(list)
    surface_region = {}
    next_loop = 1
    next_surface = 1

    for region in sorted(boundaries):
        loops = loops_from_edges(boundaries[region])
        loops.sort(key=lambda item: abs(signed_area(geo_coords, item)), reverse=True)
        if not loops:
            continue

        gmsh_loops = []
        for loop in loops:
            signed_lines = []
            for a, b in zip(loop, loop[1:] + loop[:1]):
                key = edge_key(a, b)
                if key not in line_ids:
                    line_ids[key] = next_line
                    line_defs.append((next_line, key[0], key[1]))
                    next_line += 1
                lid = line_ids[key]
                signed_lines.append(lid if (a, b) == key else -lid)
            loop_id = next_loop
            next_loop += 1
            loop_defs.append((loop_id, signed_lines))
            gmsh_loops.append((loop_id, loop))

        # Create one surface per outer loop, assigning enclosed smaller loops
        # as holes.  This handles disjoint islands conservatively.
        assigned_holes = set()
        for loop_id, loop in gmsh_loops:
            if any(
                other_id != loop_id and abs(signed_area(geo_coords, other_loop)) > abs(signed_area(geo_coords, loop))
                and polygon_contains(geo_coords, other_loop, loop)
                for other_id, other_loop in gmsh_loops
            ):
                continue
            holes = []
            for hole_id, hole in gmsh_loops:
                if hole_id == loop_id:
                    continue
                if polygon_contains(geo_coords, loop, hole):
                    holes.append(hole_id)
                    assigned_holes.add(hole_id)
            surf = next_surface
            next_surface += 1
            surface_defs.append((surf, [loop_id] + holes, region))
            surface_region[surf] = region
            physical_defs[regions.get(region, -1)].append(surf)

    profile_note = None
    if silicon_profile:
        span = material_y_span(coords, regions, triangles, silicon_material)
        if span is None:
            raise RuntimeError(f"could not find silicon material {silicon_material} in input structure")
        silicon_top, silicon_bottom = span
        top_y = silicon_top if grade_top is None else grade_top
        if active_depth is not None:
            bottom_y = top_y + active_depth
        else:
            bottom_y = silicon_bottom if grade_bottom is None else grade_bottom
        top_size = mesh_size if top_size is None else top_size
        bottom_size = mesh_size if bottom_size is None else bottom_size
        profile_note = (
            f"// silicon-profile material={silicon_material} silicon_y=[{silicon_top:.12g},{silicon_bottom:.12g}] "
            f"grade_y=[{top_y:.12g},{bottom_y:.12g}] size=[{top_size:.12g},{bottom_size:.12g}] power={grade_power:.12g}\n"
        )
    elif top_size is not None or bottom_size is not None:
        ys = [geo_coords[pid][1] for pid in used_points]
        top_y = min(ys) if grade_top is None else grade_top
        bottom_y = max(ys) if grade_bottom is None else grade_bottom
        top_size = mesh_size if top_size is None else top_size
        bottom_size = mesh_size if bottom_size is None else bottom_size
        profile_note = (
            f"// vertical-profile grade_y=[{top_y:.12g},{bottom_y:.12g}] "
            f"size=[{top_size:.12g},{bottom_size:.12g}] power={grade_power:.12g}\n"
        )
    else:
        top_y = bottom_y = None

    with open(geo_path, "w", encoding="utf-8") as handle:
        handle.write('SetFactory("Built-in");\n')
        if profile_note:
            handle.write(profile_note)
        for pid in used_points:
            x, y = geo_coords[pid]
            size = mesh_size
            if top_y is not None:
                size = graded_size(y, top_y, bottom_y, top_size, bottom_size, grade_power)
            handle.write(f"Point({pid}) = {{{x:.12g}, {y:.12g}, 0, {size:.12g}}};\n")
        for lid, a, b in line_defs:
            handle.write(f"Line({lid}) = {{{a}, {b}}};\n")
        for loop_id, signed_lines in loop_defs:
            handle.write(f"Curve Loop({loop_id}) = {{{', '.join(map(str, signed_lines))}}};\n")
        for surf, loops, _region in surface_defs:
            handle.write(f"Plane Surface({surf}) = {{{', '.join(map(str, loops))}}};\n")
        for material, surfaces in sorted(physical_defs.items()):
            handle.write(f"Physical Surface({material}) = {{{', '.join(map(str, surfaces))}}};\n")
        for surf, _loops, region in surface_defs:
            points = sorted(embedded_points_by_region.get(region, set()) - boundary_points)
            if points:
                handle.write(f"Point{{{', '.join(map(str, points))}}} In Surface{{{surf}}};\n")
        if anisotropic:
            if top_y is None:
                ys = [geo_coords[pid][1] for pid in used_points]
                top_y = min(ys)
                bottom_y = max(ys)
            hx = mesh_size if x_size is None else x_size
            if structure_top_spacing and y_top_size is None:
                depth = active_depth if active_depth is not None else 2.0
                y_top_size = measured_top_material_y_spacing(coords, regions, triangles, silicon_material, depth)
            hy0 = top_size if y_top_size is None else y_top_size
            hy1 = bottom_size if y_bottom_size is None else y_bottom_size
            if bottom_y <= top_y:
                bottom_y = top_y + 1.0
            yexpr = (
                f"({hy0:.12g} + ({hy1:.12g} - {hy0:.12g}) * "
                f"Min(Max((y - {top_y:.12g}) / ({bottom_y:.12g} - {top_y:.12g}), 0), 1))"
            )
            handle.write("Field[1] = MathEvalAniso;\n")
            handle.write(f'Field[1].M11 = "{1.0 / (hx * hx):.12g}";\n')
            handle.write('Field[1].M12 = "0";\n')
            handle.write(f'Field[1].M22 = "1 / ({yexpr} * {yexpr})";\n')
            handle.write("Background Field = 1;\n")
            handle.write("Mesh.SmoothRatio = 3;\n")
            handle.write("Mesh.AnisoMax = 1000;\n")
            handle.write("Mesh.Algorithm = 7;\n")
        handle.write("Coherence;\n")

    return len(used_points), len(line_defs), len(surface_defs), surface_region


def parse_msh(path):
    nodes = {}
    tris = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = iter(handle)
        for line in lines:
            if line.strip() == "$Nodes":
                n = int(next(lines))
                for _ in range(n):
                    parts = next(lines).split()
                    nodes[int(parts[0])] = (float(parts[1]), float(parts[2]))
            elif line.strip() == "$Elements":
                n = int(next(lines))
                for _ in range(n):
                    parts = next(lines).split()
                    eid = int(parts[0])
                    etype = int(parts[1])
                    ntags = int(parts[2])
                    if etype == 2:
                        tags = [int(x) for x in parts[3 : 3 + ntags]]
                        verts = tuple(int(x) for x in parts[3 + ntags : 3 + ntags + 3])
                        tris.append((eid, tags, verts))
    return nodes, tris


def triangle_quality(nodes, verts):
    pts = [nodes[v] for v in verts]
    lengths = []
    for i in range(3):
        ax, ay = pts[i]
        bx, by = pts[(i + 1) % 3]
        lengths.append((ax - bx) ** 2 + (ay - by) ** 2)
    area = abs(
        0.5
        * (
            (pts[1][0] - pts[0][0]) * (pts[2][1] - pts[0][1])
            - (pts[2][0] - pts[0][0]) * (pts[1][1] - pts[0][1])
        )
    )
    denom = sum(lengths)
    return 0.0 if denom <= 0 else 4.0 * math.sqrt(3.0) * area / denom


def tri_area(coords, verts):
    p0, p1, p2 = (coords[v] for v in verts)
    return abs((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])) / 2.0


def point_in_tri(point, tri_pts, eps=1.0e-10):
    x, y = point
    (x0, y0), (x1, y1), (x2, y2) = tri_pts
    den = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
    if abs(den) < eps:
        return None
    a = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) / den
    b = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) / den
    c = 1.0 - a - b
    if a >= -eps and b >= -eps and c >= -eps:
        return a, b, c
    return None


def build_old_lookup(coords, regions, triangles):
    by_material = defaultdict(list)
    for _tid, region, verts in triangles:
        mat = regions.get(region, -1)
        pts = [coords[v] for v in verts]
        box = (
            min(p[0] for p in pts),
            max(p[0] for p in pts),
            min(p[1] for p in pts),
            max(p[1] for p in pts),
        )
        by_material[mat].append((box, verts, pts))
    return by_material


def blend_solution(bary, vals, sol_ids, use_log_interp=True):
    out = []
    for i, sol_id in enumerate(sol_ids):
        if use_log_interp and sol_id in LOG_INTERPOLATED_IMPURITIES and all(v[i] > -1.0 for v in vals):
            # Match src/refine/interp.c: interpolate log(C + 1), then exp().
            out.append(
                math.exp(
                    bary[0] * math.log(vals[0][i] + 1.0)
                    + bary[1] * math.log(vals[1][i] + 1.0)
                    + bary[2] * math.log(vals[2][i] + 1.0)
                )
            )
        else:
            out.append(bary[0] * vals[0][i] + bary[1] * vals[1][i] + bary[2] * vals[2][i])
    return out


def interpolate_solution(point, material, old_lookup, old_solutions, old_coords, sol_ids, stats,
                         log_materials=None):
    use_log = True if log_materials is None else (material in log_materials)
    x, y = point
    best = None
    best_dist = float("inf")
    for box, verts, pts in old_lookup.get(material, []):
        if x < box[0] - 1.0e-9 or x > box[1] + 1.0e-9:
            continue
        if y < box[2] - 1.0e-9 or y > box[3] + 1.0e-9:
            continue
        bary = point_in_tri(point, pts)
        if bary is None:
            continue
        vals = [old_solutions.get((v, material)) for v in verts]
        if any(v is None for v in vals):
            continue
        stats["interpolated"] += 1
        return blend_solution(bary, vals, sol_ids, use_log_interp=use_log)

    for _box, verts, pts in old_lookup.get(material, []):
        for v, p in zip(verts, pts):
            vals = old_solutions.get((v, material))
            if vals is None:
                continue
            dist = (p[0] - x) ** 2 + (p[1] - y) ** 2
            if dist < best_dist:
                best_dist = dist
                best = vals
    if best is None:
        for (old_point, mat), vals in old_solutions.items():
            if mat != material or old_point not in old_coords:
                continue
            p = old_coords[old_point]
            dist = (p[0] - x) ** 2 + (p[1] - y) ** 2
            if dist < best_dist:
                best_dist = dist
                best = vals
    if best is None:
        raise RuntimeError(f"could not interpolate material {material} at {point}")
    stats["nearest_fallback"] += 1
    stats[f"nearest_fallback_material_{material}"] += 1
    return list(best)


def material_solution_integrals(coords, regions, triangles, values, sol_ids):
    integrals = defaultdict(lambda: [0.0 for _ in sol_ids])
    for item in triangles:
        if len(item) == 3:
            _tid, region, verts = item
        else:
            region, verts = item
        material = regions[region]
        node_vals = [values.get((v, material)) for v in verts]
        if any(v is None for v in node_vals):
            continue
        area = tri_area(coords, verts)
        for i in range(len(sol_ids)):
            integrals[material][i] += area * (node_vals[0][i] + node_vals[1][i] + node_vals[2][i]) / 3.0
    return integrals


def conserve_dopant_integrals(old_coords, old_regions, old_triangles, old_values,
                              new_coords, converted_tris, new_values, sol_ids, verbose):
    old_integrals = material_solution_integrals(old_coords, old_regions, old_triangles, old_values, sol_ids)
    new_integrals = material_solution_integrals(new_coords, old_regions, converted_tris, new_values, sol_ids)
    scales = {}
    for material, old_vals in old_integrals.items():
        if material == 0:
            continue
        for i, sol_id in enumerate(sol_ids):
            if sol_id not in LOG_INTERPOLATED_IMPURITIES:
                continue
            new_total = new_integrals.get(material, [0.0 for _ in sol_ids])[i]
            old_total = old_vals[i]
            if abs(new_total) <= 1.0e-30:
                continue
            scales[(material, i)] = old_total / new_total

    for (point, material), vals in new_values.items():
        for (scale_material, sol_index), scale in scales.items():
            if material == scale_material:
                vals[sol_index] *= scale

    if verbose:
        for (material, sol_index), scale in sorted(scales.items()):
            print(f"conserve material={material} sol={sol_ids[sol_index]} scale={scale:.8g}")


def _vertex_weights_from_triangles(coords, triangles, regions, material_filter):
    """Lumped-mass vertex weights for vertices of triangles whose region maps to
    a material in material_filter. weight_v = (1/3) * sum(area(T) for T touching v).
    """
    weights = defaultdict(float)
    for item in triangles:
        if len(item) == 3:
            _tid, region, verts = item
        else:
            region, verts = item
        if regions.get(region, -1) not in material_filter:
            continue
        a = tri_area(coords, verts)
        for v in verts:
            weights[v] += a / 3.0
    return weights


def _interface_vertices(triangles, regions, target_materials, neighbour_materials):
    """Return the set of vertices that lie on an edge shared between a triangle
    of `target_materials` and a triangle of `neighbour_materials`.
    """
    by_edge = defaultdict(set)
    for item in triangles:
        if len(item) == 3:
            _tid, region, verts = item
        else:
            region, verts = item
        material = regions.get(region, -1)
        for a, b in ((verts[0], verts[1]), (verts[1], verts[2]), (verts[2], verts[0])):
            by_edge[edge_key(a, b)].add(material)
    interface = set()
    for (a, b), mats in by_edge.items():
        if (mats & target_materials) and (mats & neighbour_materials):
            interface.add(a)
            interface.add(b)
    return interface


def conserve_dopant_integrals_bulk_only(old_coords, old_regions, old_triangles, old_values,
                                        new_coords, converted_tris, new_values, sol_ids,
                                        verbose,
                                        bulk_only_materials,
                                        neighbour_materials,
                                        max_scale=10.0, min_scale=0.1):
    """Conservation that holds interface vertex values fixed and absorbs all the
    correction in bulk vertex values, per material.

    For each material in `bulk_only_materials`:
      - identify interface vertices (vertices on an edge shared with a triangle
        of any material in `neighbour_materials`)
      - hold those vertex values fixed
      - scale the remaining (bulk) vertex values by a single factor per impurity
        so that the per-material integral matches the original mesh's integral
      - fall back to per-material uniform scaling if there are no bulk vertices
        or the bulk integral is too small to drive the correction

    Other materials use the existing per-material uniform scaling.
    """
    old_integrals = material_solution_integrals(old_coords, old_regions, old_triangles, old_values, sol_ids)
    interface_verts = _interface_vertices(converted_tris, old_regions, bulk_only_materials, neighbour_materials)
    weights_by_material = {
        m: _vertex_weights_from_triangles(new_coords, converted_tris, old_regions, {m})
        for m in bulk_only_materials
    }

    scales_uniform = {}
    scales_bulk = {}
    scaled_count = {"bulk": 0, "uniform": 0, "skip": 0}
    new_integrals = material_solution_integrals(new_coords, old_regions, converted_tris, new_values, sol_ids)

    for material in old_integrals:
        if material == 0:
            continue
        for i, sol_id in enumerate(sol_ids):
            if sol_id not in LOG_INTERPOLATED_IMPURITIES:
                continue
            old_total = old_integrals[material][i]
            new_total = new_integrals.get(material, [0.0 for _ in sol_ids])[i]
            if abs(new_total) <= 1.0e-30:
                scaled_count["skip"] += 1
                continue
            if material not in bulk_only_materials:
                scales_uniform[(material, i)] = old_total / new_total
                scaled_count["uniform"] += 1
                continue
            weights = weights_by_material.get(material, {})
            interface_integral = 0.0
            bulk_integral = 0.0
            for vid, w in weights.items():
                v = new_values.get((vid, material))
                if v is None:
                    continue
                if vid in interface_verts:
                    interface_integral += w * v[i]
                else:
                    bulk_integral += w * v[i]
            target_bulk = old_total - interface_integral
            if abs(bulk_integral) <= 1.0e-30 or target_bulk * bulk_integral < 0:
                scales_uniform[(material, i)] = old_total / new_total
                scaled_count["uniform"] += 1
                continue
            scale = target_bulk / bulk_integral
            if scale > max_scale or scale < min_scale:
                scales_uniform[(material, i)] = old_total / new_total
                scaled_count["uniform"] += 1
                continue
            scales_bulk[(material, i)] = scale
            scaled_count["bulk"] += 1

    for (point, material), vals in new_values.items():
        for (scale_material, sol_index), scale in scales_uniform.items():
            if material == scale_material:
                vals[sol_index] *= scale
        for (scale_material, sol_index), scale in scales_bulk.items():
            if material == scale_material and point not in interface_verts:
                vals[sol_index] *= scale

    # Post-conditions: per-material integrals match within numerical noise.
    post_integrals = material_solution_integrals(new_coords, old_regions, converted_tris, new_values, sol_ids)
    for material in old_integrals:
        if material == 0:
            continue
        for i, sol_id in enumerate(sol_ids):
            if sol_id not in LOG_INTERPOLATED_IMPURITIES:
                continue
            old_total = old_integrals[material][i]
            new_total = post_integrals.get(material, [0.0 for _ in sol_ids])[i]
            denom = max(abs(old_total), abs(new_total), 1.0e-30)
            rd = abs(new_total - old_total) / denom
            if rd > 1.0e-6 and abs(new_total) > 1.0e-20:
                raise AssertionError(
                    f"bulk-only conservation failed to match per-material integral: "
                    f"material={material} sol={sol_id} old={old_total:.6e} new={new_total:.6e} rd={rd:.3e}"
                )

    if verbose:
        for (material, sol_index), scale in sorted(scales_bulk.items()):
            print(f"conserve(bulk) material={material} sol={sol_ids[sol_index]} bulk_scale={scale:.8g}")
        for (material, sol_index), scale in sorted(scales_uniform.items()):
            print(f"conserve(uniform fallback) material={material} sol={sol_ids[sol_index]} scale={scale:.8g}")
        print(f"conservation outcome: bulk={scaled_count['bulk']} uniform={scaled_count['uniform']} "
              f"skipped={scaled_count['skip']} interface_verts={len(interface_verts)}")


def write_remeshed_str(
    path,
    nodes,
    out_tris,
    surface_region,
    old_regions,
    old_coords,
    old_triangles,
    old_solutions,
    sol_line,
    implant_line,
    model_line,
    conserve_dopants=True,
    verbose=False,
    log_materials=None,
    bulk_only_materials=None,
    bulk_only_neighbours=None,
):
    old_lookup = build_old_lookup(old_coords, old_regions, old_triangles)
    sol_ids = solution_ids(sol_line)
    stats = defaultdict(int)
    used_materials = defaultdict(set)
    converted_tris = []
    edge_count = defaultdict(int)

    for _eid, tags, verts in out_tris:
        region = surface_region.get(tags[1]) if len(tags) >= 2 else None
        if region is None:
            material = tags[0] if tags else -1
            matches = [r for r, m in old_regions.items() if m == material]
            region = matches[0] if matches else None
        if region is None:
            raise RuntimeError(f"could not map gmsh element tags {tags} to a region")
        p0, p1, p2 = (nodes[v] for v in verts)
        area2 = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])
        if area2 > 0.0:
            verts = (verts[0], verts[2], verts[1])
        material = old_regions[region]
        converted_tris.append((region, verts))
        for a, b in ((verts[1], verts[2]), (verts[2], verts[0]), (verts[1], verts[0])):
            edge_count[edge_key(a, b)] += 1
        for v in verts:
            used_materials[v].add(material)

    point_ids = {gid: i + 1 for i, gid in enumerate(sorted(nodes))}
    region_ids = {old: i + 1 for i, old in enumerate(sorted(set(r for r, _ in converted_tris)))}
    xs = [p[0] for p in nodes.values()]
    ys = [p[1] for p in nodes.values()]
    xmin, xmax = min(xs), max(xs)
    ymax = max(ys)

    def boundary_code(a, b):
        if edge_count[edge_key(a, b)] != 1:
            return -1024
        ax, ay = nodes[a]
        bx, by = nodes[b]
        tol = 1.0e-8
        if abs(ay - ymax) < tol and abs(by - ymax) < tol:
            return -1023
        if (abs(ax - xmin) < tol and abs(bx - xmin) < tol) or (
            abs(ax - xmax) < tol and abs(bx - xmax) < tol
        ):
            return -1024
        return -1022

    for old_region, verts in converted_tris:
        for a, b in ((verts[1], verts[2]), (verts[2], verts[0]), (verts[1], verts[0])):
            bc = boundary_code(a, b)
            if edge_count[edge_key(a, b)] == 1 and bc != -1024:
                used_materials[a].add(0)
                used_materials[b].add(0)

    node_values = {}
    for gid in sorted(nodes):
        for material in sorted(used_materials[gid]):
            node_values[(gid, material)] = interpolate_solution(
                nodes[gid], material, old_lookup, old_solutions, old_coords, sol_ids, stats,
                log_materials=log_materials,
            )

    if conserve_dopants:
        if bulk_only_materials:
            conserve_dopant_integrals_bulk_only(
                old_coords, old_regions, old_triangles, old_solutions,
                nodes, converted_tris, node_values, sol_ids, verbose,
                bulk_only_materials=set(bulk_only_materials),
                neighbour_materials=set(bulk_only_neighbours or ()),
            )
        else:
            conserve_dopant_integrals(
                old_coords, old_regions, old_triangles, old_solutions,
                nodes, converted_tris, node_values, sol_ids, verbose,
            )

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("v SUPREM-IV.GS B.9305\n")
        handle.write("D 2 3 3\n")
        for gid in sorted(nodes):
            pid = point_ids[gid]
            x, y = nodes[gid]
            handle.write(f"c {pid} {x:.9g} {y:.9g}  0\n")
        for old_region in sorted(region_ids):
            handle.write(f"r {region_ids[old_region]}   {old_regions[old_region]}\n")
        for tid, (old_region, verts) in enumerate(converted_tris, start=1):
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
        for gid in sorted(nodes):
            pid = point_ids[gid]
            for material in sorted(used_materials[gid]):
                vals = node_values[(gid, material)]
                # SUPREM .str files use 1-based point ids in c/t records, but
                # ig2_write writes zero-based point ids in n records.
                handle.write(f"n {pid - 1}   {material}  ")
                handle.write(" ".join(f"{v:.6e}" for v in vals))
                handle.write("\n")
        if implant_line:
            handle.write(f"{implant_line}\n")

    if verbose:
        print(f"solution interpolation: exact={stats['interpolated']} nearest={stats['nearest_fallback']}")
        for key in sorted(k for k in stats if k.startswith("nearest_fallback_material_")):
            print(f"  {key}={stats[key]}")


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--out-prefix", type=Path)
    parser.add_argument("--mesh-size", type=float, default=0.1)
    parser.add_argument("--top-mesh-size", type=float, help="boundary mesh size near the top of the structure")
    parser.add_argument("--bottom-mesh-size", type=float, help="boundary mesh size near the bottom of the structure")
    parser.add_argument("--grade-top-y", type=float, help="y coordinate where vertical mesh grading starts")
    parser.add_argument("--grade-bottom-y", type=float, help="y coordinate where vertical mesh grading reaches bottom size")
    parser.add_argument("--grade-power", type=float, default=1.0, help="vertical grading exponent; larger keeps the fine mesh deeper")
    parser.add_argument("--silicon-profile", action="store_true", help="auto-grade from the detected top of silicon to deep silicon")
    parser.add_argument("--silicon-material", type=int, default=3, help="SUPREM material id for silicon")
    parser.add_argument("--active-depth", type=float, help="distance below top silicon where the mesh reaches bottom size")
    parser.add_argument("--preserve-top-silicon", action="store_true", help="embed existing top-silicon mesh points so deck-defined near-surface spacing is preserved")
    parser.add_argument("--preserve-depth", type=float, help="depth below top silicon for --preserve-top-silicon; defaults to --active-depth or 2.0")
    parser.add_argument("--anisotropic", action="store_true", help="use Gmsh BAMG anisotropic metric instead of scalar size only")
    parser.add_argument("--x-mesh-size", type=float, help="target lateral size for --anisotropic")
    parser.add_argument("--y-top-mesh-size", type=float, help="target vertical size at top silicon for --anisotropic")
    parser.add_argument("--y-bottom-mesh-size", type=float, help="target vertical size at bottom/deep silicon for --anisotropic")
    parser.add_argument("--structure-top-spacing", action="store_true", help="derive top vertical anisotropic spacing from existing top-silicon mesh")
    parser.add_argument("--embed-top-grid", action="store_true", help="embed a deck-like anisotropic grid in the top silicon band")
    parser.add_argument("--angle", type=float, default=28.0)
    parser.add_argument("--gmsh", default="gmsh")
    parser.add_argument("--gmsh-timeout", type=float, default=60.0, help="seconds before aborting Gmsh")
    parser.add_argument("--write-str", type=Path, help="write a remeshed SUPREM .str")
    parser.add_argument("--no-conserve-dopants", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    out_prefix = args.out_prefix or args.input.with_suffix("")
    geo_path = out_prefix.with_suffix(".geo")
    msh_path = out_prefix.with_suffix(".msh")

    coords, regions, triangles, node_solutions, sol_line, implant_line, model_line, _header = parse_str(args.input)
    npts, nlines, nsurf, surface_region = write_geo(
        coords,
        regions,
        triangles,
        geo_path,
        args.mesh_size,
        args.top_mesh_size,
        args.bottom_mesh_size,
        args.grade_top_y,
        args.grade_bottom_y,
        args.grade_power,
        args.silicon_profile,
        args.silicon_material,
        args.active_depth,
        args.preserve_top_silicon,
        args.preserve_depth,
        args.anisotropic,
        args.x_mesh_size,
        args.y_top_mesh_size,
        args.y_bottom_mesh_size,
        args.structure_top_spacing,
        args.embed_top_grid,
    )
    cmd = [
        args.gmsh,
        "-2",
        "-format",
        "msh2",
        "-algo",
        "bamg" if args.anisotropic else "del2d",
        "-clscale",
        "1",
        "-o",
        str(msh_path),
        str(geo_path),
    ]
    try:
        proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=args.gmsh_timeout)
    except subprocess.TimeoutExpired as exc:
        print(exc.stdout or "")
        print(f"gmsh timed out after {args.gmsh_timeout:g} seconds")
        return 124
    if proc.returncode != 0:
        print(proc.stdout)
        return proc.returncode

    nodes, out_tris = parse_msh(msh_path)
    qualities = [triangle_quality(nodes, verts) for _eid, _tags, verts in out_tris]
    print(f"geo: {geo_path}")
    print(f"msh: {msh_path}")
    print(f"input boundary points={npts} curves={nlines} surfaces={nsurf}")
    print(f"output nodes={len(nodes)} triangles={len(out_tris)}")
    if qualities:
        print(f"quality min={min(qualities):.6g} avg={sum(qualities)/len(qualities):.6g}")
    if args.write_str:
        write_remeshed_str(
            args.write_str,
            nodes,
            out_tris,
            surface_region,
            regions,
            coords,
            triangles,
            node_solutions,
            sol_line,
            implant_line,
            model_line,
            not args.no_conserve_dopants,
            args.verbose,
        )
        print(f"str: {args.write_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
