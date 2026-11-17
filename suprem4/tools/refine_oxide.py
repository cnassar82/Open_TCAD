#!/usr/bin/env python3
"""Densify oxide-region meshes in a SUPREM .str without invoking gmsh.

Motivating problem: SUPREM's oxide-growth code adds new vertical layers only
at the silicon/oxide interface (controlled by `method grid.oxide=X`). After
several oxidation steps you can end up with a 0.66 um oxide that has only ONE
vertical node line per x position. When a Pearson implant follows, the peak
(typical ΔRp ~ 60 nm, in the middle of the oxide) lands in mesh cells that
have no nodes to deposit on, and dose is effectively lost (per-x cutline dose
drops to ~2.5 % of expected).

The existing tools/layered_remesh.py uses gmsh constrained Delaunay which
hangs on thin high-aspect-ratio oxide regions when asked for >1 layer. This
tool sidesteps gmsh: scipy.spatial.Delaunay on each oxide region's
boundary+interior point cloud, then writes the new structure using the
existing layered_remesh.write_combined_str.

The new interior nodes inherit per-material solution values from the nearest
pre-existing oxide node at the same x.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay


THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

import layered_remesh as lr  # noqa: E402
from gmsh_remesh_str import edge_key, parse_str  # noqa: E402


SiO2 = 1
SiNi = 2


def collect_region_verts(triangles, region_id):
    vert_ids = set()
    for _tid, region, verts in triangles:
        if region == region_id:
            vert_ids.update(int(v) for v in verts)
    return vert_ids


def column_y_values(coords, vert_ids):
    """Return {round(x): sorted [y, ...]} for boundary verts of the region."""
    ys_per_x = defaultdict(list)
    for v in vert_ids:
        x, y = coords[v]
        ys_per_x[round(x, 6)].append(y)
    return {x: sorted(ys) for x, ys in ys_per_x.items()}


def column_pids(coords, vert_ids):
    """{round(x): {round(y): pid, ...}} for boundary verts."""
    out = defaultdict(dict)
    for v in vert_ids:
        x, y = coords[v]
        out[round(x, 6)][round(y, 6)] = v
    return out


def triangulate_columnwise(coords, vert_ids, n_layers, boundary_buffer, next_pid):
    """Generate a structured triangulation for a rectangular oxide region by
    placing N interior y-layers at each existing boundary x-position, then
    connecting adjacent columns with consistent diagonals.

    Returns: (new_points, new_triangles, next_pid_used)
      new_points: list of (pid, x, y)
      new_triangles: list of (v1, v2, v3) global-pid triples
    """
    cols_y = column_y_values(coords, vert_ids)
    cols_pid = column_pids(coords, vert_ids)
    sorted_xs = sorted(cols_y)
    if len(sorted_xs) < 2:
        return [], [], next_pid

    # Build per-column ordered y-pid list. Boundary verts are kept as-is at
    # column top/bottom; interior layers are spaced uniformly between them.
    # IMPORTANT: We re-use the existing boundary vert pids at the top (oxide
    # surface) and bottom (oxide-Si interface) so the oxide triangles share
    # those vertices with the silicon-side triangles already in the mesh.
    cols = {}  # x -> list of (y, pid) sorted by y
    new_points = []
    for x in sorted_xs:
        col_ys = cols_y[x]
        ymin, ymax = col_ys[0], col_ys[-1]
        thickness = ymax - ymin
        if thickness < 3 * boundary_buffer:
            # too thin; keep only existing boundary verts in this column
            entries = [(y, cols_pid[x][round(y, 6)]) for y in col_ys]
        else:
            entries = [(ymin, cols_pid[x][round(ymin, 6)])]
            y0 = ymin + boundary_buffer
            y1 = ymax - boundary_buffer
            for y in np.linspace(y0, y1, n_layers + 2)[1:-1]:
                new_points.append((next_pid, x, float(y)))
                entries.append((float(y), next_pid))
                next_pid += 1
            entries.append((ymax, cols_pid[x][round(ymax, 6)]))
        cols[x] = entries

    # Connect adjacent columns. For each (x_i, x_{i+1}) pair, walk down both
    # columns using a merge-like algorithm: at each step we have two "rows"
    # (one from each column) at known y. We emit a triangle connecting the
    # next-lower pair. Diagonal choice: top-left to bottom-right (consistent).
    new_tris = []
    for i in range(len(sorted_xs) - 1):
        x_left = sorted_xs[i]
        x_right = sorted_xs[i + 1]
        left = cols[x_left]
        right = cols[x_right]
        # We need columns of the same length to do a clean structured mesh.
        # If lengths differ (e.g. step in oxide bottom at a window edge),
        # fall back to a Delaunay on this 2-column strip.
        if len(left) != len(right):
            strip_pts = []
            strip_pids = []
            for y, pid in left:
                strip_pts.append((x_left, y)); strip_pids.append(pid)
            for y, pid in right:
                strip_pts.append((x_right, y)); strip_pids.append(pid)
            arr = np.array(strip_pts, dtype=float)
            try:
                tri = Delaunay(arr)
                for s in tri.simplices:
                    new_tris.append(tuple(strip_pids[k] for k in s))
            except Exception:
                pass
            continue
        for j in range(len(left) - 1):
            _, p_lb = left[j]
            _, p_lt = left[j + 1]
            _, p_rb = right[j]
            _, p_rt = right[j + 1]
            # Two triangles per cell: (lb, rb, lt) and (rb, rt, lt)
            new_tris.append((p_lb, p_rb, p_lt))
            new_tris.append((p_rb, p_rt, p_lt))
    return new_points, new_tris, next_pid


def main(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--target-material", type=int, action="append", default=None,
                   help="material id to refine (may repeat). Default: SiO2 (1) + SiNi (2)")
    p.add_argument("--layers", type=int, default=8,
                   help="interior y layers to insert per x position (default 8)")
    p.add_argument("--boundary-buffer", type=float, default=0.01,
                   help="um buffer from oxide top/bottom for the first/last interior layer")
    args = p.parse_args(argv)

    target_materials = args.target_material or [SiO2, SiNi]

    # parse_str returns coords, regions, triangles, solutions, sol_line, implant_line, model_line, header_lines
    coords, regions, triangles, solutions, sol_line, implant_line, model_line, header = parse_str(args.input)
    # triangles: list of (tid, region, [v1, v2, v3])
    # solutions: dict {(pid, material): [val per sol_id]}
    print(f"input: {args.input}")
    print(f"  pts={len(coords)} tris={len(triangles)} regions={dict(regions)}")
    print(f"  refining materials={target_materials} with {args.layers} interior y layers per x")

    coords_dict = {int(pid): (float(x), float(y)) for pid, (x, y) in coords.items()}
    target_region_ids = [r for r, m in regions.items() if m in target_materials]
    print(f"  matching regions: {target_region_ids}")
    if not target_region_ids:
        print("  no matching regions; nothing to do")
        args.output.write_bytes(args.input.read_bytes())
        return 0

    # Build new structure
    new_coords = dict(coords_dict)
    new_tris = []  # list of (region_id, (v1, v2, v3))
    for tid, region, verts in triangles:
        if region in target_region_ids:
            continue
        new_tris.append((region, tuple(int(v) for v in verts)))

    next_pid = max(coords_dict) + 1
    new_solutions = dict(solutions)

    for r in target_region_ids:
        vert_ids = collect_region_verts(triangles, r)
        if not vert_ids:
            continue
        new_pts, local_tris, next_pid = triangulate_columnwise(
            coords_dict, vert_ids, args.layers, args.boundary_buffer, next_pid)
        if not local_tris:
            print(f"  region {r} (mat {regions[r]}): no triangulation; keeping original")
            for tid, region, verts in triangles:
                if region == r:
                    new_tris.append((region, tuple(int(v) for v in verts)))
            continue
        for pid, x, y in new_pts:
            new_coords[pid] = (x, y)
        old_count = sum(1 for tid_, region, _ in triangles if region == r)
        print(f"  region {r} (mat {regions[r]}): {len(vert_ids)} verts -> "
              f"{len(vert_ids) + len(new_pts)} verts; {old_count} -> {len(local_tris)} tris")
        for v1, v2, v3 in local_tris:
            new_tris.append((r, (v1, v2, v3)))

    # Copy solutions to new interior nodes from the same-x oxide neighbor.
    # Build {x: [(pid, material, vals)]} for old oxide nodes.
    old_oxide_at_x = defaultdict(list)
    for (pid, material), vals in solutions.items():
        if material in target_materials and pid in coords_dict:
            x, _ = coords_dict[pid]
            old_oxide_at_x[round(x, 4)].append((pid, material, vals))
    inserted = 0
    for new_pid, (nx, ny) in list(new_coords.items()):
        if new_pid in coords_dict:
            continue  # not a new node
        nearby = old_oxide_at_x.get(round(nx, 4), [])
        seen = set()
        for opid, material, vals in nearby:
            if material in seen:
                continue
            new_solutions[(new_pid, material)] = list(vals)
            seen.add(material)
        inserted += 1
    print(f"  inserted {inserted} new interior nodes; "
          f"final pts={len(new_coords)} tris={len(new_tris)}")

    # Drop coords/solutions for any pid that no triangle references — these
    # are vertices of removed oxide triangles whose role is now taken by the
    # structured-column boundary verts. SUPREM's self-test panics on orphan
    # points ("point X has no associated node").
    referenced = set()
    for _r, verts in new_tris:
        referenced.update(verts)
    orphans = set(new_coords) - referenced
    for pid in orphans:
        new_coords.pop(pid, None)
    if orphans:
        print(f"  dropped {len(orphans)} orphan vertices not used by any triangle")

    # Write via the existing layered_remesh writer (handles SUPREM-specific edge codes)
    lr.write_combined_str(
        str(args.output),
        new_coords,
        new_tris,
        regions,
        new_solutions,
        sol_line,
        implant_line,
        model_line,
    )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
