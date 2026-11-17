#!/usr/bin/env python3
"""Overlay SUPREM and Silvaco boron-vs-y profiles at a vertical cutline.

Both codes write `.str` files with a similar layout:
  c <id> <x> <y> <z>
  r <region> <material>
  t <tid> <region> <v1> <v2> <v3> ...
  s <count> <sol_id...>           (Silvaco prepends 513 to the sol list)
  n <pid> <material> <values...>  (one value per sol_id, same order as s)

This script parses both, extracts boron concentration at points along a vertical
cutline at x=X (default 35 um, a source/drain centerline for rit_pmos), and
draws both profiles on one figure for direct comparison. Also prints the peak
depth, peak concentration, and junction depth (Boron crossing 3e14 phosphorus
background) for each.
"""

import argparse
import sys
from pathlib import Path

import numpy as np


SUPREM_BORON_SOL_ID = 5
SUPREM_PHOSPHORUS_SOL_ID = 3


def _parse_str(path: Path, boron_sol_id=5, format_kind="auto"):
    """Return (coords, triangles, tri_regions, regions, boron_per_node).

    Both SUPREM and Silvaco use the same record syntax:
      c <id> <x> <y> <z>          (1-based id)
      r <region_id> <material>    (1-based region id, material 3 = silicon)
      t <tid> <region> <v1> <v2> <v3> ...   (1-based)
      s <count> <id1> <id2> ...
      n <pid> <material> <values...>        (0-based pid)
    """
    coords = {}
    triangles = []
    tri_regions = []
    regions = {}
    sol_ids = []
    boron_index = None
    n_records = {}

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            tag = parts[0]
            if tag == "c" and len(parts) >= 4:
                coords[int(parts[1]) - 1] = (float(parts[2]), float(parts[3]))
            elif tag == "t" and len(parts) >= 6:
                triangles.append((int(parts[3]) - 1, int(parts[4]) - 1, int(parts[5]) - 1))
                tri_regions.append(int(parts[2]) - 1)
            elif tag == "r" and len(parts) >= 3:
                regions[int(parts[1]) - 1] = int(parts[2])
            elif tag == "s":
                # Both SUPREM and Silvaco: `s <count> <id1> <id2> ...`
                count = int(parts[1])
                sol_ids = [int(x) for x in parts[2:2 + count]]
                if boron_sol_id in sol_ids:
                    boron_index = sol_ids.index(boron_sol_id)
            elif tag == "n" and sol_ids:
                pid = int(parts[1])
                material = int(parts[2])
                values = parts[3:3 + len(sol_ids)]
                if boron_index is not None and boron_index < len(values):
                    try:
                        n_records[(pid, material)] = float(values[boron_index])
                    except ValueError:
                        pass

    if not coords or not triangles:
        raise ValueError(f"{path} has no coordinate or triangle records")

    n_pts = max(coords) + 1
    coord_array = np.full((n_pts, 2), np.nan, dtype=float)
    for pid, (x, y) in coords.items():
        coord_array[pid] = (x, y)
    tri_array = np.array(triangles, dtype=int)
    tri_region_array = np.array(tri_regions, dtype=int)
    return coord_array, tri_array, tri_region_array, regions, n_records


def integrate_boron_in_material(coords, triangles, tri_regions, regions, n_records,
                                target_material=3):
    """Sum (triangle_area * mean_boron_at_vertices) over triangles whose region
    is `target_material`. Returns the integral in [um^2 * cm^-3] units; multiply
    by 1e-8 to get atoms per cm of depth into the page.
    """
    total = 0.0
    triangle_count = 0
    for verts, region in zip(triangles, tri_regions):
        material = regions.get(region, -1)
        if material != target_material:
            continue
        a, b, c = coords[verts[0]], coords[verts[1]], coords[verts[2]]
        area_um2 = 0.5 * abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
        # Boron at each vertex IN THIS MATERIAL specifically
        vals = []
        for v in verts:
            v_in_mat = n_records.get((v, material))
            if v_in_mat is None:
                # Try any material (some codes record only once per vertex)
                v_in_mat = boron_at_node(n_records, v)
            if v_in_mat is None or not np.isfinite(v_in_mat):
                vals = None
                break
            vals.append(v_in_mat)
        if vals is None:
            continue
        mean_b = sum(vals) / 3.0
        total += area_um2 * mean_b
        triangle_count += 1
    # Convert um^2 * cm^-3 -> atoms / cm (out-of-plane unit depth)
    atoms_per_cm = total * 1.0e-8
    return atoms_per_cm, triangle_count


def boron_at_node(n_records, pid):
    """Return Boron at a point id by trying any material it appears in."""
    candidates = [v for (p, m), v in n_records.items() if p == pid]
    if not candidates:
        return np.nan
    # Pick the max — silicon-side B is the largest where the material has B.
    return max(candidates)


def sample_vertical_cutline(coords, triangles, n_records, x_cut, y0, y1, samples):
    """Return (y_samples, boron_at_each) by barycentric interpolation."""
    pt_values = np.array([boron_at_node(n_records, pid) for pid in range(len(coords))])
    ys = np.linspace(y0, y1, samples)
    result = np.full(samples, np.nan, dtype=float)

    tri_pts = coords[triangles]  # (M, 3, 2)
    tri_xmin = np.nanmin(tri_pts[:, :, 0], axis=1)
    tri_xmax = np.nanmax(tri_pts[:, :, 0], axis=1)
    tri_ymin = np.nanmin(tri_pts[:, :, 1], axis=1)
    tri_ymax = np.nanmax(tri_pts[:, :, 1], axis=1)

    for i, y in enumerate(ys):
        candidates = np.flatnonzero(
            (x_cut >= tri_xmin - 1e-9) & (x_cut <= tri_xmax + 1e-9)
            & (y >= tri_ymin - 1e-9) & (y <= tri_ymax + 1e-9)
        )
        for tri_i in candidates:
            a, b, c = tri_pts[tri_i]
            v0 = b - a
            v1 = c - a
            v2 = np.array([x_cut, y]) - a
            denom = v0[0] * v1[1] - v1[0] * v0[1]
            if abs(denom) < 1e-30:
                continue
            u = (v2[0] * v1[1] - v1[0] * v2[1]) / denom
            v = (v0[0] * v2[1] - v2[0] * v0[1]) / denom
            w = 1.0 - u - v
            if u >= -1e-10 and v >= -1e-10 and w >= -1e-10:
                vals = pt_values[triangles[tri_i]]
                if np.all(np.isfinite(vals)):
                    result[i] = w * vals[0] + u * vals[1] + v * vals[2]
                    break
    return ys, result


def summarize(name, y, b, n_background):
    finite = np.isfinite(b) & (b > 0)
    if not np.any(finite):
        print(f"{name}: no finite boron samples")
        return
    y_f = y[finite]
    b_f = b[finite]
    peak_i = np.argmax(b_f)
    peak_y = y_f[peak_i]
    peak_b = b_f[peak_i]
    # Junction: first depth past the peak where B drops below n_background
    junction_y = None
    for yi, bi in zip(y_f[peak_i:], b_f[peak_i:]):
        if bi < n_background:
            junction_y = yi
            break
    junction_str = f"{junction_y:.4f} um" if junction_y is not None else "not crossed in range"
    print(f"{name}:")
    print(f"  peak: y = {peak_y:.4f} um, B = {peak_b:.3e} cm^-3")
    print(f"  junction (B < {n_background:.0e}): y = {junction_str}")


def main(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--suprem", type=Path, required=True, help="SUPREM .str file")
    p.add_argument("--silvaco", type=Path, required=True, help="Silvaco .str file")
    p.add_argument("--x", type=float, default=35.0, help="x position (um) of vertical cutline")
    p.add_argument("--y0", type=float, default=-0.5, help="cutline top y (um)")
    p.add_argument("--y1", type=float, default=3.0, help="cutline bottom y (um)")
    p.add_argument("--samples", type=int, default=500, help="sample count along cutline")
    p.add_argument("--background", type=float, default=3.0e14,
                   help="n-type background (phosphorus) for junction depth, default 3e14")
    p.add_argument("--out", type=Path, default=Path("compare_boron.png"),
                   help="output PNG")
    args = p.parse_args(argv)

    print(f"reading SUPREM:  {args.suprem}")
    sup_coords, sup_tris, sup_tri_regions, sup_regions, sup_n = _parse_str(args.suprem)
    print(f"reading Silvaco: {args.silvaco}")
    sil_coords, sil_tris, sil_tri_regions, sil_regions, sil_n = _parse_str(args.silvaco)

    # Dose-conservation sanity: integrate boron over silicon (material 3) in both
    sup_int, sup_n_tri = integrate_boron_in_material(
        sup_coords, sup_tris, sup_tri_regions, sup_regions, sup_n, target_material=3)
    sil_int, sil_n_tri = integrate_boron_in_material(
        sil_coords, sil_tris, sil_tri_regions, sil_regions, sil_n, target_material=3)
    expected = 2.0 * 2.0e15 * 50.0e-4   # 2 windows * dose [/cm2] * width [cm]
    print(f"\nTotal boron in silicon (atoms per cm of out-of-plane depth):")
    print(f"  SUPREM:  {sup_int:.4e}  ({sup_n_tri} silicon triangles)")
    print(f"  Silvaco: {sil_int:.4e}  ({sil_n_tri} silicon triangles)")
    print(f"  Expected from deck (2 windows x 50um x 2e15/cm2): {expected:.4e}")
    if sil_int > 0:
        print(f"  SUPREM / Silvaco ratio: {sup_int / sil_int:.4f}")

    print(f"\ncutline: x={args.x} um, y={args.y0}..{args.y1} um, {args.samples} samples\n")

    ys, sup_b = sample_vertical_cutline(sup_coords, sup_tris, sup_n,
                                        args.x, args.y0, args.y1, args.samples)
    _,  sil_b = sample_vertical_cutline(sil_coords, sil_tris, sil_n,
                                        args.x, args.y0, args.y1, args.samples)

    summarize(f"SUPREM  ({args.suprem.name})", ys, sup_b, args.background)
    print()
    summarize(f"Silvaco ({args.silvaco.name})", ys, sil_b, args.background)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    finite_sup = np.isfinite(sup_b) & (sup_b > 0)
    finite_sil = np.isfinite(sil_b) & (sil_b > 0)
    ax.semilogy(ys[finite_sup], sup_b[finite_sup], label=f"SUPREM ({args.suprem.name})",
                color="tab:blue", linewidth=1.5)
    ax.semilogy(ys[finite_sil], sil_b[finite_sil], label=f"Silvaco ({args.silvaco.name})",
                color="tab:red", linewidth=1.5, linestyle="--")
    ax.axhline(args.background, color="gray", linewidth=0.7, linestyle=":",
               label=f"n-bkgd {args.background:.0e}")
    ax.set_xlabel("y [um] (depth into structure)")
    ax.set_ylabel("Boron concentration [cm^-3]")
    ax.set_title(f"Boron profile at x = {args.x:.2f} um")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out, dpi=180)
    print(f"\nsaved: {args.out}")


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
