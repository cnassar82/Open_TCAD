#!/usr/bin/env python3
"""Compare SUPREM .str solution fields by area-weighted material integrals."""

import argparse
import math
import sys
from collections import defaultdict


def parse_str(path):
    coords = {}
    regions = {}
    triangles = []
    sol_ids = []
    values = {}

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
                triangles.append((int(parts[2]), tuple(int(x) for x in parts[3:6])))
            elif parts[0] == "s" and len(parts) >= 2:
                sol_ids = [int(x) for x in parts[2:]]
            elif parts[0] == "n" and len(parts) >= 4:
                # n point_id material values...; point_id is zero-based in SUPREM output.
                point = int(parts[1]) + 1
                material = int(parts[2])
                values[(point, material)] = [float(x) for x in parts[3:]]

    return coords, regions, triangles, sol_ids, values


def tri_area(a, b, c):
    return abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2.0


def summarize(path):
    coords, regions, triangles, sol_ids, values = parse_str(path)
    out = defaultdict(lambda: {"area": 0.0, "integral": [0.0 for _ in sol_ids]})
    missing = 0

    for region, verts in triangles:
        material = regions.get(region, -1)
        if any(v not in coords for v in verts):
            missing += 1
            continue
        node_vals = [values.get((v, material)) for v in verts]
        if any(v is None for v in node_vals):
            missing += 1
            continue
        area = tri_area(coords[verts[0]], coords[verts[1]], coords[verts[2]])
        out[material]["area"] += area
        for i in range(len(sol_ids)):
            out[material]["integral"][i] += area * (
                node_vals[0][i] + node_vals[1][i] + node_vals[2][i]
            ) / 3.0

    return sol_ids, out, missing


def rel_delta(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1.0e-30)


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline")
    parser.add_argument("candidate")
    parser.add_argument("--rel-tol", type=float, default=5.0e-2)
    parser.add_argument("--summary", action="store_true")
    args = parser.parse_args(argv)

    base_ids, base, base_missing = summarize(args.baseline)
    cand_ids, cand, cand_missing = summarize(args.candidate)
    if base_ids != cand_ids:
        print(f"FAIL: solution ids differ: {base_ids} vs {cand_ids}")
        return 1

    failed = False
    for material in sorted(set(base) | set(cand)):
        b = base.get(material)
        c = cand.get(material)
        if b is None or c is None:
            print(f"FAIL: material {material} missing from one file")
            failed = True
            continue
        if args.summary:
            print(f"material {material}: area {b['area']:.8g} -> {c['area']:.8g}")
        for i, sol in enumerate(base_ids):
            bi = b["integral"][i]
            ci = c["integral"][i]
            rd = rel_delta(bi, ci)
            if args.summary:
                bavg = bi / b["area"] if b["area"] else math.nan
                cavg = ci / c["area"] if c["area"] else math.nan
                print(
                    f"  sol {sol}: integral {bi:.8e} -> {ci:.8e} "
                    f"rel={rd:.3e} avg {bavg:.8e} -> {cavg:.8e}"
                )
            if rd > args.rel_tol:
                print(
                    f"FAIL: material {material} sol {sol} integral differs by {rd:.3e} "
                    f"({bi:.8e} vs {ci:.8e})"
                )
                failed = True

    if base_missing or cand_missing:
        print(f"FAIL: missing values baseline={base_missing} candidate={cand_missing}")
        failed = True

    if failed:
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
