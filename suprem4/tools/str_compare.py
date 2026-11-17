#!/usr/bin/env python3
"""Compare SUPREM .str files by geometry-oriented metrics."""

import argparse
import math
import sys
from collections import defaultdict


def parse_str(path):
    coords = {}
    regions = {}
    triangles = []
    nodes = {}

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            rec = parts[0]
            if rec == "c" and len(parts) >= 4:
                coords[int(parts[1])] = (float(parts[2]), float(parts[3]))
            elif rec == "r" and len(parts) >= 3:
                regions[int(parts[1])] = int(parts[2])
            elif rec == "t" and len(parts) >= 6:
                triangles.append(
                    {
                        "id": int(parts[1]),
                        "region": int(parts[2]),
                        "verts": (int(parts[3]), int(parts[4]), int(parts[5])),
                    }
                )
            elif rec == "n" and len(parts) >= 3:
                nodes[int(parts[1])] = int(parts[2])

    return {
        "coords": coords,
        "regions": regions,
        "triangles": triangles,
        "nodes": nodes,
    }


def tri_area(a, b, c):
    return abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2.0


def edge_len(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def bbox_add(box, p):
    if box is None:
        return [p[0], p[0], p[1], p[1]]
    box[0] = min(box[0], p[0])
    box[1] = max(box[1], p[0])
    box[2] = min(box[2], p[1])
    box[3] = max(box[3], p[1])
    return box


def summarize(data):
    coords = data["coords"]
    regions = data["regions"]
    by_material = defaultdict(lambda: {"area": 0.0, "bbox": None, "triangles": 0})
    min_edge = math.inf
    min_area = math.inf
    missing = 0

    for tri in data["triangles"]:
        verts = tri["verts"]
        if any(v not in coords for v in verts):
            missing += 1
            continue
        pts = [coords[v] for v in verts]
        mat = regions.get(tri["region"], -1)
        area = tri_area(pts[0], pts[1], pts[2])
        min_area = min(min_area, area)
        for i in range(3):
            min_edge = min(min_edge, edge_len(pts[i], pts[(i + 1) % 3]))
        by_material[mat]["area"] += area
        by_material[mat]["triangles"] += 1
        for p in pts:
            by_material[mat]["bbox"] = bbox_add(by_material[mat]["bbox"], p)

    bbox_all = None
    for p in coords.values():
        bbox_all = bbox_add(bbox_all, p)

    return {
        "points": len(coords),
        "nodes": len(data["nodes"]),
        "triangles": len(data["triangles"]),
        "regions": len(regions),
        "bbox": bbox_all,
        "min_edge": 0.0 if min_edge == math.inf else min_edge,
        "min_area": 0.0 if min_area == math.inf else min_area,
        "missing_tri_vertices": missing,
        "materials": dict(by_material),
    }


def rel_delta(a, b):
    scale = max(abs(a), abs(b), 1.0e-30)
    return abs(a - b) / scale


def bbox_delta(a, b):
    if a is None or b is None:
        return math.inf
    return max(abs(x - y) for x, y in zip(a, b))


def compare(base, cand, area_rel_tol, bbox_abs_tol):
    problems = []
    materials = set(base["materials"]) | set(cand["materials"])
    for mat in sorted(materials):
        bm = base["materials"].get(mat)
        cm = cand["materials"].get(mat)
        if bm is None or cm is None:
            problems.append(f"material {mat}: missing from one file")
            continue
        area_rd = rel_delta(bm["area"], cm["area"])
        if area_rd > area_rel_tol:
            problems.append(
                f"material {mat}: area differs by {area_rd:.3e} "
                f"({bm['area']:.8g} vs {cm['area']:.8g})"
            )
        bd = bbox_delta(bm["bbox"], cm["bbox"])
        if bd > bbox_abs_tol:
            problems.append(
                f"material {mat}: bbox differs by {bd:.3e} "
                f"({bm['bbox']} vs {cm['bbox']})"
            )
    if cand["missing_tri_vertices"]:
        problems.append(f"candidate has {cand['missing_tri_vertices']} triangles with missing vertices")
    return problems


def print_summary(label, summary):
    print(f"{label}:")
    print(
        f"  points={summary['points']} nodes={summary['nodes']} "
        f"triangles={summary['triangles']} regions={summary['regions']}"
    )
    print(
        f"  bbox={summary['bbox']} min_edge={summary['min_edge']:.8g} "
        f"min_area={summary['min_area']:.8g}"
    )
    for mat in sorted(summary["materials"]):
        item = summary["materials"][mat]
        print(
            f"  material {mat}: area={item['area']:.8g} "
            f"triangles={item['triangles']} bbox={item['bbox']}"
        )


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline")
    parser.add_argument("candidate")
    parser.add_argument("--area-rel-tol", type=float, default=2.0e-2)
    parser.add_argument("--bbox-abs-tol", type=float, default=5.0e-2)
    parser.add_argument("--summary", action="store_true")
    args = parser.parse_args(argv)

    base = summarize(parse_str(args.baseline))
    cand = summarize(parse_str(args.candidate))
    if args.summary:
        print_summary("baseline", base)
        print_summary("candidate", cand)

    problems = compare(base, cand, args.area_rel_tol, args.bbox_abs_tol)
    if problems:
        print("FAIL:")
        for problem in problems:
            print(f"  {problem}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
