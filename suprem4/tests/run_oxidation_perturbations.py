#!/usr/bin/env python3
"""Run small PMOS mesh/geometry perturbations and audit saved topology.

This checks process completion and mesh integrity, not physics accuracy.
Outputs and logs are retained in the printed temporary directory.
"""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from run_pmos_mesh_regression import CHECKPOINTS, make_deck, run_deck


def audit(path):
    coordinates, regions, triangles = {}, set(), {}
    for line in path.read_text().splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == "c":
            coordinates[int(p[1])] = tuple(map(float, p[2:4]))
        elif p[0] == "r":
            regions.add(int(p[1]))
        elif p[0] == "n":
            if not all(math.isfinite(float(x)) for x in p[3:]):
                raise ValueError(f"nonfinite solution at node record {p[1]}")
        elif p[0] == "t":
            tid = int(p[1])
            if tid in triangles:
                raise ValueError(f"duplicate triangle ID {tid}")
            triangles[tid] = (int(p[2]), tuple(map(int, p[3:6])),
                              tuple(map(int, p[6:9])))
    if not triangles:
        raise ValueError("no triangles")
    seen, edges, areas = set(), defaultdict(list), defaultdict(float)
    for tid, (region, vertices, neighbors) in triangles.items():
        if region not in regions or len(set(vertices)) != 3 or len(neighbors) != 3:
            raise ValueError(f"invalid triangle {tid}")
        key = tuple(sorted(vertices))
        if key in seen:
            raise ValueError(f"duplicate triangle {tid}")
        seen.add(key)
        a, b, c = [coordinates[v] for v in vertices]
        if not all(math.isfinite(x) for point in (a, b, c) for x in point):
            raise ValueError(f"nonfinite coordinates at {tid}")
        area = ((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])) / 2
        if area <= 0:
            raise ValueError(f"nonpositive area at triangle {tid}: {area}")
        areas[region] += area
        for i in range(3):
            edge = tuple(sorted((vertices[(i+1) % 3], vertices[(i+2) % 3])))
            edges[edge].append((tid, neighbors[i]))
    for edge, owners in edges.items():
        if len(owners) == 2:
            (t, n), (u, m) = owners
            valid = n == u and m == t
        else:
            valid = len(owners) == 1 and owners[0][1] < 0
        if not valid:
            raise ValueError(f"invalid edge ownership/neighbors: {edge}: {owners}")
    # Distinct vertex IDs can describe overlapping exposed edges.  Topological
    # ownership and positive triangle areas do not catch that malformed surface.
    boundary_axis_overlap(coordinates, [edge for edge, owners in edges.items()
                                        if len(owners) == 1])
    return {"triangles": len(triangles), "area_by_region_um2": dict(areas)}


def boundary_axis_overlap(coordinates, boundary_edges):
    """Reject exactly collinear horizontal/vertical boundary edge overlap."""
    horizontal, vertical = defaultdict(list), defaultdict(list)
    for a, b in boundary_edges:
        x0, y0 = coordinates[a]
        x1, y1 = coordinates[b]
        if y0 == y1:
            horizontal[y0].append((min(x0, x1), max(x0, x1), a, b))
        if x0 == x1:
            vertical[x0].append((min(y0, y1), max(y0, y1), a, b))
    for axis, groups in (("horizontal", horizontal), ("vertical", vertical)):
        for ordinate, segments in groups.items():
            segments.sort()
            for i, first in enumerate(segments):
                for second in segments[i + 1:]:
                    if second[0] >= first[1]:
                        break
                    if min(first[1], second[1]) > second[0]:
                        raise ValueError(
                            f"overlapping {axis} boundary edges at {ordinate}: "
                            f"{first[2:]} and {second[2:]}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "src/suprem4")
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    source = ROOT.parent / "decks/rit_pmos"
    work = Path(tempfile.mkdtemp(prefix="suprem4_oxidation_perturbations_"))
    print(f"workdir: {work}", flush=True)
    cases = [("spacing_minus_1pct", "0.000594", 0),
             ("spacing_plus_1pct", "0.000606", 0),
             ("gate_shift_minus_0.1um", "0.0006", -0.1),
             ("gate_shift_plus_0.1um", "0.0006", 0.1)]

    def run(case):
        name, spacing, shift = case
        directory = work / name
        directory.mkdir()
        for item in ("modelrc", "sup4gs.imp", "suprem.uk"):
            if (source / item).exists():
                shutil.copy2(source / item, directory / item)
        deck = directory / "rit_pmos.in"
        make_deck(source / "rit_pmos.in", deck, spacing)
        if shift:
            text = deck.read_text()
            start = text.index("#Pattern Gate Oxide Region")
            end = text.index("#700A Oxide Growth", start)
            block = text[start:end]
            for value in (50.0, 90.0):
                old = f"x={value:.1f}"
                if block.count(old) != 2:
                    raise ValueError(f"expected two gate coordinates {old}")
                block = block.replace(old, f"x={value+shift:.1f}")
            deck.write_text(text[:start] + block + text[end:])
        code, log = run_deck(args.binary.resolve(), deck, directory, False)
        result = {"case": name, "returncode": code, "log": str(log),
                  "checkpoints": {}, "errors": []}
        for checkpoint in CHECKPOINTS:
            try:
                result["checkpoints"][checkpoint] = audit(directory / checkpoint)
            except (OSError, ValueError, KeyError) as error:
                result["errors"].append(f"{checkpoint}: {error}")
        result["passed"] = code == 0 and not result["errors"]
        print(f"{'PASS' if result['passed'] else 'FAIL'} {name}: "
              f"exit={code}, valid checkpoints={len(result['checkpoints'])}/9", flush=True)
        return result

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        results = list(pool.map(run, cases))
    (work / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0 if all(r["passed"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
