#!/usr/bin/env python3
"""Physics metrics of a poly-line re-oxidation result (WP7, acceptance criterion 2).

Reads before.str / after.str of a poly_line_*.in deck and measures, with straight cuts through the
triangles:
  - field: oxide thickness and silicon consumed far from the line (x = --field-x)
  - poly top: oxide thickness and poly consumed at the line centre
  - poly width and sidewall oxide on both sides at the poly mid-height
  - gate-edge dielectric: oxide between silicon and poly at 0.02 / 0.05 / 0.10 um inside the
    current poly bottom edge (both edges)
  - consumption ratios (consumed / grown) on silicon and poly; expected 0.44
  - dopant inventories (P1 integral, from `core_str`) before and after

Coordinates are um, y grows into the substrate.  Material codes: 1 oxide, 3 silicon, 4 poly.
"""

import argparse
import json
import math
import subprocess
import sys

OXIDE, SILICON, POLY = 1, 3, 4
NAMES = {0: "gas", 1: "oxide", 2: "nitride", 3: "silicon", 4: "poly"}


def parse(path):
    pts, reg, tri = {}, {}, []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == "c":
                pts[int(p[1])] = (float(p[2]), float(p[3]))
            elif p[0] == "r":
                reg[int(p[1])] = int(p[2])
            elif p[0] == "t":
                tri.append((int(p[2]), int(p[3]), int(p[4]), int(p[5])))
    return [(reg[r], pts[a], pts[b], pts[c]) for r, a, b, c in tri]


def cut(tris, axis, value):
    """Intervals [(lo, hi, material)] along a line axis=value (axis 0: x=value, cut along y)."""
    o = 1 - axis
    out = []
    for mat, *v in tris:
        s = []
        for i in range(3):
            a, b = v[i], v[(i + 1) % 3]
            da, db = a[axis] - value, b[axis] - value
            if da == 0 and db == 0:
                continue
            if (da <= 0 <= db) or (db <= 0 <= da):
                t = da / (da - db)
                s.append(a[o] + t * (b[o] - a[o]))
        if len(s) >= 2 and max(s) - min(s) > 1e-12:
            out.append((min(s), max(s), mat))
    out.sort()
    merged = []
    for lo, hi, m in out:
        if merged and merged[-1][2] == m and lo <= merged[-1][1] + 1e-9:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi, m])
    return merged


def runs(intervals):
    return [(lo, hi, m) for lo, hi, m in intervals]


def first(intervals, mat, after=-math.inf):
    for lo, hi, m in intervals:
        if m == mat and hi > after:
            return lo, hi
    return None


def column(tris, x):
    """Material stack along y at x, top (most negative y) first."""
    return runs(cut(tris, 0, x))


def thickness_above(col, mat_below):
    """Oxide thickness directly above the first interval of mat_below."""
    for i, (lo, hi, m) in enumerate(col):
        if m == mat_below:
            t, j = 0.0, i - 1
            while j >= 0 and col[j][2] == OXIDE and abs(col[j][1] - (col[j + 1][0] if j + 1 < len(col) else lo)) < 1e-6:
                t += col[j][1] - col[j][0]
                j -= 1
            return t, lo
    return None, None


def inventories(core_str, path):
    if not core_str:
        return {}
    r = subprocess.run([core_str, path], capture_output=True, text=True)
    inv = {}
    for line in r.stdout.splitlines():
        if line.startswith("inventory imp"):
            p = line.split()
            imp = int(p[2].rstrip(":"))
            inv[imp] = {k: float(v) for k, v in (x.split("=") for x in p[3:])}
    return inv


def poly_metrics(tb, ta, pb):
    res = {}
    pxs = [p[0] for t in pb for p in t]
    pys = [p[1] for t in pb for p in t]
    x0, x1, ytop0, ybot0 = min(pxs), max(pxs), min(pys), max(pys)
    xc = 0.5 * (x0 + x1)
    res["poly_width_before_um"] = x1 - x0
    # poly top
    col = column(ta, xc)
    poly = first(col, POLY)
    oxtop = 0.0
    for lo, hi, m in col:
        if m == OXIDE and hi <= poly[0] + 1e-9:
            oxtop += hi - lo
    res["poly_top_oxide_A"] = oxtop * 1e4
    res["poly_top_consumed_A"] = (poly[0] - ytop0) * 1e4
    res["poly_height_um"] = poly[1] - poly[0]
    # mid-height cut
    ymid = 0.5 * (poly[0] + poly[1])
    row = cut(ta, 1, ymid)
    pi = [i for i, (lo, hi, m) in enumerate(row) if m == POLY]
    if pi:
        plo, phi = row[pi[0]][0], row[pi[-1]][1]
        res["poly_width_mid_um"] = phi - plo
        left = sum(hi - lo for lo, hi, m in row if m == OXIDE and hi <= plo + 1e-9 and lo >= plo - 0.2)
        right = sum(hi - lo for lo, hi, m in row if m == OXIDE and lo >= phi - 1e-9 and hi <= phi + 0.2)
        res["sidewall_oxide_left_A"] = left * 1e4
        res["sidewall_oxide_right_A"] = right * 1e4
        res["sidewall_consumed_left_A"] = (plo - x0) * 1e4
        res["sidewall_consumed_right_A"] = (x1 - phi) * 1e4
    # gate-edge dielectric: oxide between silicon and the poly above it, d um inside the
    # original (as-etched) poly edges, as in doc/poly_oxidation/README.md section 2
    for d in (0.02, 0.05, 0.10):
        for side, x in (("left", x0 + d), ("right", x1 - d)):
            c = column(ta, x)
            p = first(c, POLY)
            s = first(c, SILICON)
            if p and s:
                res[f"gate_edge_{side}_{d:.2f}_A"] = (s[0] - p[1]) * 1e4
    c = column(ta, xc)
    res["gate_centre_A"] = (first(c, SILICON)[0] - first(c, POLY)[1]) * 1e4
    return res


def measure(before, after, field_x, core_str=None):
    tb, ta = parse(before), parse(after)
    res = {}
    # field
    cb, ca = column(tb, field_x), column(ta, field_x)
    ox0, si0 = thickness_above(cb, SILICON)
    ox1, si1 = thickness_above(ca, SILICON)
    res["field_oxide_A"] = ox1 * 1e4
    res["field_si_consumed_A"] = (si1 - si0) * 1e4
    res["field_grown_A"] = (ox1 - ox0) * 1e4
    # poly extent before
    pb = [(a, b, c) for m, a, b, c in tb if m == POLY]
    if pb:
        res.update(poly_metrics(tb, ta, pb))
    # consumption ratios
    if res["field_grown_A"] > 0:
        res["ratio_si"] = res["field_si_consumed_A"] / res["field_grown_A"]
    top_grown = res.get("poly_top_oxide_A", 0.0)
    if top_grown > 0:
        res["ratio_poly_top"] = res["poly_top_consumed_A"] / top_grown
    if "sidewall_oxide_left_A" in res:
        for s in ("left", "right"):
            ox = res[f"sidewall_oxide_{s}_A"]
            if ox > 0:
                res[f"ratio_poly_{s}"] = res[f"sidewall_consumed_{s}_A"] / ox
    ib, ia = inventories(core_str, before), inventories(core_str, after)
    for imp in sorted(set(ib) & set(ia)):
        tb_, ta_ = sum(ib[imp].values()), sum(ia[imp].values())
        if tb_ > 0:
            res[f"inventory_{imp}_rel_change"] = (ta_ - tb_) / tb_
            res[f"inventory_{imp}_after"] = ia[imp]
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("before")
    ap.add_argument("after")
    ap.add_argument("--field-x", type=float, default=0.2)
    ap.add_argument("--core-str", help="core_str binary for dopant inventories")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    r = measure(a.before, a.after, a.field_x, a.core_str)
    if a.json:
        json.dump(r, sys.stdout, indent=1)
        print()
    else:
        for k, v in r.items():
            print(f"{k:32s} {v:.4g}" if isinstance(v, float) else f"{k:32s} {v}")


if __name__ == "__main__":
    main()
