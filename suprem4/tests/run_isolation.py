#!/usr/bin/env python3
"""Exercise LOCOS/STI geometry, including denser meshes and shifted masks.

Checks completion, saved mesh topology, and basic geometric invariants.
These are reliability exercises, not calibrated process-accuracy benchmarks.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import sys
import tempfile
import time

from run_oxidation_perturbations import audit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from run_pmos_mesh_regression import run_deck
from str_compare import parse_str, summarize

CHECKPOINTS = {
    "locos": ["initial.str", "masked.str", "field_oxide.str", "locos.str"],
    "sti": ["initial.str", "trench.str", "liner.str", "filled.str", "sti.str"],
}
CHECKPOINTS["sti_rounded"] = CHECKPOINTS["sti"]


def source_deck(name):
    folder = "sti" if name == "sti_rounded" else name
    return ROOT.parent / "decks" / folder / (name + ".in")


def trench_polygon(text):
    block = text.split("etch silicon start", 1)[1].split("structure outfile=trench.str", 1)[0]
    return [tuple(map(float, pair)) for pair in
            re.findall(r"x=([-+\d.eE]+)\s+y=([-+\d.eE]+)", block)]


def trench_void_area(text):
    """Independent polygon area above the substrate's y=0 boundary."""
    polygon = trench_polygon(text)
    clipped = []
    for a, b in zip(polygon, polygon[1:] + polygon[:1]):
        if a[1] >= 0:
            clipped.append(a)
        if (a[1] < 0 <= b[1]) or (b[1] < 0 <= a[1]):
            clipped.append((a[0]+(b[0]-a[0])*(-a[1])/(b[1]-a[1]), 0))
    return abs(sum(a[0]*b[1]-b[0]*a[1]
                   for a, b in zip(clipped, clipped[1:] + clipped[:1])))/2


def silicon_surface(data, x):
    """Intersect silicon triangles with a vertical probe (micron units)."""
    ys = []
    for tri in data["triangles"]:
        if data["regions"][tri["region"]] != 3:
            continue
        pts = [data["coords"][v] for v in tri["verts"]]
        for a, b in zip(pts, pts[1:] + pts[:1]):
            if a[0] == b[0]:
                if x == a[0]:
                    ys.extend([a[1], b[1]])
            elif min(a[0], b[0]) <= x <= max(a[0], b[0]):
                ys.append(a[1] + (b[1]-a[1])*(x-a[0])/(b[0]-a[0]))
    if not ys:
        raise ValueError(f"silicon probe at x={x} misses structure")
    return min(ys)


def make_deck(name, scale, shift):
    text = source_deck(name).read_text()
    text = re.sub(r"\bspac=([\d.]+)",
                  lambda m: f"spac={float(m[1])*scale:.10g}", text)
    if name == "locos":
        text = text.replace("etch nitride left p1.x=0",
                            f"etch nitride left p1.x={shift:.10g}")
    else:
        start = text.index("etch silicon start")
        end = text.index("structure outfile=trench.str", start)
        block = re.sub(r"\bx=([-+\d.eE]+)",
                       lambda m: f"x={float(m[1])+shift:.10g}", text[start:end])
        text = text[:start] + block + text[end:]
    return text


def geometry(name, directory, shift):
    final = parse_str(directory / ("locos.str" if name == "locos" else "sti.str"))
    summary = summarize(final)
    mats = summary["materials"]
    if set(mats) != {1, 3}:
        raise ValueError(f"unexpected final materials: {sorted(mats)}")
    if name == "locos":
        field = silicon_surface(final, shift - 0.7)
        masked = silicon_surface(final, shift + 0.7)
        edge = silicon_surface(final, shift + 0.05)
        if not (field > 0.1 and masked < 0.01 and edge > masked + 0.001):
            raise ValueError(f"missing field/mask/bird's-beak contrast: {field, masked, edge}")
        summary["silicon_surface_y_um"] = {
            "field": field, "under_mask": masked, "near_mask_edge": edge}
    else:
        trench = summarize(parse_str(directory / "trench.str"))
        expected = 2 - trench_void_area((directory / (name + ".in")).read_text())
        if abs(trench["materials"][3]["area"] - expected) > 2e-5:
            raise ValueError(f"trench silicon area differs from polygon area {expected}")
        summary["expected_trench_silicon_area_um2"] = expected
        total = sum(m["area"] for m in mats.values())
        if abs(total - 2) > 2e-5:
            raise ValueError(f"filled, planarized cross-section area is {total}, expected 2")
        if abs(mats[1]["bbox"][2]) > 1e-6:
            raise ValueError("oxide does not terminate at the y=0 planar cut")
        if not 0.35 < silicon_surface(final, shift) < 0.4:
            raise ValueError("missing oxidized trench floor")
        summary["total_area_um2"] = total
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path,
                        default=ROOT / "build/validated-native64/src/suprem4")
    parser.add_argument("--cases", nargs="+", choices=list(CHECKPOINTS),
                        default=["locos", "sti", "sti_rounded"])
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    binary = args.binary.resolve()
    work = Path(tempfile.mkdtemp(prefix="suprem4_isolation_"))
    print(f"workdir: {work}", flush=True)
    variants = [("baseline", 1, 0)]
    if not args.baseline_only:
        variants += [("fine", 2**-0.5, 0), ("left", 1, -0.025), ("right", 1, 0.025)]

    def run(case):
        name, (variant, scale, shift) = case
        directory = work / f"{name}_{variant}"
        directory.mkdir()
        deck = directory / (name + ".in")
        deck.write_text(make_deck(name, scale, shift))
        started = time.monotonic()
        code, log = run_deck(binary, deck, directory, False)
        result = {"case": name, "variant": variant, "binary": str(binary),
                  "spacing_scale": scale, "mask_shift_um": shift,
                  "returncode": code, "elapsed_seconds": time.monotonic()-started,
                  "log": str(log), "checkpoints": {}, "errors": []}
        if re.search(r"errors detected on command input|suprem4 panic:|No mesh defined", log.read_text()):
            result["errors"].append("simulator reported a command or mesh error")
        for checkpoint in CHECKPOINTS[name]:
            try:
                result["checkpoints"][checkpoint] = audit(directory / checkpoint)
            except (OSError, ValueError, KeyError) as error:
                result["errors"].append(f"{checkpoint}: {error}")
        try:
            result["geometry"] = geometry(name, directory, shift)
        except (OSError, ValueError, KeyError) as error:
            result["errors"].append(str(error))
        result["passed"] = code == 0 and not result["errors"]
        (directory / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        print(f"{'PASS' if result['passed'] else 'FAIL'} {directory.name}: "
              f"{result['errors']}", flush=True)
        return result

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        results = list(pool.map(run, [(name, variant) for name in args.cases
                                     for variant in variants]))
    (work / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0 if all(r["passed"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
