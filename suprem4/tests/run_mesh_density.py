#!/usr/bin/env python3
"""Calibrate initial point-count ratios, then run the full PMOS deck.

Uniformly scale all initial x/y line spacings. Physical dimensions and
process parameters stay fixed; adaptive processing changes later counts.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import shutil
import tempfile
import time

from run_oxidation_perturbations import ROOT, CHECKPOINTS, audit, run_deck


def point_count(path):
    return sum(line.startswith("c ") for line in path.read_text().splitlines())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "src/suprem4")
    parser.add_argument("--ratios", type=float, nargs="+", default=[0.5, 2.0])
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    if any(not 0 < ratio <= 16 for ratio in args.ratios):
        parser.error("ratios must be positive and at most 16")
    if len(set(args.ratios)) != len(args.ratios):
        parser.error("ratios must be unique")
    source = ROOT.parent / "decks/rit_pmos"
    original = (source / "rit_pmos.in").read_text().replace(
        "line y loc=0.0  spac=0.0005 tag=surf",
        "line y loc=0.0  spac=0.0006 tag=surf", 1)
    work = Path(tempfile.mkdtemp(prefix="suprem4_mesh_density_"))
    print(f"workdir: {work}", flush=True)

    def prepare(name, scale, initial_only=False):
        directory = work / name
        directory.mkdir()
        for item in ("modelrc", "sup4gs.imp", "suprem.uk"):
            if (source / item).exists():
                shutil.copy2(source / item, directory / item)
        text, changes = re.subn(
            r"(^line [xy] .*?spac=)([0-9.eE+-]+)",
            lambda m: m[1] + format(float(m[2]) * scale, ".12g"),
            original, flags=re.MULTILINE)
        if changes != 5:
            raise ValueError(f"expected 5 initial spacings, found {changes}")
        marker = "initialize phosphorus conc=3e14"
        if initial_only:
            text = text[:text.index(marker) + len(marker)]
        text = text.replace(marker, marker + "\nstructure outfile=initial_mesh.str", 1)
        deck = directory / "rit_pmos.in"
        deck.write_text(text + "\n")
        return directory, deck

    def initialize(name, scale):
        directory, deck = prepare(name, scale, True)
        code, _ = run_deck(args.binary.resolve(), deck, directory, False)
        if code:
            raise RuntimeError(f"initialization failed: {directory}, exit {code}")
        audit(directory / "initial_mesh.str")
        count = point_count(directory / "initial_mesh.str")
        print(f"{name}: spacing scale={scale:.8g}, initial points={count}", flush=True)
        return count

    baseline = initialize("baseline_initial", 1.0)
    cases = []
    for ratio in args.ratios:
        name = {0.5: "half_points", 1.0: "baseline_points", 2.0: "double_points"}.get(
            ratio, f"points_{ratio:g}x")
        target = baseline * ratio
        low, high = (1.0, 16.0) if ratio < 1 else (0.01, 1.0)
        scale = ratio ** -0.5
        for attempt in range(10):
            count = initialize(f"calibrate_{name}_{attempt}", scale)
            if abs(count / target - 1) <= 0.02:
                break
            if count > target:
                low = scale
            else:
                high = scale
            scale = (low + high) / 2
        else:
            raise RuntimeError(f"could not reach {name} within 2%")
        cases.append((name, scale, count))

    def run(case):
        name, scale, expected = case
        directory, deck = prepare(name, scale)
        started = time.monotonic()
        code, log = run_deck(args.binary.resolve(), deck, directory, False)
        result = dict(case=name, spacing_scale=scale, baseline_points=baseline,
                      binary=str(args.binary.resolve()), elapsed_seconds=time.monotonic()-started,
                      initial_points=point_count(directory / "initial_mesh.str"),
                      returncode=code, log=str(log), checkpoints={}, errors=[])
        if result["initial_points"] != expected:
            result["errors"].append("initial point count differs from calibration")
        for name in ["initial_mesh.str", *CHECKPOINTS]:
            try:
                result["checkpoints"][name] = audit(directory / name)
                result["checkpoints"][name]["points"] = point_count(directory / name)
            except (OSError, ValueError, KeyError) as error:
                result["errors"].append(f"{name}: {error}")
        result["passed"] = code == 0 and not result["errors"]
        print(f"{'PASS' if result['passed'] else 'FAIL'} {result['case']}: "
              f"exit={code}, valid checkpoints={len(result['checkpoints'])}/10", flush=True)
        return result

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        results = list(pool.map(run, cases))
    (work / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0 if all(r["passed"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
