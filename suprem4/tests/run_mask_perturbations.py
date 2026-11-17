#!/usr/bin/env python3
"""Run complete PMOS processes with translated mask edges; audit every checkpoint.

This checks completion and mesh topology, not oxidation model accuracy.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import shutil
import tempfile

from run_oxidation_perturbations import audit
from run_pmos_mesh_regression import CHECKPOINTS, make_deck, run_deck

ROOT = Path(__file__).resolve().parents[1]
GROUPS = {
    "diffusion": ("#Pattern Level #1 Diffusion", "#Implant Boron", 8),
    "gate": ("#Pattern Gate Oxide Region", "#700A Oxide Growth", 4),
    "contact": ("#Contact Etch", "#Deposit Metal", 8),
    "metal": ("#Etch Metal", None, 8),
}


def shift_masks(text, groups, shift):
    changes = []
    for group in groups:
        start_marker, end_marker, expected = GROUPS[group]
        start = text.index(start_marker)
        end = text.index(end_marker, start) if end_marker else len(text)
        block = text[start:end]
        edits = []

        def translate(match):
            old = float(match.group(1))
            new = old + shift
            edits.append({"from_um": old, "to_um": new})
            return f"x={new:.9g}"

        block, count = re.subn(r"\bx=([-+\d.eE]+)", translate, block)
        if count != expected:
            raise ValueError(f"{group}: expected {expected} mask coordinates, got {count}")
        text = text[:start] + block + text[end:]
        changes.append({"group": group, "coordinates": edits})
    return text, changes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "src/suprem4")
    parser.add_argument("--spacing", default="0.0006")
    parser.add_argument("--shift", type=float, default=0.1, help="shift magnitude in micrometers")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--groups", nargs="+", choices=["all", *GROUPS],
                        default=["all", *GROUPS], help="mask groups to translate")
    args = parser.parse_args()
    if not 0 < args.shift < 1:
        parser.error("--shift must be between 0 and 1 micrometer")
    source = ROOT.parent / "decks/rit_pmos"
    work = Path(tempfile.mkdtemp(prefix="suprem4_mask_perturbations_"))
    print(f"workdir: {work}", flush=True)
    cases = [("baseline", [], 0.0)]
    for group in args.groups:
        for sign, direction in [(-1, "left"), (1, "right")]:
            cases.append((f"{group}_{direction}", list(GROUPS) if group == "all" else [group],
                          sign * args.shift))

    def run(case):
        name, groups, shift = case
        directory = work / name
        directory.mkdir()
        for item in ("modelrc", "sup4gs.imp", "suprem.uk"):
            if (source / item).exists():
                shutil.copy2(source / item, directory / item)
        deck = directory / "rit_pmos.in"
        make_deck(source / "rit_pmos.in", deck, args.spacing)
        text, edits = shift_masks(deck.read_text(), groups, shift)
        deck.write_text(text)
        code, log = run_deck(args.binary.resolve(), deck, directory, False)
        result = {"case": name, "shift_um": shift, "mask_edits": edits,
                  "returncode": code, "log": str(log), "checkpoints": {}, "errors": []}
        for checkpoint in CHECKPOINTS:
            try:
                result["checkpoints"][checkpoint] = audit(directory / checkpoint)
            except (OSError, ValueError, KeyError) as error:
                result["errors"].append(f"{checkpoint}: {error}")
        result["passed"] = code == 0 and not result["errors"]
        (directory / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        print(f"{'PASS' if result['passed'] else 'FAIL'} {name}: "
              f"exit={code}, valid checkpoints={len(result['checkpoints'])}/9", flush=True)
        return result

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        results = list(pool.map(run, cases))
    (work / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0 if all(r["passed"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
