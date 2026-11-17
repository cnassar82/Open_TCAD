#!/usr/bin/env python3
"""Smoke test for tools/layered_remesh.py via gmsh_remesh_deck.py.

For each (spacing, --before-marker) pair, runs the deck with the layered
remesh inserted, then compares the final pmos.str against the baseline
checkpoints in the deck directory. Reports geometry and solution drift.
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


DEFAULT_CHECKPOINTS = [
    "700_ox.str",
    "cont_etch.str",
    "dep_alumin.str",
    "etch_alumin.str",
    "pmos.str",
]


def run_pipeline(repo, deck_dir, spacing, before_marker, work_dir):
    deck_src = deck_dir / "rit_pmos.in"
    deck_text = deck_src.read_text(encoding="utf-8")
    old = "line y loc=0.0  spac=0.0005 tag=surf"
    new = f"line y loc=0.0  spac={spacing} tag=surf"
    if old not in deck_text:
        raise RuntimeError(f"could not find expected surface spacing line in {deck_src}")
    deck_path = work_dir / "rit_pmos.in"
    deck_path.write_text(deck_text.replace(old, new, 1), encoding="utf-8")
    for item in ("modelrc", "sup4gs.imp", "suprem.uk"):
        src = deck_dir / item
        if src.exists():
            shutil.copy2(src, work_dir / item)
    cmd = [
        sys.executable,
        str(repo / "tools" / "gmsh_remesh_deck.py"),
        str(deck_path),
        "--before", before_marker,
        "--layered",
        "--work-dir", str(work_dir),
    ]
    proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (work_dir / "pipeline.log").write_text(proc.stdout, encoding="utf-8")
    return proc.returncode


def compare(repo, baseline, candidate, area_tol, bbox_tol, sol_tol):
    geo = subprocess.run(
        [sys.executable, str(repo / "tools" / "str_compare.py"),
         str(baseline), str(candidate),
         "--area-rel-tol", str(area_tol),
         "--bbox-abs-tol", str(bbox_tol)],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    sol = subprocess.run(
        [sys.executable, str(repo / "tools" / "str_solution_compare.py"),
         str(baseline), str(candidate),
         "--rel-tol", str(sol_tol)],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    return geo.returncode, geo.stdout.strip(), sol.returncode, sol.stdout.strip()


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deck-dir", type=Path, default=Path("/home/chris/Open_TCAD/decks/rit_pmos"))
    parser.add_argument("--spacings", nargs="+", default=["0.0005", "0.0006"])
    parser.add_argument("--before", default="#Pattern Gate Oxide Region",
                        help="comment/marker in the deck before which the layered remesh is inserted")
    parser.add_argument("--checkpoints", nargs="+", default=DEFAULT_CHECKPOINTS)
    parser.add_argument("--area-rel-tol", type=float, default=5.0e-2)
    parser.add_argument("--bbox-abs-tol", type=float, default=1.0e-1)
    parser.add_argument("--sol-rel-tol", type=float, default=2.5e-1)
    parser.add_argument("--keep", action="store_true")
    args = parser.parse_args(argv)

    repo = Path(__file__).resolve().parents[1]
    temp_root = Path(tempfile.mkdtemp(prefix="layered_remesh_smoke_"))
    failed = False
    try:
        print(f"workdir: {temp_root}")
        for spacing in args.spacings:
            run_dir = temp_root / f"spac_{spacing}"
            run_dir.mkdir()
            rc = run_pipeline(repo, args.deck_dir, spacing, args.before, run_dir)
            if rc != 0:
                print(f"FAIL spac={spacing}: pipeline exited {rc} (see {run_dir/'pipeline.log'})")
                failed = True
                continue
            print(f"PASS spac={spacing} pipeline")
            for cp in args.checkpoints:
                base = args.deck_dir / cp
                cand = run_dir / cp
                if not cand.exists():
                    print(f"FAIL spac={spacing} {cp}: missing {cand}")
                    failed = True
                    continue
                grc, gout, src, sout = compare(repo, base, cand,
                                               args.area_rel_tol, args.bbox_abs_tol, args.sol_rel_tol)
                if grc != 0:
                    print(f"FAIL spac={spacing} {cp} geometry:\n{gout}")
                    failed = True
                else:
                    print(f"PASS spac={spacing} {cp} geometry")
                if src != 0:
                    print(f"WARN spac={spacing} {cp} solution drift:\n{sout}")
                else:
                    print(f"PASS spac={spacing} {cp} solution")
        return 1 if failed else 0
    finally:
        if args.keep:
            print(f"kept: {temp_root}")
        else:
            shutil.rmtree(temp_root)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
