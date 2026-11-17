#!/usr/bin/env python3
"""Run RIT PMOS mesh regression for selected initial surface spacings."""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


CHECKPOINTS = [
    "inital_ox.str",
    "diffusion_photo_etch.str",
    "diffusion_implant.str",
    "4000_ox.str",
    "700_ox.str",
    "cont_etch.str",
    "dep_alumin.str",
    "etch_alumin.str",
    "pmos.str",
]


def make_deck(src, dst, spacing, remesh_before=None):
    text = src.read_text(encoding="utf-8")
    old = "line y loc=0.0  spac=0.0005 tag=surf"
    new = f"line y loc=0.0  spac={spacing} tag=surf"
    if old not in text:
        raise RuntimeError(f"could not find expected surface spacing line in {src}")
    text = text.replace(old, new, 1)
    if remesh_before:
        if remesh_before not in text:
            raise RuntimeError(f"could not find remesh insertion point {remesh_before!r}")
        insert = (
            "structure outfile=before_remesh.str\n"
            "remesh all\n"
            "structure outfile=after_remesh.str\n"
        )
        text = text.replace(remesh_before, insert + remesh_before, 1)
    dst.write_text(text, encoding="utf-8")


def run_deck(binary, deck, cwd, mesh_diag):
    env = os.environ.copy()
    if mesh_diag:
        env["SUPREM4_MESH_DIAG"] = "1"
    data_dir = binary.parent.parent / "data"
    if data_dir.is_dir():
        env.setdefault("SUP4KEYFILE",  str(data_dir / "suprem.uk"))
        env.setdefault("SUP4MODELRC",  str(data_dir / "modelrc"))
        env.setdefault("SUP4IMPDATA",  str(data_dir / "sup4gs.imp"))
    log_path = cwd / "suprem4.log"
    with deck.open("rb") as stdin, log_path.open("wb") as log:
        proc = subprocess.run([str(binary)], stdin=stdin, stdout=log, stderr=subprocess.STDOUT, cwd=cwd, env=env)
    return proc.returncode, log_path


def compare_checkpoint(compare_tool, baseline_dir, candidate_dir, name, area_rel_tol, bbox_abs_tol):
    baseline = baseline_dir / name
    candidate = candidate_dir / name
    return compare_files(compare_tool, baseline, candidate, area_rel_tol, bbox_abs_tol)


def compare_files(compare_tool, baseline, candidate, area_rel_tol, bbox_abs_tol):
    if not candidate.exists():
        return False, f"missing output {candidate}"
    proc = subprocess.run(
        [
            sys.executable,
            str(compare_tool),
            str(baseline),
            str(candidate),
            "--area-rel-tol",
            str(area_rel_tol),
            "--bbox-abs-tol",
            str(bbox_abs_tol),
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return proc.returncode == 0, proc.stdout.strip()


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deck-dir", type=Path, default=Path("/home/chris/Open_TCAD/decks/rit_pmos"))
    parser.add_argument("--binary", type=Path, default=Path("/home/chris/Open_TCAD/suprem4/bin/suprem4"))
    parser.add_argument("--spacings", nargs="+", default=["0.0005", "0.0006"])
    parser.add_argument("--keep", action="store_true", help="keep temporary output directory")
    parser.add_argument("--mesh-diag", action="store_true", help="enable SUPREM4_MESH_DIAG")
    parser.add_argument(
        "--remesh-before",
        help="insert structure/remesh-all/structure before the first matching deck line",
    )
    parser.add_argument("--area-rel-tol", type=float, default=2.0e-2)
    parser.add_argument("--bbox-abs-tol", type=float, default=5.0e-2)
    args = parser.parse_args(argv)

    repo = Path(__file__).resolve().parents[1]
    compare_tool = repo / "tools" / "str_compare.py"
    source_deck = args.deck_dir / "rit_pmos.in"
    temp_root = Path(tempfile.mkdtemp(prefix="suprem4_pmos_regression_"))
    failed = False

    try:
        print(f"workdir: {temp_root}")
        for spacing in args.spacings:
            run_dir = temp_root / f"spac_{spacing}"
            run_dir.mkdir()
            for item in ["modelrc", "sup4gs.imp", "suprem.uk"]:
                src = args.deck_dir / item
                if src.exists():
                    shutil.copy2(src, run_dir / item)
            deck = run_dir / "rit_pmos.in"
            make_deck(source_deck, deck, spacing, args.remesh_before)
            code, log_path = run_deck(args.binary, deck, run_dir, args.mesh_diag)
            if code != 0:
                print(f"FAIL spac={spacing}: suprem4 exited {code}; log={log_path}")
                failed = True
                continue
            if args.remesh_before:
                ok, output = compare_checkpoint(
                    compare_tool,
                    run_dir / "before_remesh.str",
                    run_dir / "after_remesh.str",
                    1.0e-6,
                    1.0e-10,
                )
                if ok:
                    print(f"PASS spac={spacing} remesh geometry")
                else:
                    print(f"FAIL spac={spacing} remesh geometry:")
                    print(output)
                    failed = True
            for checkpoint in CHECKPOINTS:
                ok, output = compare_checkpoint(
                    compare_tool,
                    args.deck_dir,
                    run_dir,
                    checkpoint,
                    args.area_rel_tol,
                    args.bbox_abs_tol,
                )
                if not ok:
                    print(f"FAIL spac={spacing} {checkpoint}:")
                    print(output)
                    failed = True
                else:
                    print(f"PASS spac={spacing} {checkpoint}")
        return 1 if failed else 0
    finally:
        if args.keep:
            print(f"kept: {temp_root}")
        else:
            shutil.rmtree(temp_root)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
