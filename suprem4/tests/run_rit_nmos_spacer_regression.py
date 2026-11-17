#!/usr/bin/env python3
"""Replay the first malformed RIT NMOS spacer boundary from a saved LDD mesh."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

from run_oxidation_perturbations import audit
from partition_audit import audit_partition

FIXTURE = Path(__file__).resolve().parent / "fixtures/rit_nmos_spacer"
COMMANDS = (
    ("spacer_deposit", "deposit oxide thick=0.5 divisions=10 space=0.04"),
    ("spacer_etch", "etch oxide dry thick=0.5"),
    ("resist_deposit", "deposit photoresist file=resist_top.dat divisions=8"),
)


def check_resist_coverage(path, mesh_audit):
    """The fixture has a 10 x 8 um domain and a flat resist top at -2 um.

    Together with the nonoverlap audit, the 100 um^2 rectangular-domain area
    checks for partial surface walks and missing material after coating.
    """
    coordinates, regions, resist_points = {}, {}, set()
    for line in path.read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "c":
            coordinates[int(parts[1])] = tuple(map(float, parts[2:4]))
        elif parts[0] == "r":
            regions[int(parts[1])] = int(parts[2])
        elif parts[0] == "t" and regions[int(parts[2])] == 7:
            resist_points.update(map(int, parts[3:6]))
    if not resist_points:
        raise ValueError("missing deposited resist")
    xmin = min(coordinates[p][0] for p in resist_points)
    xmax = max(coordinates[p][0] for p in resist_points)
    top = min(coordinates[p][1] for p in resist_points)
    area = sum(mesh_audit["area_by_region_um2"].values())
    if abs(xmin) > 1e-5 or abs(xmax - 10) > 1e-5 or abs(top + 2) > 1e-5:
        raise ValueError(f"partial resist coverage: x={xmin},{xmax}, top={top}")
    if abs(area - 100) > 1e-4:
        raise ValueError(f"coated material area {area} um^2, expected 100")
    return {"xmin_um": xmin, "xmax_um": xmax, "top_um": top,
            "coated_area_um2": area}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--device", choices=["nmos", "pmos"], default="nmos")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--timeout", type=float, default=60)
    args = parser.parse_args()
    fixture = FIXTURE.parent / f"rit_{args.device}_spacer"
    binary = args.binary.resolve()
    work = (args.output or Path(tempfile.mkdtemp(prefix="rit_nmos_spacer_"))).resolve()
    if work.exists() and any(work.iterdir()):
        parser.error(f"output folder is not empty: {work}")
    work.mkdir(parents=True, exist_ok=True)
    mesh = gzip.decompress((fixture / "ldd.str.gz").read_bytes())
    expected = json.loads((fixture / "fixture.json").read_text())["sha256"]
    if hashlib.sha256(mesh).hexdigest() != expected:
        parser.error("fixture checksum mismatch")
    (work / "ldd.str").write_bytes(mesh)
    (work / "resist_top.dat").write_bytes((fixture / "resist_top.dat").read_bytes())
    lines = ["structure infile=ldd.str"]
    for name, command in COMMANDS:
        lines += [command, f"structure outfile={name}.str"]
    lines += ["quit"]
    deck = "\n".join(lines) + "\n"
    (work / "replay.in").write_text(deck)
    data = binary.parent.parent / "data"
    env = os.environ.copy()
    for key, name in (("SUP4KEYFILE", "suprem.uk"), ("KEYFILE", "suprem.uk"),
                      ("SUP4MODELRC", "modelrc"), ("SUP4IMPDATA", "sup4gs.imp"),
                      ("IMPFILE", "sup4gs.imp")):
        env[key] = str(data / name)
    report = {"device": args.device, "binary": str(binary), "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
              "fixture_sha256": expected, "deck_sha256": hashlib.sha256(deck.encode()).hexdigest(),
              "data_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in data.iterdir() if p.name in ("suprem.uk", "modelrc", "sup4gs.imp")},
              "workdir": str(work), "stages": {}}
    build_manifest = binary.parent.parent / "build.json"
    if build_manifest.is_file():
        report["build_manifest_sha256"] = hashlib.sha256(build_manifest.read_bytes()).hexdigest()
    with (work / "replay.in").open("rb") as source, (work / "replay.log").open("wb") as log:
        try:
            process = subprocess.run([str(binary)], stdin=source, stdout=log,
                                     stderr=subprocess.STDOUT, cwd=work, env=env,
                                     timeout=args.timeout)
            report["returncode"] = process.returncode
        except subprocess.TimeoutExpired:
            report["returncode"] = None
            report["timed_out"] = True
    log_text = (work / "replay.log").read_text(errors="replace")
    report["simulator_error"] = bool(re.search(
        r"suprem4 panic:|errors detected on command input|ERROR:|runtime error:|Newton convergence failed",
        log_text, re.I))
    for name, _ in COMMANDS:
        try:
            saved = work / f"{name}.str"
            report["stages"][name] = {"status": "PASS", "audit": audit(saved),
                                     "partition": audit_partition(saved)}
            if name == "resist_deposit":
                report["stages"][name]["coverage"] = check_resist_coverage(
                    saved, report["stages"][name]["audit"])
        except (OSError, ValueError, KeyError) as error:
            report["stages"][name] = {"status": "FAIL", "error": str(error)}
    report["passed"] = report["returncode"] == 0 and not report["simulator_error"] and all(
        stage["status"] == "PASS" for stage in report["stages"].values())
    (work / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"{'PASS' if report['passed'] else 'FAIL'} RIT {args.device.upper()} spacer replay: {work}")
    for name, stage in report["stages"].items():
        print(f"  {name}: {stage['status']} {stage.get('error', '')}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
