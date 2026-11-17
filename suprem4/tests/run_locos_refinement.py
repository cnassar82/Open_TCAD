#!/usr/bin/env python3
"""Test LOCOS mesh-density sensitivity without changing geometry or process.

Substrate-only cases isolate initial mesh changes. Full-mesh cases also refine
the deposited pad/nitride layers and the oxide insertion spacing. Approximate
density multipliers are labels; measured point counts are recorded separately.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import tempfile
import time

from run_isolation import CHECKPOINTS, audit, geometry, make_deck

ROOT = Path(__file__).resolve().parents[1]
CASES = [("baseline", 1, False), ("substrate_2x", 2, False),
         ("substrate_4x", 4, False), ("full_2x", 2, True),
         ("full_4x", 4, True)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--cases", nargs="+", choices=[c[0] for c in CASES])
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    binary = args.binary.resolve()
    work = Path(tempfile.mkdtemp(prefix="suprem4_locos_refinement_"))
    print(f"workdir: {work}", flush=True)
    cases = [case for case in CASES if not args.cases or case[0] in args.cases]

    def run(case):
        name, density, layers = case
        scale = density**-0.5
        directory = work / name
        directory.mkdir()
        text = make_deck("locos", scale, 0)
        if layers:
            text = re.sub(r"\bdivisions=(\d+)",
                          lambda m: f"divisions={math.ceil(int(m[1])/scale)}", text)
            text = text.replace("grid.oxide=0.03", f"grid.oxide={0.03*scale:.10g}")
        deck = directory / "locos.in"
        deck.write_text(text)
        env = os.environ.copy()
        data = binary.parent.parent / "data"
        for variable, filename in [("SUP4KEYFILE", "suprem.uk"),
                                   ("KEYFILE", "suprem.uk"),
                                   ("SUP4MODELRC", "modelrc"),
                                   ("SUP4IMPDATA", "sup4gs.imp")]:
            env[variable] = str(data / filename)
        result = {"case": name, "nominal_initial_density": density,
                  "spacing_scale": scale, "refine_layers": layers,
                  "binary": str(binary), "directory": str(directory),
                  "timeout_seconds": args.timeout, "timed_out": False,
                  "checkpoints": {}, "errors": []}
        result["expected_process_seconds"] = 60*float(
            re.search(r"diffuse time=([\d.eE+-]+)", text)[1])
        started = time.monotonic()
        with deck.open("rb") as stdin, (directory / "suprem4.log").open("wb") as log:
            proc = subprocess.Popen([str(binary)], stdin=stdin, stdout=log,
                                    stderr=subprocess.STDOUT, cwd=directory,
                                    env=env, start_new_session=True)
            try:
                proc.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
                result["timed_out"] = True
                result["errors"].append("simulation exceeded time limit")
        result["elapsed_seconds"] = time.monotonic()-started
        result["returncode"] = proc.returncode
        log = (directory / "suprem4.log").read_text(errors="replace")
        panics = re.findall(r"suprem4 panic: ([^\n]+)", log)
        if panics:
            result["errors"].extend(panics)
        if re.search(r"errors detected on command input|No mesh defined", log):
            result["errors"].append("simulator command error")
        times = re.findall(r"Solving\s+([\d.eE+-]+)\s+\+\s+([\d.eE+-]+)"
                           r"\s+=\s+([\d.eE+-]+)", log)
        if times:
            result["last_logged_step_start_seconds"] = float(times[-1][0])
            result["last_logged_step_end_seconds"] = float(times[-1][2])
        if proc.returncode == 0 and result.get("last_logged_step_end_seconds", 0) < (
                result["expected_process_seconds"]-0.01):
            result["errors"].append("oxidation did not reach its requested duration")
        for checkpoint in CHECKPOINTS["locos"]:
            path = directory / checkpoint
            try:
                result["checkpoints"][checkpoint] = audit(path)
                result["checkpoints"][checkpoint]["points"] = sum(
                    line.startswith("c ") for line in path.read_text().splitlines())
            except (OSError, ValueError, KeyError) as error:
                result["errors"].append(f"{checkpoint}: {error}")
        if (directory / "panic.str").exists():
            try:
                result["failure_mesh_audit"] = audit(directory / "panic.str")
            except (ValueError, KeyError) as error:
                result["failure_mesh_audit_error"] = str(error)
        if (directory / "locos.str").exists():
            try:
                result["geometry"] = geometry("locos", directory, 0)
            except (OSError, ValueError, KeyError) as error:
                result["errors"].append(str(error))
        result["passed"] = proc.returncode == 0 and not result["errors"]
        (directory / "result.json").write_text(json.dumps(result, indent=2)+"\n")
        points = result["checkpoints"].get("initial.str", {}).get("points")
        print(f"{'PASS' if result['passed'] else 'FAIL'} {name}: initial={points}, "
              f"exit={proc.returncode}, elapsed={result['elapsed_seconds']:.1f}s, "
              f"{result['errors'][:1]}", flush=True)
        return result

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        results = list(pool.map(run, cases))
    baseline = next((r for r in results if r["case"] == "baseline"), None)
    if baseline and baseline["passed"]:
        for result in results:
            initial = result["checkpoints"].get("initial.str", {}).get("points")
            if initial:
                result["measured_initial_density_ratio"] = (
                    initial / baseline["checkpoints"]["initial.str"]["points"])
            if result["passed"]:
                result["oxide_area_change_percent"] = 100*(
                    result["geometry"]["materials"][1]["area"] /
                    baseline["geometry"]["materials"][1]["area"]-1)
            (Path(result["directory"]) / "result.json").write_text(
                json.dumps(result, indent=2)+"\n")
    (work / "results.json").write_text(json.dumps(results, indent=2)+"\n")
    return int(not all(r["passed"] for r in results))


if __name__ == "__main__":
    raise SystemExit(main())
