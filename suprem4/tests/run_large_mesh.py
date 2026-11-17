#!/usr/bin/env python3
"""Exercise large-region connectivity through a real SUPREM initialize/write.

Checks completion, connectivity, finite coordinates and conserved domain area.
This is a mesh construction test, not a diffusion/oxidation accuracy test.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=root / "src/suprem4")
    parser.add_argument("--keep", action="store_true")
    parser.add_argument("--cells", type=int, nargs="+", default=[160],
                        help="cells along each axis of the 160 um square")
    parser.add_argument("--timeout", type=float, default=300)
    args = parser.parse_args()
    if any(cells < 130 for cells in args.cells):
        parser.error("--cells must be at least 130 to exercise large lists")
    if len(set(args.cells)) != len(args.cells):
        parser.error("--cells must be unique")
    failed = False
    for cells in args.cells:
        if not run_case(root, args, cells):
            failed = True
    return int(failed)


def run_case(root, args, cells):
    work = Path(tempfile.mkdtemp(prefix="suprem4_large_mesh_"))
    try:
        env = os.environ.copy()
        for variable, name in (("SUP4KEYFILE", "suprem.uk"),
                               ("SUP4MODELRC", "modelrc"),
                               ("SUP4IMPDATA", "sup4gs.imp")):
            env[variable] = str(root / "data" / name)
        deck = work / "large_mesh.in"
        text = (root / "tests/large_mesh.in").read_text()
        deck.write_text(text.replace("spac=1 ", f"spac={160 / cells:.15g} "))
        started = time.monotonic()
        command = [str(args.binary.resolve()), str(deck)]
        metrics = work / "resource.json"
        if Path("/usr/bin/time").exists():
            command = ["/usr/bin/time", "-f", '{"max_rss_kib":%M}',
                       "-o", str(metrics), *command]
        with (work / "suprem4.log").open("w") as log:
            result = subprocess.Popen(
                command,
                cwd=work, env=env, stdin=subprocess.DEVNULL,
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
            )
            try:
                result.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(result.pid, signal.SIGKILL)
                result.wait()
                raise
        if result.returncode:
            raise RuntimeError(f"SUPREM exited {result.returncode}")
        elapsed = time.monotonic() - started
        # Reuse the full topology and field audit as well as domain-area checks.
        from run_oxidation_perturbations import audit
        topology = audit(work / "large_mesh.str")
        points, triangles = {}, []
        for line in (work / "large_mesh.str").read_text().splitlines():
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "c":
                point = tuple(map(float, fields[2:4]))
                if not all(map(math.isfinite, point)):
                    raise RuntimeError("nonfinite coordinate")
                points[int(fields[1])] = point
            elif fields[0] == "t":
                triangles.append(tuple(map(int, fields[3:6])))
        if len(triangles) <= 32767:
            raise RuntimeError(f"mesh too small: {len(triangles)} triangles")
        seen, area = set(), 0.0
        for triangle in triangles:
            key = tuple(sorted(triangle))
            if len(set(triangle)) != 3 or key in seen:
                raise RuntimeError("degenerate or duplicate triangle")
            seen.add(key)
            a, b, c = (points[i] for i in triangle)
            element_area = abs((b[0] - a[0]) * (c[1] - a[1])
                               - (c[0] - a[0]) * (b[1] - a[1])) / 2
            if element_area <= 0:
                raise RuntimeError("zero-area triangle")
            area += element_area
        if not math.isclose(area, 160 * 160, rel_tol=1e-6):
            raise RuntimeError(f"domain area changed: {area}")
        if len(points) != (cells + 1) ** 2 or len(triangles) != 2 * cells ** 2:
            raise RuntimeError("unexpected structured mesh dimensions")
        report = dict(binary=str(args.binary.resolve()), cells=cells, points=len(points),
                      triangles=len(triangles), area_um2=area, elapsed_seconds=elapsed,
                      audit=topology)
        if metrics.exists():
            report.update(json.loads(metrics.read_text()))
        (work / "results.json").write_text(json.dumps(report, indent=2) + "\n")
        print(f"PASS: {len(points)} points, {len(triangles)} triangles, "
              f"area={area:g} um^2, {elapsed:.2f}s, "
              f"peak RSS={report.get('max_rss_kib', 'unmeasured')} KiB")
    except (OSError, RuntimeError, ValueError, KeyError, subprocess.TimeoutExpired) as exc:
        args.keep = True
        print(f"FAIL: {exc}")
        return False
    finally:
        if args.keep:
            print(f"Artifacts: {work}")
        else:
            shutil.rmtree(work)
    return True


if __name__ == "__main__":
    raise SystemExit(main())
