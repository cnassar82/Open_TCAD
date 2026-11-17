#!/usr/bin/env python3
"""Sweep SUPREM boron-diffusion and Si/SiO2 segregation + OED parameters to
find values that best match a reference 4000_ox.str (e.g. from Silvaco).

Tunable parameters (any subset; each accepts a single value or comma-separated
list, blank means "use modelrc default, don't inject an override line"):

    --dix0    boron silicon Dix.0=<value>              (neutral B diffusion preexp)
    --dip0    boron silicon Dip.0=<value>              (B+ pair diffusion preexp)
    --theta0  interstitial silicon /oxide theta.0=<v>  (OED injection fraction)
    --seg0    boron silicon /oxide Seg.0=<value>       (Si/SiO2 segregation coef)
    --trn0    boron silicon /oxide Trn.0=<value>       (Si/SiO2 transfer rate)

For each point in the Cartesian product the script:
  1. Truncates the deck after `structure outfile=4000_ox.str`.
  2. Injects the chosen override lines after `initialize phosphorus ...`.
  3. Runs SUPREM (via bin/remesh when the deck contains a `# remesh` marker,
     otherwise bin/suprem4 directly).
  4. Reads the resulting 4000_ox.str and reports the junction depth at the
     given cutline x. Scoring is junction-only: score = ((j - ref_j)/ref_j)^2.

Default reference is junction = 1.892 um at x=35 (Silvaco for the
rit_pmos_remesh deck after 4000A oxide).
"""

import argparse
import importlib.util
import itertools
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def load_compare_module():
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location("compare_silvaco", here / "compare_silvaco.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CMP = load_compare_module()


# Parameter name -> override-line template. Each template gets `.format(v=value)`.
PARAM_OVERRIDES = {
    "dix0":      "boron silicon Dix.0={v:.6g}\n",
    "dip0":      "boron silicon Dip.0={v:.6g}\n",
    "theta0":    "interstitial silicon /oxide theta.0={v:.6g}\n",
    "seg0":      "boron silicon /oxide Seg.0={v:.6g}\n",
    "trn0":      "boron silicon /oxide Trn.0={v:.6g}\n",
    "oxgastrn0": "boron oxide /gas Trn.0={v:.6g}\n",
    "oxgastrne": "boron oxide /gas Trn.E={v:.6g}\n",
}


def parse_value_list(s, dtype=float):
    if s is None or s.strip() == "":
        return [None]  # sentinel: use modelrc default, don't inject
    return [dtype(x.strip()) for x in s.split(",") if x.strip()]


def make_truncated_deck(src_deck, dst_deck, point, stop_after="4000_ox.str"):
    """Write a truncated copy of src_deck with parameter overrides for `point`
    (dict of param_name -> value or None) injected just before the `implant`
    line (so they take effect during the thermal step of interest, and survive
    being on the right side of any `# remesh` marker that splits the deck for
    bin/remesh), and the deck cut after `structure outfile=<stop_after>`.

    Falls back to inserting after the `initialize` line if no `implant` is
    found.
    """
    text = src_deck.read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)

    impl_idx = next((i for i, ln in enumerate(lines)
                     if re.match(r"\s*implant\s+", ln, re.IGNORECASE)), None)
    init_idx = next((i for i, ln in enumerate(lines)
                     if re.match(r"\s*initialize\s+", ln, re.IGNORECASE)), None)
    insert_idx = impl_idx if impl_idx is not None else (init_idx + 1 if init_idx is not None else None)
    if insert_idx is None:
        raise RuntimeError("could not find an `implant` or `initialize` line in the deck")
    stop_idx = next((i for i, ln in enumerate(lines)
                     if f"outfile={stop_after}" in ln), None)
    if stop_idx is None:
        raise RuntimeError(f"could not find `structure outfile={stop_after}` in the deck")

    body = ["\n# === parameter overrides (sweep_boron_params.py) ===\n"]
    for name, value in point.items():
        if value is None:
            continue
        body.append(PARAM_OVERRIDES[name].format(v=value))
    body.append("# ====================================================\n\n")

    new_lines = lines[: insert_idx] + body + lines[insert_idx : stop_idx + 1]
    new_lines.append("# truncated by sweep_boron_params.py\n")
    dst_deck.write_text("".join(new_lines), encoding="utf-8")


def suprem_env(binary):
    env = os.environ.copy()
    data_dir = binary.parent.parent / "data"
    if data_dir.is_dir():
        env.setdefault("SUP4KEYFILE", str(data_dir / "suprem.uk"))
        env.setdefault("SUP4MODELRC", str(data_dir / "modelrc"))
        env.setdefault("SUP4IMPDATA", str(data_dir / "sup4gs.imp"))
    return env


def stage_deck_dir(deck, work_dir):
    for item in deck.parent.iterdir():
        if item.is_file():
            dst = work_dir / item.name
            if item.resolve() != dst.resolve():
                shutil.copy2(item, dst)


def deck_has_remesh_marker(deck_path):
    """Return True if the deck contains a `# remesh` marker line."""
    try:
        for ln in deck_path.read_text(encoding="utf-8").splitlines():
            if re.match(r"^\s*#\s*remesh\b", ln, re.IGNORECASE):
                return True
    except Exception:
        pass
    return False


def label_for_point(point):
    parts = []
    for name in ("dix0", "dip0", "theta0", "seg0", "trn0", "oxgastrn0", "oxgastrne"):
        v = point.get(name)
        if v is None:
            continue
        parts.append(f"{name}={v:.4g}")
    return "_".join(parts) if parts else "default"


def run_one(binary, src_deck, point, work_root, cutline_x, y0, y1, samples,
            background, timeout):
    """Run SUPREM on one parameter point. Returns dict of metrics."""
    label = label_for_point(point)
    work_dir = work_root / label
    work_dir.mkdir(parents=True, exist_ok=True)
    stage_deck_dir(src_deck, work_dir)
    deck = work_dir / src_deck.name
    make_truncated_deck(src_deck, deck, point)
    log_path = work_dir / "suprem.log"

    env = suprem_env(binary)
    use_remesh = deck_has_remesh_marker(deck) and binary.name != "remesh"

    if use_remesh:
        # Dispatch via bin/remesh so the deck's # remesh markers are honored.
        remesh_bin = binary.parent / "remesh"
        if remesh_bin.exists():
            cmd = [str(remesh_bin), str(deck)]
            with log_path.open("wb") as log:
                try:
                    proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                                          cwd=work_dir, env=env, timeout=timeout)
                except subprocess.TimeoutExpired:
                    return {"label": label, "point": point, "status": "timeout"}
        else:
            use_remesh = False  # fall through to direct suprem4
    if not use_remesh:
        with deck.open("rb") as stdin, log_path.open("wb") as log:
            try:
                proc = subprocess.run([str(binary)], stdin=stdin, stdout=log,
                                      stderr=subprocess.STDOUT, cwd=work_dir, env=env,
                                      timeout=timeout)
            except subprocess.TimeoutExpired:
                return {"label": label, "point": point, "status": "timeout"}

    if proc.returncode != 0:
        return {"label": label, "point": point,
                "status": f"suprem exit {proc.returncode}", "log": str(log_path)}

    str_path = work_dir / "4000_ox.str"
    if not str_path.exists():
        return {"label": label, "point": point,
                "status": "4000_ox.str missing", "log": str(log_path)}

    coords, tris, tri_regs, regs, n_rec = CMP._parse_str(str_path)
    si_int, _ = CMP.integrate_boron_in_material(coords, tris, tri_regs, regs, n_rec,
                                                target_material=3)
    ys, b = CMP.sample_vertical_cutline(coords, tris, n_rec, cutline_x, y0, y1, samples)
    import numpy as np
    finite = np.isfinite(b) & (b > 0)
    if not np.any(finite):
        return {"label": label, "point": point,
                "status": "no finite boron at cutline", "log": str(log_path)}
    y_f = ys[finite]
    b_f = b[finite]
    peak_i = int(np.argmax(b_f))
    peak_y = float(y_f[peak_i])
    junction_y = None
    for yi, bi in zip(y_f[peak_i:], b_f[peak_i:]):
        if bi < background:
            junction_y = float(yi)
            break
    return {"label": label, "point": point, "status": "ok",
            "si_int": si_int, "peak_y": peak_y,
            "junction_y": junction_y if junction_y is not None else float("nan")}


def main(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--deck", type=Path, required=True,
                   help="SUPREM deck (.in) to base the sweep on")
    p.add_argument("--binary", type=Path,
                   default=Path("/home/chris/Open_TCAD/suprem4/bin/suprem4"),
                   help="SUPREM binary; the script will auto-prefer bin/remesh "
                        "when the deck contains a `# remesh` marker")
    p.add_argument("--dix0", default="",
                   help="comma-separated Dix.0 values (cm^2/s); blank=modelrc default")
    p.add_argument("--dip0", default="",
                   help="comma-separated Dip.0 values (cm^2/s) for B+ pair diffusivity; "
                        "this is the DOMINANT diffusion knob in p-type silicon "
                        "(modelrc default 0.76 = ~20x bigger than Dix.0)")
    p.add_argument("--theta0", default="",
                   help="comma-separated interstitial /oxide theta.0 values; "
                        "blank=modelrc default (1.970)")
    p.add_argument("--seg0", default="",
                   help="comma-separated boron /oxide Seg.0 values; "
                        "blank=modelrc default (1126)")
    p.add_argument("--trn0", default="",
                   help="comma-separated boron /oxide Trn.0 values; "
                        "blank=modelrc default (1.66e-7)")
    p.add_argument("--oxgas-trn0", default="", dest="oxgastrn0",
                   help="comma-separated `boron oxide /gas Trn.0` values "
                        "(cm/s pre-exponential of B-loss velocity through "
                        "the top oxide surface). Modelrc default 27.9.")
    p.add_argument("--oxgas-trne", default="", dest="oxgastrne",
                   help="comma-separated `boron oxide /gas Trn.E` values "
                        "(eV activation for B-loss velocity). Modelrc default "
                        "2.48 — this is what makes the surface essentially "
                        "closed at 1100 C; lowering it opens the surface as "
                        "a sink.")
    p.add_argument("--ref-junction", type=float, default=1.892,
                   help="reference junction depth in um at the cutline "
                        "(Silvaco default for rit_pmos_remesh: 1.892)")
    p.add_argument("--ref-si-int", type=float, default=None,
                   help="optional reference Si-side boron integral (atoms/cm). "
                        "If set, score becomes "
                        "((junction-err)^2 + w_si*(si_int-err)^2)/(1+w_si). "
                        "Silvaco default for rit_pmos_remesh post-4000A: 9.748e12")
    p.add_argument("--w-si", type=float, default=1.0,
                   help="weight for Si-int term relative to junction term in score")
    p.add_argument("--x", type=float, default=35.0, help="cutline x position (um)")
    p.add_argument("--y0", type=float, default=-0.5)
    p.add_argument("--y1", type=float, default=3.0)
    p.add_argument("--samples", type=int, default=600)
    p.add_argument("--background", type=float, default=3.0e14)
    p.add_argument("--timeout", type=float, default=600.0,
                   help="per-run SUPREM timeout in seconds")
    p.add_argument("--work-dir", type=Path,
                   help="where to put sweep working dirs (default: tempdir)")
    p.add_argument("--csv", type=Path, help="optional CSV of all results")
    args = p.parse_args(argv)

    axes = {
        "dix0":      parse_value_list(args.dix0),
        "dip0":      parse_value_list(args.dip0),
        "theta0":    parse_value_list(args.theta0),
        "seg0":      parse_value_list(args.seg0),
        "trn0":      parse_value_list(args.trn0),
        "oxgastrn0": parse_value_list(args.oxgastrn0),
        "oxgastrne": parse_value_list(args.oxgastrne),
    }
    grid = [dict(zip(axes.keys(), vals))
            for vals in itertools.product(*axes.values())]
    n_total = len(grid)
    nontrivial = [k for k, v in axes.items() if v != [None]]
    print(f"sweep axes: {', '.join(nontrivial) if nontrivial else '(none, single default-only run)'}")
    print(f"total runs: {n_total}")

    work_root = args.work_dir if args.work_dir else Path(tempfile.mkdtemp(prefix="suprem_boron_sweep_"))
    work_root.mkdir(parents=True, exist_ok=True)
    print(f"work root: {work_root}")

    results = []
    for k, point in enumerate(grid, start=1):
        label = label_for_point(point)
        print(f"  [{k}/{n_total}] {label:<40}", end=" ... ", flush=True)
        result = run_one(args.binary, args.deck, point, work_root,
                         args.x, args.y0, args.y1, args.samples,
                         args.background, args.timeout)
        if result["status"] == "ok":
            print(f"si_int={result['si_int']:.3e}  peak={result['peak_y']:.3f}  "
                  f"junction={result['junction_y']:.3f}")
        else:
            print(f"FAILED: {result['status']}")
        results.append(result)

    # Score and sort. If --ref-si-int is given, multi-objective:
    #   score = ((junction-err)^2 + w_si * (si-err)^2) / (1 + w_si)
    ok = [r for r in results if r["status"] == "ok"]
    for r in ok:
        ju = r["junction_y"]
        if ju != ju:
            r["err_jun"] = float("inf")
            r["err_si"] = float("inf")
            r["score"] = float("inf")
            continue
        r["err_jun"] = (ju - args.ref_junction) / args.ref_junction
        if args.ref_si_int is not None:
            r["err_si"] = (r["si_int"] - args.ref_si_int) / args.ref_si_int
            r["score"] = (r["err_jun"] ** 2 + args.w_si * r["err_si"] ** 2) / (1.0 + args.w_si)
        else:
            r["err_si"] = float("nan")
            r["score"] = r["err_jun"] ** 2
    ok.sort(key=lambda r: r["score"])

    print(f"\nreference: junction={args.ref_junction:.3f} um at x={args.x:.2f}")
    obj = "multi-objective (junction+Si)" if args.ref_si_int is not None else "junction-only"
    print(f"\nResults sorted by {obj} score (lower = better fit):")
    hdr = " ".join(f"{n:>10}" for n in ("dix0", "dip0", "theta0", "seg0", "trn0", "oxgastrn0", "oxgastrne"))
    has_si = args.ref_si_int is not None
    si_col = f" {'errSi':>9}" if has_si else ""
    print(f"{hdr} {'si_int':>11} {'peak':>7} {'junction':>9} {'errJun':>9}{si_col} {'score':>10}")
    for r in ok:
        pt = r["point"]
        cells = " ".join(f"{(pt.get(n) if pt.get(n) is not None else '-'):>10.4g}"
                         if pt.get(n) is not None else f"{'-':>10}"
                         for n in ("dix0", "dip0", "theta0", "seg0", "trn0", "oxgastrn0", "oxgastrne"))
        si_cell = f" {r['err_si']:>+9.2%}" if (has_si and r["err_si"] == r["err_si"]) else ""
        print(f"{cells} {r['si_int']:>11.3e} {r['peak_y']:>7.3f} "
              f"{r['junction_y']:>9.3f} {r['err_jun']:>+9.2%}{si_cell} {r['score']:>10.5f}")
    if any(r["status"] != "ok" for r in results):
        print(f"\nFailed runs:")
        for r in results:
            if r["status"] != "ok":
                print(f"  {label_for_point(r['point'])}: {r['status']}")

    if args.csv:
        import csv
        with args.csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["dix0", "dip0", "theta0", "seg0", "trn0",
                        "oxgastrn0", "oxgastrne", "status",
                        "si_int", "peak_y", "junction_y", "err_jun", "score"])
            for r in results:
                pt = r["point"]
                w.writerow([pt.get("dix0"), pt.get("dip0"), pt.get("theta0"),
                            pt.get("seg0"), pt.get("trn0"),
                            pt.get("oxgastrn0"), pt.get("oxgastrne"),
                            r.get("status"),
                            r.get("si_int", ""), r.get("peak_y", ""),
                            r.get("junction_y", ""),
                            r.get("err_jun", ""), r.get("score", "")])
        print(f"\ncsv: {args.csv}")

    if ok:
        best = ok[0]
        pt = best["point"]
        active = {k: v for k, v in pt.items() if v is not None}
        print(f"\nBest: {active}")
        print(f"  junction={best['junction_y']:.4f} um  (ref {args.ref_junction:.4f}, "
              f"err {best['err_jun']:+.2%})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
