#!/usr/bin/env python3
"""WP7 geometry-core test suite: SUPREM oxidation with the core on, pass/fail report.

Cases
  matrix   the six poly-line decks of doc/poly_oxidation/test_decks
  random   N randomized poly-line variants (gate width, poly and gate-oxide thickness,
           ambient, temperature, time, mesh spacing), reproducible from --seed
  locos    LOCOS field oxidation (SubCmos2014 NMOS locos_mask.str), viscous and compress,
           each also run with the core off for comparison

Every poly case is accompanied by the same recipe on a flat wafer with the core off
(plain SUPREM, about a minute) as the reference for field oxide and dose drift.

Checks per poly case (acceptance criteria 1-4 of doc/plan_m0_geometry_core.md, section 5)
  complete     after.str written within the timeout
  invariants   core_str: mesh invariants, gas nodes, B-rep
  suprem       core_check_suprem: SUPREM's own validity rules (WP6a contract)
  field        field oxide within 2 % of the flat plain-SUPREM run
  ratio        consumed / grown on silicon, poly top and both poly sidewalls in 0.44 +- 5 %,
               on faces at least 4x longer than their oxide is thick (others: 2D, reported)
  area         oxide area gained / silicon + poly area consumed within 3 % of 2.2
  symmetry     left / right sidewall oxide and gate-edge dielectric within 2 %
  dose         boron drift (box measure) |d| < 1e-4, or no worse than 1.5 x the flat
               plain-SUPREM drift of the same recipe (reported as "plain-limited")
LOCOS: complete, invariants, suprem, maximum field oxide within 1 % of the core-off run,
dose drift within 1e-4 of the core-off run.

usage
  run_geometry_suite.py --build build/native64_XXXX --core-build build/core_wp6a \
      --out OUT [--matrix] [--random N --seed S] [--locos] [--jobs J] [--timeout T]
  run_geometry_suite.py --report OUT      (re-evaluate an existing run)
"""

import argparse
import concurrent.futures as cf
import json
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUPREM = HERE.parent
sys.path.insert(0, str(SUPREM / "core" / "tests"))
import poly_metrics  # noqa: E402

DECKS = SUPREM / "doc" / "poly_oxidation" / "test_decks"
MATRIX = ["base", "coarse", "dry1000_120", "narrow025", "wet900_20", "wide080"]
LOCOS_MASK = SUPREM.parent / "decks/subcmos2014/results/m0_reference_20261122/nmos/locos_mask.str"

LOCOS_RECIPE = """\
oxide Vc=300 Vr=30 Vd=25 stress.dep=true
method {model} oxide.rel=1e-2 grid.oxide=0.04
diffuse time=12 temp=800 nitrogen
diffuse time=15 temp=800 nitrogen continue
diffuse time=5 temp=825 nitrogen continue
diffuse time=5 temp=875 nitrogen continue
diffuse time=5 temp=925 nitrogen continue
diffuse time=5 temp=975 nitrogen continue
diffuse time=5 temp=1025 nitrogen continue
diffuse time=5 temp=1075 nitrogen continue
diffuse time=5 temp=1100 dryo2
diffuse time=65 temp=1100 weto2 continue
diffuse time=5 temp=1100 nitrogen continue
diffuse time=9.166666667 temp=1075 nitrogen continue
diffuse time=9.166666667 temp=1025 nitrogen continue
diffuse time=9.166666667 temp=975 nitrogen continue
diffuse time=9.166666667 temp=925 nitrogen continue
diffuse time=9.166666667 temp=875 nitrogen continue
diffuse time=9.166666667 temp=825 nitrogen continue
diffuse time=15 temp=800 nitrogen continue
etch nitride all
"""

RATIO, RATIO_TOL = 0.44, 0.05
EXPANSION = 2.2   # SUPREM's oxide / consumed volume ratio
FLAT = 4.0   # a face is planar (1D) when its length >= FLAT x its oxide thickness


# ---------------------------------------------------------------- decks

def poly_deck(p, flat=False):
    """Poly-line deck from parameters (as poly_line_base.in)."""
    w, tp, tox, s = p["width"], p["poly"], p["gox"], p["spac"]
    xl, xr = 1.0 - w / 2, 1.0 + w / 2
    lines = [
        "option quiet",
        "line x loc=0 spac=0.05 tag=left",
        f"line x loc={xl - 0.1:.4f} spac={s}",
        f"line x loc=1.0 spac={2 * s}",
        f"line x loc={xr + 0.1:.4f} spac={s}",
        "line x loc=2.0 spac=0.05 tag=right",
        "line y loc=0 spac=0.005 tag=top",
        "line y loc=0.3 spac=0.03",
        "line y loc=2 spac=0.3 tag=bottom",
        "region silicon xlo=left xhi=right ylo=top yhi=bottom",
        "boundary exposed xlo=left xhi=right ylo=top yhi=top",
        "boundary backside xlo=left xhi=right ylo=bottom yhi=bottom",
        "initialize boron conc=1e15 orient=100",
        f"deposit oxide thick={tox:.4f} divisions=2",
    ]
    if not flat:
        lines += [
            f"deposit poly thick={tp:.3f} divisions={max(4, round(tp / 0.04))} space={s}",
            f"etch poly left p1.x={xl:.4f}",
            f"etch poly right p1.x={xr:.4f}",
        ]
    lines += [
        "structure outfile=before.str",
        "method compress grid.oxide=0.01",
        f"diffuse time={p['time']:g} temp={p['temp']:g} {p['ambient']}",
        "structure outfile=after.str",
    ]
    return "\n".join(lines) + "\n"


def matrix_params(name):
    """Parameters of a test deck, read back from the deck text."""
    p = dict(width=0.4, poly=0.4, gox=0.01, spac=0.01, time=51, temp=1000, ambient="dryo2")
    text = (DECKS / f"poly_line_{name}.in").read_text()
    for line in text.splitlines():
        t = line.split()
        if line.startswith("diffuse"):
            kv = dict(x.split("=") for x in t[1:] if "=" in x)
            p.update(time=float(kv["time"]), temp=float(kv["temp"]), ambient=t[-1])
        if line.startswith("etch poly left"):
            xl = float(line.split("p1.x=")[1])
            p["width"] = 2 * (1.0 - xl)
        if line.startswith("deposit poly") and "space=0.02" in line:
            p["spac"] = 0.02
    return p


def random_params(rng):
    wet = rng.random() < 0.35
    return dict(
        width=round(rng.uniform(0.2, 0.8), 3),
        poly=round(rng.uniform(0.25, 0.5), 3),
        gox=round(rng.uniform(0.005, 0.02), 4),
        spac=rng.choice([0.01, 0.01, 0.02]),
        ambient="weto2" if wet else "dryo2",
        temp=rng.choice([850, 875, 900, 925]) if wet else rng.choice([900, 950, 1000, 1050]),
        time=round(rng.uniform(5, 25) if wet else rng.uniform(10, 120), 1),
    )


# ---------------------------------------------------------------- running

def suprem_env(build, core):
    d = build / "data"
    e = dict(os.environ)
    e.update(SUP4KEYFILE=f"{d}/suprem.uk", KEYFILE=f"{d}/suprem.uk", SUP4MODELRC=f"{d}/modelrc",
             SUP4IMPDATA=f"{d}/sup4gs.imp", IMPFILE=f"{d}/sup4gs.imp")
    for k in [k for k in e if k.startswith("SUPREM4_GEOMETRY_CORE")]:
        del e[k]
    if core:
        e["SUPREM4_GEOMETRY_CORE"] = "1"
        shm = Path("/dev/shm")
        e.setdefault("SUPREM4_CORE_TMP", str(shm) if shm.is_dir() else ".")
    return e


def run_suprem(build, workdir, deck_text, core, timeout):
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "deck.in").write_text(deck_text)
    for f in ("after.str", "log"):
        (workdir / f).unlink(missing_ok=True)
    t0 = time.time()
    with open(workdir / "deck.in") as fin, open(workdir / "log", "w") as fout:
        try:
            subprocess.run([str(build / "src" / "suprem4")], stdin=fin, stdout=fout, stderr=subprocess.STDOUT,
                           cwd=workdir, env=suprem_env(build, core), timeout=timeout,
                           preexec_fn=lambda: __import__("resource").setrlimit(
                               __import__("resource").RLIMIT_AS, (6 << 30, 6 << 30)))
            status = "done"
        except subprocess.TimeoutExpired:
            status = "timeout"
    wall = time.time() - t0
    log = (workdir / "log").read_text(errors="replace")
    info = dict(status=status, wall_s=round(wall, 1),
                complete=(workdir / "after.str").exists(),
                core_calls=log.count("[geometry-core] call"),
                regenerations=log.count("regenerated"),
                local_repairs=log.count("events local"))
    if not info["complete"]:
        solving = [ln for ln in log.splitlines() if ln.startswith("Solving")]
        info["reached"] = solving[-1] if solving else ""
        bad = [ln for ln in log.splitlines() if "panic" in ln.lower() or "error" in ln.lower()]
        info["failure"] = bad[0][:160] if bad else status
    (workdir / "run.json").write_text(json.dumps(info, indent=1))
    return info


def tool_checks(core_build, path):
    out = {}
    r = subprocess.run([str(core_build / "core_str"), str(path)], capture_output=True, text=True)
    out["invariants"] = ("invariants ok" in r.stdout and "brep ok" in r.stdout)
    chk = core_build / "core_check_suprem"
    if chk.exists():
        r = subprocess.run([str(chk), str(path)], capture_output=True, text=True)
        out["suprem"] = r.returncode == 0
        if r.returncode != 0:
            out["suprem_msg"] = (r.stdout + r.stderr).strip().splitlines()[-1:]
    return out


def deck_texts(name, params):
    """(poly deck, flat deck): the test deck file for matrix cases, else generated."""
    f = DECKS / f"poly_line_{name}.in"
    if f.exists():
        text = f.read_text()
        flat = "".join(ln for ln in text.splitlines(True)
                       if not (ln.startswith("deposit poly") or ln.startswith("etch poly")))
        return text, flat
    return poly_deck(params), poly_deck(params, flat=True)


def poly_case(args, name, params):
    d = args.out / name
    res = dict(name=name, params=params)
    deck, flat_deck = deck_texts(name, params)
    flat = run_suprem(args.build, d / "flat", flat_deck, False, args.timeout)
    run = run_suprem(args.build, d / "core", deck, True, args.timeout)
    res["run"], res["flat_run"] = run, flat
    return evaluate_poly(args, d, res)


def evaluate_poly(args, d, res):
    run = res["run"]
    if run["complete"]:
        res.update(tool_checks(args.core_build, d / "core" / "after.str"))
        res["metrics"] = poly_metrics.measure(str(d / "core" / "before.str"), str(d / "core" / "after.str"),
                                              0.2, str(args.core_build / "core_str"))
    if res["flat_run"]["complete"]:
        res["flat_metrics"] = poly_metrics.measure(str(d / "flat" / "before.str"), str(d / "flat" / "after.str"),
                                                   0.2, str(args.core_build / "core_str"))
    res["checks"] = poly_checks(res)
    res["pass"] = all(v for v in res["checks"].values() if v is not None)
    (d / "result.json").write_text(json.dumps(res, indent=1, default=str))
    return res


def poly_checks(res):
    c = dict(complete=res["run"]["complete"], invariants=None, suprem=None, field=None, ratio=None,
             area=None, symmetry=None, dose=None)
    if not c["complete"]:
        return c
    m, f = res.get("metrics", {}), res.get("flat_metrics")
    c["invariants"] = res.get("invariants")
    c["suprem"] = res.get("suprem")
    if f:
        c["field"] = abs(m["field_oxide_A"] / f["field_oxide_A"] - 1) < 0.02
    # consumed = 0.44 x grown is a 1D relation: checked only on faces that stay
    # planar on the scale of the oxide (face length >= FLAT x oxide thickness);
    # narrower faces (thick oxide on a thin line) are 2D and only reported.
    faces = [("ratio_si", None, None),
             ("ratio_poly_top", "poly_width_mid_um", "poly_top_oxide_A"),
             ("ratio_poly_left", "poly_height_um", "sidewall_oxide_left_A"),
             ("ratio_poly_right", "poly_height_um", "sidewall_oxide_right_A")]
    ok, two_d = True, []
    for k, length, ox in faces:
        r = m.get(k)
        if length and m.get(length) is not None and m.get(ox) and m[length] * 1e4 < FLAT * m[ox]:
            two_d.append(k)
            continue
        ok = ok and r is not None and abs(r / RATIO - 1) <= RATIO_TOL
    c["ratio"] = ok
    if m.get("area_expansion") is not None:
        c["area"] = abs(m["area_expansion"] / EXPANSION - 1) <= 0.03
    res["ratio_2d_faces"] = two_d
    pairs = [("sidewall_oxide_left_A", "sidewall_oxide_right_A")] + [
        (f"gate_edge_left_{x}_A", f"gate_edge_right_{x}_A") for x in ("0.02", "0.05", "0.10")]
    sym = True
    for a, b in pairs:
        if a in m and b in m and max(m[a], m[b]) > 0:
            sym = sym and abs(m[a] - m[b]) <= 0.02 * max(m[a], m[b])
    c["symmetry"] = sym
    drift = m.get("inventory_5_rel_change")
    plain = f.get("inventory_5_rel_change") if f else None
    if drift is not None:
        c["dose"] = abs(drift) < 1e-4 or (plain is not None and abs(drift) <= 1.5 * abs(plain))
    return c


def locos_case(args, model):
    d = args.out / f"locos_{model}"
    d.mkdir(parents=True, exist_ok=True)
    deck = f"option quiet\ninitialize infile={LOCOS_MASK}\n" + LOCOS_RECIPE.format(model=model) + \
        "structure outfile=after.str\n"
    res = dict(name=f"locos_{model}")
    with cf.ThreadPoolExecutor(2) as ex:
        f_off = ex.submit(run_suprem, args.build, d / "off", deck, False, args.timeout)
        f_on = ex.submit(run_suprem, args.build, d / "core", deck, True, args.timeout)
        res["off_run"], res["run"] = f_off.result(), f_on.result()
    return evaluate_locos(args, d, res)


def max_oxide(path):
    tris = poly_metrics.parse(str(path))
    best = 0.0
    xs = sorted({round(p[0], 3) for _, *v in tris for p in v})
    for x in xs[:: max(1, len(xs) // 200)]:
        col = poly_metrics.column(tris, x)
        best = max(best, sum(hi - lo for lo, hi, m in col if m == poly_metrics.OXIDE))
    return best * 1e4


def boron_total(core_build, path):
    inv = poly_metrics.inventories(str(core_build / "core_str"), str(path))
    return sum(inv.get(5, {}).values())


def evaluate_locos(args, d, res):
    c = dict(complete=res["run"]["complete"] and res["off_run"]["complete"], invariants=None, suprem=None,
             field=None, dose=None)
    if c["complete"]:
        res.update(tool_checks(args.core_build, d / "core" / "after.str"))
        c["invariants"], c["suprem"] = res.get("invariants"), res.get("suprem")
        on, off = max_oxide(d / "core" / "after.str"), max_oxide(d / "off" / "after.str")
        res["field_oxide_A"] = dict(core=on, off=off)
        c["field"] = abs(on / off - 1) < 0.01
        b0 = boron_total(args.core_build, LOCOS_MASK)
        bon, boff = boron_total(args.core_build, d / "core" / "after.str"), boron_total(args.core_build, d / "off" / "after.str")
        res["boron_drift"] = dict(core=bon / b0 - 1, off=boff / b0 - 1)
        c["dose"] = abs((bon - boff) / b0) < 1e-4
    res["checks"] = c
    res["pass"] = all(v for v in c.values() if v is not None)
    (d / "result.json").write_text(json.dumps(res, indent=1, default=str))
    return res


# ---------------------------------------------------------------- report

def report(results, out):
    lines = ["# Geometry-core suite", "",
             "| case | pass | wall s | calls | regen | field A | Si/poly ratio | sidewall A L/R | dose drift (core / plain) | failed checks |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        run, m = r["run"], r.get("metrics", {})
        failed = [k for k, v in r["checks"].items() if v is False]
        if "locos" in r["name"]:
            fo = r.get("field_oxide_A", {})
            bd = r.get("boron_drift", {})
            lines.append(f"| {r['name']} | {'PASS' if r['pass'] else 'FAIL'} | {run['wall_s']:.0f} (off {r['off_run']['wall_s']:.0f}) | "
                         f"{run['core_calls']} | {run['regenerations']} | {fo.get('core', 0):.0f} (off {fo.get('off', 0):.0f}) | | | "
                         f"{bd.get('core', 0):+.2e} / {bd.get('off', 0):+.2e} | {', '.join(failed) or run.get('failure', '')} |")
            continue
        fm = r.get("flat_metrics", {})
        two_d = r.get("ratio_2d_faces", [])
        ratio = "/".join(f"{m[k]:.3f}" + ("*" if k in two_d else "")
                         for k in ("ratio_si", "ratio_poly_top", "ratio_poly_left") if k in m)
        side = f"{m.get('sidewall_oxide_left_A', 0):.0f}/{m.get('sidewall_oxide_right_A', 0):.0f}" if m else ""
        dose = (f"{m.get('inventory_5_rel_change', 0):+.2e} / {fm.get('inventory_5_rel_change', 0):+.2e}" if m else "")
        field = f"{m.get('field_oxide_A', 0):.0f} ({fm.get('field_oxide_A', 0):.0f})" if m else ""
        lines.append(f"| {r['name']} | {'PASS' if r['pass'] else 'FAIL'} | {run['wall_s']:.0f} | {run['core_calls']} | "
                     f"{run['regenerations']} | {field} | {ratio} | {side} | {dose} | "
                     f"{', '.join(failed) or ('' if run['complete'] else run.get('failure', ''))} |")
    n = len(results)
    ok = sum(r["pass"] for r in results)
    lines += ["", f"**{ok}/{n} passed.**  Field oxide in brackets: flat wafer, plain SUPREM (core off).  '*': 2D face "
              "(shorter than 4x its oxide thickness), ratio reported, not checked.", ""]
    text = "\n".join(lines)
    (out / "report.md").write_text(text)
    (out / "results.json").write_text(json.dumps(results, indent=1, default=str))
    return text, ok, n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", type=Path, help="native64 build linked with --geometry-core")
    ap.add_argument("--core-build", type=Path, help="core build directory (core_str, core_check_suprem)")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--matrix", action="store_true")
    ap.add_argument("--random", type=int, default=0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--locos", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--timeout", type=float, default=5400)
    ap.add_argument("--report", type=Path, help="re-evaluate the cases in this output directory")
    a = ap.parse_args()

    if a.report:
        a.out = a.report
        meta = json.loads((a.out / "suite.json").read_text())
        a.build, a.core_build = Path(meta["build"]), Path(meta["core_build"])
        results = []
        for name in meta["cases"]:
            d = a.out / name
            if not (d / "result.json").exists():
                continue
            r = json.loads((d / "result.json").read_text())
            r["run"] = json.loads((d / "core" / "run.json").read_text())
            if name.startswith("locos"):
                r["off_run"] = json.loads((d / "off" / "run.json").read_text())
                results.append(evaluate_locos(a, d, r))
            else:
                r["flat_run"] = json.loads((d / "flat" / "run.json").read_text())
                results.append(evaluate_poly(a, d, r))
        text, ok, n = report(results, a.out)
        print(text)
        return 0 if ok == n else 1

    if not (a.build and a.core_build and a.out):
        ap.error("--build, --core-build and --out are required")
    a.build, a.core_build, a.out = a.build.resolve(), a.core_build.resolve(), a.out.resolve()
    a.out.mkdir(parents=True, exist_ok=True)
    cases = []
    if a.matrix:
        cases += [(n, matrix_params(n)) for n in MATRIX]
    rng = random.Random(a.seed)
    cases += [(f"rand{a.seed}_{i:02d}", random_params(rng)) for i in range(a.random)]
    names = [n for n, _ in cases] + ([f"locos_{m}" for m in ("viscous", "compress")] if a.locos else [])
    (a.out / "suite.json").write_text(json.dumps(dict(build=str(a.build), core_build=str(a.core_build),
                                                      seed=a.seed, cases=names,
                                                      params={n: p for n, p in cases}), indent=1))
    results = []
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(poly_case, a, n, p): n for n, p in cases}
        if a.locos:
            if LOCOS_MASK.exists():
                for m in ("viscous", "compress"):
                    futs[ex.submit(locos_case, a, m)] = f"locos_{m}"
            else:
                print(f"LOCOS skipped: {LOCOS_MASK} missing", file=sys.stderr)
        for f in cf.as_completed(futs):
            r = f.result()
            results.append(r)
            print(f"{r['name']}: {'PASS' if r['pass'] else 'FAIL'} ({r['run']['wall_s']:.0f} s)", flush=True)
    order = {n: i for i, n in enumerate(names)}
    results.sort(key=lambda r: order.get(r["name"], 0))
    text, ok, n = report(results, a.out)
    print(text)
    return 0 if ok == n else 1


if __name__ == "__main__":
    sys.exit(main())
