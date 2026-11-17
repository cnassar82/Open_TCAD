#!/usr/bin/env python3
"""Run a SUPREM deck with one Gmsh remesh inserted before a chosen command."""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def split_deck(deck, marker):
    lines = deck.read_text(encoding="utf-8").splitlines(keepends=True)
    for index, line in enumerate(lines):
        if marker in line:
            return lines[:index], lines[index:]
    raise RuntimeError(f"could not find marker in deck: {marker!r}")


def y_lines_from_deck(deck):
    out = []
    pattern = re.compile(r"\bline\s+y\b", re.IGNORECASE)
    value = re.compile(r"\b(loc|spac|tag)\s*=\s*([^\s#]+)", re.IGNORECASE)
    for line in deck.read_text(encoding="utf-8").splitlines():
        if not pattern.search(line):
            continue
        vals = {key.lower(): raw for key, raw in value.findall(line)}
        if "loc" not in vals:
            continue
        try:
            loc = float(vals["loc"])
            spac = float(vals["spac"]) if "spac" in vals else None
        except ValueError:
            continue
        out.append({"loc": loc, "spac": spac, "tag": vals.get("tag")})
    out.sort(key=lambda item: item["loc"])
    return out


def deck_top_silicon_profile(deck):
    y_lines = y_lines_from_deck(deck)
    if not y_lines:
        return None, None
    surface = None
    for item in y_lines:
        if item.get("tag") == "surf":
            surface = item
            break
    if surface is None:
        surface = y_lines[0]
    depth = None
    for item in y_lines:
        if item["loc"] > surface["loc"]:
            depth = item["loc"] - surface["loc"]
            break
    return surface.get("spac"), depth


def deck_lateral_spacing(deck):
    pattern = re.compile(r"\bline\s+x\b", re.IGNORECASE)
    value = re.compile(r"\bspac\s*=\s*([^\s#]+)", re.IGNORECASE)
    spacings = []
    for line in deck.read_text(encoding="utf-8").splitlines():
        if not pattern.search(line):
            continue
        match = value.search(line)
        if not match:
            continue
        try:
            spacings.append(float(match.group(1)))
        except ValueError:
            pass
    return min(spacings) if spacings else None


def suprem_env(binary):
    env = os.environ.copy()
    data_dir = binary.parent.parent / "data"
    if data_dir.is_dir():
        env.setdefault("SUP4KEYFILE", str(data_dir / "suprem.uk"))
        env.setdefault("SUP4MODELRC", str(data_dir / "modelrc"))
        env.setdefault("SUP4IMPDATA", str(data_dir / "sup4gs.imp"))
    return env


def run_suprem(binary, deck_path, cwd, log_path):
    env = suprem_env(binary)
    with deck_path.open("rb") as stdin, log_path.open("wb") as log:
        proc = subprocess.run([str(binary)], stdin=stdin, stdout=log, stderr=subprocess.STDOUT, cwd=cwd, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"suprem4 failed with exit {proc.returncode}; log={log_path}")


def stage_deck_dir(deck, work_dir):
    for item in deck.parent.iterdir():
        if item.is_file():
            dst = work_dir / item.name
            if item.resolve() != dst.resolve():
                shutil.copy2(item, dst)


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("deck", type=Path)
    parser.add_argument("--before", required=True, help="insert remesh before the first deck line containing this text")
    parser.add_argument("--binary", type=Path, default=Path("/home/chris/Open_TCAD/suprem4/bin/suprem4"))
    parser.add_argument("--work-dir", type=Path, help="directory for generated decks, structures, and logs")
    parser.add_argument("--gmsh", default="gmsh")
    parser.add_argument("--gmsh-timeout", type=float, default=60.0)
    parser.add_argument("--mesh-size", type=float, default=1.0)
    parser.add_argument("--top-silicon-size", type=float, help="override detected deck surface y spacing")
    parser.add_argument("--bottom-silicon-size", type=float, default=2.5)
    parser.add_argument("--active-depth", type=float, default=None)
    parser.add_argument("--grade-power", type=float, default=1.0)
    parser.add_argument("--embed-top-grid", action="store_true", help="experimental: embed a deck-like top silicon grid")
    parser.add_argument("--preserve-top-silicon", action="store_true", help="also embed existing top-silicon points; can create a large mesh")
    parser.add_argument("--anisotropic", action="store_true", help="use experimental Gmsh BAMG anisotropic metric")
    parser.add_argument("--layered", action="store_true", help="use layered_remesh.py (deck-driven x/y profile) instead of the scalar/anisotropic gmsh_remesh_str.py path")
    parser.add_argument("--oxide-mode", choices=("boundary-only", "interface-layers", "preserve-interior"),
                        default="boundary-only",
                        help="layered mode: how to mesh oxide regions (see layered_remesh.py --help)")
    parser.add_argument("--conservation-scope", choices=("material", "bulk-only"), default="material",
                        help="layered mode: scope of dopant-integral conservation scaling")
    parser.add_argument("--oxide-interface-layers", type=int, default=None,
                        help="layered mode: number of interface layers when --oxide-mode=interface-layers")
    parser.add_argument("--oxide-interface-distance", type=float, default=None,
                        help="layered mode: perpendicular distance (um) for first interface layer")
    parser.add_argument("--target-material", type=int, action="append", dest="target_materials",
                        help="layered mode: material id to receive deck-driven interior points; defaults to silicon")
    parser.add_argument("--no-conserve-dopants", action="store_true")
    parser.add_argument("--keep", action="store_true", help="keep temporary work directory when --work-dir is omitted")
    args = parser.parse_args(argv)

    repo = Path(__file__).resolve().parents[1]
    if args.layered:
        remesh_tool = repo / "tools" / "layered_remesh.py"
    else:
        remesh_tool = repo / "tools" / "gmsh_remesh_str.py"
    work_dir = args.work_dir
    made_temp = False
    if work_dir is None:
        work_dir = Path(tempfile.mkdtemp(prefix="suprem4_gmsh_remesh_deck_"))
        made_temp = True
    work_dir.mkdir(parents=True, exist_ok=True)
    stage_deck_dir(args.deck, work_dir)

    before_lines, after_lines = split_deck(args.deck, args.before)
    deck_top_size, deck_active_depth = deck_top_silicon_profile(args.deck)
    deck_x_size = deck_lateral_spacing(args.deck)
    top_silicon_size = args.top_silicon_size
    if top_silicon_size is None:
        top_silicon_size = args.mesh_size
    active_depth = args.active_depth
    if active_depth is None:
        active_depth = deck_active_depth
    before_str = work_dir / "before_gmsh_remesh.str"
    remesh_str = work_dir / "after_gmsh_remesh.str"
    prefix_deck = work_dir / "prefix.in"
    suffix_deck = work_dir / "suffix.in"

    try:
        prefix_deck.write_text("".join(before_lines) + f"\nstructure outfile={before_str}\n", encoding="utf-8")
        suffix_deck.write_text(f"structure infile={remesh_str}\n" + "".join(after_lines), encoding="utf-8")

        run_suprem(args.binary, prefix_deck, work_dir, work_dir / "prefix.log")

        if args.layered:
            cmd = [
                sys.executable,
                str(remesh_tool),
                str(before_str),
                "--deck", str(args.deck),
                "--out-prefix", str(work_dir / "gmsh_remesh"),
                "--gmsh", args.gmsh,
                "--gmsh-timeout", str(args.gmsh_timeout),
                "--write-str", str(remesh_str),
                "--verbose",
                "--oxide-mode", args.oxide_mode,
                "--conservation-scope", args.conservation_scope,
            ]
            if args.oxide_interface_layers:
                cmd.extend(["--oxide-interface-layers", str(args.oxide_interface_layers)])
            if args.oxide_interface_distance is not None:
                cmd.extend(["--oxide-interface-distance", str(args.oxide_interface_distance)])
            if args.no_conserve_dopants:
                cmd.append("--no-conserve-dopants")
            for mat in (args.target_materials or []):
                cmd.extend(["--target-material", str(mat)])
        else:
            cmd = [
                sys.executable,
                str(remesh_tool),
                str(before_str),
                "--out-prefix",
                str(work_dir / "gmsh_remesh"),
                "--mesh-size",
                str(args.mesh_size),
                "--silicon-profile",
                "--top-mesh-size",
                str(top_silicon_size),
                "--bottom-mesh-size",
                str(args.bottom_silicon_size),
                "--grade-power",
                str(args.grade_power),
                "--gmsh",
                args.gmsh,
                "--gmsh-timeout",
                str(args.gmsh_timeout),
                "--write-str",
                str(remesh_str),
                "--verbose",
            ]
            if deck_x_size is not None:
                cmd.extend(["--x-mesh-size", str(deck_x_size)])
            cmd.extend(["--y-bottom-mesh-size", str(args.bottom_silicon_size)])
            if args.embed_top_grid:
                cmd.append("--embed-top-grid")
            if args.anisotropic:
                cmd.append("--anisotropic")
                cmd.append("--structure-top-spacing")
            if args.active_depth is not None:
                cmd.extend(["--active-depth", str(active_depth)])
            elif active_depth is not None:
                cmd.extend(["--active-depth", str(active_depth)])
            if args.preserve_top_silicon:
                cmd.append("--preserve-top-silicon")
                if active_depth is not None:
                    cmd.extend(["--preserve-depth", str(active_depth)])
            if args.no_conserve_dopants:
                cmd.append("--no-conserve-dopants")
        remesh_log = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (work_dir / "gmsh_remesh.log").write_text(remesh_log.stdout, encoding="utf-8")
        if remesh_log.returncode != 0:
            raise RuntimeError(f"gmsh remesh failed with exit {remesh_log.returncode}; log={work_dir / 'gmsh_remesh.log'}")

        run_suprem(args.binary, suffix_deck, work_dir, work_dir / "suffix.log")

        print(f"work_dir={work_dir}")
        print(f"before={before_str}")
        print(f"after={remesh_str}")
        print(f"prefix_log={work_dir / 'prefix.log'}")
        print(f"remesh_log={work_dir / 'gmsh_remesh.log'}")
        print(f"suffix_log={work_dir / 'suffix.log'}")
        return 0
    finally:
        if made_temp and not args.keep:
            shutil.rmtree(work_dir)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
