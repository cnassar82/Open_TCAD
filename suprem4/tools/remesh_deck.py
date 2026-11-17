#!/usr/bin/env python3
"""Run a SUPREM deck, performing a remesh wherever a `# remesh` comment line
appears. Splits the deck at each marker, runs SUPREM up to the marker,
invokes the chosen remesh engine, then continues with the remeshed structure.

Deck syntax (anywhere in the .in file, treated as a comment by SUPREM):

    # remesh                                    # use defaults (gmsh-based)
    # remesh oxide-mode=preserve-interior       # override per insertion
    # remesh oxide-mode=interface-layers conservation-scope=bulk-only

  --- OR, oxide-only densification via the structured column triangulator:

    # remesh mode=refine-oxide layers=8         # 8 interior y layers per column

The refine-oxide mode bypasses gmsh entirely (which hangs on thin oxide
regions with concentrated interior points). It places N interior y nodes at
each existing boundary x position and stitches them with a structured
column-by-column triangulation. Designed for the implant-after-strip case
where SUPREM otherwise has only 1 oxide node line per x and loses Pearson
dose. Use this BEFORE the implant; downstream oxidation/etch tolerate it.

Defaults: --oxide-mode interface-layers --conservation-scope bulk-only
(when no mode= is given).
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


MARKER_RE = re.compile(r"^\s*#\s*remesh\b(.*)$", re.IGNORECASE)
DEFAULT_REMESH_ARGS = {
    "oxide-mode": "interface-layers",
    "conservation-scope": "bulk-only",
}


def split_on_markers(text):
    """Return [(chunk_text, marker_args_dict_or_None), ...]."""
    chunks = []
    current = []
    for line in text.splitlines(keepends=True):
        m = MARKER_RE.match(line)
        if m:
            args = dict(DEFAULT_REMESH_ARGS)
            for kv in m.group(1).split():
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    args[k.strip()] = v.strip()
            chunks.append(("".join(current), args))
            current = []
        else:
            current.append(line)
    chunks.append(("".join(current), None))
    return chunks


def suprem_env(binary):
    env = os.environ.copy()
    data_dir = binary.parent.parent / "data"
    if data_dir.is_dir():
        env.setdefault("SUP4KEYFILE", str(data_dir / "suprem.uk"))
        env.setdefault("SUP4MODELRC", str(data_dir / "modelrc"))
        env.setdefault("SUP4IMPDATA", str(data_dir / "sup4gs.imp"))
    return env


def _tee_subprocess(cmd, log_path, stdin=None, cwd=None, env=None):
    """Run cmd, write stdout (merged stderr) to both terminal and log_path in real time.
    Returns the exit code.
    """
    with log_path.open("wb") as log:
        proc = subprocess.Popen(cmd, stdin=stdin, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, cwd=cwd, env=env, bufsize=0)
        try:
            for line in proc.stdout:
                sys.stdout.buffer.write(line)
                sys.stdout.buffer.flush()
                log.write(line)
        finally:
            proc.stdout.close()
        return proc.wait()


def run_suprem(binary, deck_path, cwd, log_path):
    env = suprem_env(binary)
    with deck_path.open("rb") as stdin:
        rc = _tee_subprocess([str(binary)], log_path, stdin=stdin, cwd=cwd, env=env)
    if rc != 0:
        raise RuntimeError(f"suprem4 failed with exit {rc}; log={log_path}")


def stage_deck_dir(deck, work_dir):
    for item in deck.parent.iterdir():
        if item.is_file():
            dst = work_dir / item.name
            if item.resolve() != dst.resolve():
                shutil.copy2(item, dst)


def run_layered_remesh(tool, input_str, deck, output_str, work_dir, log_path, marker_args):
    cmd = [
        sys.executable, str(tool),
        str(input_str),
        "--deck", str(deck),
        "--out-prefix", str(work_dir / output_str.stem),
        "--write-str", str(output_str),
        "--verbose",
    ]
    for k, v in marker_args.items():
        if v is None or v == "":
            continue
        cmd.extend([f"--{k}", str(v)])
    rc = _tee_subprocess(cmd, log_path)
    if rc != 0:
        raise RuntimeError(f"layered_remesh failed with exit {rc}; log={log_path}")


def run_refine_oxide(tool, input_str, output_str, log_path, marker_args):
    """Dispatch to refine_oxide.py (gmsh-free structured oxide refinement).
    Recognised marker keys (besides mode=refine-oxide): layers, boundary-buffer,
    target-material (may repeat as comma-separated list).
    """
    cmd = [sys.executable, str(tool), str(input_str), str(output_str)]
    if "layers" in marker_args:
        cmd.extend(["--layers", str(marker_args["layers"])])
    if "boundary-buffer" in marker_args:
        cmd.extend(["--boundary-buffer", str(marker_args["boundary-buffer"])])
    if "target-material" in marker_args:
        for m in str(marker_args["target-material"]).split(","):
            m = m.strip()
            if m:
                cmd.extend(["--target-material", m])
    rc = _tee_subprocess(cmd, log_path)
    if rc != 0:
        raise RuntimeError(f"refine_oxide failed with exit {rc}; log={log_path}")


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("deck", type=Path, help="SUPREM deck (.in) with optional `# remesh` markers")
    parser.add_argument("--binary", type=Path, default=Path("/home/chris/Open_TCAD/suprem4/bin/suprem4"))
    parser.add_argument("--work-dir", type=Path,
                        help="directory for generated decks, structures, and logs (defaults to cwd)")
    args = parser.parse_args(argv)

    repo = Path(__file__).resolve().parents[1]
    remesh_tool = repo / "tools" / "layered_remesh.py"
    refine_oxide_tool = repo / "tools" / "refine_oxide.py"

    work_dir = args.work_dir if args.work_dir is not None else Path.cwd()
    work_dir.mkdir(parents=True, exist_ok=True)
    stage_deck_dir(args.deck, work_dir)

    deck_text = args.deck.read_text(encoding="utf-8")
    chunks = split_on_markers(deck_text)
    n_marks = sum(1 for _, m in chunks if m is not None)
    if n_marks == 0:
        print("no `# remesh` markers in deck; running SUPREM straight through")
        run_suprem(args.binary, args.deck, work_dir, work_dir / "suprem.log")
        return 0

    try:
        print(f"workdir: {work_dir}")
        print(f"found {n_marks} `# remesh` marker(s); deck splits into {len(chunks)} chunks")

        state_str = None
        for i, (chunk_text, marker_args) in enumerate(chunks):
            chunk_deck = work_dir / f"chunk_{i:02d}.in"
            chunk_log = work_dir / f"chunk_{i:02d}.log"
            if state_str is not None:
                preamble = f"structure infile={state_str}\n"
            else:
                preamble = ""
            if marker_args is not None:
                # Not the last chunk: append a structure save so the next remesh
                # has a file to read.
                save_str = work_dir / f"state_{i:02d}.str"
                chunk_deck.write_text(preamble + chunk_text +
                                      f"\nstructure outfile={save_str}\n",
                                      encoding="utf-8")
                run_suprem(args.binary, chunk_deck, work_dir, chunk_log)
                remeshed = work_dir / f"remeshed_{i:02d}.str"
                remesh_log = work_dir / f"remesh_{i:02d}.log"
                mode = marker_args.get("mode", "").lower()
                if mode == "refine-oxide":
                    # Strip 'mode' from args (consumed); rest go to refine_oxide
                    inner = {k: v for k, v in marker_args.items()
                             if k not in {"mode", "oxide-mode", "conservation-scope"}}
                    run_refine_oxide(refine_oxide_tool, save_str, remeshed,
                                     remesh_log, inner)
                    engine = "refine-oxide"
                else:
                    run_layered_remesh(remesh_tool, save_str, args.deck, remeshed,
                                       work_dir, remesh_log, marker_args)
                    engine = "layered-remesh"
                state_str = remeshed
                print(f"  chunk {i:02d}: ran SUPREM ({chunk_log.name}) -> "
                      f"{engine} ({remesh_log.name}, args={marker_args})")
            else:
                chunk_deck.write_text(preamble + chunk_text, encoding="utf-8")
                run_suprem(args.binary, chunk_deck, work_dir, chunk_log)
                print(f"  chunk {i:02d}: ran final SUPREM ({chunk_log.name})")

        print(f"done. outputs in {work_dir}")
        return 0
    except Exception:
        raise


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
