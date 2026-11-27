"""Deck rendering, local simulator jobs, and result discovery."""

from __future__ import annotations

import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time
from datetime import datetime

from electrical import collect_measurements
from deck_report import NOTE, REPORT_FILES

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "tmp" / "workbench"
GENIUS_RUNTIME = ROOT / "suprem4/experiments/rit_electrical_20261116/genius_runtime"
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
NAME = r"[A-Za-z_][A-Za-z_0-9]*"
FACTOR = re.compile(r"\{\{(" + NAME + r")=(" + NUMBER + r")\}\}")
FAILURE = re.compile(
    r"errors detected on command input|suprem4 panic:|No mesh defined|"
    r"not a valid command|\bERROR:|bad initial grid|Newton convergence failed|"
    r"runtime error:|AddressSanitizer:|Giving up|Failed in the first step|:\s*command not found\s*$", re.I
)
MPI_WARNING = re.compile(r"\bUCX\s+WARN\b")
ARTIFACTS = {".str", ".vtu", ".csv", ".dat", ".tsv", ".cgns", ".log", ".sol", ".json", ".txt",
             ".png", ".jpg", ".jpeg"}
MAX_CASES = 128


def infer_engine(path: Path, deck: str) -> str:
    genius_cards = r"^\s*(?:GLOBAL|X\.MESH|Y\.MESH|BOUNDARY|SOLVE|IMPORT|EXPORT|PROFILE)\b"
    return "Genius" if path.suffix.lower() == ".inp" or re.search(genius_cards, deck, re.M) else "SUPREM4"


def local_tif_inputs(path: Path, deck: str) -> list[Path]:
    """Stage local Genius meshes referenced by a one-line IMPORT card."""
    inputs = []
    for line in deck.splitlines():
        if not re.match(r"\s*IMPORT\s", line, re.I):
            continue
        matched = re.search(r"""\bTIFFile\s*=\s*(?:"([^"]+)"|'([^']+)'|([^\s#]+))""", line, re.I)
        if not matched:
            continue
        relative = Path(next(value for value in matched.groups() if value))
        candidate = path.parent / relative
        if not relative.is_absolute() and ".." not in relative.parts and candidate.is_file() and candidate not in inputs:
            inputs.append(candidate)
    return inputs


def factors(deck: str) -> dict[str, str]:
    deck = "".join(line for line in deck.splitlines(keepends=True) if not line.lstrip().startswith(NOTE))
    defaults = {}
    for match in FACTOR.finditer(deck):
        name, value = match.groups()
        if not math.isfinite(float(value)):
            raise ValueError(f"{name} must have a finite default.")
        if name in defaults and defaults[name] != value:
            raise ValueError(f"{name} has conflicting defaults.")
        defaults[name] = value
    remaining = FACTOR.sub("", deck)
    if "{{" in remaining or "}}" in remaining:
        raise ValueError("Use {{name=number}} for each sweep factor.")
    return defaults


def design(deck: str, levels: dict[str, str]) -> list[dict[str, str]]:
    defaults = factors(deck)
    columns = []
    for name, default in defaults.items():
        values = dict.fromkeys(re.split(r"[\s,]+", levels.get(name, default).strip()))
        if any(not re.fullmatch(NUMBER, value) or not math.isfinite(float(value)) for value in values):
            raise ValueError(f"Enter finite numbers for {name}, separated by commas or spaces.")
        columns.append(list(values))
    size = math.prod(map(len, columns))
    if size > MAX_CASES:
        raise ValueError(f"This design has {size} cases. The limit is {MAX_CASES} per batch.")
    return [dict(zip(defaults, values)) for values in itertools.product(*columns)]


def render_deck(deck, values):
    return "".join(line if line.lstrip().startswith(NOTE) else
                   FACTOR.sub(lambda match: values[match[1]], line)
                   for line in deck.splitlines(keepends=True))


def simulator(engine: str) -> tuple[Path, Path | None]:
    if engine == "SUPREM4":
        return Path(os.environ.get("TCAD_SUPREM4", ROOT / "suprem4/bin/suprem4")).absolute(), None
    runtime = GENIUS_RUNTIME if (GENIUS_RUNTIME / "bin/genius.LINUX").is_file() else ROOT / "Genius-TCAD-Open"
    binary = Path(os.environ.get("TCAD_GENIUS", runtime / "bin/genius.LINUX")).absolute()
    home = Path(os.environ.get("TCAD_GENIUS_DIR", binary.parent.parent)).absolute()
    return binary, home


def save_batch(folder: Path, batch: dict) -> None:
    temporary = folder / "batch.tmp"
    temporary.write_text(json.dumps(batch, indent=2) + "\n", encoding="utf-8")
    temporary.replace(folder / "batch.json")


def prepare_batch(path: Path, deck: str, engine: str, cases: list[dict], inputs: list[Path],
                  *, destination: Path | None = None) -> tuple[Path, dict]:
    binary, home = simulator(engine)
    if not os.access(binary, os.X_OK):
        raise ValueError(f"Simulator executable is unavailable at {binary}")
    inputs = [(file, file.relative_to(path.parent)) for file in inputs]
    if any(relative.parts[0] in {path.name, "run.log"} for _, relative in inputs) or path.name == "run.log":
        raise ValueError("The deck and input files must have distinct names from the run log and deck snapshot.")
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    folder = destination or RUNS / f"{stamp}-{path.stem}"
    folder.mkdir(parents=True)
    batch = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "source": str(path), "engine": engine, "binary": str(binary),
        "genius_home": str(home) if home else None,
        "inputs": [str(relative) for _, relative in inputs],
        "deck": path.name, "cases": [],
    }
    try:
        (folder / "template.deck").write_text(deck, encoding="utf-8")
        for index, values in enumerate(cases, start=1):
            name = f"case-{index:03}"
            case_dir = folder / name
            case_dir.mkdir()
            rendered = render_deck(deck, values)
            (case_dir / path.name).write_text(rendered.rstrip("\n") + "\n", encoding="utf-8")
            for file, relative in inputs:
                target = case_dir / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(file, target)
            batch["cases"].append({"name": name, "values": values, "status": "Queued"})
        save_batch(folder, batch)
    except (OSError, KeyError):
        shutil.rmtree(folder)
        raise
    return folder, batch


def artifact_paths(case_dir: Path, batch: dict) -> list[Path]:
    excluded = {batch["deck"], "run.log", *batch["inputs"]}
    return sorted(path.relative_to(case_dir) for path in case_dir.rglob("*")
                  if path.is_file() and (path.suffix.lower() in ARTIFACTS or path.name in REPORT_FILES)
                  and str(path.relative_to(case_dir)) not in excluded)


def has_failure(log: Path) -> bool:
    with log.open(encoding="utf-8", errors="replace") as stream:
        return any(FAILURE.search(line) and not MPI_WARNING.search(line) for line in stream)


class Job:
    """Run one case at a time, with status persisted after each transition."""

    def __init__(self, folder: Path, batch: dict):
        self.folder = folder
        self.batch = batch
        self.index = 0
        self.process = None
        self.started = None
        self.stop_time = None

    @property
    def done(self):
        return self.index == len(self.batch["cases"])

    def tick(self):
        if self.done:
            return False
        case = self.batch["cases"][self.index]
        directory = self.folder / case["name"]
        log = directory / "run.log"

        if self.process is not None:
            code = self.process.poll()
            if code is None:
                if self.stop_time is not None and time.monotonic() - self.stop_time >= 2:
                    self.terminate(signal.SIGKILL)
                return False
            case.update(returncode=code, elapsed=round(time.monotonic() - self.started, 2))
            if self.stop_time is not None:
                case["status"] = "Canceled"
            else:
                try:
                    case["status"] = "Failed" if code or has_failure(log) else "Completed"
                except OSError as error:
                    case.update(status="Failed", error=str(error))
                if self.batch["engine"] == "Genius":
                    try:
                        warnings = collect_measurements(directory, directory / self.batch["deck"], case["status"])
                        if warnings:
                            case["measurement_warnings"] = warnings
                    except OSError as error:
                        case["measurement_warnings"] = [str(error)]
            self.process = None
            self.index = len(self.batch["cases"]) if self.stop_time is not None else self.index + 1
        else:
            binary = self.batch["binary"]
            stdbuf = shutil.which("stdbuf")
            command = [stdbuf, "-oL", "-eL", binary] if stdbuf else [binary]
            env = os.environ.copy()
            if self.batch["engine"] == "Genius":
                command.extend(["-i", self.batch["deck"]])
                env["GENIUS_DIR"] = self.batch["genius_home"]
                env.pop("DISPLAY", None)
                env.setdefault("HWLOC_COMPONENTS", "-gl")
                env.setdefault("OMP_NUM_THREADS", "1")
                env["OMPI_MCA_orte_tmpdir_base"] = str(directory)
            try:
                with log.open("wb") as output, (directory / self.batch["deck"]).open("rb") as deck:
                    self.process = subprocess.Popen(
                        command, cwd=directory, env=env, start_new_session=True,
                        stdin=deck if self.batch["engine"] == "SUPREM4" else subprocess.DEVNULL,
                        stdout=output, stderr=subprocess.STDOUT,
                    )
                case["status"] = "Running"
                self.started = time.monotonic()
            except OSError as error:
                case.update(status="Failed", error=str(error))
                self.index += 1
        save_batch(self.folder, self.batch)
        return True

    def terminate(self, sig):
        try:
            os.killpg(self.process.pid, sig)
        except ProcessLookupError:
            pass

    def cancel(self):
        if self.done or self.stop_time is not None:
            return
        self.stop_time = time.monotonic()
        if self.process is not None:
            self.terminate(signal.SIGTERM)
        for case in self.batch["cases"][self.index + (self.process is not None):]:
            case["status"] = "Canceled"
        if self.process is None:
            self.index = len(self.batch["cases"])
        save_batch(self.folder, self.batch)
