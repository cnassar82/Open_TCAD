"""Saved implant-dose studies with a process-to-device run for each case."""

import csv
from datetime import datetime
from hashlib import sha256
from random import Random
import json
from math import hypot, isfinite
from pathlib import Path
import re
import shutil
from tempfile import TemporaryDirectory

from deck_report import REPORT_FILES
from contact_export import build_edges, contact_mesh, export_device
from electrical import MEASUREMENTS, PREFIX, measurement_cards
from make_genius_contacts import ContactSpec, parse_suprem
from measurements import evaluate_response, response_recipe
from runner import (MAX_CASES, Job, artifact_paths, design, factors, prepare_batch, render_deck,
                    save_batch, simulator)

MARKER = "@checkpoint@"
INPUT = "study-start.str"
OUTPUT = "study-device.str"
RESULTS = "vt-versus-dose.csv"
FLOW_RESULTS = "doe-responses.csv"
SMALL = 1e-7


def point_in_triangle(nodes, tri, x, y):
    a, b, c = (nodes[i] for i in tri.nodes)
    v1 = b.x - a.x, b.y - a.y
    v2 = c.x - a.x, c.y - a.y
    ux, uy = x - a.x, y - a.y
    det = v1[0] * v2[1] - v1[1] * v2[0]
    if abs(det) < 1e-25:
        return False
    u = (ux * v2[1] - uy * v2[0]) / det
    v = (v1[0] * uy - v1[1] * ux) / det
    return u >= -1e-9 and v >= -1e-9 and 1 - u - v >= -1e-9


def make_contact_map(nodes, triangles, regions, contacts):
    contact_mesh(triangles, regions, contacts)
    edges, numbers, _ = build_edges(triangles)
    by_number = {numbers[e]: e for e in edges}
    maps = []
    for contact in contacts:
        item = dict(name=contact.name, kind=contact.kind, target=contact.target,
                    work_function=contact.work_function)
        if contact.target == "region":
            same = [t for t in triangles if t.region == contact.region]

            def area(t):
                a, b, c = (nodes[i] for i in t.nodes)
                return abs((b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x))

            # Large triangles give anchors that are less sensitive to local remeshing.
            anchor_triangles = sorted(same, key=area, reverse=True)[:8]
            item["material"] = regions[contact.region].mat_id
            item["anchors"] = [
                [sum(getattr(nodes[n], axis) for n in tri.nodes) / 3 for axis in ("x", "y")]
                for tri in anchor_triangles
            ]
        else:
            item["segments"] = [[[nodes[n].x, nodes[n].y] for n in by_number[e]] for e in contact.edges]
        maps.append(item)
    return {"version": 1, "contacts": maps}


def segment_distance(x, y, a, b):
    dx, dy = b[0]-a[0], b[1]-a[1]
    scale = dx*dx + dy*dy
    t = min(1.0, max(0.0, ((x-a[0])*dx + (y-a[1])*dy) / scale)) if scale else 0.0
    return hypot(x-a[0]-t*dx, y-a[1]-t*dy)


def resolve_contacts(nodes, triangles, regions, recipe):
    ids = recipe["contacts"]
    edges, numbers, owners = build_edges(triangles)
    result = []
    for item in ids:
        spec = ContactSpec(item["name"], item["kind"], item["target"], work_function=float(item["work_function"]))
        if spec.target == "region":
            candidates = [t for t in triangles if regions[t.region].mat_id == item["material"]]
            hits = {t.region for x, y in item["anchors"] for t in candidates
                    if point_in_triangle(nodes, t, x, y)}
            if len(hits) != 1:
                raise ValueError(f"Contact {spec.name} matched {len(hits)} regions on this mesh.")
            spec.region = hits.pop()
        else:
            spec.edges = []
            matched = []
            for edge in edges:
                if len(owners[edge]) != 1 or regions[next(iter(owners[edge]))].mat_id != 3:
                    continue
                a, b = (nodes[n] for n in edge)
                mx, my = (a.x + b.x) / 2, (a.y + b.y) / 2
                if all(any(segment_distance(x, y, start, end) <= SMALL for start, end in item["segments"])
                       for x, y in ((a.x, a.y), (mx, my), (b.x, b.y))):
                    spec.edges.append(numbers[edge])
                    matched.append(((a.x, a.y), (b.x, b.y)))
            if not matched or any(
                not any(segment_distance((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, left, right) <= SMALL
                        for left, right in matched)
                for a, b in item["segments"]
            ):
                raise ValueError(f"Could not match the full exposed silicon contact {spec.name}.")
        result.append(spec)
    return result


def test_entries(project):
    if project["version"] == 1:
        return [{"kind": MEASUREMENTS[1], "values": project["test"]}]
    return project["tests"]


def factor_source(project):
    return "\n".join([project["process"], project.get("device_cards", ""),
                      *(value for test in test_entries(project) for value in test["values"].values())])


def flow_cases(project):
    source = factor_source(project)
    defaults = factors(source)
    unknown = project["levels"].keys() - defaults.keys()
    if unknown:
        raise ValueError("Mark these factors in the process or device settings: " + ", ".join(sorted(unknown)))
    if project.get("design", "Full factorial") == "Full factorial":
        return design(source, project["levels"])
    if project["design"] != "Latin hypercube":
        raise ValueError("Choose full factorial or Latin hypercube.")
    count = int(project["samples"])
    if not 2 <= count <= MAX_CASES:
        raise ValueError(f"Enter 2–{MAX_CASES} Latin hypercube cases.")
    random = Random(int(project.get("seed", 0)))
    columns = []
    varied = False
    for name, default in defaults.items():
        tokens = re.split(r"[\s,]+", project["levels"].get(name, default).strip())
        limits = sorted(set(map(float, tokens)))
        if not 1 <= len(limits) <= 2 or not all(map(isfinite, limits)):
            raise ValueError(f"Enter a fixed value or two finite endpoints for {name}.")
        if len(limits) == 1:
            columns.append([format(limits[0], ".17g")] * count)
            continue
        varied = True
        lo, hi = limits
        samples = []
        for i in range(count):
            t = (i + random.random()) / count
            samples.append(format((1 - t) * lo + t * hi, ".17g"))
        random.shuffle(samples)
        columns.append(samples)
    if not varied:
        raise ValueError("Give at least one factor two range endpoints for Latin hypercube sampling.")
    return [dict(zip(defaults, row)) for row in zip(*columns)]


def validate_flow_tests(project, cases):
    entries = test_entries(project)
    if not 1 <= len(entries) <= 12:
        raise ValueError("Add between 1 and 12 electrical tests.")
    prefixes = [test["values"]["prefix"] for test in entries]
    if any(not PREFIX.fullmatch(prefix) for prefix in prefixes) or len(set(map(str.casefold, prefixes))) != len(prefixes):
        raise ValueError("Use distinct, fixed output prefixes for the electrical tests.")
    names = {item["name"].casefold() for item in project["contact_map"]["contacts"]}
    for case in cases:
        for entry in entries:
            if entry["kind"] not in MEASUREMENTS:
                raise ValueError("Choose a supported electrical test.")
            values = {name: render_deck(value, case) for name, value in entry["values"].items()}
            for role in ("drain", "gate", "source"):
                if values[role].casefold() not in names:
                    raise ValueError(f"The {role} electrode in {values['prefix']} is missing from the contact map.")
            measurement_cards(entry["kind"], values)


def validate_project(project, *, cases=None):
    if project.get("version") not in {1, 2}:
        raise ValueError("Choose a process-to-device study project.")
    start = project.get("checkpoint", "")
    cards = project["process"]
    lines = cards.splitlines()
    markers = [index for index, line in enumerate(lines) if line.strip() == MARKER]
    if start:
        checkpoint = Path(start)
        if checkpoint.suffix.lower() != ".str" or not checkpoint.is_file():
            raise ValueError("Choose an existing SUPREM4 checkpoint (.str).")
        if len(markers) != 1:
            raise ValueError(f"Use one {MARKER} line where SUPREM4 should load the starting structure.")
    elif project["version"] == 1 or markers or not re.search(r"^\s*(?:initialize|init)\b", cards, re.I | re.M):
        raise ValueError(f"Choose a starting checkpoint with {MARKER}, or supply a deck with initialize.")
    if (re.search(r"^\s*(?:quit|exit)\b|^\s*structure\s+[^#\n]*\binfile\s*=", cards, re.I | re.M)
            or (start and re.search(r"^\s*(?:initialize|init)\b", cards, re.I | re.M))):
        raise ValueError("The runner manages checkpoint loading and the final quit. Remove conflicting commands.")
    if re.search(r"^\s*structure\s+[^#\n]*\boutfile\s*=\s*['\"]?(?:study-start|study-device)\.str\b",
                 cards, re.I | re.M):
        raise ValueError("Choose checkpoint output names other than study-start.str and study-device.str.")
    if cases is None:
        if project["version"] == 1:
            if set(factors(cards)) != {"dose"} or not any(
                re.match(r"\s*implant\s+[^#]*\bdose\s*=\s*\{\{dose=", line, re.I)
                for line in lines[markers[0] + 1:]
            ):
                raise ValueError("Use {{dose=number}} on the implant card after the checkpoint.")
            cases = design(cards, {"dose": project["levels"]})
            doses = [float(case["dose"]) for case in cases]
            if len(doses) < 2 or len(set(doses)) != len(doses) or min(doses) <= 0:
                raise ValueError("Enter at least two distinct, positive doses.")
        else:
            cases = flow_cases(project)
    contacts = project["contact_map"]
    if contacts.get("version") != 1 or not contacts.get("contacts"):
        raise ValueError("Choose a saved contact map.")
    if re.search(r"^\s*(?:IMPORT|SOLVE|END|EXPORT)\b", project.get("device_cards", ""), re.I | re.M):
        raise ValueError("Genius options are inserted before equilibrium. The study supplies import, solve, and export cards.")
    validate_flow_tests(project, cases)
    responses = [response_recipe(entry)[0].casefold() for entry in project.get("measurements", [])]
    if len(set(responses)) != len(responses):
        raise ValueError("Use distinct names for DOE measurement recipes.")
    return cases


def start_study(project, project_path, runs, *, retry=None):
    selected, previous = retry if retry else ([], None)
    cases = validate_project(project, cases=[case["values"] for case in selected] if retry else None)
    if not simulator("Genius")[0].is_file():
        raise ValueError("The configured Genius executable is missing.")
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    title = project_path.stem.removesuffix(".study")
    target = runs / f"{stamp}-{title}-study"
    cards = "\n".join(f"structure infile={INPUT}" if line.strip() == MARKER else line
                      for line in project["process"].splitlines())
    cards += f"\nstructure outfile={OUTPUT}\nquit\n"
    if project["version"] == 2:
        cards = "set echo\n" + cards
    with TemporaryDirectory(prefix="study-inputs-") as temp:
        staging = Path(temp)
        inputs = []
        if project["checkpoint"]:
            shutil.copy2(project["checkpoint"], staging / INPUT)
            inputs.append(staging / INPUT)
        for entry in project.get("inputs", []):
            relative = Path(entry["relative"])
            if relative.is_absolute() or ".." in relative.parts or relative in {Path(INPUT), Path(OUTPUT)}:
                raise ValueError(f"Invalid continuation input path {relative}.")
            dest = staging / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(entry["source"], dest)
            inputs.append(dest)
        folder, batch = prepare_batch(staging / "process.in", cards, "SUPREM4", cases, inputs,
                                      destination=target)
    try:
        (folder / "study-project.json").write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
        flow = project["version"] == 2
        batch.update(study=True, flow=flow, source=str(project_path), title=title)
        if retry:
            batch["retry_of"] = previous.name
        for index, case in enumerate(batch["cases"]):
            case.update(phase="process", vt_status="Waiting" if any(
                entry["kind"] == MEASUREMENTS[1] for entry in test_entries(project)) else "")
            if flow:
                case["responses"] = {}
            if retry:
                case["retry_of"] = f"{previous.name}/{selected[index]['name']}"
        save_batch(folder, batch)
        job = StudyJob(folder, batch, project)
        if retry:
            for case, old in zip(batch["cases"], selected):
                old_dir = previous / old["name"]
                if old.get("process_completed") and (old_dir / OUTPUT).is_file():
                    job.process_sources[job.process_key(case)] = old_dir
        return job
    except OSError:
        shutil.rmtree(folder)
        raise


class StudyJob:
    """Sequence the process, contact handoff, and independent Genius test jobs."""

    def __init__(self, folder, batch, project):
        self.folder, self.batch, self.project = folder, batch, project
        self.process_job = Job(folder, batch)
        self.process_sources = {}
        self.device_job = None
        self.device_case = None
        self.device_queue = []
        self.canceling = False

    @property
    def done(self):
        return self.process_job.done and self.device_case is None

    def process_key(self, case):
        deck = self.folder / case["name"] / self.batch["deck"]
        return sha256(deck.read_bytes()).digest()

    def reuse_process(self, case, source):
        directory = self.folder / case["name"]
        paths = [path for path in artifact_paths(source, self.batch)
                 if path.parts[0] != "electrical" and path.name not in {"handoff.log", "doe-measurements.json", *REPORT_FILES}
                 and not path.parts[0].startswith("study-device-device-")]
        for relative in [*paths, Path("run.log")]:
            origin = source / relative
            if not origin.is_file():
                continue
            target = directory / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(origin, target)
        case.update(status="Completed", elapsed=0,
                    process_reused_from=f"{source.parent.name}/{source.name}")

    def setup_device(self, case):
        current = self.folder / case["name"]
        checkpoint = current / OUTPUT
        mesh = parse_suprem(checkpoint, auto_name_electrodes=False)
        contacts = resolve_contacts(*mesh[:3], self.project["contact_map"])
        deck, tif = export_device(current, checkpoint, *mesh, contacts)
        template = deck.read_text(encoding="utf-8")
        head, separator, tail = template.partition("SOLVE Type=EQUILIBRIUM")
        if not separator:
            raise ValueError("The exported device deck has no equilibrium stage.")
        options = render_deck(self.project.get("device_cards", ""), case["values"]).strip()
        queue = []
        for entry in test_entries(self.project):
            values = {name: render_deck(value, case["values"]) for name, value in entry["values"].items()}
            test = measurement_cards(entry["kind"], values)
            # Keep the equilibrium export a true zero-bias artifact for new flow projects.
            if self.batch["flow"]:
                body, _, _ = tail.rpartition("END")
                test = test.replace("SOLVE Type=Equilibrium", "SOLVE Type=Equilibrium" + body, 1)
                content = f"{head}{options}\n{test}\nEND\n"
                destination = current / "electrical" / values["prefix"]
            else:
                content = f"{head}{options}\n{test}\n{tail}"
                destination = current / "electrical"
            folder, job_batch = prepare_batch(deck, content, "Genius", [{}], [tif],
                                              destination=destination)
            queue.append((entry["kind"], values["prefix"], Job(folder, job_batch)))
        self.device_case = case
        self.device_queue = queue
        self.next_test()

    def next_test(self):
        self.active_kind, self.active_prefix, self.device_job = self.device_queue.pop(0)
        case = self.device_case
        relative = self.device_job.folder.relative_to(self.folder / case["name"])
        case.update(status="Running", phase="device", device_log=str(relative / "case-001/run.log"))
        if self.active_kind == MEASUREMENTS[1] and "VT_V" not in case:
            case["vt_status"] = "Solving"

    def finish_test(self):
        case, job = self.device_case, self.device_job
        device = job.batch["cases"][0]
        case["elapsed"] = round(case.get("elapsed", 0) + device.get("elapsed", 0), 2)
        result = {"status": device["status"], "run_status": device["status"]}
        path = job.folder / device["name"] / f"{self.active_prefix}.measurement.json"
        if device["status"] == "Completed":
            try:
                report = json.loads(path.read_text(encoding="utf-8"))
                result["curves"] = len(report["curves"])
                result["status"] = ("curves_ready" if result["curves"] == report["expected_curves"]
                                    else "missing_curves")
                result.update(report.get("result", {}))
            except (OSError, ValueError, KeyError, TypeError):
                result["status"] = "Measurement report unavailable"
        if self.active_kind == MEASUREMENTS[1] and (not self.batch["flow"] or "VT_V" not in case):
            case["vt_status"] = result["status"]
            voltage = result.get("VT_V")
            if isinstance(voltage, (float, int)) and isfinite(voltage):
                case["VT_V"] = voltage
        case.setdefault("responses", {})[self.active_prefix] = result
        if self.device_queue and not self.canceling:
            self.next_test()
            return
        statuses = [response["run_status"] for response in case["responses"].values()]
        case.update(status="Canceled" if self.canceling or "Canceled" in statuses else
                    "Failed" if "Failed" in statuses else "Completed", phase="done")
        if self.batch["flow"] and not self.canceling and self.project.get("measurements"):
            self.evaluate_responses(case)
        self.device_job = self.device_case = None
        self.device_queue = []

    def evaluate_responses(self, case):
        directory = self.folder / case["name"]
        results = {}
        for entry in self.project["measurements"]:
            try:
                name, result = evaluate_response(entry, directory)
                results[name] = result
            except (OSError, ValueError, TypeError, KeyError, ImportError, RuntimeError) as error:
                results[entry["name"]] = {"status": str(error)}
        case["measurements"] = {
            name: {key: value for key, value in item.items() if key in {"value", "status"}}
            for name, item in results.items()
        }
        try:
            (directory / "doe-measurements.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n",
                                                             encoding="utf-8")
        except OSError as error:
            case["measurement_warnings"] = [str(error)]

    def write_results(self):
        flow = self.batch["flow"]
        target = self.folder / (FLOW_RESULTS if flow else RESULTS)
        temporary = target.with_suffix(".tmp")
        entries = test_entries(self.project)
        metrics = {MEASUREMENTS[0]: "curves", MEASUREMENTS[1]: "VT_V", MEASUREMENTS[2]: "BVDSS_V"}
        columns = [(entry["values"]["prefix"], metrics[entry["kind"]]) for entry in entries]
        measurements = [entry["name"] for entry in self.project.get("measurements", [])]
        names = list(factors(factor_source(self.project))) if flow else ["dose"]
        with temporary.open("w", newline="", encoding="utf-8") as output:
            table = csv.writer(output)
            if flow:
                responses = [heading for prefix, metric in columns
                             for heading in (f"{prefix}.{metric}", f"{prefix}.status")]
                measurement_columns = [heading for name in measurements
                                       for heading in (f"measure.{name}.value", f"measure.{name}.status")]
                table.writerow((*[f"factor.{name}" for name in names], "case", "run_status", "phase",
                                "process_reused_from", *responses, *measurement_columns))
                for case in self.batch["cases"]:
                    cells = []
                    for prefix, metric in columns:
                        result = case["responses"].get(prefix, {})
                        cells.extend((result.get(metric, ""), result.get("status", "")))
                    for name in measurements:
                        result = case.get("measurements", {}).get(name, {})
                        cells.extend((result.get("value", ""), result.get("status", "")))
                    table.writerow((*[case["values"][name] for name in names], case["name"], case["status"],
                                    case["phase"], case.get("process_reused_from", ""), *cells))
            else:
                table.writerow(("dose_cm-2", "VT_V", "case", "run_status", "VT_status"))
                for case in self.batch["cases"]:
                    table.writerow((case["values"]["dose"], case.get("VT_V", ""), case["name"],
                                    case["status"], case["vt_status"]))
        temporary.replace(target)
        save_batch(self.folder, self.batch)

    def tick(self):
        if self.done:
            return False
        if self.device_job is not None:
            changed = self.device_job.tick()
            if self.device_job.done:
                self.finish_test()
                self.write_results()
                return True
            return changed

        index = self.process_job.index
        case = self.batch["cases"][index]
        process_key = None
        if self.batch["flow"] and self.process_job.process is None:
            process_key = self.process_key(case)
        if process_key in self.process_sources:
            source = self.process_sources[process_key]
            try:
                self.reuse_process(case, source)
            except OSError as error:
                case.update(status="Failed", error=f"Could not reuse the process checkpoint. {error}")
            self.process_job.index += 1
            changed = True
        else:
            changed = self.process_job.tick()
        if self.process_job.index != index:
            if case["status"] == "Completed":
                current = self.folder / case["name"]
                try:
                    if not (current / OUTPUT).is_file():
                        raise ValueError("The process did not write its final checkpoint.")
                    case["process_completed"] = True
                    if self.batch["flow"]:
                        self.process_sources.setdefault(self.process_key(case), current)
                    self.setup_device(case)
                except (OSError, RuntimeError, ValueError, KeyError, IndexError) as error:
                    case.update(status="Failed", phase="contacts" if case.get("process_completed") else "process",
                                vt_status="Handoff failed" if case.get("process_completed") else "Process failed",
                                error=str(error))
                    (current / "handoff.log").write_text(str(error) + "\n", encoding="utf-8")
                    if case.get("process_completed") and self.project.get("measurements"):
                        self.evaluate_responses(case)
            else:
                case["vt_status"] = "Canceled" if self.canceling else "Process failed"
            self.write_results()
        return changed

    def cancel(self):
        if self.done or self.canceling:
            return
        self.canceling = True
        if self.device_job is not None:
            self.device_job.cancel()
        self.process_job.cancel()
        for case in self.batch["cases"]:
            if case["status"] == "Canceled":
                case["vt_status"] = "Canceled"
        self.write_results()


def retry_study(folder, batch, runs=None):
    if not batch.get("flow"):
        raise ValueError("Select a process-to-device DOE run.")
    selected = [case for case in batch["cases"] if case["status"] in {"Failed", "Canceled"}]
    if not selected:
        raise ValueError("This DOE has no failed or canceled cases to rerun.")
    project = json.loads((folder / "study-project.json").read_text(encoding="utf-8"))
    if project["checkpoint"]:
        project["checkpoint"] = str(folder / batch["cases"][0]["name"] / INPUT)
    for entry in project.get("inputs", []):
        entry["source"] = str(folder / batch["cases"][0]["name"] / entry["relative"])
    return start_study(project, Path(batch["source"]), runs or folder.parent, retry=(selected, folder))
