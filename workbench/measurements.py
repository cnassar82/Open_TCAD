"""Reusable process and electrical measurements for comparison and reports."""

import csv
from dataclasses import dataclass
import json
from math import isfinite
from pathlib import Path
import shlex
import tkinter as tk
from tkinter import filedialog, ttk

KINDS = ("Doping point", "Doping cutline", "Junction depth", "Film thickness", "VT", "BVDSS")


@dataclass(frozen=True)
class Measurement:
    kind: str
    settings: tuple[str, ...]

    @classmethod
    def from_dict(cls, value):
        kind = value["kind"]
        settings = tuple(str(item) for item in value["settings"])
        if kind not in KINDS or len(settings) != len(FIELDS[kind]):
            raise ValueError("Invalid measurement type or settings.")
        spec = cls(kind, settings)
        spec.validate()
        return spec

    def as_dict(self):
        return {"kind": self.kind, "settings": list(self.settings)}

    def validate(self):
        settings = self.settings
        try:
            numbers = {
                "Doping point": 2, "Doping cutline": 4, "Junction depth": 3,
                "Film thickness": 1, "VT": 1, "BVDSS": 1,
            }[self.kind]
            if not all(isfinite(float(value)) for value in settings[:numbers]):
                raise ValueError
            if self.kind in {"VT", "BVDSS"} and float(settings[0]) <= 0:
                raise ValueError
            if self.kind == "Junction depth" and float(settings[1]) >= float(settings[2]):
                raise ValueError
            if self.kind == "Doping cutline":
                if tuple(map(float, settings[:2])) == tuple(map(float, settings[2:4])):
                    raise ValueError
                if not 2 <= int(settings[5]) <= 2000 or settings[6] not in {"linear", "log"}:
                    raise ValueError
            if any(not item.strip() or any(char in item for char in "\n\r\x00") for item in settings):
                raise ValueError
        except (ValueError, OverflowError):
            raise ValueError("Enter finite coordinates, a positive current criterion, or a valid nonzero cut "
                             "(2–2000 samples and linear or log scale).") from None

    def report_line(self, artifact):
        value = json.dumps(self.as_dict(), separators=(",", ":"), ensure_ascii=False)
        return "#| @measure " + shlex.join([str(artifact), value])


# Each entry contains the field label and its initial value.
FIELDS = {
    "Doping point": (("X (µm)", "1"), ("Y (µm)", "0.1"), ("Field", "NetDoping")),
    "Doping cutline": (("X start (µm)", "1"), ("Y start (µm)", "0"), ("X end (µm)", "1"),
                       ("Y end (µm)", "1"), ("Field", "NetDoping"), ("Samples", "200"),
                       ("Scale", "linear")),
    "Junction depth": (("X (µm)", "1"), ("Depth search top Y (µm)", "0"),
                       ("Depth search bottom Y (µm)", "1")),
    "Film thickness": (("X (µm)", "1"), ("Material", "SiO2")),
    "VT": (("Threshold |ID| criterion (A)", "1e-7"),),
    "BVDSS": (("Breakdown |ID| criterion (A)", "2.5e-4"),),
}


def evaluate(spec, path, *, structure=None):
    """Return column names and rows. Include status rows for unresolved metrics."""
    kind, values = spec.kind, spec.settings
    if kind in {"VT", "BVDSS"}:
        from curves import read_table
        from electrical import MEASUREMENTS, breakdown, threshold

        meta = json.loads(path.read_text(encoding="utf-8"))
        expected = MEASUREMENTS[1] if kind == "VT" else MEASUREMENTS[2]
        if meta.get("kind") != expected or len(meta.get("curves", [])) != 1:
            raise ValueError(f"Choose a {expected} measurement JSON with one saved curve.")
        filename = meta["curves"][0]["file"]
        if Path(filename).name != filename:
            raise ValueError("Invalid curve filename in measurement metadata.")
        meta["criterion_A"] = float(values[0])
        extract = threshold if kind == "VT" else breakdown
        result = extract(*read_table(path.parent / filename), meta)
        return (f"{kind} (V)", "Criterion (A)", "Status"), [
            (result.get(kind + "_V", ""), meta["criterion_A"], result["status"])
        ]

    from doping import load_structure, sample_doping
    if structure is None:
        structure = load_structure(path)
    from suprem_viewer import material_thickness_rows, sample_line_profile
    if kind == "Doping point":
        x, y, field = values
        field = "Net doping" if field.casefold().replace(" ", "") == "netdoping" else field
        return ("X (µm)", "Y (µm)", "Material", field, "Status"), sample_doping(
            structure, [(float(x), float(y))], field
        )
    if kind == "Doping cutline":
        import numpy as np

        field = "NetDoping" if values[4].casefold().replace(" ", "") == "netdoping" else values[4]
        if field not in structure.fields or field == "Material":
            raise ValueError(f"{field} is not a saved scalar field.")
        distance, samples, _, _, label = sample_line_profile(
            structure, field, values[6] == "log", *map(float, values[:4]), int(values[5])
        )
        if not np.isfinite(samples).any():
            raise ValueError("The cutline has no finite samples at these coordinates.")
        return ("Distance (µm)", label), list(zip(distance.tolist(), samples.tolist()))

    x = float(values[0])
    layers = material_thickness_rows(structure, [x])
    if kind == "Film thickness":
        material = values[1]
        matches = [(row["material"], row["top_y_um"], row["bottom_y_um"], row["thickness_um"], "OK")
                   for row in layers if row["material"].casefold() == material.casefold()]
        return ("Material", "Top Y (µm)", "Bottom Y (µm)", "Thickness (µm)", "Status"), (
            matches or [(material, "", "", "", "No matching film at this X")]
        )

    import numpy as np

    # Follow crystalline-silicon intervals with the material-aware point sampler.
    search_top, search_bottom = map(float, values[1:3])
    crossings = []
    for layer in sorted(layers, key=lambda row: row["top_y_um"]):
        if layer["material_id"] != 3:
            continue
        surface = layer["top_y_um"]
        lo, hi = max(surface, search_top), min(layer["bottom_y_um"], search_bottom)
        if hi <= lo:
            continue

        def sample(points):
            rows = sample_doping(structure, [(x, float(y)) for y in points], "Net doping")
            return {row[1]: row[3] for row in rows if row[2] == layer["material"] and row[3] is not None}

        ys = np.linspace(lo, hi, 201).tolist()
        dopants = sample(ys)
        previous = None
        for y in ys:
            net = dopants.get(y)
            if net is None:
                previous = None
            elif net != 0:
                if previous and (previous[1] < 0) != (net < 0):
                    previous_y, previous_net = previous
                    crossing = previous_y - previous_net * (y - previous_y) / (net - previous_net)
                    crossings.append((crossing, crossing - surface, "Estimated from a sampled sign change"))
                previous = y, net
    return ("Junction Y (µm)", "Depth below silicon surface (µm)", "Status"), (
        crossings or [("", "", "No silicon sign change resolved in the sampled interval")]
    )


def matches_artifact(spec, path):
    return path.name.endswith(".measurement.json") if spec.kind in {"VT", "BVDSS"} else (
        path.suffix.lower() in {".str", ".vtu"}
    )


def cell(value):
    if isinstance(value, float):
        return f"{value:.8g}" if isfinite(value) else ""
    return str(value)


SCALAR_INDEX = {"Doping point": 3, "Film thickness": 3, "Junction depth": 1, "VT": 0, "BVDSS": 0}


def response_recipe(entry):
    from electrical import PREFIX

    name = entry["name"]
    artifact = Path(entry["artifact"])
    spec = Measurement.from_dict(entry["recipe"])
    if not PREFIX.fullmatch(name):
        raise ValueError("Use a response name starting with a letter or underscore.")
    if not artifact.parts or artifact.is_absolute() or ".." in artifact.parts or not matches_artifact(spec, artifact):
        raise ValueError("Use a case-relative structure or electrical measurement JSON that matches the recipe.")
    if spec.kind not in SCALAR_INDEX:
        raise ValueError("Use the checkpoint measurement dialog to compare cutlines. DOE responses are scalar.")
    return name, artifact, spec


def evaluate_response(entry, directory):
    name, artifact, spec = response_recipe(entry)
    path = (directory / artifact).resolve()
    if not path.is_relative_to(directory.resolve()):
        raise ValueError("The response artifact must be inside the case directory.")
    columns, rows = evaluate(spec, path)
    values = [row[SCALAR_INDEX[spec.kind]] for row in rows if
              isinstance(row[SCALAR_INDEX[spec.kind]], (int, float)) and isfinite(row[SCALAR_INDEX[spec.kind]])]
    result = {"columns": columns, "rows": rows}
    if len(values) == 1:
        result.update(value=values[0], status="measured")
    else:
        result["status"] = "multiple_values" if values else (rows[0][-1] if rows else "no_data")
    return name, result


class MeasurementDialog:
    def __init__(self, app, sources, peers, artifact=None):
        from ui_theme import style_table

        self.app, self.sources, self.peers, self.artifact = app, sources, peers, artifact
        self.columns, self.rows, self.profiles = (), [], []
        self.window = tk.Toplevel(app.root)
        self.window.title("Reusable measurements")
        self.window.geometry("920x570")
        body = ttk.Frame(self.window, padding=12)
        body.pack(fill="both", expand=True)
        initial_kind = KINDS[0]
        if artifact and artifact.name.endswith(".measurement.json"):
            initial_kind = "VT"
            try:
                if json.loads(sources[0][1].read_text(encoding="utf-8")).get("kind") == "BVDSS":
                    initial_kind = "BVDSS"
            except (OSError, ValueError, AttributeError):
                pass
        self.kind = tk.StringVar(value=initial_kind)
        tools = ttk.Frame(body)
        tools.pack(fill="x")
        combo = ttk.Combobox(tools, textvariable=self.kind, values=KINDS, state="readonly", width=19)
        combo.pack(side="left")
        combo.bind("<<ComboboxSelected>>", lambda _: self.build_fields())
        for label, operation in (
            ("Load recipe", self.load), ("Save recipe", self.save), ("Open files", self.choose_files),
            ("Add to report", self.add_report),
        ):
            ttk.Button(tools, text=label, command=operation).pack(side="left", padx=(6, 0))
        self.all_cases = tk.BooleanVar()
        if peers:
            ttk.Checkbutton(tools, text="Across DOE cases", variable=self.all_cases).pack(side="right")
        self.fields = ttk.Frame(body)
        self.fields.pack(fill="x", pady=(9, 3))
        self.build_fields()
        self.info = tk.StringVar(value=self.source_label())
        ttk.Label(body, textvariable=self.info, style="Muted.TLabel", wraplength=850).pack(anchor="w", pady=5)
        table_frame = ttk.Frame(body)
        table_frame.pack(fill="both", expand=True)
        self.results = ttk.Treeview(table_frame, show="headings")
        style_table(self.results)
        self.results.grid(row=0, column=0, sticky="nsew")
        vertical = ttk.Scrollbar(table_frame, command=self.results.yview)
        vertical.grid(row=0, column=1, sticky="ns")
        horizontal = ttk.Scrollbar(table_frame, orient="horizontal", command=self.results.xview)
        horizontal.grid(row=1, column=0, sticky="ew")
        self.results.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)
        footer = ttk.Frame(self.window, padding=(12, 0, 12, 12))
        footer.pack(fill="x")
        ttk.Button(footer, text="Evaluate", style="Primary.TButton", command=self.run).pack(side="right")
        ttk.Button(footer, text="Plot cutlines", command=self.plot).pack(side="right", padx=6)
        ttk.Button(footer, text="Export CSV", command=self.export).pack(side="right")

    def source_label(self):
        return "Files: " + ", ".join(path.name for _, path in self.sources) if self.sources else (
            "Open structures or electrical measurement JSON files, or launch from a run artifact."
        )

    def build_fields(self, spec=None):
        for child in self.fields.winfo_children():
            child.destroy()
        self.vars = []
        for i, (label, default) in enumerate(FIELDS[self.kind.get()]):
            variable = tk.StringVar(value=spec.settings[i] if spec else default)
            self.vars.append(variable)
            row, section = divmod(i, 3)
            ttk.Label(self.fields, text=label).grid(row=row, column=section * 2, sticky="w", pady=3)
            ttk.Entry(self.fields, textvariable=variable, width=16).grid(
                row=row, column=section * 2 + 1, sticky="w", padx=(4, 12)
            )

    def spec(self):
        return Measurement.from_dict({"kind": self.kind.get(), "settings": [var.get().strip() for var in self.vars]})

    def load(self):
        chosen = filedialog.askopenfilename(parent=self.window, filetypes=[("Measurement recipe", "*.json")])
        if chosen:
            try:
                value = json.loads(Path(chosen).read_text(encoding="utf-8"))
                spec = Measurement.from_dict(value)
                self.kind.set(spec.kind)
                self.build_fields(spec)
            except (OSError, ValueError, KeyError, TypeError) as error:
                self.app.error(error, self.window)

    def save(self):
        try:
            spec = self.spec()
            name = filedialog.asksaveasfilename(parent=self.window, defaultextension=".json",
                                                initialfile="measurement.recipe.json")
            if name:
                Path(name).write_text(json.dumps(spec.as_dict(), indent=2) + "\n", encoding="utf-8")
        except (OSError, ValueError) as error:
            self.app.error(error, self.window)

    def choose_files(self):
        initial = self.sources[0][1].parent if self.sources else Path.cwd()
        chosen = filedialog.askopenfilenames(parent=self.window, initialdir=initial, filetypes=[
            ("Structures and electrical measurements", "*.str *.vtu *.measurement.json"), ("All files", "*")
        ])
        if chosen:
            self.sources = [(str(Path(name).parent.name + "/" + Path(name).name), Path(name)) for name in chosen]
            self.artifact = None
            self.all_cases.set(False)
            self.info.set(self.source_label())

    def add_report(self):
        try:
            if self.artifact is None:
                raise ValueError("Open the dialog from a case artifact to target its case-relative name in the report.")
            self.app.insert_command(self.spec().report_line(self.artifact))
            self.info.set("Recipe inserted in the editor. Build an editor report with a matching case.")
        except ValueError as error:
            self.app.error(error, self.window)

    def run(self):
        try:
            spec = self.spec()
            sources = self.peers if self.artifact and self.all_cases.get() else self.sources
            if not sources:
                raise ValueError("Choose one or more files.")
            samples, failures = [], []
            columns = ()
            for label, path in sources:
                try:
                    if not matches_artifact(spec, path):
                        raise ValueError("The artifact type does not match this measurement.")
                    columns, rows = evaluate(spec, path)
                    samples.append((label, rows))
                except (OSError, ValueError, TypeError, KeyError, RuntimeError) as error:
                    failures.append((label, str(error)))
            self.columns = ("Source", *columns, "Message")
            self.rows = [(label, *row, "") for label, rows in samples for row in rows]
            self.rows.extend((label, *("" for _ in columns), error) for label, error in failures)
            self.profiles = [(label, columns, tuple(zip(*rows))) for label, rows in samples if rows
                             and spec.kind == "Doping cutline"]
            self.results.delete(*self.results.get_children())
            self.results.configure(columns=tuple(str(i) for i in range(len(self.columns))))
            for index, heading in enumerate(self.columns):
                self.results.heading(index, text=heading, anchor="w")
                self.results.column(index, width=150 if index == 0 else 113, stretch=True)
            for index, row in enumerate(self.rows[:2500]):
                self.results.insert("", "end", values=tuple(map(cell, row)),
                                    tags=("alternate",) if index % 2 else ())
            self.info.set(f"{len(samples)} evaluated · {len(failures)} unavailable · showing "
                          f"{min(2500, len(self.rows))} of {len(self.rows)} rows")
        except (OSError, ValueError, ImportError) as error:
            self.app.error(error, self.window)

    def plot(self):
        if self.profiles:
            self.app.plot_curves("Saved cutline comparison", [], tables=self.profiles, extraction=False)
        else:
            self.info.set("Evaluate a doping cutline before plotting.")

    def export(self):
        if not self.rows:
            self.info.set("Evaluate a measurement before exporting.")
            return
        name = filedialog.asksaveasfilename(parent=self.window, defaultextension=".csv", initialfile="measurements.csv")
        if name:
            try:
                with Path(name).open("w", newline="", encoding="utf-8") as output:
                    table = csv.writer(output)
                    table.writerow(self.columns)
                    table.writerows(tuple(map(cell, row)) for row in self.rows)
            except OSError as error:
                self.app.error(error, self.window)
