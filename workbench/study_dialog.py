"""Build and save a process-to-device DOE project."""

import json
from pathlib import Path
import re
import tkinter as tk
from tkinter import filedialog, simpledialog, ttk
from types import SimpleNamespace

from electrical import MEASUREMENTS
from electrical_wizard import ElectricalWizard
from implant_study import MARKER, OUTPUT, factor_source, flow_cases, test_entries, validate_project
from measurements import response_recipe
from runner import factors
from ui_theme import code_widget, style_table

LOAD_CARD = re.compile(r"""^[ \t]*structure[ \t]+infile[ \t]*=[ \t]*(?:"([^"]+)"|'([^']+)'|([^\s#]+))[ \t]*(?:#.*)?$""",
                       re.I | re.M)


class StudyDialog:
    def __init__(self, app, path=None):
        self.app, self.path = app, path
        project = json.loads(path.read_text(encoding="utf-8")) if path else {}
        if path and (not isinstance(project, dict) or project.get("version") not in {1, 2}):
            raise ValueError("Choose a process-to-device DOE or implant study project.")
        source = app.deck_text() if app.engine.get() == "SUPREM4" else ""
        match = LOAD_CARD.search(source)
        self.inputs = project.get("inputs", []) if project else [
            {"relative": str(item.relative_to(app.path.parent)), "source": str(item)}
            for item in app.inputs
        ] if app.path and source else []
        checkpoint = ""
        if match:
            candidate = Path(next(value for value in match.groups() if value is not None))
            if candidate.is_absolute() or app.path:
                checkpoint = str(candidate if candidate.is_absolute() else app.path.parent / candidate)
            source = LOAD_CARD.sub(MARKER, source, count=1)
        if not source:
            source = MARKER + "\n# Add process steps through the contact-ready structure.\n"
        source = re.sub(r"^[ \t]*(?:quit|exit)[ \t]*(?:#.*)?$", "", source, flags=re.I | re.M)
        self.contact_map = project.get("contact_map", {})
        self.tests = test_entries(project) if project else []
        self.measurements = project.get("measurements", [])
        levels = project.get("levels", {})
        if isinstance(levels, str):
            levels = {"dose": levels}

        win = self.window = tk.Toplevel(app.root)
        win.title("Process-to-device DOE")
        win.geometry("910x680")
        win.minsize(705, 520)
        win.transient(app.root)
        win.grab_set()
        body = ttk.Frame(win, padding=14)
        body.pack(fill="both", expand=True)
        ttk.Label(body, text="Process-to-device DOE", style="Title.TLabel").pack(anchor="w")
        tabs = ttk.Notebook(body)
        tabs.pack(fill="both", expand=True, pady=(9, 8))
        setup = ttk.Frame(tabs, padding=10)
        tabs.add(setup, text="1  Process inputs")
        self.checkpoint = tk.StringVar(win, value=project.get("checkpoint", checkpoint))
        self.contact_label = tk.StringVar(win)
        self.update_contacts()
        for i, (label, var, action) in enumerate((
            ("Starting .str (optional)", self.checkpoint, self.choose_checkpoint),
            ("Contact map", self.contact_label, self.choose_contacts),
        )):
            ttk.Label(setup, text=label).grid(row=i, column=0, sticky="w", pady=6)
            ttk.Entry(setup, textvariable=var, state="readonly").grid(row=i, column=1, sticky="ew", padx=8)
            ttk.Button(setup, text="Choose", command=action).grid(row=i, column=2)
        setup.columnconfigure(1, weight=1)
        ttk.Button(setup, text="Clear checkpoint", command=lambda: self.checkpoint.set("")).grid(
            row=2, column=2, pady=4
        )
        ttk.Label(setup, text="Load a checkpoint with @checkpoint@, or run a full process deck with initialize. "
                  "Save the contact map from a contact-ready structure using Assign contacts. "
                  "The same map is applied to each final mesh.", wraplength=640,
                  style="Muted.TLabel").grid(row=3, column=0, columnspan=3, sticky="w", pady=10)
        ttk.Label(setup, text=f"{len(self.inputs)} staged input file(s) carried into this project.",
                  style="Muted.TLabel").grid(row=4, column=0, columnspan=3, sticky="w")
        self.process = self.text_tab(tabs, "2  Process deck",
                                     "Mark numeric values with {{name=default}}. The runner writes the final "
                                     "structure and appends quit.",
                                     project.get("process", source))
        self.device = self.text_tab(tabs, "3  Device options",
                                    "Optional Genius MODEL, PMI, interface, or METHOD overrides. Use the "
                                    "same factor name and default wherever a setting should vary together.",
                                    project.get("device_cards", ""))
        electrical = ttk.Frame(tabs, padding=10)
        tabs.add(electrical, text="4  Electrical tests")
        actions = ttk.Frame(electrical)
        actions.pack(fill="x", pady=(0, 8))
        for kind in MEASUREMENTS:
            ttk.Button(actions, text="Add " + kind, command=lambda value=kind: self.edit_test(value)).pack(
                side="left", padx=(0, 5)
            )
        ttk.Button(actions, text="Remove", command=self.remove_test).pack(side="right")
        ttk.Button(actions, text="Edit", command=self.change_test).pack(side="right", padx=6)
        self.test_table = ttk.Treeview(electrical, columns=("kind", "prefix", "settings"),
                                       show="headings", selectmode="browse")
        style_table(self.test_table)
        for col, label, width in (("kind", "Test", 145), ("prefix", "Output prefix", 120),
                                  ("settings", "Test settings", 490)):
            self.test_table.heading(col, text=label, anchor="w")
            self.test_table.column(col, width=width, stretch=col == "settings")
        self.test_table.pack(fill="both", expand=True)
        self.test_table.bind("<Double-1>", lambda _: self.change_test())
        self.refresh_tests()
        doe = ttk.Frame(tabs, padding=10)
        tabs.add(doe, text="5  DOE design")
        toolbar = ttk.Frame(doe)
        toolbar.pack(fill="x")
        self.mode = tk.StringVar(win, value=project.get("design", "Full factorial"))
        ttk.Combobox(toolbar, textvariable=self.mode, values=("Full factorial", "Latin hypercube"),
                     state="readonly", width=17).pack(side="left")
        self.samples = tk.StringVar(win, value=project.get("samples", "8"))
        self.seed = tk.StringVar(win, value=project.get("seed", "0"))
        for label, var in (("LHS cases", self.samples), ("Seed", self.seed)):
            ttk.Label(toolbar, text=label).pack(side="left", padx=(10, 4))
            ttk.Entry(toolbar, textvariable=var, width=7).pack(side="left")
        ttk.Button(toolbar, text="Preview case count", command=self.preview).pack(side="right")
        ttk.Button(toolbar, text="Find factors", command=self.update_factors).pack(side="right", padx=6)
        ttk.Label(doe, text="One factor per line, for example dose=2e12,3e12,4e12. "
                  "For Latin hypercube, enter one fixed value or two range endpoints.",
                  style="Muted.TLabel", wraplength=690).pack(anchor="w", pady=8)
        self.levels = code_widget(doe, wrap="none", height=8, undo=True, padx=8, pady=6)
        self.levels.pack(fill="both", expand=True)
        self.levels.insert("1.0", "\n".join(f"{name}={value}" for name, value in levels.items()))
        responses = ttk.Frame(tabs, padding=10)
        tabs.add(responses, text="6  Measurements")
        recipe_actions = ttk.Frame(responses)
        recipe_actions.pack(fill="x", pady=(0, 8))
        ttk.Button(recipe_actions, text="Add saved recipe", command=self.add_measurement).pack(side="left")
        ttk.Button(recipe_actions, text="Remove", command=self.remove_measurement).pack(side="right")
        ttk.Label(responses, text="Run optional scalar measurements on each case. Use study-device.str for the "
                  "final process structure, or specify another case-relative artifact.",
                  style="Muted.TLabel", wraplength=690).pack(anchor="w", pady=(0, 8))
        self.measurement_table = ttk.Treeview(responses, columns=("name", "kind", "artifact"),
                                              show="headings", selectmode="browse")
        style_table(self.measurement_table)
        for column, heading, width in (("name", "Response name", 175), ("kind", "Measurement", 160),
                                       ("artifact", "Case-relative artifact", 425)):
            self.measurement_table.heading(column, text=heading, anchor="w")
            self.measurement_table.column(column, width=width, stretch=column == "artifact")
        self.measurement_table.pack(fill="both", expand=True)
        self.refresh_measurements()
        footer = ttk.Frame(body)
        footer.pack(fill="x")
        self.status = tk.StringVar(win, value=str(path) if path else "")
        ttk.Label(footer, textvariable=self.status, style="Muted.TLabel", wraplength=490).pack(side="left")
        ttk.Button(footer, text="Run DOE", style="Primary.TButton", command=self.run).pack(side="right")
        ttk.Button(footer, text="Save project", command=self.save).pack(side="right", padx=6)
        win.bind("<Escape>", lambda _: win.destroy())

    @staticmethod
    def text_tab(tabs, label, hint, content):
        panel = ttk.Frame(tabs, padding=10)
        tabs.add(panel, text=label)
        ttk.Label(panel, text=hint, style="Muted.TLabel", wraplength=685).pack(anchor="w", pady=(0, 7))
        text = code_widget(panel, wrap="none", undo=True, padx=8, pady=6)
        scroll = ttk.Scrollbar(panel, command=text.yview)
        text.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        text.pack(fill="both", expand=True)
        text.insert("1.0", content)
        return text

    def choose_checkpoint(self):
        name = filedialog.askopenfilename(parent=self.window, filetypes=[("SUPREM4 structure", "*.str")])
        if name:
            self.checkpoint.set(str(Path(name).absolute()))
            if not self.process.get("1.0", "end-1c").strip():
                self.process.insert("1.0", MARKER + "\n")

    def choose_contacts(self):
        name = filedialog.askopenfilename(parent=self.window, filetypes=[("Contact map", "*.json")])
        if not name:
            return
        try:
            value = json.loads(Path(name).read_text(encoding="utf-8"))
            if value.get("version") != 1 or not value.get("contacts"):
                raise ValueError("Choose a JSON contact map saved by the contact picker.")
            self.contact_map = value
            self.update_contacts()
        except (OSError, ValueError, TypeError, AttributeError, KeyError) as error:
            self.app.error(error, self.window)

    def update_contacts(self):
        self.contact_label.set(", ".join(item["name"] for item in self.contact_map.get("contacts", [])))

    def refresh_tests(self):
        self.test_table.delete(*self.test_table.get_children())
        for index, entry in enumerate(self.tests):
            values = entry["values"]
            description = ", ".join(f"{name}={value}" for name, value in values.items()
                                    if name not in {"source", "gate", "drain", "prefix"})
            self.test_table.insert("", "end", iid=str(index), values=(entry["kind"], values["prefix"], description),
                                   tags=("alternate",) if index % 2 else ())

    def edit_test(self, kind, index=None):
        names = [item["name"] for item in self.contact_map.get("contacts", [])]
        stub = SimpleNamespace(root=self.window, set_text=self.app.set_text,
                               deck_text=lambda: "\n".join(f"ID={name}" for name in names))

        def save(values):
            entry = {"kind": kind, "values": values}
            if index is None:
                self.tests.append(entry)
            else:
                self.tests[index] = entry
            self.refresh_tests()

        initial = self.tests[index]["values"] if index is not None else None
        wizard = ElectricalWizard(stub, kind, initial=initial, on_save=save)

        def restore(event):
            if event.widget is wizard.window and self.window.winfo_exists():
                self.window.after_idle(self.window.grab_set)

        wizard.window.bind("<Destroy>", restore, add="+")

    def change_test(self):
        if selected := self.test_table.selection():
            index = int(selected[0])
            self.edit_test(self.tests[index]["kind"], index)

    def remove_test(self):
        if selected := self.test_table.selection():
            del self.tests[int(selected[0])]
            self.refresh_tests()

    def refresh_measurements(self):
        self.measurement_table.delete(*self.measurement_table.get_children())
        for index, entry in enumerate(self.measurements):
            self.measurement_table.insert("", "end", iid=str(index),
                                          values=(entry["name"], entry["recipe"]["kind"], entry["artifact"]),
                                          tags=("alternate",) if index % 2 else ())

    def add_measurement(self):
        chosen = filedialog.askopenfilename(parent=self.window, filetypes=[("Measurement recipe", "*.json")])
        if not chosen:
            return
        try:
            path = Path(chosen)
            recipe = json.loads(path.read_text(encoding="utf-8"))
            name = simpledialog.askstring("DOE measurement", "Response name", parent=self.window,
                                          initialvalue=path.stem.removesuffix(".recipe").replace("-", "_"))
            if name is None:
                return
            initial = "" if recipe.get("kind") in {"VT", "BVDSS"} else OUTPUT
            artifact = simpledialog.askstring("DOE measurement", "Case-relative artifact path",
                                              parent=self.window, initialvalue=initial)
            if artifact is None:
                return
            entry = dict(name=name.strip(), artifact=artifact.strip(), recipe=recipe)
            response_recipe(entry)
            if entry["name"].casefold() in {item["name"].casefold() for item in self.measurements}:
                raise ValueError("A measurement with this name already exists.")
            self.measurements.append(entry)
            self.refresh_measurements()
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
            self.app.error(error, self.window)

    def remove_measurement(self):
        if selected := self.measurement_table.selection():
            del self.measurements[int(selected[0])]
            self.refresh_measurements()

    def project(self):
        levels = {}
        for line in self.levels.get("1.0", "end-1c").splitlines():
            if not line.strip():
                continue
            name, separator, value = map(str.strip, line.partition("="))
            if not separator or not name or not value or name in levels:
                raise ValueError("Use one unique factor=value(s) entry per line.")
            levels[name] = value
        return {
            "version": 2, "checkpoint": self.checkpoint.get(), "levels": levels,
            "process": self.process.get("1.0", "end-1c"), "device_cards": self.device.get("1.0", "end-1c"),
            "contact_map": self.contact_map, "inputs": self.inputs, "tests": self.tests,
            "measurements": self.measurements,
            "design": self.mode.get(), "samples": self.samples.get(), "seed": self.seed.get(),
        }

    def update_factors(self):
        try:
            project = self.project()
            defaults = factors(factor_source(project))
            self.levels.delete("1.0", "end")
            self.levels.insert("1.0", "\n".join(f"{name}={project['levels'].get(name, value)}"
                                                  for name, value in defaults.items()))
            self.status.set(f"Found {len(defaults)} factor(s).")
        except ValueError as error:
            self.app.error(error, self.window)

    def preview(self):
        try:
            cases = flow_cases(self.project())
            self.status.set(f"{len(cases)} DOE case(s) · {len(self.tests)} test(s) per case")
        except ValueError as error:
            self.app.error(error, self.window)

    def save(self):
        try:
            project = self.project()
            if self.path is None:
                name = filedialog.asksaveasfilename(
                    parent=self.window, title="Save DOE project", initialfile="process-device.study.json",
                    defaultextension=".json", filetypes=[("Study project", "*.json")],
                )
                if not name:
                    return None
                self.path = Path(name).absolute()
            temporary = self.path.with_suffix(".tmp")
            temporary.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
            temporary.replace(self.path)
            self.status.set(f"Saved {self.path.name}")
            return project
        except (OSError, ValueError) as error:
            self.app.error(error, self.window)
            return None

    def run(self):
        try:
            validate_project(self.project())
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.app.error(error, self.window)
            return
        if (project := self.save()) is not None and self.app.submit_study(project, self.path, self.window):
            self.window.destroy()
