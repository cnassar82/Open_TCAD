"""Material-aware point sampling of doping in a saved TCAD structure."""

import csv
from pathlib import Path
import re
import sys
from tkinter import filedialog, ttk
import tkinter as tk

from runner import ROOT
from setup_wizard import SetupDialog, number
from ui_theme import code_widget, style_table

DOPING_COLUMNS = {
    "Donor": ("Nd[cm-3]", "Donor", "donor"),
    "Acceptor": ("Na[cm-3]", "Acceptor", "acceptor"),
    "Net doping": ("net_doping[cm-3]", "NetDoping", "net_doping"),
}
DONORS = ("Arsenic", "Phosphorus", "Antimony", "Germanium", "Silicon", "Selenium", "Tin")
ACCEPTORS = ("Boron", "Aluminum", "Beryllium", "Carbon", "Magnesium", "Zinc")


def load_structure(path):
    viewer_path = str(ROOT / "suprem_viewer")
    if viewer_path not in sys.path:
        sys.path.insert(0, viewer_path)
    from suprem_viewer import read_structure

    return read_structure(path)


def sample_doping(structure, points, field):
    import numpy as np

    triangles = structure.triangles
    coords = structure.points[triangles]
    start = coords[:, 0]
    leg1 = coords[:, 1] - start
    leg2 = coords[:, 2] - start
    denominator = leg1[:, 0] * leg2[:, 1] - leg1[:, 1] * leg2[:, 0]
    valid = np.abs(denominator) > 1e-30
    point_fields = structure.point_fields or {}
    cell_fields = structure.cell_fields or {}
    columns = DOPING_COLUMNS.get(field, (field,))
    species = {name: i for i, name in enumerate(structure.solution_names)}

    def species_indices(names):
        return [species.get(name + "Active", species.get(name)) for name in names
                if name + "Active" in species or name in species]

    donor = species_indices(DONORS)
    acceptor = species_indices(ACCEPTORS)

    def weighted(values, weights):
        return float(np.dot(values, weights)) if np.all(np.isfinite(values)) else None

    def concentration(index, region, weights):
        for name in columns:
            if name in point_fields:
                return weighted(point_fields[name][triangles[index]], weights)
            if name in cell_fields:
                value = float(cell_fields[name][index])
                return value if np.isfinite(value) else None
        if field == "Net doping" and "Nd[cm-3]" in point_fields and "Na[cm-3]" in point_fields:
            nd, na = (point_fields[name][triangles[index]] for name in ("Nd[cm-3]", "Na[cm-3]"))
            return weighted(nd - na, weights)
        if region is None:
            return None
        nodes = [structure.solutions.get((int(node), region.material_id)) for node in triangles[index]]
        if any(value is None for value in nodes):
            return None
        signs = []
        if field in {"Donor", "Net doping"}:
            signs.extend((i, 1) for i in donor)
        if field in {"Acceptor", "Net doping"}:
            signs.extend((i, -1 if field == "Net doping" else 1) for i in acceptor)
        if field not in DOPING_COLUMNS and field in species:
            signs = [(species[field], 1)]
        if not signs and (field not in DOPING_COLUMNS or not (donor or acceptor)):
            return None
        return weighted([sum(node[i] * sign for i, sign in signs) for node in nodes], weights)

    rows = []
    for x, y in points:
        delta = np.array((x, y)) - start
        with np.errstate(divide="ignore", invalid="ignore"):
            u = (delta[:, 0] * leg2[:, 1] - delta[:, 1] * leg2[:, 0]) / denominator
            v = (leg1[:, 0] * delta[:, 1] - leg1[:, 1] * delta[:, 0]) / denominator
        w = 1 - u - v
        hits = np.flatnonzero(valid & (u >= -1e-10) & (v >= -1e-10) & (w >= -1e-10))
        if not len(hits):
            rows.append((x, y, "", None, "Outside mesh"))
            continue
        best = {}
        for index in hits:
            region = structure.regions.get(int(structure.triangle_regions[index]))
            label = region.material if region else "Unknown"
            value = concentration(index, region, [w[index], u[index], v[index]])
            if label not in best or best[label] is None:
                best[label] = value
        for material, value in best.items():
            rows.append((x, y, material, value, "OK" if value is not None else "Field unavailable"))
    return rows


class DopingDialog(SetupDialog):
    def __init__(self, app, chosen, cases):
        self.chosen = chosen
        self.cases = cases
        self.structure = load_structure(chosen)
        self.rows = []
        super().__init__(app, "Measure doping at points")
        ttk.Label(self.body, text=chosen.name, style="Muted.TLabel").pack(anchor="w")
        controls = ttk.Frame(self.body)
        controls.pack(fill="x", pady=8)
        ttk.Label(controls, text="Doping field").pack(side="left")
        species = set(DONORS + ACCEPTORS)
        choices = list(DOPING_COLUMNS) + [
            name for name in self.structure.solution_names if name.removesuffix("Active") in species
        ]
        self.field = ttk.Combobox(controls, state="readonly", values=choices, width=18)
        self.field.set("Net doping")
        self.field.pack(side="left", padx=8)
        self.apply_all = tk.BooleanVar(self.window)
        ttk.Checkbutton(controls, text="Sample matching files across the batch",
                        variable=self.apply_all, state="normal" if len(cases) > 1 else "disabled").pack(side="left")
        ttk.Label(self.body, text="Enter one X, Y pair per line in µm.").pack(anchor="w")
        self.points = code_widget(self.body, height=5, padx=10, pady=8)
        self.points.pack(fill="x", pady=(4, 8))
        ttk.Button(self.body, text="Measure points", command=self.measure).pack(anchor="w")
        table = ttk.Frame(self.body)
        table.pack(fill="both", expand=True, pady=8)
        headings = ("Case", "X (µm)", "Y (µm)", "Material", "Doping (cm⁻³)", "Status")
        self.grid = ttk.Treeview(table, columns=headings, show="headings")
        style_table(self.grid)
        for name, width in zip(headings, (138, 75, 75, 100, 112, 116)):
            self.grid.heading(name, text=name, anchor="w")
            self.grid.column(name, width=width)
        sb = ttk.Scrollbar(table, command=self.grid.yview)
        self.grid.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        self.grid.pack(side="left", fill="both", expand=True)
        self.footer("Export CSV", self.export)
        self.commit.state(["disabled"])
        self.notice.configure(text="Signed net doping is donor minus acceptor.")

    def measure(self):
        self.grid.delete(*self.grid.get_children())
        self.rows = []
        self.commit.state(["disabled"])
        try:
            rows = [re.split(r"[\s,]+", line.strip()) for line in self.points.get("1.0", "end-1c").splitlines()
                    if line.strip()]
            if not rows or any(len(row) != 2 for row in rows):
                raise ValueError("Enter one X, Y pair per line.")
            points = [(number("X", x), number("Y", y)) for x, y in rows]
            self.measured_field = self.field.get()
            selected = self.cases if self.apply_all.get() else [(self.chosen.parent.name, self.chosen)]
            for label, path in selected:
                try:
                    structure = self.structure if path == self.chosen else load_structure(path)
                    samples = sample_doping(structure, points, self.measured_field)
                except (OSError, ValueError, RuntimeError, KeyError, IndexError) as error:
                    samples = [(x, y, "", None, str(error)) for x, y in points]
                for row in samples:
                    self.rows.append((label, *row))
            for index, (case, x, y, material, value, status) in enumerate(self.rows):
                self.grid.insert("", "end", values=(case, f"{x:g}", f"{y:g}", material,
                                                   f"{value:.7g}" if value is not None else "", status),
                                 tags=("alternate",) if index % 2 else ())
            self.commit.state(["!disabled"])
            self.notice.configure(text=f"{len(self.rows)} result rows.")
        except ValueError as error:
            self.notice.configure(text=str(error))

    def export(self):
        name = filedialog.asksaveasfilename(parent=self.window, defaultextension=".csv",
                                            initialfile="doping_points.csv",
                                            filetypes=[("CSV", "*.csv")])
        if not name:
            return
        try:
            with Path(name).open("w", newline="", encoding="utf-8") as output:
                rows = csv.writer(output)
                rows.writerow(("case", "x_um", "y_um", "material", "field", "value_cm-3", "status"))
                rows.writerows((case, x, y, material, self.measured_field, value, status)
                               for case, x, y, material, value, status in self.rows)
            self.notice.configure(text=f"Saved {name}")
        except OSError as error:
            self.notice.configure(text=str(error))
