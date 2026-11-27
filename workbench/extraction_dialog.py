"""Extract constant-current voltages from the axes of an open curve plot."""

import csv
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, ttk

from electrical import MEASUREMENTS, breakdown, column_index, current_crossing, threshold
from ui_theme import style_table

VT = "VT"
BVDSS = "BVDSS"
COLUMNS = ("Curve", "Voltage (V)", "Last sampled (V)", "Peak |ID| (A)", "Status")


class ExtractionDialog:
    def __init__(self, parent, tables, x, y, measurement=None):
        self.tables, self.x, self.y = tables, x, y
        self.measurement = measurement
        self.results = []
        self.window = tk.Toplevel(parent)
        self.window.title("Extract VT or BVDSS")
        self.window.geometry("850x390")
        self.window.minsize(660, 320)
        self.window.transient(parent)
        self.window.grab_set()
        self.window.bind("<Escape>", lambda _: self.window.destroy())
        body = ttk.Frame(self.window, padding=15)
        body.pack(fill="both", expand=True)
        labels = tables[0][1]
        ttk.Label(body, text="Extract VT or BVDSS", style="Title.TLabel").pack(anchor="w")
        ttk.Label(body, text=f"Voltage: {labels[x]}     Drain current: {labels[y]}",
                  style="Muted.TLabel").pack(anchor="w", pady=(2, 10))

        controls = ttk.Frame(body)
        controls.pack(fill="x")
        self.method = ttk.Combobox(controls, values=(VT, BVDSS), state="readonly", width=10)
        self.method.pack(side="left")
        self.method.set(BVDSS if measurement and measurement["kind"] == BVDSS else VT)
        ttk.Label(controls, text="Drain current criterion (A)").pack(side="left", padx=(12, 6))
        self.criterion = ttk.Entry(controls, width=13)
        self.criterion.pack(side="left")
        if measurement and "criterion_A" in measurement:
            self.criterion.insert(0, f"{measurement['criterion_A']:g}")
        ttk.Button(controls, text="Extract", style="Primary.TButton", command=self.extract).pack(side="left", padx=8)
        self.export_button = ttk.Button(controls, text="Export CSV", command=self.export, state="disabled")
        self.export_button.pack(side="right")
        self.hint = ttk.Label(body, style="Muted.TLabel", wraplength=710)
        self.hint.pack(anchor="w", pady=(8, 6))
        self.method.bind("<<ComboboxSelected>>", self.update_hint)
        self.criterion.bind("<Return>", lambda _: self.extract())
        self.update_hint()

        frame = ttk.Frame(body)
        frame.pack(fill="both", expand=True)
        self.grid = ttk.Treeview(frame, columns=COLUMNS, show="headings")
        style_table(self.grid)
        for name, width in zip(COLUMNS, (180, 105, 130, 120, 190)):
            self.grid.heading(name, text=name, anchor="w")
            self.grid.column(name, width=width, minwidth=75)
        bar = ttk.Scrollbar(frame, command=self.grid.yview)
        self.grid.configure(yscrollcommand=bar.set)
        bar.pack(side="right", fill="y")
        self.grid.pack(side="left", fill="both", expand=True)
        self.notice = ttk.Label(body, style="Muted.TLabel")
        self.notice.pack(anchor="w", pady=(8, 0))
        self.criterion.focus_set()

    def update_hint(self, _=None):
        if self.method.get() == VT:
            hint = "Use VGS on X and ID on Y from a transfer curve at fixed, nonzero VDS."
        else:
            hint = "Use VDS on X and ID on Y from an off-state sweep with VGS = 0 V."
        self.hint.configure(text=hint + " The criterion uses |ID| and the result keeps the voltage sign.")

    def extract(self):
        self.grid.delete(*self.grid.get_children())
        self.export_button.state(["disabled"])
        self.results = []
        try:
            criterion = float(self.criterion.get())
        except ValueError:
            self.notice.configure(text="Enter a positive current criterion in amperes.")
            return
        method = self.method.get()
        labels = self.tables[0][1]
        spec = self.measurement
        native = False
        if spec is not None and spec["kind"] == (BVDSS if method == BVDSS else MEASUREMENTS[1]):
            try:
                native = (self.x, self.y) == (column_index(labels, spec["x"]), column_index(labels, spec["y"]))
            except ValueError:
                pass
        for label, columns, data in self.tables:
            try:
                if native:
                    extractor = breakdown if method == BVDSS else threshold
                    result = extractor(columns, data, {**spec, "criterion_A": criterion})
                else:
                    result = current_crossing(data[self.x], data[self.y], criterion, low_drive_first=True)
            except ValueError as error:
                self.notice.configure(text=str(error))
                self.results = []
                self.grid.delete(*self.grid.get_children())
                return
            voltage = result.get("voltage_V", result.get(f"{method}_V"))
            last = result.get("last_voltage_V", result.get("last_Vds_V"))
            row = (label, voltage, last, result.get("max_abs_Id_A"), result["status"])
            self.results.append(row)
            rendered = [f"{value:.7g}" if value is not None else "" for value in row[1:4]]
            self.grid.insert("", "end", values=(label, *rendered, row[4].replace("_", " ")),
                             tags=("alternate",) if len(self.results) % 2 == 0 else ())
        self.export_button.state(["!disabled"])
        self.result_metadata = (method, criterion, labels[self.x], labels[self.y])
        reached = sum(row[1] is not None for row in self.results)
        self.notice.configure(text=f"{reached} of {len(self.results)} curves reached the criterion.")

    def export(self):
        name = filedialog.asksaveasfilename(parent=self.window, initialfile="curve_extractions.csv",
                                           defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not name:
            return
        try:
            with Path(name).open("w", encoding="utf-8", newline="") as output:
                writer = csv.writer(output)
                writer.writerow(("method", "criterion_A", "voltage_column", "current_column", "curve",
                                 "voltage_V", "last_voltage_V", "max_abs_Id_A", "status"))
                writer.writerows((*self.result_metadata, *row) for row in self.results)
            self.notice.configure(text=f"Saved {Path(name).name}")
        except OSError as error:
            self.notice.configure(text=str(error))
