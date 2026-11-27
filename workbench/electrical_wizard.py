"""Genius measurement forms."""

import re

from electrical import MEASUREMENTS, measurement_cards
from process_wizard import CardForm, ProcessWizard
from setup_wizard import SetupDialog
from runner import factors, render_deck
from tkinter import ttk


class ElectricalWizard(SetupDialog):
    def __init__(self, app, kind, *, initial=None, on_save=None):
        self.kind, self.on_save = kind, on_save
        super().__init__(app, kind)
        hint = ("Numeric fields can use {{name=number}} to define DOE factors." if on_save else
                "Insert after the device model and METHOD cards, before END.")
        ttk.Label(self.body, text=hint,
                  style="Muted.TLabel").pack(anchor="w")
        region = ttk.Frame(self.body)
        region.pack(fill="both", expand=True, pady=(6, 0))
        self.form = CardForm(ProcessWizard.scroll_form(region), self.refresh)
        deck = app.deck_text()
        names = sorted(set(re.findall(r"\b(?:ID|Label)\s*=\s*['\"]?([\w.:-]+)", deck, re.I)))

        def contact(name):
            matches = [label for label in names if name.lower() in label.lower()]
            return matches[0] if len(matches) == 1 else name

        for name in ("Drain", "Gate", "Source"):
            self.form.field(name.lower(), f"{name} electrode", default=contact(name))
        defaults = ("id_vds", "id_vgs", "bvdss")
        self.form.field("prefix", "Output prefix", default=defaults[MEASUREMENTS.index(kind)])
        if kind == MEASUREMENTS[0]:
            self.form.field("levels", "Gate biases (V)", "Separate fixed VGS values with commas or spaces.",
                            default="0, 0.5, 1, 1.5, 2")
            self.form.field("stop", "VDS stop (V)", "Sweep starts at 0 V.", default="2")
        elif kind == MEASUREMENTS[1]:
            self.form.field("drain_bias", "Fixed VDS (V)", default="0.05")
            self.form.field("start", "VGS start (V)", default="0")
            self.form.field("stop", "VGS stop (V)", default="2")
            self.form.field("criterion", "VT current criterion (A)", "Optional. Use the total drain current "
                            "specified for your device and fixed VDS.", default="")
        else:
            self.form.field("stop", "VDS stop (V)", "Apply the sign for your device; VGS is held at 0 V.", default="20")
            self.form.field("criterion", "Drain current criterion (A)", "Magnitude of total drain current in the run.",
                            default="2.5e-4")
            self.form.field("impact", "Impact ionization", choices=("Existing model", "Local", "Nonlocal"),
                            default="Existing model")
            self.form.field("region", "Semiconductor region", "Required when setting an impact model.",
                            default="Silicon")
        self.form.field("sweep_step", "Sweep step magnitude (V)", default="0.1")
        if kind != "BVDSS":
            self.form.field("ramp_step", "Bias ramp step magnitude (V)", default="0.1")
        self.preview_box()
        self.footer("Save test" if self.on_save else "Insert into deck", self.insert)
        for key, value in (initial or {}).items():
            if key in self.form.values:
                self.form.values[key].set(value)
        self.refresh()

    def get_values(self):
        return {name: var.get().strip() for name, var in self.form.values.items()}

    def insert(self):
        if not self.on_save:
            super().insert()
            return
        values = self.get_values()
        self.window.destroy()
        self.on_save(values)

    def refresh(self, *_):
        if not hasattr(self, "preview"):
            return
        try:
            args = self.get_values()
            if self.on_save:
                defaults = factors("\n".join(args.values()))
                args = {name: render_deck(value, defaults) for name, value in args.items()}
            self.show_preview(measurement_cards(self.kind, args),
                              "Set the source and body biases to 0 V in the device deck.")
        except ValueError as error:
            self.show_preview("", str(error))
