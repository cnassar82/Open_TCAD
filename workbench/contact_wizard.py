"""Embed the existing graphical contact picker in the workbench."""

from pathlib import Path
import json
import tkinter as tk
from tkinter import filedialog, ttk

from contact_export import TERMINALS, export_device
from genius_contact_gui import ContactDialog, ContactGui
from implant_study import make_contact_map
from runner import ROOT


class PresetContactDialog(ContactDialog):
    def __init__(self, parent, options, name, kind):
        self.options = dict(options)
        super().__init__(parent, "Assign contact", name, kind)

    def body(self, master):
        picker = ttk.Frame(master)
        picker.pack(fill="x", pady=(2, 4))
        ttk.Label(picker, text="Terminal").pack(side="left", padx=4)
        choice = ttk.Combobox(picker, state="readonly", values=(*self.options, "Custom"), width=16)
        choice.set(self.default_name if self.default_name in self.options else "Custom")
        choice.pack(side="left", padx=8)
        fields = ttk.Frame(master)
        fields.pack(fill="both", expand=True)
        entry = super().body(fields)

        def select(_):
            name = choice.get()
            if name in self.options:
                self.name_var.set(name)
                self.kind_var.set(self.options[name])
            else:
                self.name_var.set("")
                entry.focus_set()

        choice.bind("<<ComboboxSelected>>", select)
        return entry


class ContactWizard(ContactGui):
    def __init__(self, app, structure):
        self.app = app
        window = tk.Toplevel(app.root)
        window.transient(app.root)
        window.geometry("1160x730")
        window.minsize(790, 510)
        self.preset = tk.StringVar(window, value="MOSFET")
        self.next_terminal = tk.StringVar(window)
        try:
            super().__init__(window, structure, suggest=False)
        except (OSError, RuntimeError, ValueError, KeyError, IndexError):
            window.destroy()
            raise
        window.title(f"Assign contacts | {structure.name}")

    def _build_widgets(self):
        bar = ttk.Frame(self.root, style="Toolbar.TFrame", padding=(12, 8))
        bar.pack(fill="x", padx=8, pady=(8, 0))
        ttk.Label(bar, text="Device type", style="Toolbar.TLabel").pack(side="left", padx=(0, 7))
        selector = ttk.Combobox(bar, values=tuple(TERMINALS), textvariable=self.preset,
                                state="readonly", width=12)
        selector.pack(side="left")
        ttk.Label(bar, textvariable=self.next_terminal, style="Toolbar.TLabel").pack(side="left", padx=10)
        selector.bind("<<ComboboxSelected>>", lambda _: self.update_terminals())
        super()._build_widgets()
        self.save_button.configure(text="Generate Genius deck", style="Primary.TButton")
        ttk.Button(self.save_button.master, text="Save contact map", command=self.save_map).pack(
            fill="x", pady=(5, 0))

    def pending(self):
        used = {contact.name.casefold() for contact in self.contacts}
        return [(name, kind) for name, kind in TERMINALS[self.preset.get()] if name.casefold() not in used]

    def update_terminals(self):
        pending = self.pending()
        text = "To assign: " + ", ".join(name for name, _ in pending) if pending else "Custom contact names accepted"
        self.next_terminal.set(text)

    def refresh_contacts(self):
        super().refresh_contacts()
        self.update_terminals()

    def ask_contact(self, default_name, default_kind):
        pending = self.pending()
        if pending:
            default_name, default_kind = pending[0]
        return PresetContactDialog(self.root, pending, default_name, default_kind)

    def save_outputs(self):
        # Leave the current editor untouched when the user cancels saving it or the export.
        if not self.app.confirm_save():
            return
        destination = filedialog.askdirectory(parent=self.root, title="Choose a folder for the Genius device",
                                               initialdir=ROOT / "tmp")
        if not destination:
            return
        try:
            deck, tif = export_device(Path(destination), self.path, self.nodes, self.tris, self.regions,
                                      self.sol_codes, self.sol_rows, self.contacts)
            self.app.load_generated_device(deck, tif)
        except (OSError, RuntimeError, ValueError, KeyError) as error:
            self.app.error(error, self.root)
            return
        self.root.destroy()

    def save_map(self):
        try:
            recipe = make_contact_map(self.nodes, self.tris, self.regions, self.contacts)
            chosen = filedialog.asksaveasfilename(
                parent=self.root, title="Save contact map", initialdir=self.path.parent,
                initialfile=f"{self.path.stem}.contacts.json", defaultextension=".json",
                filetypes=[("Contact map", "*.json")],
            )
            if chosen:
                Path(chosen).write_text(json.dumps(recipe, indent=2) + "\n", encoding="utf-8")
                self.status.set(f"Saved {Path(chosen).name}")
        except (OSError, ValueError, KeyError) as error:
            self.app.error(error, self.root)
