"""Edit deck annotations and insert references to report data."""

from pathlib import Path
from shlex import quote
import tkinter as tk
from tkinter import ttk, filedialog

from deck_report import markdown_lines
from ui_theme import code_widget

KINDS = {
    "Structure image": ("structure", (("file", "Structure file", "final.str"),
                                      ("field", "Field", "Material"), ("scale", "Scale", "linear"))),
    "Doping profile": ("profile", (("file", "Structure file", "final.str"),
                                   ("x0", "Start X (µm)", "1"), ("y0", "Start Y (µm)", "0"),
                                   ("x1", "End X (µm)", "1"), ("y1", "End Y (µm)", "1"),
                                   ("field", "Field", "NetDoping"), ("scale", "Scale", "linear"))),
    "Point measurement": ("point", (("file", "Structure file", "final.str"),
                                    ("x", "X (µm)", "1"), ("y", "Y (µm)", "0.1"),
                                    ("field", "Field", "NetDoping"))),
    "Extracted value": ("value", (("file", "Measurement JSON", "id_vgs.measurement.json"),
                                  ("key", "Value path", "result.VT_V"), ("unit", "Unit", "V"))),
    "Image file": ("image", (("file", "Image file", ""),)),
    "DOE factor": ("factor", (("name", "Factor name", "dose"), ("unit", "Unit", ""))),
}
OPTIONS = {"field", "scale", "unit", "title"}


class ReportDialog:
    def __init__(self, app, *, markdown=False):
        self.app = app
        self.markdown = markdown
        win = self.window = tk.Toplevel(app.root)
        win.title("Insert Markdown" if markdown else "Insert report data")
        win.transient(app.root)
        win.grab_set()
        win.bind("<Escape>", lambda _: win.destroy())
        self.panel = ttk.Frame(win, padding=16)
        self.panel.pack(fill="both", expand=True)
        ttk.Label(self.panel, text=win.title(), style="Title.TLabel").pack(anchor="w")
        if markdown:
            win.geometry("680x420")
            ttk.Label(self.panel, text="Use Markdown headings, lists, tables, or fenced code blocks.",
                      style="Muted.TLabel").pack(anchor="w", pady=(4, 9))
            self.notes = code_widget(self.panel, wrap="word", undo=True, padx=9, pady=7)
            self.notes.pack(fill="both", expand=True)
            self.notes.focus_set()
            self.new_slide = tk.BooleanVar(win)
            ttk.Checkbutton(self.panel, text="Start a new slide here", variable=self.new_slide).pack(anchor="w", pady=(7, 3))
        else:
            win.geometry("670x470")
            chooser = ttk.Frame(self.panel)
            chooser.pack(fill="x", pady=10)
            ttk.Label(chooser, text="Include").pack(side="left")
            self.kind = tk.StringVar(win, value="Structure image")
            types = ttk.Combobox(chooser, textvariable=self.kind, values=tuple(KINDS), state="readonly", width=22)
            types.pack(side="left", padx=8)
            types.bind("<<ComboboxSelected>>", lambda _: self.fields())
            self.form = ttk.Frame(self.panel)
            self.form.pack(fill="both", expand=True)
            self.form.columnconfigure(1, weight=1)
            self.fields()

        actions = ttk.Frame(self.panel)
        actions.pack(fill="x", pady=(9, 0))
        self.hint = tk.StringVar(win)
        ttk.Label(actions, textvariable=self.hint, style="Muted.TLabel", wraplength=415).pack(side="left")
        ttk.Button(actions, text="Insert into deck", style="Primary.TButton", command=self.insert).pack(side="right")
        ttk.Button(actions, text="Cancel", command=win.destroy).pack(side="right", padx=7)

    def fields(self):
        for child in self.form.winfo_children():
            child.destroy()
        self.values = {}
        kind, fields = KINDS[self.kind.get()]
        selected = self.app.artifacts.selection()
        for row, (key, label, default) in enumerate((*fields, ("title", "Caption (optional)", ""))):
            if key == "file" and selected:
                suffix = Path(selected[0]).suffix.lower()
                if ((kind in {"structure", "profile", "point"} and suffix in {".str", ".vtu"}) or
                        (kind == "value" and suffix == ".json") or (kind == "image" and suffix in {".png", ".jpg", ".jpeg"})):
                    default = selected[0]
            value = self.values[key] = tk.StringVar(self.window, value=default)
            ttk.Label(self.form, text=label).grid(row=row, column=0, sticky="w", pady=4)
            if key == "scale":
                field = ttk.Combobox(self.form, textvariable=value, values=("linear", "log"), state="readonly")
            else:
                field = ttk.Entry(self.form, textvariable=value)
            field.grid(row=row, column=1, sticky="ew", padx=9, pady=4)
            if key == "file":
                ttk.Button(self.form, text="Choose", command=self.choose_file).grid(row=row, column=2)
        ttk.Label(self.form, text="Use case-relative paths for simulated outputs. Local images can also be beside the deck.",
                  style="Muted.TLabel", wraplength=545).grid(row=len(fields) + 1, column=0, columnspan=3,
                                                            sticky="w", pady=(12, 0))

    def choose_file(self):
        kind = KINDS[self.kind.get()][0]
        case_dir = self.app.directory()
        deck_dir = self.app.path.parent if self.app.path else None
        roots = [case_dir] if kind != "image" else [deck_dir, case_dir]
        roots = [path.resolve() for path in roots if path is not None]
        chosen = filedialog.askopenfilename(parent=self.window, initialdir=roots[0] if roots else None)
        if not chosen:
            return
        chosen = Path(chosen).resolve()
        for root in roots:
            if chosen.is_relative_to(root):
                self.values["file"].set(chosen.relative_to(root).as_posix())
                self.hint.set("")
                return
        self.hint.set("Choose a file in the selected case or enter its relative name.")

    def insert(self):
        if self.markdown:
            text = self.notes.get("1.0", "end-1c")
            if self.new_slide.get():
                text = "@slide\n" + text
            if not text.strip():
                self.hint.set("Enter some Markdown first.")
                return
        else:
            kind = KINDS[self.kind.get()][0]
            arguments = [f"@{kind}"]
            for name, field in self.values.items():
                value = field.get().strip()
                if name in OPTIONS:
                    if value:
                        arguments.append(f"{name}={quote(value)}")
                elif not value or "\n" in value:
                    self.hint.set("Fill in the required fields.")
                    return
                else:
                    arguments.append(quote(value))
            text = " ".join(arguments)
        self.app.workspace.select(0)
        self.app.insert_command(markdown_lines(text))
        self.window.destroy()
