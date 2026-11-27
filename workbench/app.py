"""Native deck editor, design of experiments, and results workbench."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from tempfile import mkdtemp
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk

from curves import read_table
from deck_report import EDITOR_FILENAME, FILENAME as REPORT_FILENAME, REPORT_FILES, write_report
from implant_study import FLOW_RESULTS, RESULTS as STUDY_RESULTS, retry_study, start_study
from study_dialog import StudyDialog
from doping import DopingDialog
from measurements import MeasurementDialog
from electrical import MEASUREMENTS, column_index
from electrical_wizard import ElectricalWizard
from extraction_dialog import ExtractionDialog
from report_dialog import ReportDialog
from process_inspector import ProcessInspector
from process_wizard import OPERATIONS, ProcessWizard
from runner import (
    FACTOR, NAME, NUMBER, ROOT, RUNS, Job, artifact_paths, design, factors,
    infer_engine, local_tif_inputs, prepare_batch,
)
from suprem_cards import STRING, Parameter, parameter_token
from setup_wizard import MaskWizard, StructureWizard

from ui_theme import (
    BG, CODE_GUTTER, CODE_MUTED, CURVES, EDGE, ERROR, FG, MUTED, PANEL,
    SELECTION, SUCCESS, SYNTAX, TINT, WARNING, apply_theme, code_widget, menu_widget,
    style_table,
)

DECKS = ROOT / "decks"
SAMPLE = DECKS / "diode/my_diode.in"
VIEWER = ROOT / "suprem_viewer/suprem_viewer.py"
PLOT_FILES = {".csv", ".tsv", ".dat"}
MESH_FILES = {".str", ".vtu"}
CARDS = {
    "SUPREM4": {
        "Mesh line": "line x loc=0.0 spac=0.1 tag=left",
        "Silicon region": "region silicon xlo=left xhi=right ylo=surf yhi=back",
        "Initialize wafer": "initialize boron conc=1e15",
    },
    "Genius": {
        "Import mesh": "IMPORT TIFFile=device.tif",
        "Method": "METHOD Type=DDML1 NS=Basic LS=MUMPS MaxIt=30",
        "Solve equilibrium": "SOLVE Type=Equilibrium",
        "Voltage sweep": (
            "SOLVE Type=DC Vscan=Anode Vstart=0 Vstep=0.05 Vstop=1.2 out.prefix=forward_iv"
        ),
        "Export structure": "EXPORT VTKFile=final.vtu",
        "End deck": "END",
    },
}
TOKEN = re.compile(
    r"(?P<report>^[ \t]*#\|.*)|(?P<comment>[#!].*)|(?P<factor>" + FACTOR.pattern + r")|"
    r"(?P<command>^\s*[\w.]+)|(?P<key>[\w.]+(?=\s*=))|"
    r"(?P<number>(?<![\w.])" + NUMBER + r"(?![\w.]))",
    re.M,
)


def scrollable(parent, widget, horizontal=False):
    parent.rowconfigure(0, weight=1)
    parent.columnconfigure(0, weight=1)
    widget.grid(row=0, column=0, sticky="nsew")
    vertical = ttk.Scrollbar(parent, orient="vertical", command=widget.yview)
    vertical.grid(row=0, column=1, sticky="ns")
    widget.configure(yscrollcommand=vertical.set)
    if horizontal:
        bar = ttk.Scrollbar(parent, orient="horizontal", command=widget.xview)
        bar.grid(row=1, column=0, sticky="ew")
        widget.configure(xscrollcommand=bar.set)


class Workbench:
    def __init__(self, root):
        self.root = root
        self.root.title("Open TCAD Workbench")
        self.root.geometry("1370x870")
        self.root.minsize(960, 580)
        apply_theme(root)
        self.path = None
        self.saved = ""
        self.inputs = []
        self.jobs = {}
        self.batches = {}
        self.deck_files = {}
        self.highlight_pending = None
        self.log_position = None
        self.viewed_case = None
        self.closing = False
        self.engine = tk.StringVar(value="SUPREM4")
        self.status = tk.StringVar(value="Ready")
        self.deck_label = tk.StringVar()
        self.input_label = tk.StringVar()
        self.cursor_label = tk.StringVar()

        self.build()
        self.load_history()
        self.populate_library()
        self.set_deck(SAMPLE, SAMPLE.read_text(encoding="utf-8"))
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.root.bind("<Control-o>", lambda _: self.open_deck())
        self.root.bind("<Control-s>", lambda _: self.save_deck())
        self.root.bind("<F5>", lambda _: self.submit())
        self.root.after(200, self.tick)

    def build(self):
        top = ttk.Frame(self.root, style="Header.TFrame", padding=(24, 15, 22, 15))
        top.pack(fill="x")
        ttk.Label(top, text="Open TCAD Workbench", style="Header.TLabel").pack(side="left")
        ttk.Button(top, text="Open result", style="Header.TButton",
                   command=self.open_result).pack(side="right")
        study = ttk.Menubutton(top, text="Process-to-device DOE")
        study.pack(side="right", padx=(8, 0))
        study_menu = menu_widget(study)
        study_menu.add_command(label="New DOE project", command=self.open_study)
        study_menu.add_command(label="Open saved project", command=lambda: self.open_study(reopen=True))
        study.configure(menu=study_menu)
        report = ttk.Menubutton(top, text="Report")
        report.pack(side="right", padx=(8, 0))
        report_menu = menu_widget(report)
        report_menu.add_command(label="Insert Markdown", command=lambda: self.insert_report_notes(markdown=True))
        report_menu.add_command(label="Insert report data", command=self.insert_report_notes)
        report_menu.add_command(label="Reusable measurements", command=self.open_measurements)
        report_menu.add_separator()
        report_menu.add_command(label="Preview editor layout", command=self.preview_report)
        report_menu.add_command(label="Build report from selected case", command=self.build_case_report)
        report_menu.add_command(label="Build editor report with selected case",
                                command=lambda: self.build_case_report(editor=True))
        report.configure(menu=report_menu)
        ttk.Combobox(top, textvariable=self.engine, values=tuple(CARDS), state="readonly",
                     width=10).pack(side="right", padx=12)
        ttk.Label(top, text="Simulator", style="HeaderMeta.TLabel").pack(side="right")

        footer = ttk.Frame(self.root, style="Surface.TFrame", padding=(20, 8))
        footer.pack(side="bottom", fill="x")
        ttk.Separator(self.root).pack(side="bottom", fill="x")
        ttk.Label(footer, textvariable=self.status, style="PanelMuted.TLabel").pack(side="left")
        ttk.Label(footer, text="Ctrl+S  Save   /   F5  Run deck",
                  style="PanelMuted.TLabel").pack(side="right")

        body = ttk.Panedwindow(self.root, orient="horizontal")
        body.pack(fill="both", expand=True, padx=18, pady=(14, 10))
        sidebar = ttk.Frame(body, style="Panel.TFrame", padding=(12, 12, 12, 8), width=230)
        body.add(sidebar, weight=0)
        center = ttk.Frame(body)
        body.add(center, weight=4)
        sidebar.pack_propagate(False)

        ttk.Label(sidebar, text="Deck library", style="PanelHeading.TLabel").pack(anchor="w", pady=(0, 10))
        library_actions = ttk.Frame(sidebar, style="Surface.TFrame")
        library_actions.pack(fill="x", pady=(0, 10))
        ttk.Button(library_actions, text="New", command=self.new_deck).pack(side="left")
        ttk.Button(library_actions, text="Open", command=self.open_deck).pack(side="left", padx=5)
        ttk.Separator(sidebar).pack(fill="x", pady=(0, 6))
        library = ttk.Frame(sidebar, style="Surface.TFrame")
        library.pack(fill="both", expand=True)
        self.library = ttk.Treeview(library, show="tree", selectmode="browse")
        self.library.column("#0", width=195)
        scrollable(library, self.library, horizontal=True)
        self.library.bind("<Double-1>", self.open_library)
        self.library.bind("<Return>", self.open_library)
        ttk.Label(sidebar, text="Double-click a deck to open it.",
                  style="PanelMuted.TLabel").pack(anchor="w", pady=(9, 0))

        self.workspace = ttk.Notebook(center)
        self.workspace.pack(fill="both", expand=True)
        editor = ttk.Frame(self.workspace, padding=12)
        self.workspace.add(editor, text="Deck editor")
        self.build_editor(editor)
        self.results = ttk.Frame(self.workspace, padding=12)
        self.workspace.add(self.results, text="Run results")
        self.build_results()

    def build_editor(self, parent):
        actions = ttk.Frame(parent)
        actions.pack(fill="x", pady=(0, 10))
        ttk.Label(actions, textvariable=self.deck_label, style="Section.TLabel").pack(side="left")
        ttk.Button(actions, text="Run deck", style="Primary.TButton",
                   command=self.submit).pack(side="right")
        ttk.Button(actions, text="Design DOE", command=self.open_doe).pack(side="right", padx=6)
        ttk.Button(actions, text="Save as", command=lambda: self.save_deck(True)).pack(side="right")
        ttk.Button(actions, text="Save", command=self.save_deck).pack(side="right", padx=6)

        tools = ttk.Frame(parent, style="Toolbar.TFrame", padding=(8, 6))
        tools.pack(fill="x", pady=(0, 6))
        insert = ttk.Menubutton(tools, text="Insert command", direction="below")
        insert.pack(side="left")
        menu = menu_widget(insert)
        menu.configure(postcommand=lambda: self.populate_commands(menu))
        insert.configure(menu=menu)
        process = ttk.Menubutton(tools, text="Process wizards", direction="below")
        process.pack(side="left", padx=(6, 0))
        process_menu = menu_widget(process)
        process_menu.add_command(label="Set up starting structure", command=self.open_structure_wizard)
        process_menu.add_command(label="Import GDS mask edges", command=self.open_mask_wizard)
        process_menu.add_command(label="Assign contacts to checkpoint", command=self.choose_contact_structure)
        process_menu.add_command(label="Inspect process", command=self.inspect_process)
        process_menu.add_command(label="Inspect saved process folder", command=self.inspect_saved_process)
        process_menu.add_command(label="Run with step log", command=lambda: self.submit(trace=True))
        process_menu.add_separator()
        for operation in OPERATIONS:
            process_menu.add_command(label=f"{operation}",
                                     command=lambda value=operation: self.open_process_wizard(value))
        process.configure(menu=process_menu)
        electrical_menu = menu_widget(process)
        for kind in MEASUREMENTS:
            electrical_menu.add_command(label=f"{kind}", command=lambda value=kind: ElectricalWizard(self, value))
        electrical_menu.add_separator()
        electrical_menu.add_command(label="Assign contacts to checkpoint", command=self.choose_contact_structure)

        def update_process_menu(*_):
            if self.engine.get() == "SUPREM4":
                process.configure(text="Process wizards", menu=process_menu)
            else:
                process.configure(text="Electrical wizards", menu=electrical_menu)

        self.engine.trace_add("write", update_process_menu)
        update_process_menu()
        ttk.Button(tools, text="Mark factor", command=self.mark_factor).pack(side="left", padx=6)
        ttk.Button(tools, text="Run inputs", command=self.choose_inputs).pack(side="left")
        ttk.Label(tools, textvariable=self.input_label, style="Toolbar.TLabel").pack(side="left", padx=8)
        ttk.Label(tools, textvariable=self.cursor_label, style="Toolbar.TLabel").pack(side="right", padx=(6, 0))

        frame = ttk.Frame(parent)
        frame.pack(fill="both", expand=True)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        self.gutter = tk.Canvas(frame, width=52, bg=CODE_GUTTER, highlightthickness=0)
        self.gutter.grid(row=0, column=0, sticky="ns")
        self.editor = code_widget(frame, undo=True, padx=12, pady=8)
        self.editor.grid(row=0, column=1, sticky="nsew")
        vertical = ttk.Scrollbar(frame, command=self.editor.yview)
        vertical.grid(row=0, column=2, sticky="ns")
        horizontal = ttk.Scrollbar(frame, orient="horizontal", command=self.editor.xview)
        horizontal.grid(row=1, column=1, sticky="ew")

        def yscroll(first, last):
            vertical.set(first, last)
            self.draw_lines()
        self.editor.configure(yscrollcommand=yscroll, xscrollcommand=horizontal.set)
        for tag, color in SYNTAX.items():
            self.editor.tag_configure(tag, foreground=color)
        self.editor.tag_raise("sel")
        self.editor.bind("<<Modified>>", self.modified)
        self.editor.bind("<Configure>", lambda _: self.draw_lines())
        self.editor.bind("<KeyRelease>", self.show_cursor)
        self.editor.bind("<ButtonRelease-1>", self.show_cursor)

    def build_results(self):
        head = ttk.Frame(self.results)
        head.pack(fill="x", pady=(0, 9))
        ttk.Label(head, text="Run history", style="Section.TLabel").pack(side="left")
        ttk.Button(head, text="Open case folder", command=self.open_case_folder).pack(side="right")
        ttk.Button(head, text="Inspect process", command=lambda: self.inspect_process(result=True)).pack(side="right")
        ttk.Button(head, text="Cancel batch", command=self.cancel_batch).pack(side="right", padx=6)
        cases = ttk.Menubutton(head, text="DOE results")
        cases.pack(side="right", padx=6)
        case_menu = menu_widget(cases)
        case_menu.add_command(label="Open response table", command=self.view_flow_results)
        case_menu.add_command(label="Export response CSV", command=self.export_study)
        case_menu.add_command(label="Plot VT versus dose (implant study)", command=self.plot_study)
        case_menu.add_command(label="Rerun failed or canceled cases", command=self.retry_flow)
        cases.configure(menu=case_menu)

        table = ttk.Frame(self.results)
        table.pack(fill="both", expand=True)
        self.runs = ttk.Treeview(table, columns=("state", "settings", "vt", "time"), height=7)
        for col, label, width in (
            ("#0", "Batch / case", 285), ("state", "Status", 110),
            ("settings", "Engine / factor values", 245), ("vt", "VT (V) / status", 176), ("time", "Elapsed", 75),
        ):
            self.runs.heading(col, text=label, anchor="w")
            self.runs.column(col, width=width, stretch=col in {"#0", "settings"})
        style_table(self.runs)
        self.runs.tag_configure("batch", background=TINT)
        for state, color in (("Completed", SUCCESS), ("Running", WARNING), ("Failed", ERROR)):
            self.runs.tag_configure(state, foreground=color)
        scrollable(table, self.runs, horizontal=True)
        self.runs.bind("<<TreeviewSelect>>", self.case_changed)

        details = ttk.Notebook(self.results)
        details.pack(fill="both", expand=True, pady=(10, 0))
        self.details = details
        self.log_tab = ttk.Frame(details, padding=(3, 4))
        details.add(self.log_tab, text="Simulator log")
        self.log = code_widget(self.log_tab, state="disabled", height=9, padx=9, pady=7)
        scrollable(self.log_tab, self.log, horizontal=True)

        self.artifact_tab = ttk.Frame(details, padding=(3, 4))
        details.add(self.artifact_tab, text="Artifacts")
        controls = ttk.Frame(self.artifact_tab, style="Toolbar.TFrame", padding=6)
        controls.pack(fill="x", pady=(0, 5))
        ttk.Button(controls, text="Open selected", command=self.open_artifact).pack(side="left")
        ttk.Button(controls, text="Compare DOE curves", command=self.compare_curves).pack(side="left", padx=7)
        ttk.Button(controls, text="New deck from checkpoint", command=self.resume_checkpoint).pack(side="left")
        measurement = ttk.Menubutton(controls, text="Measurements")
        measurement.pack(side="left", padx=7)
        measure_menu = menu_widget(measurement)
        measure_menu.add_command(label="Measure doping points", command=self.measure_doping)
        measure_menu.add_command(label="Reusable measurements", command=self.open_measurements)
        measurement.configure(menu=measure_menu)
        ttk.Button(controls, text="Assign contacts", command=self.assign_artifact_contacts).pack(side="left")
        file_table = ttk.Frame(self.artifact_tab)
        file_table.pack(fill="both", expand=True)
        self.artifacts = ttk.Treeview(file_table, columns=("type",), show="tree headings", height=6)
        style_table(self.artifacts)
        self.artifacts.heading("#0", text="File", anchor="w")
        self.artifacts.heading("type", text="Format", anchor="w")
        self.artifacts.column("#0", width=380)
        self.artifacts.column("type", width=100, stretch=False)
        scrollable(file_table, self.artifacts)
        self.artifacts.bind("<Double-1>", lambda _: self.open_artifact())

    def populate_library(self):
        folders = {}
        for path in sorted(DECKS.rglob("*")):
            if path.suffix.lower() not in {".in", ".inp"} or not path.is_file():
                continue
            parent = ""
            relative = path.relative_to(DECKS)
            for i in range(1, len(relative.parts)):
                directory = relative.parts[:i]
                if directory not in folders:
                    folders[directory] = self.library.insert(parent, "end", text=directory[-1])
                parent = folders[directory]
            iid = self.library.insert(parent, "end", text=path.name)
            self.deck_files[iid] = path
            if path == SAMPLE:
                self.library.item(parent, open=True)
                self.library.selection_set(iid)
                self.library.see(iid)

    def populate_commands(self, menu):
        menu.delete(0, "end")
        for label, command in CARDS[self.engine.get()].items():
            menu.add_command(label=label, command=lambda value=command: self.insert_command(value))
        if self.engine.get() == "SUPREM4":
            menu.add_separator()
            menu.add_command(label="Save checkpoint", command=self.save_checkpoint)
            menu.add_command(label="Save checkpoint and stop", command=lambda: self.save_checkpoint(stop=True))
            menu.add_command(label="Load checkpoint", command=self.load_checkpoint)

    @staticmethod
    def checkpoint_card(field, path):
        return "structure " + parameter_token(Parameter(field, STRING, ""), str(path))

    def save_checkpoint(self, stop=False):
        value = simpledialog.askstring("Save checkpoint", "Output structure filename",
                                       initialvalue="checkpoint.str", parent=self.root)
        if value is None:
            return
        name = value.strip()
        if not name or Path(name).name != name or name in {".", ".."}:
            self.error("Enter a filename within the run folder.")
            return
        if not name.lower().endswith(".str"):
            name += ".str"
        try:
            card = self.checkpoint_card("outfile", name)
        except ValueError as error:
            self.error(error)
            return
        self.insert_command(card + ("\nquit" if stop else ""))

    def load_checkpoint(self):
        chosen = filedialog.askopenfilename(
            title="Load SUPREM4 checkpoint", parent=self.root,
            initialdir=self.path.parent if self.path else DECKS,
            filetypes=[("SUPREM4 structures", "*.str"), ("All files", "*")],
        )
        if chosen:
            try:
                self.insert_command(self.checkpoint_card("infile", Path(chosen).absolute()))
            except ValueError as error:
                self.error(error)

    def resume_checkpoint(self):
        selected = self.artifacts.selection()
        key, case = self.selected_case()
        if (not selected or case is None or self.batches[key]["engine"] != "SUPREM4"
                or Path(selected[0]).suffix.lower() != ".str"):
            self.status.set("Select a SUPREM4 structure (.str) in the artifact list.")
            return
        if case["status"] in {"Queued", "Running"}:
            self.status.set("Wait for this case to finish before opening its checkpoint.")
            return
        self.branch_checkpoint(RUNS / key / case["name"] / selected[0])

    def branch_checkpoint(self, path):
        if not path.is_file():
            self.error(f"Checkpoint file is missing: {path}")
            return
        try:
            card = self.checkpoint_card("infile", path)
        except ValueError as error:
            self.error(error)
            return
        if self.confirm_save():
            self.set_deck(None, f"# Reapply custom model and method cards before loading.\n{card}\n\n",
                          engine="SUPREM4")
            self.editor.mark_set("insert", "end-1c")
            self.show_cursor()
            self.status.set(f"Add the next process step after loading {path.name}.")

    def inspect_process(self, result=False):
        try:
            if not result:
                ProcessInspector(self)
                return
            key, case = self.selected_case()
            if case is None or self.batches[key]["engine"] != "SUPREM4":
                self.status.set("Select a SUPREM4 run case to inspect.")
                return
            batch = self.batches[key]
            directory = RUNS / key / case["name"]
            if case["status"] in {"Queued", "Running"} and case.get("phase") != "device":
                self.status.set("Wait for this case to finish writing process checkpoints.")
                return
            ProcessInspector(self, directory=directory, snapshot=directory / batch["deck"],
                             origin=Path(batch["source"]))
        except (OSError, ValueError) as error:
            self.error(error)

    def inspect_saved_process(self):
        chosen = filedialog.askopenfilename(parent=self.root, title="Select the process deck in a saved run folder",
                                            initialdir=self.path.parent if self.path else DECKS,
                                            filetypes=[("SUPREM4 deck", "*.in *.deck"), ("All files", "*")])
        if chosen:
            path = Path(chosen)
            try:
                ProcessInspector(self, directory=path.parent, snapshot=path, origin=path)
            except (OSError, ValueError) as error:
                self.error(error)

    def focus_editor_line(self, line):
        self.workspace.select(0)
        self.editor.mark_set("insert", f"{line}.0")
        self.editor.see("insert")
        self.editor.focus_set()
        self.show_cursor()

    def open_process_wizard(self, operation):
        try:
            ProcessWizard(self, operation)
        except (OSError, ValueError) as error:
            self.error(error)

    def open_study(self, reopen=False):
        path = None
        if reopen:
            chosen = filedialog.askopenfilename(parent=self.root, filetypes=[("Study project", "*.json")])
            if not chosen:
                return
            path = Path(chosen)
        try:
            StudyDialog(self, path)
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
            self.error(error)

    def submit_study(self, project, path, parent=None):
        if self.closing:
            return False
        try:
            job = start_study(project, path, RUNS)
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.error(error, parent)
            return False
        self.track_batch(job.folder, job.batch, job)
        return True

    def choose_contact_structure(self):
        chosen = filedialog.askopenfilename(
            title="Choose a SUPREM4 checkpoint for contacts", parent=self.root,
            initialdir=self.path.parent if self.path else DECKS,
            filetypes=[("SUPREM4 structures", "*.str"), ("All files", "*")],
        )
        if chosen:
            self.open_contact_wizard(Path(chosen))

    def assign_artifact_contacts(self):
        selected = self.artifacts.selection()
        key, case = self.selected_case()
        if (not selected or case is None or self.batches[key]["engine"] != "SUPREM4"
                or Path(selected[0]).suffix.lower() != ".str"):
            self.status.set("Select a SUPREM4 .str artifact to assign contacts.")
            return
        if case["status"] in {"Queued", "Running"}:
            self.status.set("Wait for this case to finish before assigning contacts.")
            return
        self.open_contact_wizard(RUNS / key / case["name"] / selected[0])

    def open_contact_wizard(self, path):
        try:
            from contact_wizard import ContactWizard

            ContactWizard(self, path)
        except (ImportError, OSError, RuntimeError, ValueError, KeyError, IndexError) as error:
            self.error(error)

    def load_generated_device(self, deck, tif):
        self.set_deck(deck, deck.read_text(encoding="utf-8"), engine="Genius")
        self.inputs = [tif]
        self.refresh_editor()
        self.status.set(f"Created {deck.name} with a staged mesh. Run deck to check equilibrium.")

    def open_structure_wizard(self):
        try:
            StructureWizard(self)
        except (OSError, ValueError, KeyError, StopIteration) as error:
            self.error(f"Could not open the starting-structure wizard. {error}")

    def open_mask_wizard(self):
        chosen = filedialog.askopenfilename(
            title="Open GDS layout", parent=self.root,
            initialdir=self.path.parent if self.path else DECKS,
            filetypes=[("GDSII layouts", "*.gds *.gdsii *.strm"), ("All files", "*")],
        )
        if chosen:
            try:
                MaskWizard(self, chosen)
            except ImportError:
                self.error("GDS import requires gdstk. See the workbench guide for the launch command.")
            except (OSError, ValueError, RuntimeError) as error:
                self.error(error)

    def insert_command(self, command):
        line = self.editor.get("insert linestart", "insert lineend")
        pos = "insert lineend" if line.strip() else "insert"
        self.editor.insert(pos, ("\n" if line.strip() else "") + command + "\n")
        self.editor.focus_set()

    def deck_text(self):
        return self.editor.get("1.0", "end-1c")

    def set_deck(self, path, body, engine=None):
        self.path = path
        self.saved = body if path else ""
        self.inputs = []
        self.engine.set(engine or infer_engine(path, body))
        self.editor.delete("1.0", "end")
        self.editor.insert("1.0", body)
        self.editor.edit_reset()
        self.editor.mark_set("insert", "1.0")
        self.editor.xview_moveto(0)
        self.editor.yview_moveto(0)
        self.refresh_editor()
        self.workspace.select(0)
        self.editor.focus_set()

    def refresh_editor(self):
        self.deck_label.set((self.path.name if self.path else "Untitled") + ("  •" if self.dirty else ""))
        self.input_label.set(f"{len(self.inputs)} staged" if self.inputs else "")
        self.show_cursor()

    @property
    def dirty(self):
        return self.deck_text() != self.saved

    def modified(self, _):
        if self.editor.edit_modified():
            self.editor.edit_modified(False)
            self.refresh_editor()
            if self.highlight_pending:
                self.root.after_cancel(self.highlight_pending)
            self.highlight_pending = self.root.after(95, self.highlight)

    def highlight(self):
        self.highlight_pending = None
        for tag in SYNTAX:
            self.editor.tag_remove(tag, "1.0", "end")
        for match in TOKEN.finditer(self.deck_text()):
            self.editor.tag_add(match.lastgroup, f"1.0+{match.start()}c", f"1.0+{match.end()}c")
        self.draw_lines()

    def draw_lines(self):
        self.gutter.delete("all")
        index = self.editor.index("@0,0")
        while (info := self.editor.dlineinfo(index)) is not None:
            self.gutter.create_text(41, info[1] + info[3] / 2, text=index.split(".")[0],
                                    anchor="e", fill=CODE_MUTED, font="TkFixedFont")
            following = self.editor.index(f"{index}+1line")
            if self.editor.compare(following, ">=", "end"):
                break
            index = following

    def show_cursor(self, _=None):
        line, col = self.editor.index("insert").split(".")
        self.cursor_label.set(f"Ln {line}   Col {int(col) + 1}")

    def confirm_save(self):
        if not self.dirty:
            return True
        answer = messagebox.askyesnocancel(
            "Unsaved deck", "Save changes to the current deck?", parent=self.root,
        )
        return answer is False or (answer is True and self.save_deck())

    def new_deck(self):
        if self.confirm_save():
            engine = self.engine.get()
            if engine == "SUPREM4":
                starter = SAMPLE.read_text(encoding="utf-8").replace("time=60", "time={{anneal=60}}")
            else:
                source = (DECKS / "diode/diode_iv.inp").read_text(encoding="utf-8")
                forward, _, _ = source.partition("EXPORT  VTKFile=diode_forward.vtu")
                starter = "\n".join(line for line in forward.splitlines()
                                    if "reverse_iv" not in line and "diode_reverse" not in line)
                starter = starter.replace("Vstop=1.2", "Vstop={{stop_voltage=1.2}}")
                starter += "\nEXPORT  VTKFile=diode_forward.vtu\n\nEND\n"
            self.set_deck(None, starter, engine)

    def open_library(self, _):
        selection = self.library.selection()
        if selection and selection[0] in self.deck_files:
            self.open_deck(self.deck_files[selection[0]])

    def open_deck(self, path=None):
        if not self.confirm_save():
            return
        chosen = path or filedialog.askopenfilename(
            parent=self.root, initialdir=self.path.parent if self.path else DECKS,
            filetypes=[("Simulator decks", "*.in *.inp *.deck"), ("All files", "*")],
        )
        if chosen:
            path = Path(chosen).absolute()
            try:
                body = path.read_text(encoding="utf-8")
                self.set_deck(path, body)
                if self.engine.get() == "Genius":
                    self.inputs = local_tif_inputs(path, body)
                    self.refresh_editor()
                self.status.set(str(path))
            except (OSError, UnicodeError) as error:
                self.error(error)

    def save_deck(self, save_as=False):
        path = self.path
        if save_as or path is None:
            name = filedialog.asksaveasfilename(
                parent=self.root, initialdir=path.parent if path else DECKS,
                initialfile=path.name if path else "",
                defaultextension=".inp" if self.engine.get() == "Genius" else ".in",
            )
            if not name:
                return False
            path = Path(name).absolute()
        try:
            body = self.deck_text()
            path.write_text(body, encoding="utf-8")
        except OSError as error:
            self.error(error)
            return False
        if self.path != path and (self.path is None or self.path.parent != path.parent):
            self.inputs = []
        self.path, self.saved = path, body
        self.refresh_editor()
        self.status.set(f"Saved {path}")
        return True

    def choose_inputs(self):
        if self.path is None:
            self.status.set("Save the deck to choose files relative to its folder.")
            return
        dialog = tk.Toplevel(self.root)
        dialog.title("Run inputs")
        dialog.transient(self.root)
        dialog.grab_set()
        body = ttk.Frame(dialog, padding=14)
        body.pack(fill="both", expand=True)
        ttk.Label(body, text="Files are copied into every case with their paths relative to the deck.").pack(anchor="w")
        files = tk.Listbox(body, bg=PANEL, fg=FG, selectbackground=SELECTION, width=60, height=7)
        files.pack(fill="both", expand=True, pady=8)
        selected = list(self.inputs)

        def refresh():
            files.delete(0, "end")
            for file in selected:
                files.insert("end", str(file.relative_to(self.path.parent)))

        def add():
            for name in filedialog.askopenfilenames(parent=dialog, initialdir=self.path.parent):
                file = Path(name).absolute()
                if file == self.path or not file.is_relative_to(self.path.parent):
                    messagebox.showerror("Run inputs", "Choose a separate file within the deck folder.", parent=dialog)
                    continue
                if file not in selected:
                    selected.append(file)
            refresh()

        def remove():
            for index in reversed(files.curselection()):
                selected.pop(index)
            refresh()

        def apply():
            self.inputs = selected
            self.refresh_editor()
            dialog.destroy()

        actions = ttk.Frame(body)
        actions.pack(fill="x")
        ttk.Button(actions, text="Add files", command=add).pack(side="left")
        ttk.Button(actions, text="Remove", command=remove).pack(side="left", padx=6)
        ttk.Button(actions, text="Apply", style="Primary.TButton", command=apply).pack(side="right")
        refresh()

    def mark_factor(self):
        ranges = self.editor.tag_ranges("sel")
        value = self.editor.get(*ranges) if ranges else ""
        if not re.fullmatch(NUMBER, value):
            self.status.set("Select a numeric value in the editor, then click Mark factor.")
            return
        name = simpledialog.askstring("Mark factor", f"Name for the selected value ({value})",
                                      parent=self.root)
        if name is None:
            return
        if not re.fullmatch(NAME, name):
            self.error("Start the factor name with a letter or underscore. Use letters, digits, and underscores.")
            return
        self.editor.replace(*ranges, "{{" + name + "=" + value + "}}")

    def open_doe(self):
        deck = self.deck_text()
        try:
            defaults = factors(deck)
            if not defaults:
                raise ValueError("Select a numeric value in the editor and click Mark factor to start a DOE.")
        except ValueError as error:
            self.error(error)
            return
        dialog = tk.Toplevel(self.root)
        dialog.title("Design of experiments")
        dialog.transient(self.root)
        dialog.grab_set()
        body = ttk.Frame(dialog, padding=16)
        body.pack(fill="both", expand=True)
        ttk.Label(body, text="Full-factorial design", style="Title.TLabel").grid(
            row=0, column=0, columnspan=2, sticky="w")
        ttk.Label(body, text="Enter the levels for each factor, separated by commas or spaces.",
                  style="Muted.TLabel").grid(row=1, column=0, columnspan=2, sticky="w", pady=(2, 10))
        levels = {}
        for row, (name, default) in enumerate(defaults.items(), start=2):
            ttk.Label(body, text=f"{name}  (default {default})").grid(row=row, column=0, sticky="w", pady=3)
            entry = ttk.Entry(body, width=29)
            entry.insert(0, default)
            entry.grid(row=row, column=1, sticky="ew", padx=(10, 0), pady=3)
            levels[name] = entry
        status = ttk.Label(body, style="Muted.TLabel")
        status.grid(row=len(levels) + 2, column=0, columnspan=2, sticky="w", pady=(10, 5))
        body.columnconfigure(1, weight=1)

        def preview(_=None):
            try:
                cases = design(deck, {name: entry.get() for name, entry in levels.items()})
                status.configure(text=f"{len(cases)} case{'s' if len(cases) != 1 else ''}")
                run.state(["!disabled"])
                return cases
            except ValueError as error:
                status.configure(text=str(error))
                run.state(["disabled"])
                return None

        def launch():
            if (cases := preview()) and self.submit(cases, parent=dialog):
                dialog.destroy()

        run = ttk.Button(body, text="Run DOE", style="Primary.TButton", command=launch)
        run.grid(row=len(levels) + 3, column=1, sticky="e", pady=(6, 0))
        for field in levels.values():
            field.bind("<KeyRelease>", preview)
        next(iter(levels.values())).focus_set()
        preview()

    def submit(self, cases=None, parent=None, *, trace=False, deck_override=None):
        if self.closing:
            return False
        deck = self.deck_text() if deck_override is None else deck_override
        if not deck.strip():
            self.error("Enter a deck before running the simulator.", parent)
            return False
        try:
            if cases is None:
                cases = [factors(deck)]
            path = self.path or DECKS / ("untitled.inp" if self.engine.get() == "Genius" else "untitled.in")
            if trace and self.engine.get() == "SUPREM4":
                deck = "set echo\n" + deck
            folder, batch = prepare_batch(path, deck, self.engine.get(), cases, self.inputs)
        except (OSError, ValueError) as error:
            self.error(error, parent)
            return False
        self.track_batch(folder, batch, Job(folder, batch))
        return True

    def track_batch(self, folder, batch, job):
        self.jobs[folder.name] = job
        self.add_batch(folder.name, batch)
        first = self.case_id(folder.name, batch["cases"][0])
        self.runs.selection_set(first)
        self.runs.focus(first)
        self.runs.see(first)
        self.workspace.select(self.results)
        self.details.select(self.log_tab)
        count = len(batch["cases"])
        self.status.set(f"Queued {count} case{'s' if count != 1 else ''} in {folder.name}")

    def load_history(self):
        for path in sorted(RUNS.glob("*/batch.json")):
            try:
                batch = json.loads(path.read_text(encoding="utf-8"))
                self.add_batch(path.parent.name, batch, expanded=False)
            except (OSError, ValueError, KeyError, TypeError) as error:
                self.status.set(f"Unable to load {path.parent.name}: {error}")

    def add_batch(self, key, batch, expanded=True):
        self.batches[key] = batch
        created = batch["created"].replace("T", " ")[:16]
        label = batch.get("title", Path(batch["source"]).stem.removesuffix(".study")) if batch.get("study") else batch["deck"]
        self.runs.insert("", 0, iid=key, text=f"{label}  ·  {created}", open=expanded,
                         tags=("batch",))
        for case in batch["cases"]:
            self.runs.insert(key, "end", iid=self.case_id(key, case), text=case["name"])
        self.update_batch(key)

    @staticmethod
    def case_id(key, case):
        return f"{key}/{case['name']}"

    @staticmethod
    def values_label(case):
        return ",  ".join(f"{name}={value}" for name, value in case["values"].items())

    def update_batch(self, key):
        batch = self.batches[key]
        finished = 0
        for index, case in enumerate(batch["cases"]):
            status = case["status"]
            interrupted_study = batch.get("study") and status == "Completed" and case.get("phase") == "process"
            if key not in self.jobs and (status in {"Queued", "Running"} or interrupted_study):
                status = "Untracked"
            elif status in {"Completed", "Failed", "Canceled"}:
                finished += 1
            elapsed = f"{case['elapsed']:.2f} s" if "elapsed" in case else ""
            tags = (status, "alternate") if index % 2 else (status,)
            display_status = case.get("phase", status).capitalize() if batch.get("study") and status == "Running" else status
            vt = f"{case['VT_V']:.6g}" if "VT_V" in case else case.get("vt_status", "").replace("_", " ")
            self.runs.item(self.case_id(key, case),
                           values=(display_status, self.values_label(case), vt, elapsed), tags=tags)
        engine = "SUPREM4 → Genius" if batch.get("study") else batch["engine"]
        self.runs.item(key, values=(f"{finished} / {len(batch['cases'])}", engine, "", ""))

    def selected_case(self):
        selected = self.runs.selection()
        if not selected:
            return None, None
        key, _, name = selected[0].partition("/")
        batch = self.batches[key]
        case = next((case for case in batch["cases"] if case["name"] == name), batch["cases"][0])
        return key, case

    def directory(self):
        key, case = self.selected_case()
        return RUNS / key / case["name"] if case else None

    def case_changed(self, _=None):
        directory = self.directory()
        if directory != self.viewed_case:
            self.viewed_case = directory
            self.log_position = None
            self.set_text(self.log, "")
        self.refresh_log()
        self.refresh_artifacts()

    def refresh_log(self):
        if self.viewed_case is None:
            return
        path = self.viewed_case / "run.log"
        _, case = self.selected_case()
        error = case.get("error", "") if case else ""
        if case:
            device_log = self.viewed_case / case.get("device_log", "electrical/case-001/run.log")
            handoff_log = self.viewed_case / "handoff.log"
            if case.get("phase") in {"device", "done"} and device_log.is_file():
                path = device_log
            elif case.get("phase") == "contacts" and handoff_log.is_file():
                path = handoff_log
        try:
            size = path.stat().st_size
            if (path, size, error) == self.log_position:
                return
            with path.open("rb") as stream:
                stream.seek(max(0, size - 160_000))
                log = stream.read(160_000).decode("utf-8", errors="replace")
        except OSError:
            log, size = "", 0
        if (path, size, error) == self.log_position:
            return
        previously_at_end = self.log.yview()[1] > 0.99 or not self.log_position or self.log_position[0] != path
        self.log_position = path, size, error
        self.set_text(self.log, (log + "\n" + error).strip("\n"), follow=previously_at_end)

    @staticmethod
    def set_text(widget, value, follow=False):
        top = widget.yview()[0]
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", value)
        widget.configure(state="disabled")
        if follow:
            widget.see("end")
        else:
            widget.yview_moveto(top)

    def refresh_artifacts(self):
        selected = self.artifacts.selection()
        self.artifacts.delete(*self.artifacts.get_children())
        key, case = self.selected_case()
        if case is None:
            return
        batch = self.batches[key]
        directory = RUNS / key / case["name"]
        paths = artifact_paths(directory, batch)
        if batch.get("study") and (directory / "run.log").is_file():
            paths.insert(0, Path("run.log"))
        for index, path in enumerate(paths):
            label = str(path)
            if batch.get("study") and path == Path("run.log"):
                label = "run.log (SUPREM4)"
            self.artifacts.insert("", "end", iid=str(path), text=label,
                                  values=(path.suffix[1:].upper(),), tags=("alternate",) if index % 2 else ())
        if selected and self.artifacts.exists(selected[0]):
            self.artifacts.selection_set(selected[0])

    def tick(self):
        for key, job in list(self.jobs.items()):
            try:
                if job.tick():
                    self.update_batch(key)
                    if self.selected_case()[0] == key:
                        self.refresh_artifacts()
            except OSError as error:
                self.update_batch(key)
                self.error(f"Could not save run history for {key}.\n{error}")
            if job.done:
                del self.jobs[key]
                self.update_batch(key)
                statuses = [case["status"] for case in job.batch["cases"]]
                parts = [f"{statuses.count(state)} {state.lower()}"
                         for state in ("Completed", "Failed", "Canceled") if state in statuses]
                self.status.set(f"{job.batch['deck']}  ·  " + ", ".join(parts))
        if self.closing and not self.jobs:
            self.root.destroy()
            return
        self.refresh_log()
        self.root.after(200, self.tick)

    def cancel_batch(self):
        key, _ = self.selected_case()
        if key not in self.jobs:
            self.status.set("Select an active batch to cancel its running and queued cases.")
            return
        try:
            self.jobs[key].cancel()
            self.update_batch(key)
            self.status.set("Canceling selected batch...")
        except OSError as error:
            self.error(error)

    def open_case_folder(self):
        directory = self.directory()
        if directory is None:
            self.status.set("Select a case to open its folder.")
            return
        self.launch(["xdg-open", str(directory)])

    def insert_report_notes(self, markdown=False):
        ReportDialog(self, markdown=markdown)

    def preview_report(self):
        name = self.path.name if self.path else "Untitled deck"
        output = Path(mkdtemp(prefix="tcad-report-")) / REPORT_FILENAME
        self.show_report(output, self.deck_text(), name=name,
                         deck_dir=self.path.parent if self.path else None,
                         metadata=("Editor preview",))

    def build_case_report(self, editor=False):
        key, case = self.selected_case()
        if case is None:
            self.status.set("Select a run case to provide the report data.")
            return
        if case["status"] in {"Queued", "Running"}:
            self.status.set("Wait for the case to finish before building its report.")
            return
        batch = self.batches[key]
        directory = RUNS / key / case["name"]
        try:
            deck = self.deck_text() if editor else (directory / batch["deck"]).read_text(encoding="utf-8")
            origin = self.path.parent if editor and self.path else Path(batch["source"]).parent
            name = self.path.name if editor and self.path else batch.get("title", batch["deck"])
            metadata = (("Current editor notes" if editor else name), case["name"], case["status"],
                        ", ".join(f"{sample}={value}" for sample, value in case["values"].items()))
            output = directory / (EDITOR_FILENAME if editor else REPORT_FILENAME)
            if self.show_report(output, deck, name=name, deck_dir=origin, case_dir=directory,
                                factors=case["values"], metadata=metadata):
                self.refresh_artifacts()
        except (OSError, UnicodeError) as error:
            self.error(error)

    def show_report(self, output, deck, *, name, **options):
        try:
            warnings = write_report(output, deck, title=f"{name} · Workbench report", **options)
        except ImportError:
            self.error("Report rendering needs Mistune. Install mistune and restart the workbench.")
            return False
        except (OSError, ValueError) as error:
            self.error(error)
            return False
        self.launch(["xdg-open", str(output)])
        detail = f" · {len(warnings)} source issue(s) shown in the report" if warnings else ""
        self.status.set(f"Report saved to {output}{detail}")
        return True

    def open_result(self):
        name = filedialog.askopenfilename(
            parent=self.root, initialdir=self.path.parent if self.path else DECKS,
            filetypes=[
                ("TCAD results", "*.str *.vtu *.csv *.tsv *.dat *.log *.json *.txt *.sol"),
                ("All files", "*"),
            ],
        )
        if name:
            self.show_artifact(Path(name))

    def open_artifact(self):
        chosen = self.artifacts.selection()
        directory = self.directory()
        if chosen and directory is not None:
            self.show_artifact(directory / chosen[0])
        else:
            self.status.set("Select an artifact to open it.")

    def open_measurements(self):
        selection = self.artifacts.selection()
        key, case = self.selected_case()
        sources, peers, artifact = [], [], None
        if case and selection:
            candidate = Path(selection[0])
            if candidate.suffix.lower() in MESH_FILES or candidate.name.endswith(".measurement.json"):
                artifact = candidate
                folder = RUNS / key
                sources = [(case["name"], folder / case["name"] / artifact)]
                peers = [(item["name"] + (" · " + self.values_label(item) if item["values"] else ""),
                          folder / item["name"] / artifact) for item in self.batches[key]["cases"]
                         if item["status"] == "Completed"]
        MeasurementDialog(self, sources, peers, artifact)

    def measure_doping(self):
        chosen = self.artifacts.selection()
        key, case = self.selected_case()
        if not chosen or case is None or Path(chosen[0]).suffix.lower() not in MESH_FILES:
            self.status.set("Select a saved .str or .vtu structure to measure doping.")
            return
        if case["status"] in {"Queued", "Running"}:
            self.status.set("Wait for this case to finish before sampling its structure.")
            return
        folder = RUNS / key
        target = folder / case["name"] / chosen[0]
        cases = []
        for entry in self.batches[key]["cases"]:
            path = folder / entry["name"] / chosen[0]
            if entry["status"] == "Completed" and path.is_file():
                label = self.values_label(entry)
                cases.append((f"{entry['name']} • {label}" if label else entry["name"], path))
        try:
            DopingDialog(self, target, cases)
        except (ImportError, OSError, ValueError, RuntimeError, KeyError, IndexError) as error:
            self.error(error)

    def show_artifact(self, path):
        suffix = path.suffix.lower()
        if path.name.endswith(".measurement.json") and self.open_measurement(path):
            return
        if path.name in REPORT_FILES or suffix in {".png", ".jpg", ".jpeg"}:
            self.launch(["xdg-open", str(path)])
        elif suffix in MESH_FILES:
            self.launch([sys.executable, str(VIEWER), str(path)])
        elif suffix in PLOT_FILES:
            self.plot_curves(path.name, [(path.name, path)])
        elif suffix == ".cgns":
            self.launch(["xdg-open", str(path)])
        else:
            try:
                with path.open("rb") as stream:
                    value = stream.read(200_000).decode("utf-8", errors="replace")
            except OSError as error:
                self.error(error)
                return
            window = tk.Toplevel(self.root)
            window.title(path.name)
            window.geometry("850x460")
            text = code_widget(window, state="disabled", padx=10, pady=8)
            scrollable(window, text, horizontal=True)
            self.set_text(text, value)

    def open_measurement(self, path):
        try:
            spec = json.loads(path.read_text(encoding="utf-8"))
            sources = [(curve["name"], path.parent / curve["file"]) for curve in spec["curves"]
                       if Path(curve["file"]).name == curve["file"]]
            if not sources:
                return False
            title = spec["kind"]
            if result := spec.get("result"):
                voltage_key = "BVDSS_V" if spec["kind"] == "BVDSS" else "VT_V"
                detail = result["status"].replace("_", " ")
                if voltage_key in result:
                    detail = (f"{voltage_key.removesuffix('_V')} = {result[voltage_key]:.6g} V "
                              f"at |ID| = {result['criterion_A']:.5g} A")
                title += f" • {detail}"
            if spec.get("simulator_status") == "Failed" or len(sources) < spec.get("expected_curves", len(sources)):
                title += " • partial simulation"
            self.plot_curves(title, sources, axes=(spec["x"], spec["y"]), measurement=spec)
            return True
        except (OSError, ValueError, TypeError, KeyError) as error:
            self.error(error)
            return True

    def compare_curves(self):
        chosen = self.artifacts.selection()
        key, _ = self.selected_case()
        if not chosen or key is None or Path(chosen[0]).suffix.lower() not in PLOT_FILES:
            self.status.set("Select a CSV, TSV, or DAT artifact to compare the matching file across cases.")
            return
        curves = []
        for case in self.batches[key]["cases"]:
            path = RUNS / key / case["name"] / chosen[0]
            if case["status"] == "Completed" and path.is_file():
                label = self.values_label(case) or case["name"]
                curves.append((label, path))
        if not curves:
            self.status.set("No completed case has this artifact yet.")
            return
        self.plot_curves(chosen[0], curves)

    def plot_study(self):
        key, _ = self.selected_case()
        if key is None or not self.batches[key].get("study"):
            self.status.set("Select an implant-study batch to plot VT versus dose.")
            return
        cases = self.batches[key]["cases"]
        if not all("dose" in case["values"] for case in cases):
            self.status.set("This DOE has no factor named dose. Open its response table to compare cases.")
            return
        ordered = sorted(cases, key=lambda case: float(case["values"]["dose"]))
        samples = [(float(case["values"]["dose"]), case.get("VT_V", float("nan"))) for case in ordered]
        if not any("VT_V" in case for case in cases):
            self.status.set("This study has no extracted VT values yet.")
            return
        def open_point(index):
            row = self.case_id(key, ordered[index])
            self.runs.item(key, open=True)
            self.runs.selection_set(row)
            self.runs.see(row)
            self.workspace.select(self.results)
            self.details.select(self.artifact_tab)
            self.root.lift()

        self.plot_curves("VT versus implant dose", [],
                         tables=[("VT", ("Dose (cm⁻²)", "VT (V)"), tuple(zip(*samples)))],
                         extraction=False, pick_point=open_point)
        self.status.set("Select a point to inspect its case files.")

    def retry_flow(self):
        key, _ = self.selected_case()
        if key is None or not self.batches[key].get("flow") or key in self.jobs:
            self.status.set("Select a finished process-to-device DOE batch to retry.")
            return
        try:
            job = retry_study(RUNS / key, self.batches[key])
            self.track_batch(job.folder, job.batch, job)
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.error(error)

    def view_flow_results(self):
        key, _ = self.selected_case()
        if key is None or not self.batches[key].get("study"):
            self.status.set("Select a process-to-device DOE or implant study.")
            return
        batch = self.batches[key]
        source = RUNS / key / (FLOW_RESULTS if batch.get("flow") else STUDY_RESULTS)
        if not source.is_file():
            self.status.set("The DOE will publish its table after the first case advances.")
            return
        try:
            with source.open(encoding="utf-8", newline="") as stream:
                headers, *responses = csv.reader(stream)
        except (OSError, ValueError) as error:
            self.error(error)
            return
        window = tk.Toplevel(self.root)
        window.title("DOE response table")
        window.geometry("1100x510")
        frame = ttk.Frame(window, padding=12)
        frame.pack(fill="both", expand=True)
        table = ttk.Treeview(frame, columns=tuple(map(str, range(len(headers)))), show="headings")
        style_table(table)
        for index, heading in enumerate(headers):
            table.heading(str(index), text=heading, anchor="w")
            table.column(str(index), width=max(105, min(245, 9 * len(heading))), stretch=False)
        for i, row in enumerate(responses):
            table.insert("", "end", iid=row[headers.index("case")], values=row,
                         tags=("alternate",) if i % 2 else ())
        scrollable(frame, table, horizontal=True)

        def inspect(_=None):
            if selection := table.selection():
                case_id = f"{key}/{selection[0]}"
                self.runs.item(key, open=True)
                self.runs.selection_set(case_id)
                self.runs.see(case_id)
                self.workspace.select(self.results)
                self.details.select(self.artifact_tab)
                self.root.lift()

        table.bind("<Double-1>", inspect)
        table.bind("<Return>", inspect)
        ttk.Label(window, padding=(12, 0, 12, 8), text="Double-click a case to inspect its artifacts.",
                  style="Muted.TLabel").pack(anchor="w")

    def export_study(self):
        key, _ = self.selected_case()
        if key is None or not self.batches[key].get("study"):
            self.status.set("Select a DOE or implant-study batch to export its responses.")
            return
        filename = FLOW_RESULTS if self.batches[key].get("flow") else STUDY_RESULTS
        source = RUNS / key / filename
        if not source.is_file():
            self.status.set("The study has not produced results yet.")
            return
        target = filedialog.asksaveasfilename(parent=self.root, initialfile=filename, defaultextension=".csv",
                                              filetypes=[("CSV table", "*.csv")])
        if target:
            try:
                if Path(target).resolve() != source.resolve():
                    shutil.copyfile(source, target)
                self.status.set(f"Saved {target}")
            except OSError as error:
                self.error(error)

    def plot_curves(self, title, sources, axes=None, measurement=None, *,
                    tables=None, extraction=True, pick_point=None):
        try:
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
            from matplotlib.figure import Figure

            if tables is None:
                tables = [(label, *read_table(path)) for label, path in sources]
            labels = tables[0][1]
            if any(columns != labels for _, columns, _ in tables):
                raise ValueError("The selected case files have different columns or column order.")
        except (ImportError, OSError, ValueError) as error:
            self.error(error)
            return

        window = tk.Toplevel(self.root)
        window.title(f"Curves • {title}")
        window.geometry("960x640")
        controls = ttk.Frame(window, padding=10)
        controls.pack(fill="x")
        options = [f"{i}  {name}" for i, name in enumerate(labels, 1)]
        x_axis = ttk.Combobox(controls, values=options, state="readonly", width=24)
        y_axis = ttk.Combobox(controls, values=options, state="readonly", width=24)
        for label, axis in (("X", x_axis), ("Y", y_axis)):
            ttk.Label(controls, text=label).pack(side="left", padx=(0, 4))
            axis.pack(side="left", padx=(0, 10))
        x_axis.current(0)
        y_axis.current(1)
        for i, label in enumerate(labels) if axes is None else ():
            if label.startswith("Vapp(") and len(set(tables[0][2][i])) > 1:
                x_axis.current(i)
                current = "I(" + label.split("(", 1)[1].split(")", 1)[0] + ")"
                if (match := next((j for j, name in enumerate(labels) if name.startswith(current)), None)) is not None:
                    y_axis.current(match)
                break
        if axes is not None:
            for axis, label in zip((x_axis, y_axis), axes):
                try:
                    axis.current(column_index(labels, label))
                except ValueError as error:
                    self.status.set(str(error))

        fig = Figure(facecolor=BG, layout="constrained")
        ax = fig.add_subplot()
        plot = FigureCanvasTkAgg(fig, master=window)
        plot.get_tk_widget().pack(side="top", fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(plot, window, pack_toolbar=False)
        toolbar.configure(background=TINT)
        for child in toolbar.winfo_children():
            child.configure(background=TINT)
            if isinstance(child, (tk.Button, tk.Checkbutton)):
                child.configure(activebackground=SELECTION)
        toolbar.pack(fill="x")
        scale = tk.BooleanVar(window, value=False)
        if pick_point is not None:
            def select_point(event):
                if not toolbar.mode and len(event.ind):
                    pick_point(int(event.ind[0]))
            plot.mpl_connect("pick_event", select_point)

        def render(_=None):
            x, y = x_axis.current(), y_axis.current()
            ax.clear()
            ax.set_facecolor(PANEL)
            ax.tick_params(colors=MUTED)
            for spine in ax.spines.values():
                spine.set_color(EDGE)
            ax.set_xlabel(labels[x], color=FG)
            ax.set_ylabel(("|" + labels[y] + "|") if scale.get() else labels[y], color=FG)
            ax.set_title(title, loc="left", color=FG)
            ax.grid(color=EDGE)
            for i, (label, _, data) in enumerate(tables):
                yy = [abs(value) for value in data[y]] if scale.get() else data[y]
                ax.plot(data[x], yy, marker="o" if pick_point is not None else ".",
                        ms=5 if pick_point is not None else 3, lw=1.5, label=label,
                        color=CURVES[i % len(CURVES)],
                        linestyle=("-", "--", "-.", ":")[(i // len(CURVES)) % 4],
                        **({"picker": 5} if pick_point is not None else {}))
            ax.set_yscale("log" if scale.get() else "linear", **({"nonpositive": "mask"} if scale.get() else {}))
            if len(tables) > 1:
                ax.legend(facecolor=PANEL, labelcolor=FG, edgecolor=EDGE, fontsize=8,
                          ncols=min(4, (len(tables) + 7) // 8))
            plot.draw_idle()
            toolbar.update()

        ttk.Checkbutton(controls, text="Log |Y|", variable=scale, command=render).pack(side="left")
        if extraction:
            ttk.Button(
                controls, text="Extract VT / BVDSS",
                command=lambda: ExtractionDialog(window, tables, x_axis.current(), y_axis.current(),
                                                 measurement),
            ).pack(side="right")
        for axis in (x_axis, y_axis):
            axis.bind("<<ComboboxSelected>>", render)
        render()

    def launch(self, command):
        try:
            subprocess.Popen(command, start_new_session=True)
        except OSError as error:
            self.error(error)

    def error(self, error, parent=None):
        messagebox.showerror("Open TCAD Workbench", str(error), parent=parent or self.root)

    def close(self):
        if not self.confirm_save():
            return
        if self.jobs:
            if not messagebox.askyesno("Active runs", "Cancel the active batches and close the workbench?",
                                       parent=self.root):
                return
            for job in self.jobs.values():
                try:
                    job.cancel()
                except OSError as error:
                    self.error(error)
                    return
            self.closing = True
            self.status.set("Canceling active batches...")
        else:
            self.root.destroy()


def main():
    root = tk.Tk()
    Workbench(root)
    root.mainloop()


if __name__ == "__main__":
    main()
