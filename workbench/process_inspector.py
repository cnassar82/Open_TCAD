"""Inspect the process deck and its recorded SUPREM4 checkpoints."""

from dataclasses import dataclass
from pathlib import Path
import re
import shlex
import tkinter as tk
from tkinter import ttk

from deck_report import checkpoint_figure
from ui_theme import code_widget, style_table


OUTPUT = re.compile(r"""\boutfile\s*=\s*(?:"([^"]+)"|'([^']+)'|([^\s#]+))""", re.I)
COMMENT = re.compile(r"""("[^"]*"|'[^']*')|[#!].*$""")


@dataclass(frozen=True)
class Step:
    first: int
    last: int
    card: str
    output: str = ""
    capture: bool = False


def steps(deck):
    """Group physical continuation lines and identify complete capture points."""
    cards = []
    parts = []
    depth = 0
    for number, line in enumerate(deck.splitlines(), 1):
        line = COMMENT.sub(lambda match: match[1] or "", line).strip()
        if not line:
            continue
        if not parts:
            first = number
        parts.append(line.removesuffix("\\").rstrip())
        if line.endswith("\\"):
            continue
        card = " ".join(parts)
        parts.clear()
        command = card.split()[0].casefold()
        if command in {"end", "}"}:
            depth = max(0, depth - 1)
        output = OUTPUT.search(card) if command == "structure" else None
        name = next((value for value in output.groups() if value is not None), "") if output else ""
        tokens = set(re.findall(r"\b\w+\b", card.casefold()))
        capture = command in {"deposit", "diffuse", "implant", "initialize", "init", "etch"}
        if command == "etch" and {"start", "continue"} & tokens:
            capture = False
        cards.append(Step(first, number, card, name, capture and depth == 0))
        if command in {"for", "{"}:
            depth += 1
    return cards


def checkpoint_path(directory, step, cards):
    # Earlier references to an overwritten filename cannot identify its old state.
    if not directory or not step.output or any(
        other.last > step.last and other.output and Path(other.output) == Path(step.output) for other in cards
    ):
        return None
    path = (directory / step.output).resolve()
    if path.is_relative_to(directory.resolve()) and path.suffix.lower() == ".str" and path.is_file():
        return path
    return None


def log_span(log, cards, selected):
    """Extract the command's output from an echo-enabled process log."""
    from collections import deque

    if cards and cards[0].card.casefold() == "set echo":
        if selected == 0:
            return None
        cards, selected = cards[1:], selected - 1
    if isinstance(log, str):
        log = log.splitlines(keepends=True)
    card_index = 0
    captured = deque(maxlen=500)
    continuation = ""
    for line in log:
        visible = COMMENT.sub(lambda match: match[1] or "", line).strip()
        echo = continuation + visible
        if visible.endswith("\\"):
            continuation = echo[:-1].rstrip() + " "
        else:
            continuation = ""
            if card_index < len(cards) and " ".join(echo.split()) == " ".join(cards[card_index].card.split()):
                if card_index == selected + 1:
                    break
                card_index += 1
        if card_index > selected:
            captured.append(line)
    return "".join(captured)[-24000:] if captured else None


class ProcessInspector:
    def __init__(self, app, *, directory=None, snapshot=None, origin=None):
        self.app, self.directory, self.snapshot, self.origin = app, directory, snapshot, origin
        self.deck = snapshot.read_text(encoding="utf-8") if snapshot else app.deck_text()
        self.steps = steps(self.deck)
        self.canvas = None
        self.window = tk.Toplevel(app.root)
        self.window.title("Process-step inspector")
        self.window.geometry("1130x710")
        body = ttk.Frame(self.window, padding=12)
        body.pack(fill="both", expand=True)
        ttk.Label(body, text=str(snapshot or app.path or "Editor deck"),
                  style="Muted.TLabel").pack(anchor="w", pady=(0, 6))
        actions = ttk.Frame(body)
        actions.pack(fill="x", pady=(0, 7))
        for label, action in (
            ("Go to editor line", self.go_to_line), ("Insert capture", self.capture),
            ("Run through step", self.run_through),
            ("Compare checkpoints", self.compare), ("Branch from checkpoint", self.branch),
            ("Measure checkpoints", self.measure),
            ("Add views to report", self.add_report),
        ):
            ttk.Button(actions, text=label, command=action).pack(side="left", padx=(0, 5))
        panes = ttk.Panedwindow(body, orient="vertical")
        panes.pack(fill="both", expand=True)
        table = ttk.Frame(panes)
        panes.add(table, weight=1)
        self.tree = ttk.Treeview(table, columns=("line", "command", "output"), show="headings",
                                 selectmode="extended", height=8)
        for key, label, width in (("line", "Line", 65), ("command", "Process card", 590),
                                  ("output", "Checkpoint", 270)):
            self.tree.heading(key, text=label, anchor="w")
            self.tree.column(key, width=width, stretch=key == "command")
        style_table(self.tree)
        scroll = ttk.Scrollbar(table, command=self.tree.yview)
        self.tree.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        self.tree.pack(fill="both", expand=True)
        for i, step in enumerate(self.steps):
            path = checkpoint_path(directory, step, self.steps)
            ready = " (saved)" if path else " (unavailable)" if step.output and directory else ""
            self.tree.insert("", "end", iid=str(i), values=(step.first, step.card, step.output + ready),
                             tags=("alternate",) if i % 2 else ())
        self.tabs = ttk.Notebook(panes)
        panes.add(self.tabs, weight=3)
        self.preview = ttk.Frame(self.tabs)
        self.tabs.add(self.preview, text="Structure")
        log_tab = ttk.Frame(self.tabs)
        self.tabs.add(log_tab, text="Step log")
        self.log = code_widget(log_tab, state="disabled", height=12, padx=7, pady=5)
        log_scroll = ttk.Scrollbar(log_tab, command=self.log.yview)
        self.log.configure(yscrollcommand=log_scroll.set)
        log_scroll.pack(side="right", fill="y")
        self.log.pack(fill="both", expand=True)
        self.note = tk.StringVar(self.window, value="Select a step. Ctrl-click a second checkpoint to compare.")
        ttk.Label(body, textvariable=self.note, style="Muted.TLabel").pack(anchor="w", pady=(6, 0))
        self.tree.bind("<<TreeviewSelect>>", self.select)
        if self.steps:
            initial = next((str(i) for i in reversed(range(len(self.steps)))
                            if checkpoint_path(directory, self.steps[i], self.steps)), "0")
            self.tree.selection_set(initial)
            self.tree.see(initial)

    def selection(self):
        return sorted(int(i) for i in self.tree.selection())

    def paths(self):
        return [path for i in self.selection()
                if (path := checkpoint_path(self.directory, self.steps[i], self.steps))]

    def select(self, _=None):
        indices = self.selection()
        if not indices:
            return
        i = indices[-1]
        step = self.steps[i]
        path = checkpoint_path(self.directory, step, self.steps)
        previous = next((file for card in reversed(self.steps[:i])
                         if (file := checkpoint_path(self.directory, card, self.steps))), None)
        if path:
            self.note.set(f"{path.name} · Ctrl-click another saved checkpoint to compare.")
        elif previous:
            self.note.set(f"No checkpoint recorded at line {step.last}. Showing preceding checkpoint {previous.name}.")
        else:
            self.note.set(f"No saved checkpoint at or before line {step.last}.")
        self.show_structure([path or previous] if path or previous else [])
        log = next((path for name in ("run.log", "suprem4.log")
                    if self.directory and (path := self.directory / name).is_file()), None)
        if log:
            try:
                with log.open(encoding="utf-8", errors="replace") as stream:
                    message = log_span(stream, self.steps, i)
                if message is None:
                    with log.open("rb") as stream:
                        stream.seek(max(0, log.stat().st_size - 24000))
                        recent = stream.read().decode("utf-8", errors="replace")
                    message = "This step has no matched command echo. Recent case log follows.\n\n" + recent
            except OSError as error:
                message = str(error)
        else:
            message = "Open a run case to inspect its log. Use Run with step log to record command echoes."
        self.app.set_text(self.log, message)

    def show_structure(self, paths):
        for child in self.preview.winfo_children():
            child.destroy()
        self.canvas = None
        if not paths:
            ttk.Label(self.preview, text="Select a saved checkpoint to view its material cross-section.",
                      padding=20).pack(anchor="nw")
            return
        try:
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            from doping import load_structure

            figure = checkpoint_figure([load_structure(path) for path in paths])
            self.canvas = FigureCanvasTkAgg(figure, master=self.preview)
            self.canvas.draw()
            self.canvas.get_tk_widget().pack(fill="both", expand=True)
        except (ImportError, OSError, ValueError, RuntimeError) as error:
            ttk.Label(self.preview, text=str(error), wraplength=880, padding=20).pack(anchor="nw")

    def compare(self):
        if len(self.selection()) != 2 or len(paths := self.paths()) != 2:
            self.note.set("Select two saved checkpoint rows with Ctrl-click.")
            return
        self.note.set("Both cross-sections use the same X and Y limits and material colors.")
        self.tabs.select(self.preview)
        self.show_structure(paths)

    def measure(self):
        paths = self.paths()
        if paths and len(paths) == len(self.selection()):
            from measurements import MeasurementDialog

            artifact = paths[0].relative_to(self.directory.resolve()) if len(paths) == 1 else None
            MeasurementDialog(self.app, [(path.name, path) for path in paths], [], artifact)
        else:
            self.note.set("Select one or more saved checkpoints to apply a measurement.")

    def branch(self):
        if len(self.selection()) == 1 and (paths := self.paths()):
            self.app.branch_checkpoint(paths[0])
        else:
            self.note.set("Select one saved checkpoint to branch.")

    def editor_matches(self):
        if self.app.engine.get() != "SUPREM4":
            return False
        if not self.snapshot:
            return self.app.deck_text() == self.deck
        # Rendered DOE snapshots can differ numerically from the editor.
        deck = self.deck.rstrip()
        editor = self.app.deck_text().rstrip()
        return editor == deck or editor == deck.removeprefix("set echo\n")

    def go_to_line(self):
        if len(self.selection()) != 1 or not self.editor_matches():
            self.note.set("Select one step from the current editor deck. Reopen the inspector after changing the deck.")
            return
        line = self.steps[self.selection()[0]].first
        if self.snapshot and self.deck.startswith("set echo\n") and not self.app.deck_text().startswith("set echo\n"):
            line = max(1, line - 1)
        self.app.focus_editor_line(line)

    def capture_target(self):
        if len(self.selection()) != 1 or self.snapshot or not self.editor_matches():
            self.note.set("Open Inspect process from the current editor and select one complete process step.")
            return None
        step = self.steps[self.selection()[0]]
        if not step.capture:
            self.note.set("Select a completed deposit, implant, diffuse, initialize, or etch card.")
            return None
        name = f"step-{step.last:04}.str"
        used = {Path(card.output).name for card in self.steps if card.output}
        suffix = 2
        while name in used:
            name = f"step-{step.last:04}-{suffix}.str"
            suffix += 1
        return step, name

    def capture(self):
        if target := self.capture_target():
            step, name = target
            self.app.focus_editor_line(step.last)
            self.app.editor.insert(f"{step.last}.end", "\n" + self.app.checkpoint_card("outfile", name))
            self.app.status.set(f"Added {name}. Run the updated process deck to capture this structure.")
            self.window.destroy()

    def run_through(self):
        if target := self.capture_target():
            step, name = target
            prefix = "\n".join(self.deck.splitlines()[:step.last])
            deck = f"{prefix}\n{self.app.checkpoint_card('outfile', name)}\nquit\n"
            if self.app.submit(parent=self.window, trace=True, deck_override=deck):
                self.window.destroy()

    def add_report(self):
        paths = self.paths()
        if len(paths) != len(self.selection()) or not 1 <= len(paths) <= 2:
            self.note.set("Select one or two saved checkpoint rows.")
            return
        if self.app.path not in {self.origin, self.snapshot} and not self.editor_matches():
            self.note.set("Open this run’s source deck in the editor before adding views to its report.")
            return
        directive = "@compare" if len(paths) == 2 else "@structure"
        names = " ".join(shlex.quote(str(path.relative_to(self.directory.resolve()))) for path in paths)
        self.app.insert_command(f"#| {directive} {names}")
        self.note.set("Added to the editor report. Build editor report with this case to render the views.")
