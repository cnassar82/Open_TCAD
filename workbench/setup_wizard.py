"""Starting structure and layout mask dialogs for SUPREM4 decks."""

from pathlib import Path
from math import isfinite
from tkinter import filedialog, ttk
import tkinter as tk

from layout import fmt, invert_intervals, layer_pairs, mask_commands, merge_mesh_lines, read_layout, section
from process_wizard import CardForm, ProcessWizard
from suprem_cards import catalog, key_file, numeric_literal
from ui_theme import ACCENT, EDGE, FG, MUTED, PANEL, code_widget


def number(label, value, positive=False):
    parsed = numeric_literal(value)
    if parsed is None or "{{" in value or (positive and parsed <= 0):
        raise ValueError(f"{label} requires a {'positive ' if positive else ''}number.")
    return parsed


def structure_cards(values):
    def get(key, label, positive=False):
        return number(label, values[key], positive)

    left = get("origin", "Left X")
    width = get("width", "Width", True)
    depth = get("depth", "Depth", True)
    sx = get("dx", "Lateral spacing", True)
    surface = get("surface", "Surface spacing", True)
    bulk = get("bulk", "Backside spacing", True)
    ratio = get("ratio", "Mesh growth ratio", True)
    if ratio <= 1:
        raise ValueError("Mesh growth ratio must be greater than 1.")
    right = left + width
    if not isfinite(right) or right <= left:
        raise ValueError("Choose a width with a distinct, finite right edge.")

    def mesh_line(axis, location, spacing, tag=""):
        return f"line {axis} loc={fmt(location)} spac={fmt(spacing)}" + (f" tag={tag}" if tag else "")

    commands = ["mode two.dim", mesh_line("x", left, sx, "left"), mesh_line("x", right, sx, "right"),
                mesh_line("y", 0, surface, "surf"), mesh_line("y", depth, bulk, "back"),
                f"region {values['material']} xlo=left xhi=right ylo=surf yhi=back"]
    for boundary, tag in (("surface_bc", "surf"), ("back_bc", "back")):
        if selected := values[boundary]:
            commands.append(f"boundary {selected} xlo=left xhi=right ylo={tag} yhi={tag}")
    initialize = f"initialize orient={values['orientation']} interval.r={fmt(ratio)}"
    if dopant := values["dopant"]:
        concentration = get("conc", "Background concentration", True)
        initialize += f" {dopant} conc={fmt(concentration)}"
    commands.append(initialize)
    return "\n".join(commands)


class SetupDialog:
    def __init__(self, app, title):
        self.app = app
        self.window = tk.Toplevel(app.root)
        self.window.title(title)
        self.window.transient(app.root)
        self.window.geometry("850x680")
        self.window.minsize(680, 500)
        self.window.grab_set()
        self.window.bind("<Escape>", lambda _: self.window.destroy())
        self.body = ttk.Frame(self.window, padding=14)
        self.body.pack(fill="both", expand=True)
        ttk.Label(self.body, text=title, style="Title.TLabel").pack(anchor="w", pady=(0, 8))

    def preview_box(self):
        ttk.Label(self.body, text="Deck preview", style="Section.TLabel").pack(anchor="w", pady=(6, 4))
        self.preview = code_widget(self.body, height=8, state="disabled", padx=10, pady=8)
        self.preview.pack(fill="x")

    def footer(self, label, action):
        footer = ttk.Frame(self.body, padding=(0, 8, 0, 0))
        footer.pack(fill="x")
        self.notice = ttk.Label(footer, style="Muted.TLabel", wraplength=450)
        self.notice.pack(side="left")
        self.commit = ttk.Button(footer, text=label, style="Primary.TButton", command=action)
        self.commit.pack(side="right")
        ttk.Button(footer, text="Cancel", command=self.window.destroy).pack(side="right", padx=8)

    def show_preview(self, value, message=""):
        self.app.set_text(self.preview, value)
        self.commit.state(["!disabled"] if value else ["disabled"])
        self.notice.configure(text=message)

    def insert(self):
        commands = self.preview.get("1.0", "end-1c")
        self.window.destroy()
        self.app.insert_command(commands)


class StructureWizard(SetupDialog):
    def __init__(self, app):
        cards = catalog(key_file())
        materials = tuple(item.name for item in cards["region"]["mater"].children)
        dopants = tuple(item.name for item in cards["initialize"]["impurity"].children)
        boundary = tuple(item.name for item in cards["boundary"]["type"].children)
        super().__init__(app, "Set up starting structure")
        region = ttk.Frame(self.body)
        region.pack(fill="both", expand=True)
        self.form = CardForm(ProcessWizard.scroll_form(region), self.refresh)
        for name, label, default, hint in (
            ("origin", "Left X (µm)", "0", "Use the same X origin for the GDS cross-section."),
            ("width", "Width (µm)", "3", "Lateral extent of the simulation window."),
            ("depth", "Depth (µm)", "5", "Depth below the surface, which is at Y = 0."),
            ("dx", "Lateral spacing (µm)", "0.1", ""),
            ("surface", "Surface spacing (µm)", "0.02", ""),
            ("bulk", "Backside spacing (µm)", "0.25", ""),
            ("ratio", "Mesh growth ratio", "1.5", "Maximum ratio between adjacent mesh intervals."),
        ):
            self.form.field(name, label, hint, default=default)
        self.form.field("material", "Substrate material", choices=materials, default="silicon")
        self.form.field("orientation", "Crystal orientation", choices=("100", "110", "111"), default="100")
        self.form.field("dopant", "Background dopant", choices=("", *dopants), default="boron")
        self.form.field("conc", "Concentration (cm⁻³)", default="1e15")
        self.form.field("surface_bc", "Top boundary", choices=("", *boundary), default="exposed")
        self.form.field("back_bc", "Bottom boundary", choices=("", *boundary), default="backside")
        self.preview_box()
        self.footer("Start new deck", self.start_deck)
        self.refresh()

    def refresh(self, *_):
        try:
            commands = structure_cards({name: var.get().strip() for name, var in self.form.values.items()})
            self.show_preview(commands, "Rectangular 2D substrate with tagged mesh boundaries.")
        except ValueError as error:
            self.show_preview("", str(error))

    def start_deck(self):
        commands = self.preview.get("1.0", "end-1c")
        self.window.grab_release()
        if not self.app.confirm_save():
            self.window.grab_set()
            return
        self.window.destroy()
        self.app.set_deck(None, commands + "\n\n", engine="SUPREM4")
        self.app.editor.mark_set("insert", "end-1c")
        self.app.show_cursor()


class MaskWizard(SetupDialog):
    def __init__(self, app, source):
        self.library = read_layout(source)
        self.cells = {cell.name: cell for cell in self.library.cells}
        tops = sorted(cell.name for cell in self.library.top_level())
        selections = [*tops, *(name for name in sorted(self.cells) if name not in tops)]
        # Use a small fraction of the database grid to preserve one-DBU features.
        self.tolerance = max(self.library.precision * 1e6 * 1e-5, 1e-12)
        self.result = None
        materials = tuple(item.name for item in catalog(key_file())["etch"]["material"].children)
        super().__init__(app, "Import GDS mask edges")
        ttk.Label(self.body, text=Path(source).name, style="Muted.TLabel").pack(anchor="w")
        region = ttk.Frame(self.body)
        region.pack(fill="both", expand=True)
        self.form = CardForm(ProcessWizard.scroll_form(region), self.invalidate)
        self.form.field("cell", "Cell", "References, arrays, and paths are included.",
                        choices=selections, default=selections[0])
        self.layer = self.form.field("layer", "Layer / datatype", choices=())
        self.form.values["cell"].trace_add("write", self.change_cell)
        for name, label, default in (("x0", "Start X (µm)", "0"), ("y0", "Start Y (µm)", "0"),
                                     ("x1", "End X (µm)", "3"), ("y1", "End Y (µm)", "0"),
                                     ("origin", "Deck X at start (µm)", "0")):
            self.form.field(name, label, default=default)
        self.form.field("tone", "Use layout polygons as", choices=("Openings", "Protected"),
                        default="Openings")
        self.form.field("mode", "Deck output", choices=("Mask opening etches", "Mesh lines"),
                        default="Mask opening etches")
        self.form.field("material", "Material to etch", choices=materials, default="photoresist")
        self.form.field("top", "Etch top Y (µm)", "Start above the film.", default="-2")
        self.form.field("bottom", "Etch bottom Y (µm)", "End below the film to be removed.", default="0.5")
        self.form.field("spacing", "Edge mesh spacing (µm)", "Used for mesh line output.", default="0.02")
        controls = ttk.Frame(self.body)
        controls.pack(fill="x", pady=(5, 2))
        ttk.Button(controls, text="Compute cut", command=self.compute).pack(side="left")
        self.copy = ttk.Button(controls, text="Copy intervals", command=self.copy_intervals)
        self.copy.pack(side="left", padx=7)
        self.cross_section = tk.Canvas(self.body, height=42, bg=PANEL, highlightthickness=1,
                                       highlightbackground=EDGE)
        self.cross_section.pack(fill="x")
        self.cross_section.bind("<Configure>", lambda _: self.draw_cut())
        self.preview_box()
        self.footer("Insert into deck", self.insert)
        self.change_cell()
        self.invalidate()

    def change_cell(self, *_):
        cell = self.cells[self.form.get("cell")]
        pairs = [f"{layer} / {datatype}" for layer, datatype in layer_pairs(cell)]
        self.layer.configure(values=pairs)
        self.form.values["layer"].set(pairs[0] if pairs else "")
        if (bounds := cell.bounding_box()) is not None:
            (xmin, ymin), (xmax, ymax) = bounds
            for name, value in (("x0", xmin), ("y0", (ymin + ymax) / 2),
                                ("x1", xmax), ("y1", (ymin + ymax) / 2)):
                self.form.values[name].set(fmt(value))

    def invalidate(self, *_):
        if not hasattr(self, "preview"):
            return
        self.result = None
        self.copy.state(["disabled"])
        self.draw_cut()
        self.show_preview("", "Compute the cut to preview the mask intervals.")

    def compute(self):
        self.invalidate()
        try:
            form = self.form
            if not (spec := form.get("layer")):
                raise ValueError("Choose a cell with a mask layer.")
            layer, datatype = (int(part) for part in spec.split("/"))
            start, end = (tuple(number(name, form.get(name)) for name in names)
                          for names in (("x0", "y0"), ("x1", "y1")))
            origin = number("Deck X at start", form.get("origin"))
            polygons = self.cells[form.get("cell")].get_polygons(layer=layer, datatype=datatype)
            length, intervals = section(polygons, start, end, tolerance=self.tolerance)
            if not isfinite(origin + length) or origin + length <= origin:
                raise ValueError("The cut must map to a finite, nonzero span in deck X.")
            if form.get("tone") == "Protected":
                intervals = invert_intervals(intervals, length, self.tolerance)
            mode = form.get("mode")
            options = {}
            if mode == "Mesh lines":
                options["spacing"] = number("Edge mesh spacing", form.get("spacing"), True)
            else:
                options["top"] = number("Etch top", form.get("top"))
                options["bottom"] = number("Etch bottom", form.get("bottom"))
                if options["top"] >= options["bottom"]:
                    raise ValueError("Etch bottom Y must be greater than top Y.")
                options["material"] = form.get("material")
            commands = mask_commands(intervals, length, origin, mode=mode,
                                     tolerance=self.tolerance, **options)
            self.result = origin, length, intervals
            self.draw_cut()
            if intervals:
                self.copy.state(["!disabled"])
            advice = ("Edge lines merge into the starting X mesh." if mode == "Mesh lines"
                      else "Insert after depositing the material to etch.")
            if not commands:
                advice = "This cut has no interior mask edges." if intervals else "This cut has no openings."
            self.show_preview(commands, f"{len(intervals)} opening(s) over {fmt(length)} µm. {advice}")
        except (ValueError, RuntimeError, OverflowError) as error:
            self.show_preview("", str(error))

    def draw_cut(self):
        canvas = self.cross_section
        canvas.delete("all")
        if not self.result:
            return
        origin, length, intervals = self.result
        left, right = 12, max(13, canvas.winfo_width() - 12)
        canvas.create_line(left, 14, right, 14, fill=MUTED, width=2)
        for lo, hi in intervals:
            start = left + (right - left) * lo / length
            end = left + (right - left) * hi / length
            canvas.create_rectangle(start, 8, end, 20, fill=ACCENT, outline="")
        canvas.create_text(left, 27, anchor="nw", fill=FG, text=f"{fmt(origin)} µm")
        canvas.create_text(right, 27, anchor="ne", fill=FG, text=f"{fmt(origin + length)} µm")

    def copy_intervals(self):
        origin, _, intervals = self.result
        content = "left_um,right_um\n" + "".join(f"{fmt(origin + lo)},{fmt(origin + hi)}\n"
                                                  for lo, hi in intervals)
        self.app.root.clipboard_clear()
        self.app.root.clipboard_append(content)
        self.notice.configure(text="Copied the opening intervals in deck X coordinates.")

    def insert(self):
        if self.form.get("mode") != "Mesh lines":
            super().insert()
            return
        try:
            deck = self.app.deck_text()
            merged = merge_mesh_lines(deck, self.preview.get("1.0", "end-1c"))
        except ValueError as error:
            self.notice.configure(text=str(error))
            return
        if merged != deck:
            self.app.editor.edit_separator()
            self.app.editor.replace("1.0", "end-1c", merged)
            self.app.editor.edit_separator()
        self.app.status.set("Merged GDS edge lines into the starting mesh." if merged != deck
                            else "The starting mesh already has these edge locations.")
        self.window.destroy()
