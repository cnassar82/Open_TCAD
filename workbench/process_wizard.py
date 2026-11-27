"""Guided SUPREM4 process cards with a live deck preview."""

import tkinter as tk
from tkinter import ttk

from suprem_cards import BOOLEAN, CHOICE, catalog, key_file, numeric_literal, parameter_token, polygon_cards
from ui_theme import BG, code_widget

OPERATIONS = ("Deposit", "Implant", "Diffuse", "Oxidize", "Etch")
MODEL_CARDS = (
    "method", "oxide", "material", "vacancy", "interstitia", "arsenic",
    "phosphorus", "antimony", "boron", "gold", "cesium", "trap", "beryllium",
    "carbon", "germanium", "magnesium", "selenium", "isilicon", "tin", "zinc", "generic",
)
ETCH_MODES = {
    "all": "Remove all", "dry": "Directional depth", "left": "Left of a line",
    "right": "Right of a line", "physical": "Material rates", "polygon": "Polygon",
    "file": "Contour file",
}
LABELS = {
    "thick": "Thickness (µm)", "concentrati": "Dopant concentration (cm⁻³)",
    "divisions": "Mesh divisions", "space": "Surface spacing (µm)", "spacing": "Spacing (µm)",
    "square": "Square film edges", "temp": "Temperature (°C)", "time": "Time",
    "pressure": "Pressure", "ambient": "Ambient", "solid.sol": "Solid solubility",
    "gas.conc": "Gas concentration", "file": "Contour file", "defectm": "Diffusion model",
    "oxidem": "Oxidation model", "preconditio": "Preconditioner", "blkmeth": "Block solver",
    "mater": "Second material", "orient": "Crystal orientation", "dryo2": "Dry O₂",
    "weto2": "Wet O₂", "ammonia": "Ammonia (inert in this build)",
    "dose": "Dose (cm⁻²)", "energy": "Energy (keV)", "angle": "Angle from vertical (°)",
    "model": "Distribution", "range": "Projected range (µm)",
}
MODEL_LABELS = {
    "fermi": "Fermi level", "two.dim": "Full 2D", "steady": "2D steady state",
    "full.cpl": "Fully coupled", "vertical": "Vertical growth", "compress": "Compressible oxide",
    "viscous": "Viscous oxide", "erfc": "Analytic ERFC",
}


class CardForm:
    def __init__(self, parent, changed):
        self.frame = parent
        self.changed = changed
        self.values = {}
        self.row = 0

    def reset_rows(self):
        for child in self.frame.winfo_children():
            child.destroy()
        self.row = 0

    def variable(self, name, default=""):
        if name not in self.values:
            value = tk.StringVar(self.frame, value=default)
            value.trace_add("write", self.changed)
            self.values[name] = value
        return self.values[name]

    def get(self, name):
        return self.values[name].get().strip() if name in self.values else ""

    def field(self, name, label, hint="", choices=None, default=""):
        variable = self.variable(name, default)
        ttk.Label(self.frame, text=label).grid(row=self.row, column=0, sticky="w", pady=4)
        if choices is None:
            control = ttk.Entry(self.frame, textvariable=variable, width=24)
        else:
            control = ttk.Combobox(self.frame, textvariable=variable, values=choices,
                                   state="readonly", width=23)
        control.grid(row=self.row, column=1, sticky="ew", padx=10, pady=4)
        ttk.Label(self.frame, text=hint, style="Muted.TLabel", wraplength=355).grid(
            row=self.row, column=2, sticky="w", pady=4)
        self.row += 1
        return control

    def parameter(self, item, *, label=None, hint=None, options=None, default=""):
        label = label or LABELS.get(item.name, item.name)
        hint = item.description if hint is None else hint
        if item.kind == CHOICE:
            children = [child for child in item.children if options is None or child.name in options]
            menu = [""] + [self.choice_text(child.name) for child in children]
            self.field(item.name, label, hint, choices=menu, default=self.choice_text(default))
        elif item.kind == BOOLEAN:
            self.field(item.name, label, hint, choices=("", "On", "Off"), default=default)
        else:
            self.field(item.name, label, hint, default=default)

    @staticmethod
    def choice_text(name):
        label = MODEL_LABELS.get(name) or (LABELS.get(name) if name in {"dryo2", "weto2", "ammonia"} else None)
        return f"{name} · {label}" if label else name

    def choice(self, name):
        return self.get(name).partition(" · ")[0]

    def tokens(self, parameters):
        tokens = []
        for item in parameters:
            if item.kind == CHOICE:
                picked = self.choice(item.name)
                if picked:
                    tokens.append(picked)
                    tokens.extend(self.tokens(item[picked].children))
            elif value := self.get(item.name):
                tokens.append(parameter_token(item, value))
                if item.kind == BOOLEAN and value == "On":
                    tokens.extend(self.tokens(item.children))
        return tokens


class ProcessWizard:
    def __init__(self, app, operation):
        self.app = app
        self.operation = operation
        self.cards = catalog(key_file())
        required = {"deposit", "diffuse", "etch", "method", self.card_name(operation)}
        if missing := required - self.cards.keys():
            raise ValueError("SUPREM4 catalog is missing " + ", ".join(sorted(missing)))
        self.card = self.cards[self.card_name(operation)]
        self.model_lines = []
        self.model_forms = {}
        self.etch_fields = []
        self.polygon = None
        self.points = ""
        self.window = tk.Toplevel(app.root)
        self.window.title(f"{operation} wizard")
        self.window.transient(app.root)
        self.window.geometry("860x710")
        self.window.minsize(720, 530)
        self.window.grab_set()
        body = ttk.Frame(self.window, padding=16)
        body.pack(fill="both", expand=True)
        ttk.Label(body, text=operation, style="Title.TLabel").pack(anchor="w")
        self.tabs = ttk.Notebook(body)
        self.tabs.pack(fill="both", expand=True, pady=(8, 10))

        step = ttk.Frame(self.tabs, padding=8)
        self.tabs.add(step, text="Process step")
        self.step = CardForm(self.scroll_form(step), self.update_preview)
        self.build_step()

        models = ttk.Frame(self.tabs, padding=8)
        self.tabs.add(models, text="Models and coefficients")
        heading = ttk.Frame(models)
        heading.pack(fill="x", pady=(0, 8))
        names = [name for name in MODEL_CARDS if name in self.cards]
        ttk.Label(heading, text="Model card").pack(side="left", padx=(0, 8))
        self.model_choice = ttk.Combobox(heading, state="readonly", values=names, width=17)
        self.model_choice.set("method")
        self.model_choice.pack(side="left")
        self.model_choice.bind("<<ComboboxSelected>>", self.build_model)
        self.add_model = ttk.Button(heading, text="Queue model card", command=self.queue_model)
        self.add_model.pack(side="right")
        self.remove_model = ttk.Button(heading, text="Remove last queued", command=self.unqueue_model)
        self.remove_model.pack(side="right", padx=6)
        self.model_box = self.scroll_form(models)
        self.build_model()

        ttk.Label(body, text="Deck preview", style="Section.TLabel").pack(anchor="w")
        self.preview = code_widget(body, height=6, wrap="word", state="disabled", padx=10, pady=8)
        self.preview.pack(fill="x", pady=(5, 7))
        footer = ttk.Frame(body)
        footer.pack(fill="x")
        self.notice = ttk.Label(footer, style="Muted.TLabel", wraplength=530)
        self.notice.pack(side="left")
        self.insert = ttk.Button(footer, text="Insert into deck", style="Primary.TButton", command=self.finish)
        self.insert.pack(side="right")
        ttk.Button(footer, text="Cancel", command=self.window.destroy).pack(side="right", padx=7)
        self.window.bind("<Escape>", lambda _: self.window.destroy())
        self.update_preview()

    @staticmethod
    def card_name(operation):
        return "diffuse" if operation == "Oxidize" else operation.lower()

    @staticmethod
    def scroll_form(parent):
        region = ttk.Frame(parent)
        region.pack(fill="both", expand=True)
        canvas = tk.Canvas(region, bg=BG, highlightthickness=0, height=210)
        canvas.pack(side="left", fill="both", expand=True)
        scroll = ttk.Scrollbar(region, command=canvas.yview)
        scroll.pack(side="right", fill="y")
        canvas.configure(yscrollcommand=scroll.set)
        frame = ttk.Frame(canvas, padding=(6, 2))
        frame.columnconfigure(2, weight=1)
        item = canvas.create_window(0, 0, anchor="nw", window=frame)
        canvas.bind("<Configure>", lambda event: canvas.itemconfigure(item, width=event.width))
        frame.bind("<Configure>", lambda _: canvas.configure(scrollregion=canvas.bbox("all")))
        for direction, movement in ((4, -4), (5, 4)):
            canvas.bind(f"<Button-{direction}>", lambda _, n=movement: canvas.yview_scroll(n, "units"))
            frame.bind(f"<Button-{direction}>", lambda _, n=movement: canvas.yview_scroll(n, "units"))
        return frame

    def build_step(self):
        form = self.step
        if self.operation == "Etch":
            form.field("mode", "Etch mode", "Choose how SUPREM4 defines the etch boundary.",
                       choices=tuple(ETCH_MODES.values()), default=ETCH_MODES["all"])
            self.etch_body = CardForm(ttk.Frame(form.frame), self.update_preview)
            self.etch_body.frame.columnconfigure(2, weight=1)
            self.etch_body.frame.grid(row=1, column=0, columnspan=3, sticky="ew")
            form.values["mode"].trace_add("write", self.build_etch)
            self.build_etch()
            return

        if self.operation in {"Deposit", "Implant"}:
            defaults = ({"material": "oxide", "thick": "0.1", "divisions": "3"}
                        if self.operation == "Deposit" else
                        {"impurity": "boron", "dose": "3e12", "energy": "30", "model": "pearson"})
            for name, value in defaults.items():
                form.parameter(self.card[name], default=value)
            for item in self.card.children:
                if item.name not in defaults:
                    form.parameter(item)
            return

        form.parameter(self.card["temp"], hint="Furnace temperature in °C.", default="900")
        form.parameter(self.card["time"], label="Time (min)", default="10")
        oxidants = {"dryo2", "weto2"}
        options = oxidants if self.operation == "Oxidize" else {item.name for item in self.card["ambient"].children} - oxidants
        form.parameter(self.card["ambient"], options=options,
                       default="dryo2" if self.operation == "Oxidize" else "nitrogen")
        for item in self.card.children:
            if item.name not in {"temp", "time", "ambient"}:
                hint = None
                if item.name in {"ramp", "start", "stop"}:
                    hint = "Legacy card field. Use successive fixed-temperature steps to model a ramp."
                form.parameter(item, hint=hint)

    def etch_mode(self):
        return next(mode for mode, label in ETCH_MODES.items() if label == self.step.get("mode"))

    def build_etch(self, *_):
        if self.polygon is not None:
            self.points = self.polygon.get("1.0", "end-1c")
            self.polygon = None
        mode = self.etch_mode()
        form = self.etch_body
        form.reset_rows()
        self.etch_fields = []
        if mode != "physical":
            material = self.card["material"]
            options = {item.name for item in material.children}
            if mode == "all":
                options.discard("gaas")  # The validated key rejects GaAs with the all flag.
            if form.choice("material") and form.choice("material") not in options:
                form.values["material"].set("")
            form.parameter(material, options=options, default="oxide")
            self.etch_fields.append(material)
        if mode == "polygon":
            ttk.Label(form.frame, text="Vertices (µm)").grid(row=form.row, column=0, sticky="nw", pady=5)
            self.polygon = code_widget(form.frame, height=5, width=25, padx=8, pady=6)
            self.polygon.grid(row=form.row, column=1, sticky="ew", padx=10, pady=5)
            ttk.Label(form.frame, text="Enter one X, Y pair per line. SUPREM4 closes the polygon.",
                      style="Muted.TLabel", wraplength=355).grid(row=form.row, column=2, sticky="w")
            self.polygon.insert("1.0", self.points)
            self.polygon.edit_modified(False)
            self.polygon.bind("<<Modified>>", self.polygon_changed)
            form.row += 1
        modes = self.card["type_etch"]
        fields = {
            "file": ["file"], "left": ["p1.x", "p1.y", "p2.x", "p2.y", "spacing"],
            "right": ["p1.x", "p1.y", "p2.x", "p2.y", "spacing"], "polygon": ["spacing"],
        }
        for item in modes.children:
            if item.name == mode:
                self.etch_fields.extend(item.children)
        self.etch_fields.extend(self.card[name] for name in fields.get(mode, ()))
        for item in self.etch_fields:
            if item.name == "material":
                continue
            hint = item.description
            if item.name == "time":
                hint = "Duration in the time units used by the etch rates."
            elif item.name.startswith("r."):
                hint = f"{hint} (µm per unit time)"
            form.parameter(item, label="Depth (µm)" if item.name == "thick" else None, hint=hint)
        self.update_preview()

    def polygon_changed(self, _):
        if self.polygon is not None and self.polygon.edit_modified():
            self.polygon.edit_modified(False)
            self.update_preview()

    def build_model(self, _=None):
        name = self.model_choice.get()
        if name not in self.model_forms:
            self.model_forms[name] = CardForm(self.model_box, self.update_preview)
        self.model = self.model_forms[name]
        self.model.reset_rows()
        for item in self.cards[name].children:
            self.model.parameter(item)
            for child in item.children if item.kind == BOOLEAN else ():
                self.model.parameter(child)
        self.update_preview()

    def model_command(self, name=None):
        name = name or self.model_choice.get()
        tokens = self.model_forms[name].tokens(self.cards[name].children)
        return " ".join((name, *tokens)) if tokens else ""

    def queue_model(self):
        try:
            if command := self.model_command():
                self.model_lines.append(command)
                del self.model_forms[self.model_choice.get()]
                self.build_model()
        except ValueError as error:
            self.notice.configure(text=str(error))

    def unqueue_model(self):
        if self.model_lines:
            self.model_lines.pop()
            self.update_preview()

    def step_command(self):
        form = self.step
        if self.operation == "Deposit":
            if not form.choice("material"):
                raise ValueError("Choose a material to deposit.")
            if not form.get("thick") and not form.get("file"):
                raise ValueError("Enter a thickness or a contour file.")
            if not form.get("file"):
                check_positive(form.get("thick"), "Thickness")
            check_positive(form.get("divisions"), "Mesh divisions")
            params = [item for item in self.card.children if item.name != "thick" or not form.get("file")]
        elif self.operation == "Implant":
            if not form.choice("impurity") or not form.get("dose") or not form.get("energy"):
                raise ValueError("Choose an implant species and enter dose and energy.")
            check_positive(form.get("dose"), "Implant dose")
            check_positive(form.get("energy"), "Implant energy")
            params = self.card.children
        elif self.operation == "Etch":
            mode = self.etch_mode()
            form = self.etch_body
            material = form.choice("material")
            required = {"dry": "thick", "physical": "time", "left": "p1.x", "right": "p1.x", "file": "file"}
            if mode == "all" and not material:
                raise ValueError("Choose which material to remove.")
            if mode in required and not form.get(required[mode]):
                raise ValueError(f"Enter {required[mode]} for this etch.")
            if mode == "polygon":
                return polygon_cards(material, self.polygon.get("1.0", "end-1c"), form.get("spacing"))
            tokens = form.tokens(self.etch_fields)
            if mode == "physical":
                check_positive(form.get("time"), "Etch time")
                rates = [form.get(item.name) for item in self.etch_fields if item.name.startswith("r.")]
                if not any(value and numeric_literal(value) != 0 for value in rates):
                    raise ValueError("Enter at least one positive material etch rate.")
                for value in rates:
                    check_positive(value, "Etch rate", allow_zero=True)
            elif mode == "dry":
                check_positive(form.get("thick"), "Etch depth")
            return " ".join(("etch", *((mode,) if mode != "file" else ()), *tokens))
        else:
            for name in ("temp", "time"):
                if not form.get(name):
                    raise ValueError("Enter both furnace temperature and time.")
            if not form.choice("ambient"):
                raise ValueError("Choose the furnace ambient.")
            check_positive(form.get("time"), "Furnace time")
            check_positive(form.get("pressure"), "Furnace pressure")
            params = self.card.children
        return " ".join((self.card.name, *form.tokens(params)))

    def update_preview(self, *_):
        if not hasattr(self, "preview"):
            return
        self.remove_model.state(["!disabled"] if self.model_lines else ["disabled"])
        cards = list(self.model_lines)
        self.add_model.state(["disabled"])
        try:
            for name in self.model_forms:
                if draft := self.model_command(name):
                    cards.append(draft)
                    if name == self.model_choice.get():
                        self.add_model.state(["!disabled"])
            cards.append(self.step_command())
            self.insert.state(["!disabled"])
            note = "Model cards appear before the process step."
        except ValueError as error:
            self.insert.state(["disabled"])
            note = str(error)
        self.preview.configure(state="normal")
        self.preview.delete("1.0", "end")
        self.preview.insert("1.0", "\n".join(cards))
        self.preview.configure(state="disabled")
        self.notice.configure(text=note)

    def finish(self):
        value = self.preview.get("1.0", "end-1c")
        self.window.destroy()
        self.app.insert_command(value)


def check_positive(value, label, allow_zero=False):
    number = numeric_literal(value)
    if number is not None and (number < 0 or (number == 0 and not allow_zero)):
        qualifier = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{label} must be {qualifier}.")
