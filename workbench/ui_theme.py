"""Shared desktop styles for the workbench and its dialogs."""

import tkinter as tk
from tkinter import font, ttk

BG = "#f6f8fb"
PANEL = "#ffffff"
TINT = "#eef2f7"
EDGE = "#d4dce7"
FG = "#182538"
MUTED = "#53647a"
ACCENT = "#2356a8"
ACCENT_HOVER = "#194687"
SELECTION = "#dbe9fd"
HEADER = "#142338"
HEADER_FG = "#f6f9fe"
HEADER_MUTED = "#c0cde0"

CODE_BG = "#152133"
CODE_FG = "#e3ebf6"
CODE_GUTTER = "#101b2b"
CODE_MUTED = "#9eafc4"
CODE_SELECTION = "#36577f"
SYNTAX = {
    "report": "#9dd8e1",
    "comment": CODE_MUTED,
    "command": "#9cc9ff",
    "key": "#cbd6e9",
    "number": "#f3c58e",
    "factor": "#8fdec8",
}
SUCCESS = "#17624b"
WARNING = "#805218"
ERROR = "#a12b32"
CURVES = ("#2456a8", "#087e72", "#aa531d", "#7449a9", "#b43a62", "#52667a")


def apply_theme(root):
    theme = ttk.Style(root)
    theme.theme_use("clam")
    families = set(font.families(root))
    family = next((name for name in ("Inter", "Aptos", "Noto Sans", "DejaVu Sans") if name in families),
                  "TkDefaultFont")
    mono = next((name for name in ("JetBrains Mono", "Liberation Mono", "DejaVu Sans Mono") if name in families),
                "TkFixedFont")
    for named in ("TkDefaultFont", "TkTextFont", "TkMenuFont", "TkHeadingFont"):
        face = font.nametofont(named, root)
        face.configure(family=family if family != "TkDefaultFont" else face.actual("family"), size=10)
    fixed = font.nametofont("TkFixedFont", root)
    fixed.configure(family=mono if mono != "TkFixedFont" else fixed.actual("family"), size=10)
    family = font.nametofont("TkDefaultFont", root).actual("family")

    root.configure(background=BG)
    theme.configure(".", font="TkDefaultFont", background=BG, foreground=FG, troughcolor=TINT,
                    bordercolor=EDGE, lightcolor=EDGE, darkcolor=EDGE, fieldbackground=PANEL,
                    selectbackground=SELECTION, selectforeground=FG, insertcolor=FG, arrowcolor=MUTED)
    theme.configure("Header.TFrame", background=HEADER)
    theme.configure("Header.TLabel", background=HEADER, foreground=HEADER_FG,
                    font=(family, 19, "bold"))
    theme.configure("HeaderMeta.TLabel", background=HEADER, foreground=HEADER_MUTED)
    theme.configure("Panel.TFrame", background=PANEL, bordercolor=EDGE, borderwidth=1, relief="solid")
    theme.configure("Surface.TFrame", background=PANEL)
    theme.configure("Toolbar.TFrame", background=TINT, bordercolor=EDGE, borderwidth=1, relief="solid")
    theme.configure("Muted.TLabel", foreground=MUTED)
    theme.configure("PanelMuted.TLabel", background=PANEL, foreground=MUTED)
    theme.configure("Toolbar.TLabel", background=TINT, foreground=MUTED)
    theme.configure("Title.TLabel", font=(family, 17, "bold"))
    theme.configure("Section.TLabel", font=(family, 11, "bold"))
    theme.configure("PanelHeading.TLabel", background=PANEL, font=(family, 11, "bold"))

    for name in ("TButton", "TMenubutton"):
        theme.configure(name, background=PANEL, padding=(10, 6), borderwidth=1, relief="flat",
                        focusthickness=1, focuscolor=ACCENT)
        theme.map(name, background=[("disabled", TINT), ("pressed", SELECTION), ("active", TINT)],
                  foreground=[("disabled", MUTED)], bordercolor=[("focus", ACCENT)])
    theme.configure("Primary.TButton", background=ACCENT, foreground="#ffffff",
                    bordercolor=ACCENT, lightcolor=ACCENT, darkcolor=ACCENT, focuscolor="#ffffff")
    theme.map("Primary.TButton",
              background=[("disabled", TINT), ("pressed", ACCENT_HOVER), ("active", ACCENT_HOVER)],
              foreground=[("disabled", MUTED), ("!disabled", "#ffffff")])
    theme.configure("Header.TButton", background="#233851", foreground=HEADER_FG,
                    bordercolor="#4a607a", lightcolor="#4a607a", darkcolor="#4a607a",
                    focuscolor=HEADER_FG)
    theme.map("Header.TButton", background=[("pressed", "#354f6f"), ("active", "#304762")])
    theme.configure("TButton", width=0)

    theme.configure("TEntry", padding=(6, 4))
    theme.configure("TCombobox", padding=(6, 4), background=PANEL, arrowsize=12)
    for name in ("TEntry", "TCombobox"):
        theme.map(name, bordercolor=[("focus", ACCENT)], fieldbackground=[("disabled", TINT)])
    theme.map("TCombobox", foreground=[("readonly", FG)],
              fieldbackground=[("readonly", PANEL), ("disabled", TINT)])
    root.option_add("*TCombobox*Listbox.background", PANEL)
    root.option_add("*TCombobox*Listbox.foreground", FG)
    root.option_add("*TCombobox*Listbox.selectBackground", SELECTION)
    root.option_add("*TCombobox*Listbox.selectForeground", FG)

    theme.configure("TNotebook", background=BG, borderwidth=0, tabmargins=(0, 0, 0, 0))
    theme.configure("TNotebook.Tab", padding=(16, 8), background=TINT, foreground=MUTED)
    theme.map("TNotebook.Tab", background=[("selected", PANEL), ("active", SELECTION)],
              foreground=[("selected", ACCENT), ("active", FG)],
              bordercolor=[("selected", EDGE)], expand=[("selected", (0, 0, 0, 0))])
    theme.configure("Treeview", background=PANEL, fieldbackground=PANEL, rowheight=28,
                    borderwidth=0, indent=17)
    theme.configure("Treeview.Heading", background=TINT, foreground=FG, padding=(7, 7),
                    font=(family, 10, "bold"), relief="flat")
    theme.map("Treeview", background=[("selected", SELECTION)], foreground=[("selected", FG)])
    theme.map("Treeview.Heading", background=[("active", SELECTION)])
    for direction in ("Vertical", "Horizontal"):
        theme.configure(f"{direction}.TScrollbar", width=12, arrowsize=10, background="#c8d3e2",
                        troughcolor=TINT, bordercolor=EDGE, lightcolor=TINT, darkcolor=TINT)
        theme.map(f"{direction}.TScrollbar", background=[("active", "#aebfd3")])
    theme.configure("TPanedwindow", sashwidth=9)
    theme.configure("TSeparator", background=EDGE)


def menu_widget(parent):
    return tk.Menu(parent, tearoff=False, background=PANEL, foreground=FG,
                   activebackground=SELECTION, activeforeground=FG,
                   borderwidth=1, relief="solid", font="TkMenuFont")


def code_widget(parent, *, wrap="none", **kwargs):
    return tk.Text(
        parent, background=CODE_BG, foreground=CODE_FG, insertbackground=CODE_FG,
        selectbackground=CODE_SELECTION, selectforeground=CODE_FG, relief="flat",
        highlightthickness=0, font="TkFixedFont", wrap=wrap, **kwargs,
    )


def style_table(tree):
    tree.tag_configure("alternate", background="#f1f5fa")
