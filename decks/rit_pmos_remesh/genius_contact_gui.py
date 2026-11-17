#!/usr/bin/env python3
"""GUI for assigning Genius contacts to a SUPREM .str structure."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from tkinter import filedialog, messagebox, simpledialog, ttk
import tkinter as tk

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.widgets import RectangleSelector
import matplotlib.tri as mtri
import numpy as np

from make_genius_contacts import (
    ContactSpec,
    Node,
    build_edges,
    parse_suprem,
    region_extents,
    suggested_contacts,
    write_deck,
    write_tif,
)


MATERIAL_COLORS = {
    1: "#8ecae6",
    2: "#bde0fe",
    3: "#90be6d",
    4: "#adb5bd",
    6: "#f4a261",
    7: "#ffafcc",
}


class ContactDialog(simpledialog.Dialog):
    def __init__(
        self,
        parent,
        title: str,
        default_name: str = "",
        default_kind: str = "ohmic",
        default_work: float = 5.25,
    ):
        self.default_name = default_name
        self.default_kind = default_kind
        self.default_work = default_work
        self.result: tuple[str, str, float] | None = None
        super().__init__(parent, title)

    def body(self, master):
        self.name_var = tk.StringVar(value=self.default_name)
        self.kind_var = tk.StringVar(value=self.default_kind)
        self.work_var = tk.StringVar(value=f"{self.default_work:g}")

        ttk.Label(master, text="Name").grid(row=0, column=0, sticky="w", padx=4, pady=4)
        name_entry = ttk.Entry(master, textvariable=self.name_var, width=24)
        name_entry.grid(row=0, column=1, sticky="ew", padx=4, pady=4)

        ttk.Label(master, text="Type").grid(row=1, column=0, sticky="w", padx=4, pady=4)
        ttk.Radiobutton(master, text="Ohmic", variable=self.kind_var, value="ohmic").grid(
            row=1, column=1, sticky="w", padx=4, pady=2
        )
        ttk.Radiobutton(master, text="Gate", variable=self.kind_var, value="gate").grid(
            row=2, column=1, sticky="w", padx=4, pady=2
        )

        ttk.Label(master, text="Gate work function").grid(
            row=3, column=0, sticky="w", padx=4, pady=4
        )
        ttk.Entry(master, textvariable=self.work_var, width=10).grid(
            row=3, column=1, sticky="w", padx=4, pady=4
        )
        master.columnconfigure(1, weight=1)
        return name_entry

    def validate(self) -> bool:
        name = self.name_var.get().strip()
        if not name:
            messagebox.showerror("Contact", "Contact name is required.", parent=self)
            return False
        try:
            work = float(self.work_var.get())
        except ValueError:
            messagebox.showerror("Contact", "Work function must be numeric.", parent=self)
            return False
        self.result = (name, self.kind_var.get(), work)
        return True


class ContactGui:
    def __init__(self, root: tk.Tk, path: Path):
        self.root = root
        self.path = path
        self.nodes, self.tris, self.regions, self.sol_codes, self.sol_rows = parse_suprem(
            path, auto_name_electrodes=False
        )
        self.contacts: list[ContactSpec] = suggested_contacts(self.nodes, self.tris, self.regions)
        self.selected_region: int | None = None
        self.selected_edges: set[int] = set()
        self.edge_lookup = self._edge_lookup()
        self.region_extents = region_extents(self.nodes, self.tris)
        self.mode_var = tk.StringVar(value="region")
        self.mesh_var = tk.BooleanVar(value=False)
        self.y_stretch_var = tk.StringVar(value="1")
        self.status = tk.StringVar(value="")
        self.selector: RectangleSelector | None = None
        self.highlight_artists = []

        self.root.title(f"Genius Contact GUI - {path.name}")
        self._build_widgets()
        self.redraw()
        self.refresh_contacts()

    def _edge_lookup(self):
        boundary_edges, edge_number, edge_regions = build_edges(self.tris)
        return {
            edge_number[edge]: {
                "nodes": edge,
                "regions": edge_regions[edge],
                "mid": (
                    0.5 * (self.nodes[edge[0]].x + self.nodes[edge[1]].x),
                    0.5 * (self.nodes[edge[0]].y + self.nodes[edge[1]].y),
                ),
            }
            for edge in boundary_edges
        }

    def _build_widgets(self):
        top = ttk.Frame(self.root, padding=6)
        top.pack(side=tk.TOP, fill=tk.X)
        left = ttk.Frame(self.root, padding=(6, 0, 6, 6))
        left.pack(side=tk.LEFT, fill=tk.Y)
        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        ttk.Radiobutton(top, text="Region", variable=self.mode_var, value="region").pack(
            side=tk.LEFT
        )
        ttk.Radiobutton(top, text="Edge", variable=self.mode_var, value="edge").pack(
            side=tk.LEFT, padx=(6, 0)
        )
        ttk.Radiobutton(top, text="Box edges", variable=self.mode_var, value="box").pack(
            side=tk.LEFT, padx=(6, 12)
        )
        ttk.Checkbutton(top, text="mesh", variable=self.mesh_var, command=self.redraw).pack(
            side=tk.LEFT
        )
        ttk.Label(top, text="Y stretch").pack(side=tk.LEFT, padx=(12, 4))
        ttk.Entry(top, textvariable=self.y_stretch_var, width=7).pack(side=tk.LEFT)
        ttk.Button(top, text="Redraw", command=self.redraw).pack(side=tk.LEFT, padx=8)
        ttk.Button(top, text="Clear selection", command=self.clear_selection).pack(
            side=tk.LEFT
        )

        ttk.Label(left, text="Contacts").pack(anchor="w")
        columns = ("name", "type", "target")
        self.tree = ttk.Treeview(left, columns=columns, show="headings", height=14)
        for column, width in (("name", 110), ("type", 80), ("target", 120)):
            self.tree.heading(column, text=column)
            self.tree.column(column, width=width, anchor=tk.W)
        self.tree.pack(fill=tk.BOTH, expand=True, pady=(4, 6))

        ttk.Button(left, text="Add from selection", command=self.add_contact).pack(
            fill=tk.X
        )
        ttk.Button(left, text="Edit contact", command=self.edit_contact).pack(
            fill=tk.X, pady=(4, 0)
        )
        ttk.Button(left, text="Delete contact", command=self.delete_contact).pack(
            fill=tk.X, pady=(4, 0)
        )
        ttk.Button(left, text="Save TIF + deck", command=self.save_outputs).pack(
            fill=tk.X, pady=(14, 0)
        )
        ttk.Label(left, textvariable=self.status, wraplength=310, anchor="w").pack(
            fill=tk.X, pady=(10, 0)
        )

        self.fig = Figure(figsize=(10, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()
        self.canvas.mpl_connect("button_press_event", self.on_click)
        self.mode_var.trace_add("write", lambda *_: self.configure_selector())
        self.tree.bind("<Double-1>", lambda _event: self.edit_contact())

    def points_array(self):
        max_idx = max(self.nodes)
        points = np.zeros((max_idx, 2), dtype=float)
        for idx, node in self.nodes.items():
            points[idx - 1] = (node.x, node.y)
        return points

    def triangles_array(self):
        return np.array([[node - 1 for node in tri.nodes] for tri in self.tris], dtype=int)

    def material_values(self):
        material_ids = sorted({region.mat_id for region in self.regions.values()})
        mat_to_value = {mat: i for i, mat in enumerate(material_ids)}
        return np.array(
            [mat_to_value[self.regions[tri.region].mat_id] for tri in self.tris], dtype=float
        )

    def redraw(self):
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        points = self.points_array()
        triangles = self.triangles_array()
        self.triangulation = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
        values = self.material_values()
        material_count = max(len(set(values)), 1)
        cmap = matplotlib.colormaps["tab20"].resampled(material_count)
        self.ax.tripcolor(
            self.triangulation,
            facecolors=values,
            shading="flat",
            edgecolors="k" if self.mesh_var.get() else "none",
            linewidth=0.12,
            cmap=cmap,
        )
        self.ax.set_title(f"{self.path.name}: assign contacts")
        self.ax.set_xlabel("x [um]")
        self.ax.set_ylabel("y [um]")
        try:
            y_stretch = float(self.y_stretch_var.get())
            if y_stretch <= 0:
                raise ValueError
        except ValueError:
            y_stretch = 1.0
            self.y_stretch_var.set("1")
        self.ax.set_aspect(y_stretch, adjustable="datalim", anchor="C")
        self.ax.invert_yaxis()
        self.draw_contact_overlays()
        self.draw_selection()
        self.fig.tight_layout()
        self.canvas.draw_idle()
        self.configure_selector()

    def configure_selector(self):
        if self.selector is not None:
            self.selector.set_active(False)
        if self.mode_var.get() == "box":
            self.selector = RectangleSelector(
                self.ax,
                self.on_box_select,
                useblit=True,
                button=[1],
                minspanx=0,
                minspany=0,
                spancoords="data",
                interactive=False,
            )
        self.status.set("Mode: select a region, edge, or drag a box around edges.")

    def draw_contact_overlays(self):
        for contact in self.contacts:
            if contact.target == "region" and contact.region in self.region_extents:
                xmin, xmax, ymin, ymax = self.region_extents[contact.region]
                rect = Rectangle(
                    (xmin, ymin),
                    xmax - xmin,
                    ymax - ymin,
                    fill=False,
                    edgecolor="black",
                    linewidth=1.4,
                    linestyle="--",
                )
                self.ax.add_patch(rect)
                self.ax.text(
                    0.5 * (xmin + xmax),
                    0.5 * (ymin + ymax),
                    contact.name,
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="black",
                    bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
                )
            elif contact.target == "boundary":
                self.plot_edges(set(contact.edges or []), color="#d00000", linewidth=2.0)
                mids = [self.edge_lookup[e]["mid"] for e in contact.edges or [] if e in self.edge_lookup]
                if mids:
                    x = sum(p[0] for p in mids) / len(mids)
                    y = sum(p[1] for p in mids) / len(mids)
                    self.ax.text(
                        x,
                        y,
                        contact.name,
                        ha="center",
                        va="center",
                        fontsize=9,
                        color="#d00000",
                        bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
                    )

    def draw_selection(self):
        if self.selected_region is not None and self.selected_region in self.region_extents:
            xmin, xmax, ymin, ymax = self.region_extents[self.selected_region]
            self.ax.add_patch(
                Rectangle(
                    (xmin, ymin),
                    xmax - xmin,
                    ymax - ymin,
                    fill=False,
                    edgecolor="#0057ff",
                    linewidth=2.0,
                )
            )
        self.plot_edges(self.selected_edges, color="#0057ff", linewidth=2.4)

    def plot_edges(self, edges: set[int], color: str, linewidth: float):
        for edge_no in edges:
            item = self.edge_lookup.get(edge_no)
            if not item:
                continue
            a, b = item["nodes"]
            na, nb = self.nodes[a], self.nodes[b]
            self.ax.plot([na.x, nb.x], [na.y, nb.y], color=color, linewidth=linewidth)

    def on_click(self, event):
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return
        if self.mode_var.get() == "region":
            finder = self.triangulation.get_trifinder()
            tri_idx = int(finder(event.xdata, event.ydata))
            if tri_idx < 0:
                return
            self.selected_region = self.tris[tri_idx].region
            self.selected_edges.clear()
            region = self.regions[self.selected_region]
            self.status.set(
                f"Selected region {self.selected_region}: {region.name}, material {region.mat_id}"
            )
        elif self.mode_var.get() == "edge":
            edge = self.nearest_edge(event.xdata, event.ydata)
            if edge is None:
                return
            if edge in self.selected_edges:
                self.selected_edges.remove(edge)
            else:
                self.selected_edges.add(edge)
            self.selected_region = None
            self.status.set(f"Selected {len(self.selected_edges)} edge(s).")
        self.redraw()

    def nearest_edge(self, x: float, y: float) -> int | None:
        xvals = [node.x for node in self.nodes.values()]
        yvals = [node.y for node in self.nodes.values()]
        tol = 0.02 * max(max(xvals) - min(xvals), max(yvals) - min(yvals))
        best_edge = None
        best_distance = float("inf")
        for edge_no, item in self.edge_lookup.items():
            a, b = item["nodes"]
            distance = segment_distance(x, y, self.nodes[a], self.nodes[b])
            if distance < best_distance:
                best_distance = distance
                best_edge = edge_no
        return best_edge if best_distance <= tol else None

    def on_box_select(self, eclick, erelease):
        if None in (eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata):
            return
        xmin, xmax = sorted((eclick.xdata, erelease.xdata))
        ymin, ymax = sorted((eclick.ydata, erelease.ydata))
        for edge_no, item in self.edge_lookup.items():
            x, y = item["mid"]
            if xmin <= x <= xmax and ymin <= y <= ymax:
                self.selected_edges.add(edge_no)
        self.selected_region = None
        self.status.set(f"Selected {len(self.selected_edges)} edge(s).")
        self.redraw()

    def add_contact(self):
        if self.selected_region is None and not self.selected_edges:
            messagebox.showinfo("Contact", "Select a region or one or more edges first.")
            return
        default_name = ""
        default_kind = "ohmic"
        if self.selected_region is not None:
            region = self.regions[self.selected_region]
            default_name = region.name
            if region.mat_id in {4, 6}:
                default_name = ""
        dialog = ContactDialog(self.root, "Add Contact", default_name, default_kind)
        if dialog.result is None:
            return
        name, kind, work = dialog.result
        if any(contact.name == name for contact in self.contacts):
            messagebox.showerror("Contact", f"Contact name {name!r} already exists.")
            return
        if self.selected_region is not None:
            self.contacts = [
                contact
                for contact in self.contacts
                if not (contact.target == "region" and contact.region == self.selected_region)
            ]
            self.contacts.append(
                ContactSpec(
                    name=name,
                    kind=kind,
                    target="region",
                    region=self.selected_region,
                    work_function=work,
                )
            )
        else:
            used = {
                edge
                for contact in self.contacts
                if contact.target == "boundary"
                for edge in (contact.edges or [])
            }
            overlap = used & self.selected_edges
            if overlap:
                messagebox.showerror(
                    "Contact",
                    f"{len(overlap)} selected edge(s) already belong to another contact.",
                )
                return
            self.contacts.append(
                ContactSpec(
                    name=name,
                    kind=kind,
                    target="boundary",
                    edges=sorted(self.selected_edges),
                    work_function=work,
                )
            )
        self.clear_selection()
        self.refresh_contacts()

    def delete_contact(self):
        selected = self.tree.selection()
        if not selected:
            return
        indices = sorted((int(item) for item in selected), reverse=True)
        for index in indices:
            if 0 <= index < len(self.contacts):
                del self.contacts[index]
        self.refresh_contacts()
        self.redraw()

    def edit_contact(self):
        selected = self.tree.selection()
        if len(selected) != 1:
            return
        index = int(selected[0])
        if not (0 <= index < len(self.contacts)):
            return
        contact = self.contacts[index]
        dialog = ContactDialog(
            self.root,
            "Edit Contact",
            default_name=contact.name,
            default_kind=contact.kind,
            default_work=contact.work_function,
        )
        if dialog.result is None:
            return
        name, kind, work = dialog.result
        if any(i != index and item.name == name for i, item in enumerate(self.contacts)):
            messagebox.showerror("Contact", f"Contact name {name!r} already exists.")
            return
        contact.name = name
        contact.kind = kind
        contact.work_function = work
        self.refresh_contacts()
        self.redraw()

    def refresh_contacts(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        for i, contact in enumerate(self.contacts):
            if contact.target == "region":
                target = f"region {contact.region}"
            else:
                target = f"{len(contact.edges or [])} edges"
            self.tree.insert("", tk.END, iid=str(i), values=(contact.name, contact.kind, target))
        self.status.set(f"{len(self.contacts)} contact(s) defined.")

    def clear_selection(self):
        self.selected_region = None
        self.selected_edges.clear()
        self.redraw()

    def save_outputs(self):
        if not self.contacts:
            if not messagebox.askyesno("Save", "No contacts are defined. Save anyway?"):
                return
        default_tif = self.path.with_name(f"{self.path.stem}_contacts.tif")
        tif_name = filedialog.asksaveasfilename(
            title="Save TIF",
            initialfile=default_tif.name,
            initialdir=str(default_tif.parent),
            defaultextension=".tif",
            filetypes=[("TIF", "*.tif"), ("All files", "*.*")],
        )
        if not tif_name:
            return
        tif_path = Path(tif_name)
        default_deck = tif_path.with_suffix(".in")
        deck_name = filedialog.asksaveasfilename(
            title="Save Genius deck",
            initialfile=default_deck.name,
            initialdir=str(default_deck.parent),
            defaultextension=".in",
            filetypes=[("Genius input", "*.in"), ("All files", "*.*")],
        )
        if not deck_name:
            return
        try:
            counts = write_tif(
                tif_path,
                self.nodes,
                self.tris,
                self.regions,
                self.sol_codes,
                self.sol_rows,
                body_tol=1e-6,
                contacts=self.contacts,
            )
            write_deck(Path(deck_name), tif_path, self.contacts)
        except Exception as exc:
            messagebox.showerror("Save failed", str(exc))
            return
        messagebox.showinfo(
            "Saved",
            f"Wrote {tif_path.name} and {Path(deck_name).name}\n"
            f"{counts['nodes']} nodes, {counts['triangles']} triangles, "
            f"{counts['edges']} boundary/interface edges.",
        )


def segment_distance(x: float, y: float, a: Node, b: Node) -> float:
    px, py = x, y
    ax, ay = a.x, a.y
    bx, by = b.x, b.y
    dx = bx - ax
    dy = by - ay
    denom = dx * dx + dy * dy
    if denom <= 0:
        return math.hypot(px - ax, py - ay)
    t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / denom))
    qx = ax + t * dx
    qy = ay + t * dy
    return math.hypot(px - qx, py - qy)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure", nargs="?", type=Path, help="SUPREM .str file")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = tk.Tk()
    path = args.structure
    if path is None:
        selected = filedialog.askopenfilename(
            parent=root,
            title="Open SUPREM structure",
            filetypes=[("SUPREM structure", "*.str"), ("All files", "*.*")],
        )
        if not selected:
            return 0
        path = Path(selected)
    ContactGui(root, path)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
