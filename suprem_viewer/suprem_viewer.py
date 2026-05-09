#!/usr/bin/env python3
"""
View SUPREM-IV.GS .str and VTK .vtu structure files directly.

The parser handles the Stanford/SUPREM convention where c/t/r records are
one-based, while n solution records are zero-based. VTU support covers
GENIUS-style compressed appended-data unstructured grids with scalar point
and cell fields.
"""

from __future__ import annotations

import argparse
import base64
import csv
import math
import struct
import sys
import xml.etree.ElementTree as ET
import zlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np


MATERIALS = {
    0: "Vacuum",
    1: "SiO2",
    2: "Si3N4",
    3: "Si",
    4: "Poly",
    5: "OxNi",
    6: "Al",
    7: "Photoresist",
    8: "GaAs",
    9: "SiC6H",
    10: "SiC4H",
    11: "SiC3C",
}

DOPANTS = {
    0: "Vacancy",
    1: "Interstitial",
    2: "Arsenic",
    3: "Phosphorus",
    4: "Antimony",
    5: "Boron",
    6: "BF2",
    7: "Aluminum",
    8: "Nitrogen",
    10: "O2",
    11: "H2O",
    12: "Trap",
    13: "Gold",
    14: "Potential",
    15: "Sxx",
    16: "Syy",
    17: "Sxy",
    18: "Cs",
    19: "DELA",
    20: "ArsenicActive",
    21: "PhosphorusActive",
    22: "AntimonyActive",
    23: "BoronActive",
    24: "AluminumActive",
    25: "NitrogenActive",
    31: "Beryllium",
    32: "BerylliumActive",
    33: "Magnesium",
    34: "MagnesiumActive",
    35: "Selenium",
    36: "SeleniumActive",
    37: "Silicon",
    38: "SiliconActive",
    39: "Tin",
    40: "TinActive",
    41: "Germanium",
    42: "GermaniumActive",
    43: "Zinc",
    44: "ZincActive",
    45: "Carbon",
    46: "CarbonActive",
    47: "Generic",
    48: "GenericActive",
    50: "GRN",
    51: "Ga",
    60: "Psi",
    61: "Electron",
    62: "Hole",
    70: "XVEL",
    71: "YVEL",
}


@dataclass
class Region:
    index: int
    material_id: int
    material: str


@dataclass
class SupremStructure:
    path: Path
    points: np.ndarray
    triangles: np.ndarray
    triangle_regions: np.ndarray
    regions: dict[int, Region]
    solution_ids: list[int]
    solution_names: list[str]
    solutions: dict[tuple[int, int], np.ndarray]
    point_fields: dict[str, np.ndarray] | None = None
    cell_fields: dict[str, np.ndarray] | None = None

    @property
    def fields(self) -> list[str]:
        names = ["Material"]
        names.extend(self.solution_names)
        for field_names in (self.point_fields, self.cell_fields):
            if field_names:
                for name in field_names:
                    if name not in names:
                        names.append(name)
        derived = ["Donor", "Acceptor", "NetDoping"]
        for name in derived:
            if name not in names:
                names.append(name)
        return names


def read_suprem(path: str | Path) -> SupremStructure:
    path = Path(path)
    points: dict[int, tuple[float, float]] = {}
    triangles: list[tuple[int, int, int]] = []
    triangle_regions: list[int] = []
    regions: dict[int, Region] = {}
    solution_ids: list[int] = []
    solution_names: list[str] = []
    solutions: dict[tuple[int, int], np.ndarray] = {}

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if not raw.strip():
                continue
            parts = raw.split()
            tag = parts[0]
            if tag == "c":
                idx = int(parts[1]) - 1
                points[idx] = (float(parts[2]), float(parts[3]))
            elif tag == "t":
                region = int(parts[2]) - 1
                c1 = int(parts[3]) - 1
                c2 = int(parts[4]) - 1
                c3 = int(parts[5]) - 1
                triangles.append((c1, c2, c3))
                triangle_regions.append(region)
            elif tag == "r":
                idx = int(parts[1]) - 1
                material_id = int(parts[2])
                material = MATERIALS.get(material_id, f"material_{material_id}")
                regions[idx] = Region(idx, material_id, material)
            elif tag == "s":
                count = int(parts[1])
                solution_ids = [int(value) for value in parts[2 : 2 + count]]
                solution_names = [
                    DOPANTS.get(value, f"species_{value}") for value in solution_ids
                ]
            elif tag == "n":
                if not solution_ids:
                    continue
                node = int(parts[1])
                material_id = int(parts[2])
                values = np.array(
                    [float(value) for value in parts[3 : 3 + len(solution_ids)]],
                    dtype=float,
                )
                solutions[(node, material_id)] = values

    if not points:
        raise ValueError(f"{path} has no SUPREM coordinate records")
    if not triangles:
        raise ValueError(f"{path} has no SUPREM triangle records")

    point_array = np.full((max(points) + 1, 2), np.nan, dtype=float)
    for idx, point in points.items():
        point_array[idx] = point

    return SupremStructure(
        path=path,
        points=point_array,
        triangles=np.array(triangles, dtype=int),
        triangle_regions=np.array(triangle_regions, dtype=int),
        regions=regions,
        solution_ids=solution_ids,
        solution_names=solution_names,
        solutions=solutions,
    )


VTK_DTYPES = {
    "Float32": np.dtype("<f4"),
    "Float64": np.dtype("<f8"),
    "Int8": np.dtype("i1"),
    "UInt8": np.dtype("u1"),
    "Int16": np.dtype("<i2"),
    "UInt16": np.dtype("<u2"),
    "Int32": np.dtype("<i4"),
    "UInt32": np.dtype("<u4"),
    "Int64": np.dtype("<i8"),
    "UInt64": np.dtype("<u8"),
}


def _strip_namespace(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _decode_vtu_base64_segment(text: str, offset: int) -> tuple[bytes, int]:
    end = text.find("=", offset)
    if end < 0:
        chunk = text[offset:]
        padding = "=" * ((4 - len(chunk) % 4) % 4)
        return base64.b64decode(chunk + padding), len(text)
    while end + 1 < len(text) and text[end + 1] == "=":
        end += 1
    return base64.b64decode(text[offset : end + 1]), end + 1


def _decode_vtu_appended_array(text: str, offset: int, compressed: bool) -> bytes:
    if compressed:
        header, pos = _decode_vtu_base64_segment(text, offset)
        block_count, block_size, last_block_size = struct.unpack_from("<III", header, 0)
        sizes = struct.unpack_from("<" + "I" * block_count, header, 12)
        compressed_payload = bytearray()
        total_size = sum(sizes)
        while len(compressed_payload) < total_size:
            block, pos = _decode_vtu_base64_segment(text, pos)
            compressed_payload.extend(block)
        chunks = []
        payload_pos = 0
        for size in sizes:
            block = compressed_payload[payload_pos : payload_pos + size]
            payload_pos += size
            chunks.append(zlib.decompress(block))
        payload = b"".join(chunks)
        expected = block_size * max(block_count - 1, 0) + last_block_size
        return payload[:expected]
    header, pos = _decode_vtu_base64_segment(text, offset)
    byte_count = struct.unpack_from("<I", header, 0)[0]
    payload, _ = _decode_vtu_base64_segment(text, pos)
    return payload[:byte_count]


def _read_vtu_array(
    element: ET.Element, appended: str, compressed: bool
) -> np.ndarray:
    dtype_name = element.attrib["type"]
    dtype = VTK_DTYPES.get(dtype_name)
    if dtype is None:
        raise ValueError(f"unsupported VTU array type {dtype_name!r}")

    fmt = element.attrib.get("format", "ascii")
    components = int(element.attrib.get("NumberOfComponents", "1"))
    if fmt == "appended":
        payload = _decode_vtu_appended_array(appended, int(element.attrib["offset"]), compressed)
        values = np.frombuffer(payload, dtype=dtype).copy()
    elif fmt == "ascii":
        text = element.text or ""
        values = np.fromstring(text, sep=" ", dtype=dtype)
    else:
        raise ValueError(f"unsupported VTU array format {fmt!r}")

    if components > 1:
        values = values.reshape((-1, components))
    return values


def _vtu_region_map(root: ET.Element) -> dict[int, Region]:
    regions: dict[int, Region] = {}
    for element in root.iter():
        if _strip_namespace(element.tag) != "region":
            continue
        if "id" not in element.attrib:
            continue
        idx = int(element.attrib["id"])
        material = element.attrib.get("material") or element.attrib.get("name") or "region"
        material_id = next(
            (mid for mid, name in MATERIALS.items() if name == material),
            idx,
        )
        regions[idx] = Region(idx, material_id, material)
    return regions


def read_vtu(path: str | Path) -> SupremStructure:
    path = Path(path)
    raw = path.read_bytes()
    root = ET.fromstring(raw)
    compressed = "compressor" in root.attrib

    appended_text = ""
    for element in root.iter():
        if _strip_namespace(element.tag) == "AppendedData":
            appended_text = element.text or ""
            break
    appended_text = "".join(appended_text.split())
    appended = ""
    if appended_text:
        if not appended_text.startswith("_"):
            raise ValueError("VTU appended data is missing the '_' marker")
        appended = appended_text[1:]

    point_data: dict[str, np.ndarray] = {}
    cell_data: dict[str, np.ndarray] = {}
    points = None
    connectivity = None
    offsets = None
    cell_types = None

    current_section = ""
    for element in root.iter():
        tag = _strip_namespace(element.tag)
        if tag in {"PointData", "CellData", "Points", "Cells"}:
            current_section = tag
            continue
        if tag != "DataArray":
            continue

        name = element.attrib.get("Name", "")
        components = int(element.attrib.get("NumberOfComponents", "1"))
        if current_section in {"PointData", "CellData"} and components != 1:
            continue
        array = _read_vtu_array(element, appended, compressed)

        if current_section == "Points":
            points = array[:, :2]
        elif current_section == "Cells":
            if name == "connectivity":
                connectivity = array.astype(int)
            elif name == "offsets":
                offsets = array.astype(int)
            elif name == "types":
                cell_types = array.astype(int)
        elif current_section == "PointData" and name and components == 1:
            point_data[name] = array.astype(float)
        elif current_section == "CellData" and name and components == 1:
            cell_data[name] = array.astype(float)

    if points is None or connectivity is None or offsets is None:
        raise ValueError(f"{path} is missing VTU points or cells")

    triangles: list[tuple[int, int, int]] = []
    source_cell_indices: list[int] = []
    start = 0
    for cell_i, stop in enumerate(offsets):
        cell = connectivity[start:stop]
        start = int(stop)
        vtk_type = int(cell_types[cell_i]) if cell_types is not None else 5
        if vtk_type == 5 and len(cell) == 3:
            triangles.append((int(cell[0]), int(cell[1]), int(cell[2])))
            source_cell_indices.append(cell_i)
        elif vtk_type in {7, 9} and len(cell) >= 4:
            first = int(cell[0])
            for i in range(1, len(cell) - 1):
                triangles.append((first, int(cell[i]), int(cell[i + 1])))
                source_cell_indices.append(cell_i)

    if not triangles:
        raise ValueError(f"{path} has no plottable triangle cells")

    regions = _vtu_region_map(root)
    region_field = cell_data.get("region")
    if region_field is None:
        triangle_regions = np.zeros(len(triangles), dtype=int)
        if not regions:
            regions[0] = Region(0, 0, "region_0")
    else:
        triangle_regions = region_field[source_cell_indices].astype(int)
        if not regions:
            for region_id in sorted(set(int(value) for value in triangle_regions)):
                regions[region_id] = Region(region_id, region_id, f"region_{region_id}")

    expanded_cell_data = {
        name: values[source_cell_indices].astype(float)
        for name, values in cell_data.items()
        if len(values) >= max(source_cell_indices) + 1
    }

    return SupremStructure(
        path=path,
        points=points.astype(float),
        triangles=np.array(triangles, dtype=int),
        triangle_regions=triangle_regions,
        regions=regions,
        solution_ids=[],
        solution_names=[],
        solutions={},
        point_fields=point_data,
        cell_fields=expanded_cell_data,
    )


def read_structure(path: str | Path) -> SupremStructure:
    suffix = Path(path).suffix.lower()
    if suffix == ".vtu":
        return read_vtu(path)
    return read_suprem(path)


def _solution_index(structure: SupremStructure, names: list[str]) -> int | None:
    for name in names:
        if name in structure.solution_names:
            return structure.solution_names.index(name)
    return None


def material_entries(structure: SupremStructure) -> list[tuple[int, str]]:
    entries = {}
    for region in structure.regions.values():
        entries[region.material_id] = region.material
    return sorted(entries.items())


def material_values(structure: SupremStructure) -> np.ndarray:
    values = np.full(len(structure.triangles), np.nan, dtype=float)
    material_to_value = {
        material_id: i for i, (material_id, _) in enumerate(material_entries(structure))
    }
    for i, region_id in enumerate(structure.triangle_regions):
        region = structure.regions.get(int(region_id))
        if region is not None:
            values[i] = material_to_value[region.material_id]
    return values


def add_material_legend(ax, structure: SupremStructure, cmap) -> None:
    from matplotlib.patches import Patch

    handles = []
    for i, (material_id, material) in enumerate(material_entries(structure)):
        handles.append(
            Patch(
                facecolor=cmap(i),
                edgecolor="black",
                linewidth=0.5,
                label=f"{material}  id {material_id}",
            )
        )
    ax.legend(
        handles=handles,
        title="Materials",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=True,
    )


def apply_gui_layout(fig, field: str) -> None:
    if field == "Material":
        fig.subplots_adjust(left=0.08, right=0.78, bottom=0.1, top=0.9)
    else:
        fig.subplots_adjust(left=0.08, right=0.88, bottom=0.1, top=0.9)


def scalar_values(structure: SupremStructure, field: str) -> np.ndarray:
    if field == "Material":
        return material_values(structure)

    point_fields = structure.point_fields or {}
    cell_fields = structure.cell_fields or {}
    if field == "Donor":
        for name in ("Nd[cm-3]", "Donor", "donor"):
            if name in point_fields or name in cell_fields:
                return scalar_values(structure, name)
    elif field == "Acceptor":
        for name in ("Na[cm-3]", "Acceptor", "acceptor"):
            if name in point_fields or name in cell_fields:
                return scalar_values(structure, name)
    elif field == "NetDoping":
        for name in ("net_doping[cm-3]", "NetDoping", "net_doping"):
            if name in point_fields or name in cell_fields:
                return scalar_values(structure, name)

    if field in cell_fields:
        return cell_fields[field].astype(float)
    if field in point_fields:
        values = np.full(len(structure.triangles), np.nan, dtype=float)
        point_values = point_fields[field]
        for tri_i, tri in enumerate(structure.triangles):
            samples = point_values[tri]
            finite = np.isfinite(samples)
            if np.any(finite):
                values[tri_i] = float(np.mean(samples[finite]))
        return values

    if field == "Donor":
        indices = [
            _solution_index(structure, ["ArsenicActive", "Arsenic"]),
            _solution_index(structure, ["PhosphorusActive", "Phosphorus"]),
            _solution_index(structure, ["AntimonyActive", "Antimony"]),
        ]
    elif field == "Acceptor":
        indices = [
            _solution_index(structure, ["BoronActive", "Boron"]),
            _solution_index(structure, ["AluminumActive", "Aluminum"]),
        ]
    elif field == "NetDoping":
        donor = scalar_values(structure, "Donor")
        acceptor = scalar_values(structure, "Acceptor")
        return donor - acceptor
    else:
        if field not in structure.solution_names:
            raise ValueError(f"unknown field {field!r}")
        indices = [structure.solution_names.index(field)]

    values = np.full(len(structure.triangles), np.nan, dtype=float)
    for tri_i, region_id in enumerate(structure.triangle_regions):
        region = structure.regions.get(int(region_id))
        if region is None:
            continue
        material_id = region.material_id
        tri = structure.triangles[tri_i]
        samples = []
        for node in tri:
            data = structure.solutions.get((int(node), material_id))
            if data is None:
                continue
            total = 0.0
            for idx in indices:
                if idx is not None:
                    total += data[idx]
            samples.append(total)
        if samples:
            values[tri_i] = float(np.mean(samples))
    return values


def transform_values(
    values: np.ndarray, field: str, log10: bool
) -> tuple[np.ndarray, str]:
    title_field = field
    if field != "Material" and log10:
        if field == "NetDoping":
            values = np.sign(values) * np.log10(np.maximum(np.abs(values), 1.0))
            title_field = f"sign*log10(abs({field}))"
        else:
            values = np.log10(np.maximum(values, 1.0))
            title_field = f"log10({field})"
    return values, title_field


def sample_line_profile(
    structure: SupremStructure,
    field: str,
    log10: bool,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    samples: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, str]:
    if field == "Material":
        raise ValueError("line profiles are for scalar fields, not Material")
    if samples < 2:
        raise ValueError("line profile needs at least 2 samples")

    point_fields = structure.point_fields or {}
    if field in point_fields:
        point_values, title_field = transform_values(point_fields[field].astype(float), field, log10)
    else:
        face_values = scalar_values(structure, field)
        face_values, title_field = transform_values(face_values, field, log10)
        point_values = face_to_point_values(structure, face_values)

    xs = np.linspace(x0, x1, samples)
    ys = np.linspace(y0, y1, samples)
    values = interpolate_line_values(structure, point_values, xs, ys)
    distances = np.hypot(xs - x0, ys - y0)
    return distances, values, xs, ys, title_field


def interpolate_line_values(
    structure: SupremStructure, point_values: np.ndarray, xs: np.ndarray, ys: np.ndarray
) -> np.ndarray:
    coords = structure.points[structure.triangles]
    mins = np.nanmin(coords, axis=1)
    maxs = np.nanmax(coords, axis=1)
    a = coords[:, 0, :]
    b = coords[:, 1, :]
    c = coords[:, 2, :]
    v0 = b - a
    v1 = c - a
    denom = v0[:, 0] * v1[:, 1] - v1[:, 0] * v0[:, 1]
    values = np.full(len(xs), np.nan, dtype=float)

    for i, (x, y) in enumerate(zip(xs, ys)):
        candidates = (
            (x >= mins[:, 0] - 1e-12)
            & (x <= maxs[:, 0] + 1e-12)
            & (y >= mins[:, 1] - 1e-12)
            & (y <= maxs[:, 1] + 1e-12)
            & (np.abs(denom) > 1e-30)
        )
        candidate_indices = np.flatnonzero(candidates)
        if not len(candidate_indices):
            continue
        v2 = np.array([x, y], dtype=float) - a[candidate_indices]
        candidate_v0 = v0[candidate_indices]
        candidate_v1 = v1[candidate_indices]
        candidate_denom = denom[candidate_indices]
        u = (v2[:, 0] * candidate_v1[:, 1] - candidate_v1[:, 0] * v2[:, 1]) / candidate_denom
        v = (candidate_v0[:, 0] * v2[:, 1] - v2[:, 0] * candidate_v0[:, 1]) / candidate_denom
        w = 1.0 - u - v
        inside = np.flatnonzero((u >= -1e-10) & (v >= -1e-10) & (w >= -1e-10))
        if not len(inside):
            continue
        tri_index = candidate_indices[int(inside[0])]
        tri_values = point_values[structure.triangles[tri_index]]
        if np.all(np.isfinite(tri_values)):
            values[i] = (
                w[inside[0]] * tri_values[0]
                + u[inside[0]] * tri_values[1]
                + v[inside[0]] * tri_values[2]
            )
    return values


def plot_line_profile(
    structure: SupremStructure,
    field: str,
    log10: bool,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    samples: int,
    output: str | None = None,
    show: bool = True,
):
    import matplotlib.pyplot as plt

    distances, values, _, _, title_field = sample_line_profile(
        structure, field, log10, x0, y0, x1, y1, samples
    )
    finite = np.isfinite(values)
    if not np.any(finite):
        raise ValueError("line does not intersect plottable scalar data")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(distances[finite], values[finite], s=18)
    ax.plot(distances[finite], values[finite], linewidth=0.8, alpha=0.6)
    ax.set_title(f"{structure.path.name}: {title_field} line profile")
    ax.set_xlabel("distance along line [um]")
    ax.set_ylabel(title_field)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if output:
        fig.savefig(output, dpi=180)
    if show:
        plt.show()
    else:
        plt.close(fig)


def save_line_profile_csv(
    path: str | Path,
    distances: np.ndarray,
    xs: np.ndarray,
    ys: np.ndarray,
    values: np.ndarray,
    value_name: str,
) -> None:
    with Path(path).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["distance_um", "x_um", "y_um", value_name])
        for distance, x, y, value in zip(distances, xs, ys, values):
            writer.writerow([distance, x, y, value])


def plot_structure(
    structure: SupremStructure,
    field: str = "Material",
    log10: bool = False,
    mesh: bool = False,
    contours: int = 0,
    vmin: float | None = None,
    vmax: float | None = None,
    output: str | None = None,
    show: bool = True,
):
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    triangulation = mtri.Triangulation(
        structure.points[:, 0], structure.points[:, 1], structure.triangles
    )
    values = scalar_values(structure, field)
    values, title_field = transform_values(values, field, log10)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title(f"{structure.path.name}: {title_field}")
    ax.set_xlabel("x [um]")
    ax.set_ylabel("y [um]")
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()

    if field == "Material":
        cmap = plt.get_cmap("tab20", max(len(material_entries(structure)), 1))
        artist = ax.tripcolor(
            triangulation,
            facecolors=values,
            shading="flat",
            edgecolors="k" if mesh else "none",
            linewidth=0.15,
            cmap=cmap,
            vmin=-0.5,
            vmax=max(len(material_entries(structure)) - 0.5, 0.5),
        )
        add_material_legend(ax, structure, cmap)
    else:
        finite = np.isfinite(values)
        if not np.any(finite):
            raise ValueError(f"field {field!r} has no plottable values")
        artist = ax.tripcolor(
            triangulation,
            facecolors=values,
            shading="flat",
            edgecolors="k" if mesh else "none",
            linewidth=0.12,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
        )
        fig.colorbar(artist, ax=ax, label=title_field)
        if contours > 0:
            point_values = face_to_point_values(structure, values)
            finite_points = np.isfinite(point_values)
            if np.count_nonzero(finite_points) >= 3:
                ax.tricontour(
                    triangulation,
                    point_values,
                    levels=contours,
                    colors="white",
                    linewidths=0.7,
                )

    if mesh and field != "Material":
        ax.triplot(triangulation, color="black", linewidth=0.1, alpha=0.4)

    fig.tight_layout()
    if output:
        fig.savefig(output, dpi=180)
    if show:
        plt.show()
    else:
        plt.close(fig)


def face_to_point_values(structure: SupremStructure, face_values: np.ndarray) -> np.ndarray:
    sums = np.zeros(len(structure.points), dtype=float)
    counts = np.zeros(len(structure.points), dtype=float)
    for tri, value in zip(structure.triangles, face_values):
        if not np.isfinite(value):
            continue
        for node in tri:
            sums[node] += value
            counts[node] += 1.0
    point_values = np.full(len(structure.points), np.nan, dtype=float)
    valid = counts > 0
    point_values[valid] = sums[valid] / counts[valid]
    return point_values


def launch_gui(
    path: str | Path,
    initial_vmin: float | None = None,
    initial_vmax: float | None = None,
):
    import tkinter as tk
    from tkinter import filedialog, ttk

    import matplotlib

    matplotlib.use("TkAgg")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    import matplotlib.tri as mtri

    structure = read_structure(path)

    root = tk.Tk()
    root.title(f"SUPREM Structure Viewer - {structure.path.name}")

    controls = ttk.Frame(root, padding=6)
    controls.pack(side=tk.TOP, fill=tk.X)
    line_controls = ttk.Frame(root, padding=(6, 0, 6, 6))
    line_controls.pack(side=tk.TOP, fill=tk.X)

    field_var = tk.StringVar(value="Material")
    log_var = tk.BooleanVar(value=False)
    mesh_var = tk.BooleanVar(value=False)
    auto_scale_var = tk.BooleanVar(
        value=initial_vmin is None and initial_vmax is None
    )
    min_var = tk.StringVar(
        value="" if initial_vmin is None else f"{initial_vmin:.6g}"
    )
    max_var = tk.StringVar(
        value="" if initial_vmax is None else f"{initial_vmax:.6g}"
    )
    xmin = float(np.nanmin(structure.points[:, 0]))
    xmax = float(np.nanmax(structure.points[:, 0]))
    ymin = float(np.nanmin(structure.points[:, 1]))
    ymax = float(np.nanmax(structure.points[:, 1]))
    xmid = 0.5 * (xmin + xmax)
    line_x0_var = tk.StringVar(value=f"{xmid:.6g}")
    line_y0_var = tk.StringVar(value=f"{ymin:.6g}")
    line_x1_var = tk.StringVar(value=f"{xmid:.6g}")
    line_y1_var = tk.StringVar(value=f"{ymax:.6g}")
    line_samples_var = tk.StringVar(value="200")

    ttk.Label(controls, text="Field").pack(side=tk.LEFT)
    field_box = ttk.Combobox(
        controls,
        textvariable=field_var,
        values=structure.fields,
        state="readonly",
        width=24,
    )
    field_box.pack(side=tk.LEFT, padx=(4, 10))
    ttk.Checkbutton(controls, text="log10", variable=log_var).pack(side=tk.LEFT)
    ttk.Checkbutton(controls, text="mesh", variable=mesh_var).pack(side=tk.LEFT)
    ttk.Checkbutton(controls, text="auto scale", variable=auto_scale_var).pack(
        side=tk.LEFT, padx=(10, 4)
    )
    ttk.Label(controls, text="min").pack(side=tk.LEFT)
    min_entry = ttk.Entry(controls, textvariable=min_var, width=10)
    min_entry.pack(side=tk.LEFT, padx=(4, 8))
    ttk.Label(controls, text="max").pack(side=tk.LEFT)
    max_entry = ttk.Entry(controls, textvariable=max_var, width=10)
    max_entry.pack(side=tk.LEFT, padx=(4, 8))

    ttk.Label(line_controls, text="Line x0").pack(side=tk.LEFT)
    ttk.Entry(line_controls, textvariable=line_x0_var, width=9).pack(
        side=tk.LEFT, padx=(4, 8)
    )
    ttk.Label(line_controls, text="y0").pack(side=tk.LEFT)
    ttk.Entry(line_controls, textvariable=line_y0_var, width=9).pack(
        side=tk.LEFT, padx=(4, 8)
    )
    ttk.Label(line_controls, text="x1").pack(side=tk.LEFT)
    ttk.Entry(line_controls, textvariable=line_x1_var, width=9).pack(
        side=tk.LEFT, padx=(4, 8)
    )
    ttk.Label(line_controls, text="y1").pack(side=tk.LEFT)
    ttk.Entry(line_controls, textvariable=line_y1_var, width=9).pack(
        side=tk.LEFT, padx=(4, 8)
    )
    ttk.Label(line_controls, text="samples").pack(side=tk.LEFT)
    ttk.Entry(line_controls, textvariable=line_samples_var, width=7).pack(
        side=tk.LEFT, padx=(4, 8)
    )

    fig = Figure(figsize=(10, 6), dpi=100)
    ax = fig.add_subplot(111)
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    toolbar = NavigationToolbar2Tk(canvas, root)
    toolbar.update()

    status = tk.StringVar(value="")
    ttk.Label(root, textvariable=status, anchor="w").pack(side=tk.BOTTOM, fill=tk.X)

    triangulation = mtri.Triangulation(
        structure.points[:, 0], structure.points[:, 1], structure.triangles
    )

    def redraw(*_):
        nonlocal ax
        fig.clear()
        ax = fig.add_subplot(111)
        field = field_var.get()
        values = scalar_values(structure, field)
        title_field = field
        if field != "Material" and log_var.get():
            if field == "NetDoping":
                values = np.sign(values) * np.log10(np.maximum(np.abs(values), 1.0))
                title_field = f"sign*log10(abs({field}))"
            else:
                values = np.log10(np.maximum(values, 1.0))
                title_field = f"log10({field})"

        vmin = None
        vmax = None
        if field != "Material" and not auto_scale_var.get():
            try:
                vmin = float(min_var.get()) if min_var.get().strip() else None
                vmax = float(max_var.get()) if max_var.get().strip() else None
            except ValueError:
                status.set("Scale min/max must be numbers, or enable auto scale.")
                canvas.draw_idle()
                return
            if vmin is not None and vmax is not None and vmin >= vmax:
                status.set("Scale min must be less than scale max.")
                canvas.draw_idle()
                return

        ax.set_title(f"{structure.path.name}: {title_field}")
        ax.set_xlabel("x [um]")
        ax.set_ylabel("y [um]")
        ax.set_aspect("equal", adjustable="datalim", anchor="C")
        ax.invert_yaxis()

        if field == "Material":
            cmap = matplotlib.colormaps["tab20"].resampled(
                max(len(material_entries(structure)), 1)
            )
            artist = ax.tripcolor(
                triangulation,
                facecolors=values,
                shading="flat",
                edgecolors="k" if mesh_var.get() else "none",
                linewidth=0.15,
                cmap=cmap,
                vmin=-0.5,
                vmax=max(len(material_entries(structure)) - 0.5, 0.5),
            )
            add_material_legend(ax, structure, cmap)
        else:
            artist = ax.tripcolor(
                triangulation,
                facecolors=values,
                shading="flat",
                edgecolors="k" if mesh_var.get() else "none",
                linewidth=0.12,
                cmap="viridis",
                vmin=vmin,
                vmax=vmax,
            )
        if mesh_var.get() and field != "Material":
            ax.triplot(triangulation, color="black", linewidth=0.1, alpha=0.4)
        try:
            x0, y0, x1, y1, _ = read_line_inputs()
            ax.plot([x0, x1], [y0, y1], color="red", linewidth=1.4, linestyle="--")
        except ValueError:
            pass
        if field != "Material":
            cbar = fig.colorbar(artist, ax=ax)
            cbar.set_label(title_field)
        apply_gui_layout(fig, field)
        finite = values[np.isfinite(values)]
        if finite.size:
            if field != "Material" and auto_scale_var.get():
                min_var.set(f"{finite.min():.6g}")
                max_var.set(f"{finite.max():.6g}")
            status.set(
                f"{field}: min={finite.min():.6g}, max={finite.max():.6g}, "
                f"triangles={len(structure.triangles)}, points={len(structure.points)}"
            )
        canvas.draw_idle()

    def save_png():
        filename = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG image", "*.png"), ("All files", "*.*")],
        )
        if filename:
            fig.savefig(filename, dpi=180)

    def read_line_inputs() -> tuple[float, float, float, float, int]:
        x0 = float(line_x0_var.get())
        y0 = float(line_y0_var.get())
        x1 = float(line_x1_var.get())
        y1 = float(line_y1_var.get())
        samples = int(line_samples_var.get())
        if samples < 2:
            raise ValueError("line profile needs at least 2 samples")
        if x0 == x1 and y0 == y1:
            raise ValueError("line endpoints must be different")
        return x0, y0, x1, y1, samples

    def open_line_profile():
        field = field_var.get()
        if field == "Material":
            status.set("Choose a scalar field before plotting a line profile.")
            return
        try:
            x0, y0, x1, y1, samples = read_line_inputs()
            distances, values, xs, ys, title_field = sample_line_profile(
                structure, field, log_var.get(), x0, y0, x1, y1, samples
            )
        except ValueError as exc:
            status.set(str(exc))
            return
        finite = np.isfinite(values)
        if not np.any(finite):
            status.set("Line does not intersect plottable scalar data.")
            return

        profile = tk.Toplevel(root)
        profile.title(f"Line Profile - {title_field}")
        profile_controls = ttk.Frame(profile, padding=6)
        profile_controls.pack(side=tk.TOP, fill=tk.X)
        profile_fig = Figure(figsize=(8, 5), dpi=100)
        profile_ax = profile_fig.add_subplot(111)
        profile_ax.scatter(distances[finite], values[finite], s=18)
        profile_ax.plot(distances[finite], values[finite], linewidth=0.8, alpha=0.6)
        profile_ax.set_title(f"{structure.path.name}: {title_field} line profile")
        profile_ax.set_xlabel("distance along line [um]")
        profile_ax.set_ylabel(title_field)
        profile_ax.grid(True, alpha=0.3)
        profile_fig.tight_layout()
        profile_canvas = FigureCanvasTkAgg(profile_fig, master=profile)
        profile_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        profile_toolbar = NavigationToolbar2Tk(profile_canvas, profile)
        profile_toolbar.update()

        def save_profile_png():
            filename = filedialog.asksaveasfilename(
                parent=profile,
                defaultextension=".png",
                filetypes=[("PNG image", "*.png"), ("All files", "*.*")],
            )
            if filename:
                profile_fig.savefig(filename, dpi=180)

        def save_profile_csv():
            filename = filedialog.asksaveasfilename(
                parent=profile,
                defaultextension=".csv",
                filetypes=[("CSV file", "*.csv"), ("All files", "*.*")],
            )
            if filename:
                save_line_profile_csv(filename, distances, xs, ys, values, title_field)

        ttk.Button(profile_controls, text="Save PNG", command=save_profile_png).pack(
            side=tk.LEFT
        )
        ttk.Button(profile_controls, text="Save CSV", command=save_profile_csv).pack(
            side=tk.LEFT, padx=(8, 0)
        )
        status.set(
            f"Line profile: {np.count_nonzero(finite)} finite samples of {samples}."
        )
        profile_canvas.draw_idle()

    ttk.Button(controls, text="Redraw", command=redraw).pack(side=tk.LEFT, padx=8)
    ttk.Button(controls, text="Apply Scale", command=redraw).pack(side=tk.LEFT)
    ttk.Button(controls, text="Save PNG", command=save_png).pack(side=tk.LEFT)
    ttk.Button(line_controls, text="Plot Line", command=open_line_profile).pack(
        side=tk.LEFT, padx=(8, 0)
    )

    field_box.bind("<<ComboboxSelected>>", redraw)
    log_var.trace_add("write", redraw)
    mesh_var.trace_add("write", redraw)
    auto_scale_var.trace_add("write", redraw)
    redraw()
    root.mainloop()


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure", help="SUPREM .str or VTK .vtu structure file")
    parser.add_argument("--field", default="Material", help="field to plot")
    parser.add_argument("--log10", action="store_true", help="plot log10(field)")
    parser.add_argument("--mesh", action="store_true", help="overlay mesh edges")
    parser.add_argument("--contours", type=int, default=0, help="number of contour levels")
    parser.add_argument("--vmin", type=float, help="minimum scalar color scale value")
    parser.add_argument("--vmax", type=float, help="maximum scalar color scale value")
    parser.add_argument("--save", help="save PNG instead of opening a window")
    parser.add_argument(
        "--line",
        nargs=4,
        type=float,
        metavar=("X0", "Y0", "X1", "Y1"),
        help="plot a 1-D scalar profile along this line",
    )
    parser.add_argument(
        "--line-samples",
        type=int,
        default=200,
        help="number of samples for --line",
    )
    parser.add_argument("--line-save", help="save the 1-D line profile PNG")
    parser.add_argument("--line-csv", help="save the 1-D line profile CSV")
    parser.add_argument("--list", action="store_true", help="list available fields and exit")
    parser.add_argument("--gui", action="store_true", help="open the Tk viewer (default)")
    parser.add_argument(
        "--plot",
        action="store_true",
        help="open a one-shot matplotlib plot instead of the Tk viewer",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if args.vmin is not None and args.vmax is not None and args.vmin >= args.vmax:
        raise ValueError("--vmin must be less than --vmax")

    if (
        not args.save
        and not args.list
        and not args.plot
        and not args.line
        and not args.line_save
        and not args.line_csv
    ):
        launch_gui(args.structure, initial_vmin=args.vmin, initial_vmax=args.vmax)
        return 0

    structure = read_structure(args.structure)
    if args.list:
        print("Fields:")
        for field in structure.fields:
            print(f"  {field}")
        print("\nRegions:")
        for region in sorted(structure.regions.values(), key=lambda item: item.index):
            print(f"  {region.index}: {region.material} (material id {region.material_id})")
        print(f"\nPoints: {len(structure.points)}")
        print(f"Triangles: {len(structure.triangles)}")
        print(f"Solution records: {len(structure.solutions)}")
        return 0

    if args.line or args.line_save or args.line_csv:
        if args.line is None:
            raise ValueError("--line is required with --line-save or --line-csv")
        x0, y0, x1, y1 = args.line
        if args.line_csv:
            distances, values, xs, ys, title_field = sample_line_profile(
                structure,
                args.field,
                args.log10,
                x0,
                y0,
                x1,
                y1,
                args.line_samples,
            )
            save_line_profile_csv(args.line_csv, distances, xs, ys, values, title_field)
        if args.line_save or not args.line_csv:
            plot_line_profile(
                structure,
                args.field,
                args.log10,
                x0,
                y0,
                x1,
                y1,
                args.line_samples,
                output=args.line_save,
                show=args.line_save is None,
            )
        return 0

    plot_structure(
        structure,
        field=args.field,
        log10=args.log10,
        mesh=args.mesh,
        contours=args.contours,
        vmin=args.vmin,
        vmax=args.vmax,
        output=args.save,
        show=args.save is None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
