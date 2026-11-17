#!/usr/bin/env python3
"""Convert a SUPREM structure to a Genius TIF with named contacts.

The Genius SUPREM importer reads the mesh and doping, but it does not create
boundary records.  This script keeps the SUPREM mesh and solution data and adds
boundary/interface edges.  The default CLI behavior keeps the original PMOS
convenience naming: three aluminum regions become Source/Gate/Drain and the
silicon backside becomes Body.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


SUPREM_MATERIAL_TO_TIF = {
    1: "Ox",
    2: "Nit",
    3: "Si",
    4: "Elec",
    6: "Elec",
}

DOPING_CODE_TO_NAME = {
    2: "Arsenic",
    3: "Phosphorus",
    4: "Antimony",
    5: "Boron",
    7: "Aluminum",
    8: "Nitrogen",
    20: "ArsenicActive",
    21: "PhosphorusActive",
    22: "AntimonyActive",
    23: "BoronActive",
    24: "AluminumActive",
    25: "NitrogenActive",
}


@dataclass
class Node:
    idx: int
    x: float
    y: float
    h: float


@dataclass
class Tri:
    idx: int
    region: int
    nodes: tuple[int, int, int]
    neigh: tuple[int, int, int]


@dataclass
class Region:
    idx: int
    mat_id: int
    name: str

    @property
    def tif_material(self) -> str:
        return SUPREM_MATERIAL_TO_TIF.get(self.mat_id, "Vacuum")


@dataclass
class SolutionRow:
    node_index_zero_based: int
    material_id: int
    values: list[float]


@dataclass
class ContactSpec:
    name: str
    kind: str
    target: str
    region: int | None = None
    edges: list[int] | None = None
    work_function: float = 5.25

    @property
    def genius_type(self) -> str:
        return "GateContact" if self.kind == "gate" else "OhmicContact"


def parse_suprem(path: Path, auto_name_electrodes: bool = True):
    nodes: dict[int, Node] = {}
    tris: list[Tri] = []
    regions: dict[int, Region] = {}
    sol_codes: list[int] = []
    sol_rows: list[SolutionRow] = []

    for raw in path.read_text().splitlines():
        parts = raw.split()
        if not parts:
            continue
        flag = parts[0]
        if flag == "c":
            idx = int(parts[1])
            nodes[idx] = Node(idx, float(parts[2]), float(parts[3]), float(parts[4]))
        elif flag == "t":
            idx = int(parts[1])
            region = int(parts[2])
            tri_nodes = (int(parts[3]), int(parts[4]), int(parts[5]))
            neigh = (int(parts[6]), int(parts[7]), int(parts[8]))
            tris.append(Tri(idx, region, tri_nodes, neigh))
        elif flag == "r":
            idx = int(parts[1])
            mat_id = int(parts[2])
            regions[idx] = Region(idx, mat_id, f"region_{idx - 1}")
        elif flag == "s":
            sol_codes = [int(v) for v in parts[2 : 2 + int(parts[1])]]
        elif flag == "n":
            node_index_zero_based = int(parts[1])
            material_id = int(parts[2])
            values = [float(v) for v in parts[3:]]
            sol_rows.append(SolutionRow(node_index_zero_based, material_id, values))

    if not nodes or not tris or not regions:
        raise RuntimeError(f"{path} does not look like a complete SUPREM structure")

    if auto_name_electrodes:
        name_electrodes(nodes, tris, regions)
    return nodes, tris, regions, sol_codes, sol_rows


def region_extents(
    nodes: dict[int, Node], tris: list[Tri]
) -> dict[int, tuple[float, float, float, float]]:
    extents: dict[int, list[float]] = {}
    for tri in tris:
        box = extents.setdefault(
            tri.region, [float("inf"), -float("inf"), float("inf"), -float("inf")]
        )
        for node_id in tri.nodes:
            node = nodes[node_id]
            box[0] = min(box[0], node.x)
            box[1] = max(box[1], node.x)
            box[2] = min(box[2], node.y)
            box[3] = max(box[3], node.y)
    return {idx: tuple(box) for idx, box in extents.items()}


def name_electrodes(nodes: dict[int, Node], tris: list[Tri], regions: dict[int, Region]) -> None:
    extents = region_extents(nodes, tris)
    aluminum = [r for r in regions.values() if r.mat_id == 6]
    aluminum.sort(key=lambda r: 0.5 * (extents[r.idx][0] + extents[r.idx][1]))
    labels = ["Source", "Gate", "Drain"]
    for region, label in zip(aluminum, labels):
        region.name = label


def suggested_contacts(
    nodes: dict[int, Node],
    tris: list[Tri],
    regions: dict[int, Region],
    body_tol: float = 1e-6,
) -> list[ContactSpec]:
    extents = region_extents(nodes, tris)
    aluminum = [r for r in regions.values() if r.mat_id == 6]
    aluminum.sort(key=lambda r: 0.5 * (extents[r.idx][0] + extents[r.idx][1]))
    labels = ["Source", "Gate", "Drain"]
    contacts = [
        ContactSpec(
            name=label,
            kind="gate" if label == "Gate" else "ohmic",
            target="region",
            region=region.idx,
        )
        for region, label in zip(aluminum, labels)
    ]

    _boundary_edges, edge_number, edge_regions = build_edges(tris)
    body_edges = find_body_edges(nodes, regions, edge_regions, edge_number, body_tol)
    if body_edges:
        contacts.append(
            ContactSpec(name="Body", kind="ohmic", target="boundary", edges=body_edges)
        )
    return contacts


def tri_edges(tri: Tri) -> tuple[tuple[int, int], tuple[int, int], tuple[int, int]]:
    a, b, c = tri.nodes
    return ((a, b), (b, c), (c, a))


def normalized_edge(edge: tuple[int, int]) -> tuple[int, int]:
    a, b = edge
    return (a, b) if a < b else (b, a)


def build_edges(tris: list[Tri]):
    edge_to_tris: dict[tuple[int, int], list[tuple[Tri, int]]] = defaultdict(list)
    for tri in tris:
        for side, edge in enumerate(tri_edges(tri)):
            edge_to_tris[normalized_edge(edge)].append((tri, side))

    boundary_edges: list[tuple[int, int]] = []
    edge_regions: dict[tuple[int, int], set[int]] = {}
    for edge, owners in edge_to_tris.items():
        regions = {tri.region for tri, _side in owners}
        if len(owners) == 1 or len(regions) > 1:
            edge_regions[edge] = regions
            boundary_edges.append(edge)

    boundary_edges.sort()
    edge_number = {edge: i + 1 for i, edge in enumerate(boundary_edges)}
    return boundary_edges, edge_number, edge_regions


def find_body_edges(
    nodes: dict[int, Node],
    regions: dict[int, Region],
    edge_regions: dict[tuple[int, int], set[int]],
    edge_number: dict[tuple[int, int], int],
    tol: float,
) -> list[int]:
    silicon_regions = {idx for idx, region in regions.items() if region.tif_material == "Si"}
    ymax = max(node.y for node in nodes.values())
    body = []
    for edge, owners in edge_regions.items():
        if not (owners & silicon_regions):
            continue
        n1, n2 = nodes[edge[0]], nodes[edge[1]]
        if abs(n1.y - ymax) <= tol and abs(n2.y - ymax) <= tol:
            body.append(edge_number[edge])
    body.sort()
    return body


def region_boundary_map(
    edge_regions: dict[tuple[int, int], set[int]],
    edge_number: dict[tuple[int, int], int],
) -> dict[int, list[int]]:
    by_region: dict[int, list[int]] = defaultdict(list)
    for edge, owners in edge_regions.items():
        for region in owners:
            by_region[region].append(edge_number[edge])
    for edges in by_region.values():
        edges.sort()
    return by_region


def solution_name_to_index(sol_codes: list[int]) -> dict[str, int]:
    names = [DOPING_CODE_TO_NAME.get(code, f"Code{code}") for code in sol_codes]
    return {name: idx for idx, name in enumerate(names)}


def pick(row: SolutionRow, name_to_index: dict[str, int], *names: str) -> float:
    for name in names:
        idx = name_to_index.get(name)
        if idx is not None and idx < len(row.values):
            return row.values[idx]
    return 0.0


def acceptor_donor(row: SolutionRow, name_to_index: dict[str, int]) -> tuple[float, float]:
    accept = (
        pick(row, name_to_index, "BoronActive", "Boron")
        + pick(row, name_to_index, "AluminumActive", "Aluminum")
    )
    donor = (
        pick(row, name_to_index, "PhosphorusActive", "Phosphorus")
        + pick(row, name_to_index, "ArsenicActive", "Arsenic")
        + pick(row, name_to_index, "AntimonyActive", "Antimony")
        + pick(row, name_to_index, "NitrogenActive", "Nitrogen")
    )
    return accept, donor


def write_tif(
    out_path: Path,
    nodes: dict[int, Node],
    tris: list[Tri],
    regions: dict[int, Region],
    sol_codes: list[int],
    sol_rows: list[SolutionRow],
    body_tol: float,
    contacts: list[ContactSpec] | None = None,
) -> dict[str, int]:
    boundary_edges, edge_number, edge_regions = build_edges(tris)
    by_region = region_boundary_map(edge_regions, edge_number)
    if contacts is None:
        body_edges = find_body_edges(nodes, regions, edge_regions, edge_number, body_tol)
        if not body_edges:
            raise RuntimeError("No backside Body edges found. Try increasing --body-tol.")
        contacts = [
            ContactSpec(name=region.name, kind="ohmic", target="region", region=region.idx)
            for region in regions.values()
            if region.name in {"Source", "Drain"}
        ]
        contacts.extend(
            ContactSpec(name=region.name, kind="gate", target="region", region=region.idx)
            for region in regions.values()
            if region.name == "Gate"
        )
        contacts.append(
            ContactSpec(name="Body", kind="ohmic", target="boundary", edges=body_edges)
        )
    body_edge_count = sum(
        len(contact.edges or [])
        for contact in contacts
        if contact.target == "boundary" and contact.name == "Body"
    )
    region_contact_names = {
        contact.region: contact.name
        for contact in contacts
        if contact.target == "region" and contact.region is not None
    }

    name_to_index = solution_name_to_index(sol_codes)
    material_name = {
        0: "Vacuum",
        1: "Ox",
        2: "Nit",
        3: "Si",
        4: "Elec",
        6: "Elec",
    }

    with out_path.open("w") as fout:
        fout.write("h TIF V1.2.1 MEDICI generated_by_make_genius_contacts.py\n")
        fout.write("cd GEN blnk blnk blnk cart2D 1.0 0.0\n")
        fout.write("cg 300\n")

        for idx in sorted(nodes):
            node = nodes[idx]
            fout.write(f"c {idx:8d} {node.x: .12e} {node.y: .12e} {node.h: .12e}\n")

        for edge in boundary_edges:
            eidx = edge_number[edge]
            fout.write(f"e {eidx:8d} {edge[0]:8d} {edge[1]:8d} 0\n")

        for idx in sorted(regions):
            region = regions[idx]
            name = region_contact_names.get(idx, region.name)
            fout.write(f"r {idx:4d} {region.tif_material:<14s} {name}\n")
            for eidx in by_region.get(idx, []):
                fout.write(f"b {eidx:8d}\n")

        interface_index = max(regions) + 1
        for contact in contacts:
            if contact.target != "boundary":
                continue
            edges = sorted(set(contact.edges or []))
            if not edges:
                continue
            fout.write(f"i {interface_index:4d} Vacuum         {contact.name} 0\n")
            for eidx in edges:
                fout.write(f"j {eidx:8d}\n")
            interface_index += 1

        for tri in tris:
            fout.write(
                f"t {tri.idx:8d} {tri.region:8d}"
                f" {tri.nodes[0]:8d} {tri.nodes[1]:8d} {tri.nodes[2]:8d}"
                f" {tri.neigh[0]:8d} {tri.neigh[1]:8d} {tri.neigh[2]:8d}\n"
            )

        if sol_rows:
            fout.write("s 2 Accept Donor\n")
            for row in sol_rows:
                accept, donor = acceptor_donor(row, name_to_index)
                medici_node_index = row.node_index_zero_based + 1
                mat = material_name.get(row.material_id, "Vacuum")
                fout.write(
                    f"n {medici_node_index:8d} {mat:<14s}"
                    f" {accept: .12e} {donor: .12e}\n"
                )

    counts = {
        "nodes": len(nodes),
        "triangles": len(tris),
        "edges": len(boundary_edges),
        "body_edges": body_edge_count,
    }
    return counts


def write_deck(out_path: Path, tif_path: Path, contacts: list[ContactSpec] | None = None) -> None:
    if contacts is None:
        contacts = [
            ContactSpec("Source", "ohmic", "region"),
            ContactSpec("Gate", "gate", "region"),
            ContactSpec("Drain", "ohmic", "region"),
            ContactSpec("Body", "ohmic", "boundary"),
        ]
    fout = out_path.open("w")
    with fout:
        fout.write("GLOBAL TExternal=300 DopingScale=1e18 Z.Width=1.0\n\n")
        fout.write(f'IMPORT TIFFile="{tif_path.name}"\n\n')
        for contact in contacts:
            if contact.target == "region":
                if contact.kind == "gate":
                    fout.write(
                        f"CONTACT  Type=GateContact  ID={contact.name} "
                        f"WorkFunction={contact.work_function:g} Res=0 Cap=0 Ind=0\n"
                    )
                else:
                    fout.write(
                        f"CONTACT  Type=OhmicContact ID={contact.name} Res=0 Cap=0 Ind=0\n"
                    )
        for contact in contacts:
            if contact.target == "boundary":
                fout.write(
                    f"BOUNDARY Type={contact.genius_type} ID={contact.name} "
                    "Res=0 Cap=0 Ind=0\n"
                )
        fout.write("\n")
        fout.write("METHOD Type=DDML1 NS=Basic LS=MUMPS MaxIt=30\n")
        fout.write("SOLVE Type=EQUILIBRIUM\n\n")
        fout.write(f"EXPORT VTKFile={tif_path.stem}_equilibrium.vtu\n\n")
        fout.write("END\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="pmos.str", type=Path)
    parser.add_argument("--tif", default="pmos_contacts.tif", type=Path)
    parser.add_argument("--deck", default="rit_pmos_contacts.in", type=Path)
    parser.add_argument("--body-tol", default=1e-6, type=float)
    args = parser.parse_args()

    nodes, tris, regions, sol_codes, sol_rows = parse_suprem(args.input)
    counts = write_tif(args.tif, nodes, tris, regions, sol_codes, sol_rows, args.body_tol)
    write_deck(args.deck, args.tif)

    print(f"wrote {args.tif}")
    print(f"wrote {args.deck}")
    print(
        "mesh: "
        f"{counts['nodes']} nodes, {counts['triangles']} triangles, "
        f"{counts['edges']} boundary/interface edges, {counts['body_edges']} Body edges"
    )
    print("contacts: Source, Gate, Drain, Body")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
