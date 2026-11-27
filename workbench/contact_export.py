"""Prepare explicitly assigned SUPREM4 contacts for the existing Genius exporter."""

from dataclasses import replace
from math import isfinite
from pathlib import Path
import shutil
import sys
import tempfile

from electrical import IDENTIFIER
from runner import ROOT

VIEWER_DIR = str(ROOT / "suprem_viewer")
if VIEWER_DIR not in sys.path:
    sys.path.insert(0, VIEWER_DIR)

from make_genius_contacts import (  # noqa: E402
    SUPREM_MATERIAL_TO_TIF, build_edges, write_deck, write_tif,
)

TERMINALS = {
    "MOSFET": (("Source", "ohmic"), ("Drain", "ohmic"), ("Gate", "gate"), ("Body", "ohmic")),
    "Diode": (("Anode", "ohmic"), ("Cathode", "ohmic")),
    "BJT": (("Emitter", "ohmic"), ("Base", "ohmic"), ("Collector", "ohmic")),
    "Custom": (),
}
CONDUCTORS = {4, 6}


def contact_mesh(triangles, regions, contacts):
    if not contacts:
        raise ValueError("Assign at least one contact before generating a Genius deck.")
    unsupported = {region.mat_id for region in regions.values()} - SUPREM_MATERIAL_TO_TIF.keys()
    if unsupported:
        materials = ", ".join(map(str, sorted(unsupported)))
        raise ValueError(f"The Genius exporter cannot translate SUPREM4 material IDs {materials}. "
                         "Use a device checkpoint with supported materials.")
    names = [contact.name.casefold() for contact in contacts]
    if any(not IDENTIFIER.fullmatch(contact.name) for contact in contacts) or len(set(names)) != len(names):
        raise ValueError("Contact names must be unique Genius labels, without spaces.")
    edges, numbers, owners = build_edges(triangles)
    by_number = {numbers[edge]: edge for edge in edges}
    links = {index: set() for index, region in regions.items() if region.mat_id in CONDUCTORS}
    for shared in owners.values():
        conductive = shared & links.keys()
        for index in conductive:
            links[index].update(conductive - {index})

    group_of = {}
    for index in links:
        if index in group_of:
            continue
        stack = [index]
        group_of[index] = index
        while stack:
            for other in links[stack.pop()]:
                if other not in group_of:
                    group_of[other] = index
                    stack.append(other)

    assigned_groups = {}
    selected_edges = set()
    for contact in contacts:
        if contact.kind not in {"ohmic", "gate"}:
            raise ValueError(f"Unsupported contact type for {contact.name}.")
        if contact.kind == "gate" and (not isfinite(contact.work_function) or contact.work_function <= 0):
            raise ValueError(f"Enter a positive gate work function for {contact.name}.")
        if contact.target == "region":
            group = group_of.get(contact.region)
            if group is None:
                raise ValueError(f"Select a metal or polysilicon region for {contact.name}. "
                                 "Use exposed edges for a direct semiconductor contact.")
            if group in assigned_groups:
                other = assigned_groups[group]
                raise ValueError(f"{other} and {contact.name} are on the same connected conductor.")
            assigned_groups[group] = contact.name
            continue
        if contact.target != "boundary" or contact.kind == "gate":
            raise ValueError(f"Assign a gate to a metal or polysilicon region for {contact.name}.")
        if not contact.edges:
            raise ValueError(f"Select at least one exposed semiconductor edge for {contact.name}.")
        for number in contact.edges:
            adjoining = owners.get(by_number.get(number), ())
            if len(adjoining) != 1 or regions[next(iter(adjoining))].mat_id != 3:
                raise ValueError(f"{contact.name} includes an edge that is not exposed silicon.")
            if number in selected_edges:
                raise ValueError(f"An edge for {contact.name} is already assigned to another contact.")
            selected_edges.add(number)

    unassigned = sorted(set(group_of.values()) - assigned_groups.keys())
    if unassigned:
        ids = ", ".join(map(str, unassigned))
        raise ValueError(f"Assign a terminal to each conductive island before export. Unassigned regions: {ids}.")

    representatives = sorted({group_of.get(i, i) for i in regions})
    renumber = {old: i for i, old in enumerate(representatives, start=1)}

    def region_number(index):
        return renumber[group_of.get(index, index)]

    new_regions = {number: replace(regions[index], idx=number)
                   for index, number in renumber.items()}
    new_tris = [replace(tri, region=region_number(tri.region)) for tri in triangles]
    _, new_edge_numbers, _ = build_edges(new_tris)
    new_contacts = [
        replace(contact, region=region_number(contact.region)) if contact.target == "region"
        else replace(contact, edges=sorted(new_edge_numbers[by_number[edge]] for edge in contact.edges))
        for contact in contacts
    ]
    return new_tris, new_regions, new_contacts


def export_device(destination, source, nodes, triangles, regions, codes, rows, contacts):
    triangles, regions, contacts = contact_mesh(triangles, regions, contacts)
    folder = Path(tempfile.mkdtemp(prefix=f"{source.stem}-device-", dir=destination))
    tif = folder / "device.tif"
    deck = folder / "device.inp"
    try:
        write_tif(tif, nodes, triangles, regions, codes, rows, body_tol=1e-6, contacts=contacts)
        write_deck(deck, tif, contacts)
    except (OSError, ValueError, RuntimeError):
        shutil.rmtree(folder)
        raise
    return deck, tif
