# PMOS mask translation regression

The current source-built simulator completed eleven full PMOS processes:
a fresh control and ten mask translations of ±0.1 µm, using initial surface
spacing 0.0006 µm. Each translation preserves mask widths, vertical mask
coordinates, initial mesh lines, and process settings. The mask groups are
diffusion openings, gate oxide opening, contact openings, and metal etch
openings. Each group moves independently in two cases; two further cases
move all groups together.

All eleven processes returned zero. All 99 saved checkpoint meshes passed
the existing audit for finite coordinates and solution values, positive
triangle areas, valid vertex/region references, duplicate triangles, edge
ownership, and reciprocal neighbors. No simulator changes were needed for
these runs.

| Mask group shifted | Left: final oxide area change | Right: final oxide area change |
| --- | ---: | ---: |
| All groups | -0.3301% | -0.2622% |
| Diffusion | -0.2514% | -0.1636% |
| Gate | -0.0516% | -0.0789% |
| Contact | <0.0001% | <0.0001% |
| Metal | <0.0001% | +0.0102% |

Changes are relative to the fresh control's final oxide area of
50.3530564244 µm². Silicon area changes by less than 0.0073% in every case.
These measurements establish completion and the audited mesh properties
for the tested offsets. They do not establish physical accuracy or resolve
the older stored-reference oxide-area discrepancy.

Repeat with:

```sh
python3 suprem4/tests/run_mask_perturbations.py --jobs 2
```

The runner accepts `--binary`, `--spacing`, and `--shift` (positive magnitude
in micrometers), retains every generated deck, structure, and simulator log,
and writes aggregate and per-case JSON results. Coordinate replacement counts
are checked before simulation, and exact edited coordinates are recorded.

Results from this run:
`/tmp/suprem4_mask_perturbations_7lpeoase/results.json`.
Area measurements:
`/tmp/suprem4_mask_perturbations_7lpeoase/summary.json`.
