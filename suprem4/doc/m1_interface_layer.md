# M1 interface boundary layer (2026-11-25)

Origin: M1.0(b) (`m1_0_core_on_off.md`).  Near-surface phosphorus after the sacrificial (Kooi)
oxidation was still moving by 9 % between the 0.7x and 0.5x production meshes, because the pile-up
sits in the first 20 nm of silicon and the interface moves through the mesh.  Decision (user,
2026-11-25): the core keeps a boundary layer that moves with the Si/oxide interface (option A of
the first attempt: a structured zone owned by the core).

Opt-in: `UpdateOptions::interface_h` (um, 0 = off), from SUPREM via `SUPREM4_CORE_INTERFACE_H`,
`SUPREM4_CORE_INTERFACE_DEPTH` (default 0.05), `SUPREM4_CORE_INTERFACE_GROWTH` (default 1.5).
Recommended value: `SUPREM4_CORE_INTERFACE_H=0.002`.  The core only runs on vertical-model steps
with `SUPREM4_CORE_ALL_MODELS=1`.  Code: `core/src/remesh.cpp` (`class Interface`, steps 1, 3-, 3a,
the zone rounds of steps 4-5).

## Design
- **Zone**: the consumed side (silicon, poly) of every interface with SiO2.  Its point set is
  exactly a set of targets: one column below each interface vertex, along the vertex normal, at
  levels h0 uniform to `interface_depth`, then growing by `interface_growth` per row (rounded to
  h0) until the spacing reaches the material's own finest size there.  Along the interface the
  column spacing is the interface vertex spacing, so elements are anisotropic and the point count
  scales with the number of rows, not with 1/h^2.
- **Lagrangian**: an old point standing for a target (within 0.25 spacings in the anisotropic
  metric, along the interface at most 0.1 x its depth, same region) is kept with its values; other
  old points in the zone are removed; unmatched targets are inserted and interpolated.  Rows keep
  their phase: a column's levels are aligned to the existing old points under its vertex (absolute
  lattice x.n only for new columns).  So between calls only the rows the interface eats and the
  odd row at the deep end change.
- **Conservation**: unchanged, `update_fields` restores each dopant's dose per material over the
  changed nodes.  A per-patch restoration (each connected cavity of the update gets back its own
  old inventory) was tried and removed: it was needed only while new interface vertices lacked
  their column (fixed structurally, see below), the 1-D result is the same without it (P @ 10 nm
  1.7722 vs 1.7715e17), and with the layer off it moved boron next to the field-oxide interface
  enough to change compress LOCOS by 0.5 % (7185 vs 7147 A) and to make rand7_05 fail (coupling
  -3.6, then SUPREM's retry budget).

## Rules that each fixed a real failure
| rule | failure without it |
|---|---|
| Zone membership uses the same depth as the row placement | 200-225 points removed and re-inserted at the same coordinates every call |
| Interface vertices created by the encroachment repair get their column in the same call (steps 4-5 loop, up to 3 rounds) | fan triangles under the new vertex; next call interpolated a linear ramp, the per-patch restoration then scaled the bulk by 0.66 |
| Chains that end on the interface (reflecting side, Si/nitride) get vertices on the row levels, and the chain coarsening (step 1) leaves them | edge fans at the chain end drained the pile-up into the first column (3.5e18 at x = 0.05) |
| Row phase from the column's old points | the absolute lattice rotates with the interface in a bird's beak: ~500 points churned per call in LOCOS (15-40 with the rule) |
| Kept points lean at most 0.1 x depth along the interface | LOCOS: first-row points kept 3.7 nm off their column at 2 nm depth, interface edges with an opposite angle of 148 deg, couplings to -0.79 (47 edges) |
| h0 >= 2 nm | SUPREM's own grid_sub still runs with the core on and kills edges below 0.5 nm (abstol 5e-8 cm) or closing by half in one step; at h0 = 1 nm it removed ~10k zone points and the profile became laterally uneven |

## 1-D acceptance (flat wafer, Kooi recipe of the 0.5 um flow, n-well P 6.39e16)
Deck: 825-900 C dry/wet oxidation (~1000 A), N2 anneals, core on all models, `grid.oxide=0.005`,
`oxide.gdt=0.016`.  P at 10 nm below the final interface:

| run | P @ 10 nm | surface (0.3 nm) |
|---|---|---|
| zone 2 nm, grid.oxide 0.01 / 0.005 / 0.0025 | 1.713 / 1.719 / 1.722e17 | 2.32-2.34e17 |
| zone 1 nm (grid_sub kills, see above) | 1.724e17 | 2.35e17 |
| core off, 2 nm start mesh, any grid.oxide | 1.762e17 | 2.58e17 |
| zone 2 nm, no pre-oxidation N2 anneal: 10 nm / 2 nm start mesh | 1.7715 / 1.7733e17 | 2.41e17 |

- Converged to 0.5 % in the oxide spacing, 0.3 % in h0, 0.1 % in the start mesh.  Lateral
  uniformity 0.6 %, no coupling below -0.05.
- **Core-off references are biased.**  SUPREM's vertical model never adds nodes inside grown
  oxide (`grid.oxide` has no effect with the core off), so every core-off run has one element
  across ~1000 A of oxide.  That overstates the back-flux from the oxide into silicon during the
  following inert anneals: +10 % at the surface, +2.3 % at 10 nm.  The converged reference is the
  core with a refined oxide, not plain SUPREM on a fine mesh.
- With the pre-oxidation N2 anneal, a 10 nm start mesh stays 3.5 % low: that anneal runs on the
  user's mesh before any oxide exists, so the core (and the layer) is not active yet.  Not a layer
  issue; the deck needs a fine initial surface mesh, or the core needs to run on inert anneals.

## LOCOS (suite case `locos_viscous`, 2 nm layer)
| run | field oxide | boron drift | negative couplings (< -0.05) | wall time | core calls |
|---|---|---|---|---|---|
| core off | 7157 A | -5.8e-4 | 8 | 1023 s | 0 |
| core, no layer | 7161 A | 5.4e-5 | 10 (worst -0.26) | 1033 s | 71 |
| core, layer (before lean cap) | 7111 A | 2.4e-4 | 47 (worst -0.79) | 8412 s | 1024 |
| core, layer (final: lean cap, per-material restoration) | 7112 A | 9.0e-5 | 9 (8 silicon, worst -0.068; 1 oxide, -0.096) | 9445 s | 842 |

Field oxide -0.6 % vs core off (suite limit 1 %), dose drift within the strict 1e-4 criterion.
The run with the layer takes 5.5x more SUPREM time steps (995 vs 181) and each step costs about
1.3x (10.8k vs 7.2k points).  The finer mesh shortens SUPREM's own time steps.  That is an M1.2/M1.3
item: the layer is opt-in and meant for the oxidations where near-surface doping matters, not for
the field oxidation.

The worst coupling of -0.45 from call 1 onwards is in the input structure (pad oxide and nitride
of the LOCOS mask, interface edges up to 1900 A long with the oxide apex at 132 deg), with or
without the layer: an M1.4 item.  The layer's own defect was the 47 silicon edges (lean, fixed).

## Open
- Pre-oxidation inert anneals are not covered (core not active).
- Runtime: the layer multiplies SUPREM's step count in 2-D viscous LOCOS (see above).
- Electricals of the 0.5 um flow with the layer on (the PMOS buried channel is the sensitive
  device) after the suite passes.
