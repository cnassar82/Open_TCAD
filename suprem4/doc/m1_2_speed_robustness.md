# M1.2/M1.3: step count and LOCOS robustness with the geometry core (2026-11-26)

## Profile (`SUPREM4_PROFILE=1`)
One line per SUPREM time step: requested and taken step, the limit that set it (diffusion LTE
Ktd, grid motion Kto, the core's half-area limit, collision limit, SUPREM's back_pedal), and wall
time per phase (oxide flow, core write/run/read, diffusion, grid_sub).  `SUPREM4_FE_DIAG=1` adds
the viscous Newton history and, on a rejected step, the nodes with the largest free-dof residual.

| run (LOCOS unless noted) | steps | limited by the area limit | wall time split |
|---|---|---|---|
| viscous, core off | 196 | - | flow 91 % |
| viscous, core | 169 | 6 | flow 81 %, core 9 %, file exchange 1 % |
| viscous, core + 2 nm layer | 995 | 109 of 121 oxidation steps (partial run) | flow 54 %, core 32 % |
| poly `base`, core | 100 | 82 | flow 67 %, core 16 % |

Measured on the same machine load: viscous LOCOS with the core costs 875 s vs 815 s core off
(+7 %).  In-memory exchange (M1.3a) would save 1-3 %: not worth doing now.

## Step count: look-ahead removal (`UpdateOptions::step_dt`)
The core limits every SUPREM step so that no triangle loses half its area.  Where the interface
sweeps a mesh fixed in the consumed material, a point (or row, or chain vertex) close ahead of it
makes that triangle halve step after step until the core removes the point (0.82 -> 0.59 -> 0.40
-> 0.24 -> 0.14 nm at the poly corners): 4-5 steps per point.  SUPREM's own rule
(ChooseKillNodes) removes nodes that will be overrun within the step; the core now does the same,
with SUPREM's wanted step (`geometry_core_set_step`) and the node velocities (XVEL/YVEL in the
exchanged structure):
- interior points closer than near x size + 2 x their approach within the step are removed (and
  no point is inserted there);
- the first row of an interface-layer column starts below near x h0 + 2 x the interface advance;
- a chain vertex on a short segment (now or within the step) that will be within tol of its chord
  during the step is removed (a consumed corner flattening).  Long segments are never merged by
  the look-ahead: a flat chain is collinear at every vertex (tried: it stripped flat chains and
  the encroachment repair put the vertices back, 140 per call).
Poly `base`: 100 -> 66 steps.  `SUPREM4_CORE_LOOKAHEAD=0`: off.  The remaining area limits are
true corners (gate edge, interface junctions), not removable within the geometry tolerance.

## LOCOS failure (0.5 um NMOS, converged deck, core on)
SUPREM's viscous flow (stress-dependent viscosity) did not converge 60 s into the wet field
oxidation ("finite element continuation failed"); plain SUPREM passes.  Reduced to a regression
cell (`tests/fixtures/locos_nmos_cell.in`: x = 1.0-2.5 um, reflecting at the active nitride and the
window centre; 2300 points, fails in 3 min).

Root cause: needle fans in the oxide.  The core kept the oxide at its floor size (grid.oxide =
40 nm: one element through the field oxide) while the interface below was divided every 3-12 nm
(the silicon side's 2 nm grid): a coarse region on a finely divided boundary can only connect by
fans of needles (worst 2.5-3 deg; plain SUPREM: 10.7 deg).  The flow's Newton stalls on them.  Two
sources, both fixed:
1. Encroachment splits at the midpoint cascade: a silicon grid column 2 nm below a 20 nm interface
   edge encroaches one half again after each split (4 halvings, 2.6 nm segments).  Now the split
   is at the foot of the perpendicular from the encroaching apex (inner 60 % of the edge), which
   resolves it with one split (Ruppert/Shewchuk practice).
2. Flow quality: oxide/nitride triangles below 15 deg (`flow_min_angle`,
   `SUPREM4_CORE_FLOW_MIN_ANGLE`) get the midpoint of their longest interior edge with a clearance
   relative to their own size, not the material floor: longest-edge bisection with
   re-triangulation, i.e. Delaunay refinement graded from the boundary to the material's size.
   The silicon layer zone (anisotropic by design) is not touched.

Cell result: core on completes with no continuation cut-back, with or without the recovery below;
field oxide 7110 vs 7134 A core off (-0.3 %), bird's beak at the nitride edge equal within 1 A;
worst oxide angle 8.9 deg (core off: 2.4 deg); 1665 oxide triangles (core off 145); wall 525 s vs
388 s.

## Safety net: the flow never ends a run on one failed solve
- With a converged velocity field to fall back on (warm start), a continuation failure is not
  fatal: the step keeps that field (as between SUPREM's own flow re-solves, oxide_redo), is
  halved, the core rebuilds the mesh at its next call, and the flow is re-solved next step;
  9 failures in a row end the run as before (`SUPREM4_FLOW_RECOVERY=0`: off).
- With the core: a solve that needs more than 4 continuation cut-backs, or a Newton that stalls
  (3 iterations with line-search steps below 1e-3 and < 0.1 % residual decrease), gives up at
  once: a creeping continuation had taken 25 min for 2 s of process time.
Verified by replaying the old failing state with the flow quality off: gives up after 5 cut-backs
and continues with the previous field.  Core off: the paths only change runs that would have
panicked (bit-identical otherwise).

## Full-size NMOS LOCOS
The converged NMOS deck's LOCOS (`experiments/m10b_mesh_20261124/locos_fail/locos.in`: dry 0-300 s,
wet 300-4200 s, 1100 C, viscous, stress-dependent) completes with the core: 1568 steps, 1401 core
calls (1 full regeneration), 18.2k points; 5 continuation cut-backs, 38 Newton stalls accepted on
residual, no flow recovery and no solve given up.  Field oxide 7111 A (regression cell core off:
7134 A, core on 7110 A).  Oxide: 8333 triangles, worst angle 5.7 deg (5 below 10 deg).  Worst
coupling of the last call -0.93.  Wall 21.4 ks on a loaded machine (with the suite running).

## Suite
58/58 on the final build (`tools/run_geometry_suite.py --matrix --random 50 --seed 7 --locos`,
`experiments/m1_final_suite_20261126`).  LOCOS: viscous 7167 A vs core off 7157 A (+0.14 %), 801 s
vs 1116 s; compress (now at grid.oxide 0.01, see below) 7118 A vs 7099 A (+0.27 %), 3146 s vs 3086 s.

`locos_compress` at grid.oxide 0.04 had failed the 1 % check (7167 vs 7092 A).  Oxide spacing study,
field oxide in A:

| grid.oxide (um) | core off | core on |
|---|---|---|
| 0.04 | 7092 | 7167 |
| 0.02 | 7093 | 7179 |
| 0.01 | 7099 | 7118 |

The core's compress answer is not converged at 0.04 (it moved 7147-7185 A between builds this
week); at 0.01 on and off agree to 0.27 %.  The case now runs at 0.01 (`LOCOS_GRID_OXIDE`); the
tolerance is unchanged.

## Open
- The last oxidation step's interface motion is not cleaned up (the core runs at the start of a
  step): one silicon interface edge of the cell ends at coupling -1.2 (core off: -0.73).  A final
  core call after the last oxidizing step (M1.4).
- Point churn from the flow quality refinement: the interior removal (3b) measures clearance
  against the flat oxide floor and removes most of the graded points that step 4 inserts, which
  step 4 then inserts again (regression cell: 401 of 635 calls change > 100 points, ~250 each;
  offline, repeated updates of a fixed structure: 149 removed / 118 inserted per call, 2 with
  `flow_min_angle` = 0).  Likely most of the cell's 3.4x runtime vs core off.  Tried: one graded
  floor (boundary spacing + 0.5 x distance) for both, removal at 0.75 x: churn gone offline
  (worst oxide angle 15 deg) but in the run the oxide coarsened during growth, 316 cut-backs and
  36 recoveries (`experiments/m10b_mesh_20261124/locos_cell/graded_floor_try1.patch`).  Not
  merged.
  Merged (2026-11-27): an exact cavity test in 3b.  A flow-material point that 3b would remove
  for proximity is kept when the best re-triangulation of its star (max-min-angle DP over the
  star polygon, up to 24 vertices) would fall below `flow_min_angle` and below the star's own
  minimum angle: the removal must not create what step 4 then repairs.  Regression cell: wall
  1337 -> 730 s, removed points 100k -> 50k, calls changing > 100 points 401 -> 141, 0 cut-backs,
  field oxide 7103 A (core off 7134 A), worst oxide angle 13.1 deg.  Suite 58/58
  (`experiments/m1_cavity_suite_20261126`; LOCOS viscous 7171 A / 772 s, compress 7130 A /
  2974 s).  Tried on top: removing only an independent set per call (no two neighbours): 12 %
  less churn but 6 cut-backs and worst angle 9.1 deg; not merged
  (`locos_cell/cavity_plus_indep.patch`).  The remaining removals (~200 oxide points per late
  call, ~10 nm from the interface) track the growing oxide.
- A final core call after the last oxidation step was tried: no effect on the cell (worst
  coupling -0.536 with or without; the -1.2 was an earlier build).  Not merged.
- Layer runtime: LOCOS with the 2 nm layer still needs to be re-measured with the look-ahead.
