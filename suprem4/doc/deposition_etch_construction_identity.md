# Depositing and etching the same planar layer must cancel

2026-11-16 investigation. The original defect and evidence are preserved below.
An opt-in implementation now passes the bounded cancellation and retained-film
tests: see [construction repair](process_surface_construction_repair.md).
Full fresh CMOS remains under validation. A minimal planar reproduction does not involve oxidation,
Newton iteration, a mask corner, or a deforming mesh.

Follow-up research: [how other TCAD packages handle these issues](tcad_repair_methods_research_20261116.md).
It distinguishes analytic construction identity from resolution-controlled
implicit geometry, and gives an implementation and validation order.

## Reproduction

```
python3 suprem4/tests/run_deposit_etch_identity.py \
  --binary suprem4/build/native64_n6whaa5l/src/suprem4 \
  --output /tmp/unique_planar_identity_result
```

The command deposits 0.5 microns of oxide on a flat, 12-micron-wide silicon
surface and dry-etches 0.5 microns. It repeats at three vertical positions.
The input and result partitions are valid; the failure is retained material.
The test fails on any positive remaining oxide area. It does not delete a
region or apply an area threshold to manufacture a passed result.

Evidence: `experiments/conservative_connectivity_20261116/planar_identity_gate`.

| Initial surface y (microns) | Remaining oxide area (square microns) |
| --- | --- |
| 0 | 3.929017298277349e-6 |
| 0.1 | 8.731149551727442e-6 |
| -0.427 | 1.3096724327591163e-6 |

## Identified construction inconsistencies

The domain offset backend chooses a dyadic integer grid from the largest
coordinate magnitude, quantizes the source curve, offsets that curve, then
retains the original unquantized source as the inner material interface.
For the zero-height test, the grid is 1.1641532182693481e-6 microns. The stored
coating area is 6.0000040624292375 square microns rather than the analytic
6 times the stored width correction. Dry etching translates that already
rounded outer surface. The excess becomes a finite oxide film.

There is a second inconsistency: deposition converts thickness to cm in a
double, while dry etching stores the converted thickness in a float. Equal
nominal recipe parameters therefore need not produce equal translated lengths.
Binary32 storage can also lose source-coordinate information when a coating
moves a boundary to a larger coordinate magnitude. Merely changing the dry-etch
local variable to double cannot generally restore information already lost.

A legacy-offset control, using the same engine with the domain backend disabled,
passes the zero-height and +0.1 cases but leaves 4.365574775863721e-7 square
microns at -0.427. Its output is in `planar_identity_legacy_control`. Reverting
to that backend is not a general construction-identity fix.

These experiments prove the planar artifact. They do not prove that every thin
oxide fragment in RIT is artificial. A conformal coating over a sloped surface
can have a real remnant after a vertical etch; provenance must distinguish it.

## Required direction

Follow the CAD separation between authoritative geometric constructions and
an approximate simulation mesh. Preserve source facet identity and nominal
offset/translation parameters through deposition and etching. Build shared
material constraints once from those constructions. A cancellation can then
be established from its construction, rather than inferred from a small area
or a proximity tolerance after meshing.

Research clarification: this provenance-based construction is a recommendation
for our analytic operations, not the only architecture used in TCAD. Implicit
material surfaces with thickness-aware refinement are another documented route.
An approximate geometry mode requires declared error bounds and convergence;
the strict zero-remnant gate above applies to the intended analytic planar
identity and is not a universal bitwise-identity claim about level-set solvers.

Coordinate conversion must have a declared error model. If an exact or higher
precision construction cannot be embedded in legacy binary32 coordinates while
preserving the required material ordering, report that limitation or use a
higher precision representation. Do not accept a crossing, a zero-area cell,
an arbitrary enlarged etch depth or deletion of a finite fragment.

Retain the independent partition and chemical-inventory gates for every
connectivity change. Exact predicates on already-rounded points cannot recover
the intended construction; CGAL explicitly documents this distinction in its
[corefinement manual](https://doc.cgal.org/6.1.2/Polygon_mesh_processing/index.html).

## Time-step experiments are separate

`SUPREM4_POSITIVITY_STEP=1` is an experimental, default-off predictor which
re-estimates chemical depletion after connectivity changes, without resetting
potentials or defects. It derives a startup estimate from M*C and the current
loss flux. It is not a proof of nonlinear TR-BDF2 positivity and does not replace
trial rejection or error estimation. The old unchanged retry budget remains.

`astra_constraints_rezone_depletion` advances NMOS beyond the previous failure,
to 1822 seconds of the continued oxidation schedule, but both device runs time
out at 600 seconds of wall time without an annealed checkpoint. PMOS reaches
about 81 seconds of the first 825 C oxidation. Always-on prediction and a
first-order implicit comparison also time out. None is an accepted full flow.
This numerical progress does not repair the planar construction identity.

[TR-BDF2 positivity analysis](https://arxiv.org/abs/1510.04303) establishes that
L-stability does not imply unconditional positivity and analyzes hybrid implicit
variants. Use that theory when designing the integrator, rather than treating
negative intermediate concentrations as evidence for an unlimited retry tree.
