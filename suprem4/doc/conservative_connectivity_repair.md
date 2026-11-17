# Connectivity changes must transfer chemical fields

Status, 2026-11-16: the captured shared-interface refinement passes independent
partition and chemical-inventory checks in ordinary and ASan/UBSan builds.
Full RIT CMOS is still failing in source/drain activation. This is a bounded
operation repair, not a completed CMOS validation or launcher promotion.

The subsequent minimal deposit/etch failure and time-step comparisons are
recorded in [deposition_etch_construction_identity.md](deposition_etch_construction_identity.md).

## Reproduced cause

`vertical_constraints.c` refines the shared oxide/poly interface to sample the
same vertical deformation as the more finely resolved oxide free surface.
Subdivision prepares bounded chemical fields using the complete P1 supports of
the affected material nodes. A shared material node can span multiple legacy
region labels; it has only one concentration across that bookkeeping seam.

In `astra_constraint_connect_trace`, independent integration finds exactly
unchanged boron and phosphorus inventories immediately after subdivision.
`bd_connect()` then calls `geom()`, which automatically changes Delaunay
diagonals. The old `tri_flip()` changed nodal integration weights without
transferring concentrations. Polysilicon boron inventory fell by
0.001319434635627 relative and phosphorus by 0.000759856056911. Both meshes
passed the partition audit. Positive geometry alone missed the corruption.

## Implemented contract

For a fixed-coordinate rezone, `tri_flip()` prepares its field transfer before
mutating connectivity. It computes each of the four affected nodes' complete
old and new P1 supports, including unchanged incident triangles. A bounded
least-change projection preserves the inventory of each transported chemical.
Original local minimum and maximum concentrations bound the result. An affine
field already integrates to the same inventory on both triangulations and is
left unchanged. An infeasible transfer refuses the operation before mutation.
Active concentrations retain their pre-transfer active/total fractions.

`include/conservative_transfer.h` shares that projection and species semantics
with shared-interface insertion. STEADY/PSEUDO algebraic variables, trap
fractions and insulator point-defect placeholders are not transported chemical
inventory. They must not be constrained as concentration-storage equations.

This is a local lumped-mass constrained projection, not an implementation of a
supermesh Galerkin projection. It preserves the stated inventory, bounds and
affine-field invariants; it does not establish a global convergence order for
the complete adaptive process simulator.

Implicit moving-mesh trials are a separate case. `time.c` carries old/mid/new
fields with their corresponding mass terms and commits `nd[].sol` only after
success. `mesh_trial_active()` identifies this interval. A trial flip leaves
committed nodal fields untouched; static rezoning transfers them. Mixing those
two responsibilities would corrupt the retry input without transferring the
integrator's actual history arrays.

## Evidence

- `tests/run_mesh_trial.sh` exercises the real connectivity primitive under
  ASan/UBSan, checking a non-affine chemical inventory, bounds, affine-field
  preservation, committed fields during an implicit trial and rollback.
  The old flipper from `native64__63eakrl` fails the new inventory assertion.
- `tests/constraint_transfer.c` checks the shared projection, constants,
  infeasible bounds and algebraic-field semantics.
- `tests/inventory_audit.py` independently integrates saved triangles and
  chemical fields by material. The vertical-constraint replay now requires
  this check; a simulator exit of zero cannot hide a failed transfer.
- `astra_constraint_conservative_flip` (ordinary `native64_a9xbvu_4`) and
  `astra_constraint_conservative_flip_sanitized` (`native64_bxuopbdh`) both pass
  the 6-second captured replay and independent partition/inventory checks.
  All six material/species inventory differences are below output resolution.
- The 13 existing focused numerical tests pass; evidence is in
  `experiments/conservative_connectivity_20261116/focused_checks.json`.

## Remaining failure and next diagnostic

Both devices resumed from unchanged source/drain checkpoints in
`astra_conservative_constraints_activation` fail before `annealed.str`:
NMOS at 6.0304226086023833 seconds of the 825 C oxidation command and PMOS at
20.929270441409475 seconds. Retries report TR/BDF solves as well as preceding
rounded-motion rejections. NMOS's saved failed-step mesh passes the independent
partition audit. Restarting that saved state with an estimated small startup
step runs a further three seconds, so the saved geometry alone does not explain
the original failure. Capture the actual failing solve and history, rather than
attributing every failure to mesh geometry.

Some attempted constraint subdivisions are correctly refused because the
rounded point produces a nonpositive child in an extremely thin layer. Keep
that guard. Do not remove the oxide remnant, loosen conservation, expand retry
budgets or accept crossing boundaries to manufacture a completed run.

## Research basis

[Farrell et al., conservative interpolation via supermesh construction](https://www.sciencedirect.com/science/article/pii/S0045782509001315)
describes why ordinary interpolation destroys conservation during adaptation
and develops a bounded conservative transfer for linear elements.
[Alauzet et al., mesh adaptation and P1-conservative interpolation](https://pages.saclay.inria.fr/frederic.alauzet/download/Vanharen_Non-manifold%20anisotropic%20mesh%20adaptation%20application%20to%20fluid-structure%20interaction.pdf)
sets conservation, P1 exactness and bounds as requirements across insertion,
collapse and swapping. Those invariants guide the tests above.

[CGAL's corefinement documentation](https://doc.cgal.org/6.1.2/Polygon_mesh_processing/index.html)
explicitly distinguishes exact predicates from exact constructions: rounded
intersection coordinates can introduce self-intersections even with correct
predicates. [IuE's process mesh discussion](https://www.iue.tuwien.ac.at/phd/cervenka/node14.html)
requires moving-mesh algorithms to prevent overlap of thin material layers.
Neither source justifies deleting a finite material layer or accepting an
invalid rounded subdivision.
