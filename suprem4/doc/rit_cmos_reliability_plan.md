# RIT CMOS reliability implementation plan

Implementation handoff: [Sol’s mesh repair guide](sol_mesh_repair_implementation_guide.md).

Prepared 2026-11-14 for implementation by Sol. This is a plan, not a claim that
any candidate is ready for use. Preserve the current dirty workspace and the
validated launcher while carrying it out.

For the immediate spacer failure, follow
[Sol's spacer debugging guide](rit_spacer_debugging_guide.md). It gives the
captured geometry, competing hypotheses, instrumentation points and the next
bounded acceptance gate.

## Implementation progress, 2026-11-15

Bounded NMOS and PMOS spacer replays now pass strict ASan/UBSan (full-engine
leak checking remains open). PMOS resumed heavy-implant polarity also passes.
The original startup parser UB is repaired; earlier statements below and in
archived results describe historical candidates. See the current
[repair report](rit_spacer_boundary_repair.md) for evidence and scope.

Full RIT CMOS is not accepted. Activation still breaks mesh invariants, and
NMOS body phosphorus transport/polarity remains unresolved. The experimental
coarsening rewrite was withdrawn after an exposed-point/GAS-node failure,
even though its panic mesh passed independent geometry checks. The validated
launcher remains unchanged.

The immediate implementation priority is the
[mesh replacement contract](mesh_patch_replacement_contract.md). Replace the
complete multi-material removal operation with a validated atomic proposal;
do not continue adding downstream cleanup guards. Geometry, ownership and
conservative field transfer must pass independently before fresh full flows
and the broader robustness matrix can establish acceptance.

## Objective and scope

Run the separate NMOS and PMOS process cuts derived from
`rit_process/SubCmos2014_0p5um.pdf`, from substrate initialization through first
patterned metal, with valid material geometry, defensible dopant profiles,
bounded failure handling, and demonstrated mesh and mask robustness.

Use the existing process translation in `decks/subcmos2014/README.md` as the
recipe ledger. The document's 0.5 um is layout lambda; the drawn gate is 1 um.
Mask coordinates are assumed, not reconstructed mask artwork. Preserve the
specified species, energies, doses, and furnace budgets during engine debugging.
Do not alter the recipe to make a numerical failure disappear.

Two acceptance levels must be reported separately:

* Engineering reliability: correct execution and numerical convergence for the
  declared models and idealized layout through first metal.
* Process fidelity: validated oxidation, segregation and diffusion models against
  independent references or measurements. Matching nominal PDF thicknesses alone
  is not calibration. Electrical I-V validation is a later project.

The current vertical post-LOCOS oxidation and flat resist-top approximations
must remain visibly labeled. Full viscous oxidation and conformal deposition
are separate required reliability cases, not silently substituted by simpler
models when a test fails. Shared CMOS layouts, second metal, and 3D are outside
this first milestone.

## Evidence at handoff

Read result JSON and logs rather than relying on old narrative status files.
`subcmos2014_reliability.md` and its validation JSON predate recent candidates.

| Evidence | Established result | Limit |
| --- | --- | --- |
| `decks/subcmos2014/results/exposure_rebuilt/{nmos,pmos}/result.json` | Both pass through spacers on `native64_h7qfs_d8` | Resumed after LDD; not fresh full flows |
| `decks/subcmos2014/results/implant_bounds/nmos/result.json` | Abort, return code -6, missing source/drain checkpoint on `native64_1nq7dojj` | Must minimize the command/mesh failure |
| `decks/subcmos2014/results/implant_bounds/pmos/result.json` | Simulator exits 0 but drain polarity fails | B = 2.52e19, P = 1.80e20 cm^-3 at probe; not an accepted device |
| Recent implant quadrature and etch sanitizer tests | Focused tests pass | Not evidence for complete process correctness |
| Current `geom/misc.c` area determinant edit | Translated double arithmetic proposed for interpolation cancellation | Latest edit still needs a matching build and regression verification |

The installed `build/validated-native64` resolves to `native64_2r55t88h` at
handoff. Do not promote a candidate until the release gates below pass.
Temporary debugging artifacts may disappear; archive useful ones immediately.

## Execution rules

1. Work on one minimized failure and one hypothesis at a time. Record the failing
   command, source state, binary, data files, and checkpoint before editing.
2. Each fix needs a regression that fails with the relevant fix removed and
   passes with it restored. Add an independent oracle where correctness matters.
   A test that also passes the broken version is coverage, not proof of the fix.
3. Preserve user changes. Snapshot the current source/diff and build inputs before
   splitting fixes into reviewable units. Never reset this workspace wholesale.
4. Do not swallow errors, delete inconvenient small regions, widen tolerances to
   absorb invalid topology, or update golden output without explaining the change.
5. Use resumed runs for diagnosis. Only fresh runs from the same source/data
   snapshot count toward full-process acceptance.
6. Record outcomes as PASS, FAIL, or NOT RUN; include physical checks independently
   of process exit status. Keep a single manifest/report instead of inferring
   progress from the newest build-directory name.

## Work package 1 — Reproducible failures and independent validation

Primary files: `decks/subcmos2014/run.py`, `validate.py`,
`suprem4/tools/run_pmos_mesh_regression.py`, and `suprem4/tests/`.

Deliverables:

* Persistent small fixtures for the NMOS post-spacer abort, PMOS wrong drain,
  legacy PMOS photoresist/gate etch interpolation failure, thin boundary fragment,
  shared/tangent crossings, conformal spacer-topography offset failure, and
  viscous oxidation nonconvergence. Archive the predecessor structure and exact
  next command, not just screenshots or a whole multi-hour log.
* A manifest recording source snapshot/diff hash, compiler/options, executable
  hash, model/implant/key-file hashes, generated deck hash, environment, elapsed
  time, peak memory, and whether the run was fresh or resumed.
* Explicit paths for all model and implant data. Remove dependence on inherited
  environment defaults selecting different tables in different runners.
* Independent mesh audit: positive orientation, finite values, valid indices,
  reciprocal adjacency, manifold edge ownership, boundary loops, no crossings,
  no overlapping triangle interiors, and consistent material partitions.
  Positive triangle areas alone do not establish a valid mesh.
* Integrals of each chemical species by material and region, interface geometry,
  connected components, dielectric thickness profiles, and contact connectivity.
  Distinguish total chemical concentration from electrically active concentration.
* Deliberately corrupted fixtures to prove the audit detects missing faces,
  overlapping faces, wrong material labels, broken neighbors, and invalid fields.

Gate: one command reproduces each known failure with a useful stage report;
independent audit catches injected corruption. Set timeouts for every fixture.

## Work package 2 — Establish a coherent geometry contract

Primary files: `geom/misc.c`, `refine/skel.c`, `sp_edge.c`, `bfill.c`,
`triheur.c`, `triang.c`, `offset.c`, and `dbase/make_db.c`.

Order:

1. Verify the translated determinant fix with the captured point-on-edge case
   and equivalent translated/scaled geometries. Inspect all callers and return
   precision; fixing one determinant does not establish a consistent kernel.
2. Centralize orientation, point-on-segment, segment intersection and point-in-
   triangle predicates. Use double arithmetic and robust/adaptive signs where
   cancellation demands it. Keep numerical error bounds separate from intentional
   geometric snapping. Test large origins with small features.
3. Audit coordinate storage and I/O. Double intermediates cannot restore bits
   already lost in float storage. If storage resolution cannot represent the
   supported thin layers at supported extents, migrate geometry storage to double
   in a bounded change, auditing structs, prototypes, format strings and ABI uses.
4. Node all intersections before Boolean face traversal, including collinear
   overlap, tangencies, shared endpoints and multiple crossings on one edge.
   Assign shared intersection identity once; never move a junction after another
   arc has recorded it. Define behavior for point-touching material components.
5. Triangulate only the classified polygon interior; preserve its constrained
   boundary. Verify polygon-versus-triangle area and absence of crossings/overlap.
   Refinement must preserve the boundary, not progressively round it away.
6. Derive exposed interfaces from valid topology. Review the gas-node cleanup:
   it must not hide a malformed partition or remove a real exposed boundary.
7. Minimize the conformal-offset failure and test offset self-intersection
   resolution with analytic stepped and rounded profiles.

Known test weakness: the existing thin-boundary refinement fixture reportedly
passed even without the protective guard. Recover the actual graded-spacing
failure and demonstrate a failing control before accepting that guard as a fix.

Decision gate: if the existing arc walker cannot satisfy these contracts without
case-specific exceptions, isolate polygon Boolean operations behind a narrow API
and evaluate a robust planar-arrangement implementation. Require deterministic
boundary/material output and a field-transfer strategy before replacing it.
Do not rewrite unrelated diffusion code as part of this decision.

Gate: all minimized geometry cases pass independent audits and sanitizers;
translations, traversal reversal, and harmless edge subdivisions preserve the
same physical partition. No material-area discrepancy is dismissed solely
because the picture looks plausible.

## Work package 3 — Bounded and conservative implantation

Primary files: `implant/imp_qeqv.c`, `imp_zeqv.c`, `imp_qpprs.c`, `pearson.c`,
`implant.c`, `include/implant.h`, and implant table readers.

* Review the new clipped quadrature and inverse together with profile
  normalization. They must use the same domain and final partial interval.
* Define finite input checks, valid zero/full-dose behavior, monotonic inverse,
  bounded iteration, and explicit error propagation. A returned NaN must not
  become a field value or a successful process step.
* Verify material traversal, shadowing and mask penetration using planar stacks,
  partially masked implants, and thin layers. Test BF2 table identity and units.
* Test analytic constant/linear profiles, sub-grid layer thickness, tail doses,
  round trips and integrated deposited dose. Account explicitly for ions stopped
  in resist/oxide and particles leaving the represented domain.
* Preserve state on failed implantation. Add failure injection to prove it.

Gate: analytic tests pass with stated error bounds, no hangs/invalid fields,
and integrated species accounting agrees with the declared implant model.

## Work package 4 — Explain the wrong dopant profiles

Do this before claiming the transistors work. Do not tune implants to force the
expected polarity.

1. Add checkpoints immediately before/after every implant and anneal around poly
   doping, gate patterning, LDD, body tie and source/drain formation.
2. Locate the first unexpected phosphorus increment in the NMOS body and PMOS
   drain; compute integrated inventories and material-resolved profiles.
3. Check for unintended poly/silicon contact, missing oxide, wrong implant-window
   geometry, region traversal errors, field transfer errors and interface flux.
4. Reproduce poly/oxide/silicon transport in an independently checked 1D stack
   with the same coefficients, oxide thickness and thermal history. Thin body
   oxide may permit real penetration; do not assume every wrong polarity is a
   solver bug. Establish which mechanism explains it quantitatively.
5. Refine across the oxide and halve time steps independently. Compare full
   profiles, integrated transfer and junction locations, not a single probe.
6. If the declared model genuinely predicts the wrong device, report a model or
   recipe-translation discrepancy and resolve it with source evidence. Keep it
   separate from engine reliability and from experimental calibration.

Gate: the phosphorus source and transport path are explained; conservation and
convergence pass; source/drain/body polarity and junction checks pass throughout
specified sampling regions, not only at three convenient points.

## Work package 5 — Transactional process steps and oxidation

Primary files: `dbase/mesh_trial.c`, `diffuse/`, `oxide/`, `finel/`, and process
command error propagation.

* Validate before committing a process operation. On failure preserve the last
  valid structure, report the command and reason, and return failure to the runner.
* Audit rollback scope. The current mesh-trial interface explicitly forbids
  adding, deleting or renumbering entities between save and restore; it is not a
  general transaction for etch/deposition/remeshing. Use a complete snapshot or
  build-and-swap operation for topology changes.
* Reject failed Newton solves or inverted moving meshes. Retry with bounded
  step reduction, a minimum step and an attempt limit; never continue silently
  after 'Giving up'. Restore concentrations, coordinates, connectivity, solver
  history and elapsed physical time consistently.
* Validate planar oxidation against the same model's 1D/reference solution, then
  LOCOS corners. Track consumed silicon and grown oxide with the model's volume
  relation; total domain area is not conserved during oxidation.
* Add local interface refinement based on geometry and concentration gradients.
  Specify supported minimum feature resolution; reject unrepresentable layers
  explicitly. Demonstrate time-step and mesh convergence for full viscous growth.

Gate: forced failures restore valid state and return nonzero; successful retries
agree with a deliberately smaller-step run. Both LOCOS and thin-oxide viscous
cases complete without changing models mid-run.

## Work package 6 — Full flows, perturbations, and scaling

First run fresh nominal NMOS and PMOS through all 16 existing checkpoints to
`metal1`. Then run the following staged matrix; avoid a full Cartesian explosion
until the inexpensive tests pass.

| Dimension | Required checks |
| --- | --- |
| Mesh density | Approximately 0.5x, 1x, 2x and 8x baseline point counts; report actual counts at each stage. Changing spacing by two is not doubling a 2D mesh. |
| Critical interfaces | Independently refine gate/body oxide, spacer edges, LOCOS corners and junction gradients; include deeper substrate where required by gradients. |
| Masks | Shift relevant masks left/right by 10 and 25 nm, individually; include a few combined worst-case overlaps. These are physical changes, so require smooth plausible changes, not identical outputs. |
| Coordinate origin | Translate complete structures and masks together without changing the physical problem; require invariant geometry/fields after translation is removed. |
| Small geometry changes | Perturb film thicknesses by +/-1%; test rounded/sharp corners and near-coincident boundaries within the supported resolution. |
| Time steps | Halve maximum steps for oxidation and critical anneals independently of mesh refinement. |
| Restart | Compare fresh and restarted runs at wells, LOCOS, gate and spacers on the same build and data. |
| Existing capabilities | Legacy RIT PMOS, LOCOS, rounded STI, mask-shift cases and focused sanitizer suites. |
| Resource growth | Record peak memory, wall time, mesh counts and sparse-system size; test allocation failure and checked index/size arithmetic. |

For changed meshes, compare fields on common physical probes/sections, not node
IDs. Compare material areas, dielectric profiles, LOCOS encroachment, junction
depths and integrated dopants. Write tolerances before evaluating the candidate:
use analytic/error bounds for geometry, discretization studies for PDE quantities,
and the represented boundary/source flux for mass balance. As an initial
engineering target, seek <=2% change between the two finest runs for nonzero
bulk metrics, with separately specified absolute bounds for tiny layers and
junction displacements. This is a proposed convergence target, not a proven
accuracy or calibration claim. If it fails, investigate rather than widen it.

Gate: nominal fresh flows pass, required variants have no crashes/hangs/invalid
meshes, and convergence trends support the stated accuracy. A coarse mesh may
fail a resolution criterion explicitly; it must not produce an accepted but
incorrect structure. Document the minimum supported resolution.

## Release and handoff

Produce one machine-readable results manifest and one concise human report with:

* exact accepted source/build/data identifiers;
* fresh/resumed provenance and PASS/FAIL/NOT RUN for every required case;
* structural, conservation, convergence and polarity results;
* resource scaling and remaining physical-model limitations;
* reproducible build/test commands and locations of final NMOS/PMOS structures.

Only then switch the launcher, verify its resolved binary/data, and run a short
smoke case through the launcher itself. Preserve the previous validated build
for rollback. Do not equate first-metal process acceptance with calibrated
transistor I-V predictions.

## Next assignment to Sol

The captured spacer replay and independent geometry gate now pass; see the
repair report before revisiting that diagnosis. Start work package 4 by saving
chemical-species inventories immediately before/after the poly implant and each
activation segment, then locate the first unexpected phosphorus transfer into
the NMOS body and PMOS active silicon. Use a planar poly/oxide/silicon reference
to distinguish real model transport from interface/field-transfer errors.
Preserve the recipe. The earlier area-determinant fix still needs its separate
captured failing-control test; startup parser UBSan and the other reliability
packages remain tracked work, not accepted results.
