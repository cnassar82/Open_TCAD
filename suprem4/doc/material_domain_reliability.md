# Material domains and coupled-solver reliability

Updated 2026-11-15. Candidate implementation; full RIT CMOS acceptance remains open.

The latest follow-up is [open-source TCAD research and implementation](open_source_tcad_repair_research.md).
The Clipper2 domain backend passes captured ILD/spacer geometry, ordinary and
sanitized, and completes the archived PMOS annealed-to-metal continuation with
16/16 checks and independent final partition audits. Fresh PMOS still fails
on a real nitride crossing during LOCOS; fresh NMOS fails at a reversing
spacer-side boundary. See the current table in the
[repair report](rit_spacer_boundary_repair.md). The launcher is not promoted.

The following investigation is retained as historical evidence.

The stable-weight candidate completes the unchanged PMOS
activation schedule from the archived source/drain checkpoint, in ordinary and
strict ASan/UBSan builds. The final ordinary partition audit passes 12,310
triangles and 23,679 candidate overlap pairs. Subsequent fresh-flow testing
exposed a straight-side removal gap, now implemented with ordinary/sanitized
captured replay passes. The preceding combined ordinary candidate also passes
the archived full activation schedule; the newer Kooi build has not rerun it.
Promoting coordinate operands before subtraction in `back_pedal` removes the
captured fresh-PMOS Kooi assertion. The new ordinary and sanitized candidates
pass the reduced Kooi replay, and an ordinary continuation from `locos.str`
reaches LDD. It then fails while clipping the spacer oxide boundary; fresh
NMOS body polarity and post-activation ILD remain open. Full CMOS is not
accepted. Historical failures below are retained as causal controls.

## Root cause and evidence

The spacer sequence deposits 0.5 um oxide and dry-etches 0.5 um. Its output
contains extremely thin oxide remnants recorded as separate regions from the
underlying oxide. These regions share material nodes and use the same current
material-global physics. Their region-number boundaries nevertheless constrain
the mesh as if they were physical interfaces.

In `decks/subcmos2014/results/live_interface_fresh/pmos/spacers.str`, region 15
has 61 triangles, area approximately 4.105e-7 um², and minimum dimensionless
quality approximately 1.64e-7. These slivers precede implantation and activation.
In a captured later state, improving flips for three bad oxide triangles cross
only same-material region boundaries; same-region alternatives cannot resolve
the bad geometry. Deleting or straightening a physical boundary is not needed.

Two additional numerical defects aggravated this representation problem:

* `ck_flip` used float slope/intercept arithmetic to decide admissibility.
  Exact orientation of the stored coordinates accepts a captured valid quad
  that the old arithmetic rejected.
* `moving_bound` temporarily moved float coordinates backward and forward to
  assemble flux. Those round trips changed geometry. Flux assembly now uses
  local double vectors without writing coordinates.

Even continuous motion that preserves orientation can invert a binary32-stored
sliver. The first captured example has determinant 7.47877e-19 cm² before
motion, 7.47882e-19 under continuous motion, and -6.59892e-19 after rounding.
`mesh_trial_move` therefore checks the actual rounded candidate before mutation.
This is necessary validity checking, not a cure for artificial constraints.
Trace artifacts are in `decks/subcmos2014/results/sol_motion_trace/`.

## Applying CAD and meshing principles

Open CASCADE's [same-domain unification](https://www.occt3d.com/dev/doc/refman/html/class_shape_upgrade___unify_same_domain.html)
provides the architectural analogy: model-equivalent pieces may be unified,
while protected shapes and modification history remain explicit. Gmsh's
[conformal Boolean fragments](https://www.gmsh.info/doc/texinfo/)
reinforce the need for shared interfaces rather than independently generated
coincident boundaries. These are design references; no CAD library is linked.

The new `material_domain_plan` constructs a quotient of oriented triangle
incidence. It joins edge-adjacent regions only when material and model key agree
and neither is protected. Coordinates, triangles and fields are not edited by
this operation. It rejects invalid/nonmanifold edge incidence. Because the
legacy skeleton requires disk-shaped regions, a union with multiple boundary
loops or a pinched boundary retains its partition. This is a compatibility
bridge, not full boundary-representation support for holes.

The legacy adapter is explicitly selected with `SUPREM4_MATERIAL_DOMAINS=oxide`.
Only oxide is eligible under the current material-global model. A same-material
seam represented by distinct nodal identities protects both regions. Other
materials remain protected. Future region-specific oxide models must supply
different model keys or protection; material name alone is insufficient.

The adapter runs at diffusion command boundaries and after native oxide
creation, outside Boolean operand iteration. It prepares replacement region
records before changing labels, preserves representative non-derived metadata,
and logs the old-to-new region map. Geometry and field arrays are unchanged by
the relabel. The subsequent legacy derived-connectivity rebuild can perform
edge flips: this report does **not** establish conservative field projection
through all such legacy flips or rollback for failures inside that rebuild.

The independent point-removal operation still preserves its input partition
per region. Domain unification is a separate, explicit model-equivalence
decision, not permission for a coarsener to erase arbitrary interfaces.

CGAL's [polygon mesh processing guidance](https://doc.cgal.org/6.1.1/Polygon_mesh_processing/index.html)
distinguishes robust predicates from valid constructed coordinates. Its
[Mesh_2 documentation](https://doc.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_2/Chapter_main.html)
also explains the limitations imposed by small input angles. This supports
fixing domain constraints before expecting refinement alone to cure slivers.

## Controlled comparison and current result

The diagnostic script `suprem4/experiments/material_domain_20261115/coalesce.py`
changes only region labels/records in the archived source/drain structure. It
asserts that material geometry, connectivity and every field token remain
unchanged. The report is
`decks/subcmos2014/results/sol_material_domain_input_report.json`.
This changed-input experiment is causal evidence, not production acceptance.

| Run | Result |
| --- | --- |
| Checked motion, original region partition | Retry limit at 72.2058 s; repeated rounded-motion rejection |
| Diagnostic same-domain input | Reaches 1717.9189 s; no mesh-motion rejections; solver retry limit |
| Integrated domain operation, unchanged input | Reaches 1716.1226 s; no mesh-motion rejections; solver retry limit |

The initial integrated candidate is `suprem4/build/native64__faneiky/src/suprem4`.
It maps 29 input regions to 3 domains. The PMOS run reaches the 1000 C hold but
does not produce `annealed.str`: **full activation FAILS**. The reported failure
is at retry depth 11, trial 129, dt 0.0019293034 s. A bounded retry failure does
not prove the equations are unsolvable; it identifies unfinished solver work.

An independent audit of its `step_failure.str` passes: 12,348 positive triangles,
23,762 candidate overlap pairs checked, 940 boundary edges. This verifies the
stored partition against inversion, overlap and boundary-intersection checks;
it does not prove physics, absence of missing material, or correct time history.
The audit is saved beside the failed checkpoint.

The initial ordinary and strict ASan/UBSan candidates pass the captured tiny PMOS
activation reproducer. Ordinary NMOS and PMOS spacer tests and the LOCOS strip
test also pass. Unit tests cover topology/model protection, hole retention,
allocation-failure atomicity, motion rejection/rollback, and read-only flux
assembly. See `material_domain_reliability_results.json` for artifact paths and
hashes. The validated launcher remains `native64_2r55t88h`.

## Follow-up: the apparent nonlinear failure originated in geometry weights

A diagnostic replay isolates the first nonfinite residual at time
1714.4981057308769 s, TR dt 2.3145671853107026 s. Boron and phosphorus residuals
are infinite at silicon node 6985 because its assembled control-volume area is
negative infinity; the initial concentrations are finite. At the next failing
BDF stage another nodal area is positive infinity. This is not evidence that
Newton damping is the remedy. See
`decks/subcmos2014/results/sol_domain_residual_trace/pmos/suprem4.log`.

The old `chp_fix` derives angle denominators by subtracting squared side lengths.
In valid high-aspect-ratio triangle 6448 from the failed checkpoint, cancellation
rounds an acute-angle denominator to zero. Division then creates infinite
couplings and nodal areas. This is a defect in evaluation of the discretization,
even when the mesh itself has positive orientation.

The candidate `geometry.c` now computes the equivalent denominators with vertex
dot products, uses double intermediates, and avoids reloading rounded side
lengths during area calculation. `tri_geom` computes the same circumcentric
area terms as edge-squared times absolute half-cotangent divided by four;
this removes a second cancellation between squared circumradii near right
angles. Neither change edits geometry, relaxes a tolerance, clips a weight,
nor changes a material boundary. The legacy float database is retained.

`tests/geometry_weights.c` fails with the old code on the captured triangle and
passes with the new code. It checks finite nonnegative weights and triangle
area conservation for the corrected obtuse and acute/right cases, including
cyclic permutations, scales, and 200 additional shapes (219 cases total).
This does not establish global conservation for every legacy obtuse-edge
policy, mesh adaptation, or moving-boundary step; those remain acceptance work.

Separately, the trace confirms ignored comparison-solve failures in the old
time controller, including a failed estimate at 1297.7098446557916 s followed by
an increased timestep recommendation. `time.c` now rejects failed/nonfinite
comparison estimates before publishing the stage. `test_time_stage.py` checks
first TR, BDF and comparison failure, nonfinite comparison output, successful
integration and a fresh retry. The old source fails the last two rejection
cases; the fixed source passes all five scenarios under ASan/UBSan. Correcting
this alone still fails the full replay at 1718.9151 s, so it did not resolve the
geometry-weight root cause.

The combined ordinary candidate is `suprem4/build/native64_m23qw5jt/src/suprem4`;
the sanitizer build is `native64_zhq5f5_p`. Their results are recorded in the
machine-readable manifest. Do not substitute a passing unit test for full-flow
acceptance.

## Follow-up: straight side-boundary subdivisions

Fresh NMOS and PMOS share the pad-oxide prefix. With interior-only removal,
141 points are removed but a redundant silicon point remains on each vertical
side. Oxidation then approaches an almost zero-height silicon triangle. A
diagnostic checkpoint at 1906.44 s contains triangle 1 at these coordinates
(um): `(0, .010042106168839382)`, `(.18658874978427775,
.010042106168839382)`, `(0, .01004210730570776)`. The third vertex subdivides a
straight side; it is not the moving silicon/oxide junction.

The proposal now supports a single open fan only when its two boundary rays
are exactly collinear and the removed point lies between their endpoints.
It substitutes the identical directed segment, checks positive triangles,
attachments and quality, and preserves the region. The adapter requires two
matching non-exposed boundary codes, one material node, and the same conservative
field-transfer and prepare/commit guarantees as interior removal. Exposed,
backside, multi-material and noncollinear boundary points remain protected.

The preceding ordinary candidate is `native64_n0paar46`; sanitizer candidate
`native64_k8ofsow3`. The six-second oxidation replay in
`tests/fixtures/rit_pad_side_boundary/` passes, removes two side points, and
passes independent partition checks. The C commit test checks boundary-code
preservation, species inventory, refusal for differing codes, and unchanged
semantic state under each of 25 injected allocation failures for both interior
and side-boundary operations. Kernel tests independently check exact area and
the directed boundary after combining the collinear subdivisions.

The removal inventory tests use the piecewise-linear triangular field integral
(triangle area/3 per vertex). Diffusion assembles its time terms with the legacy
`earea` control-volume weights. Equivalence of these inventory measures through
adaptation has not been established. A physically conservative remapper must be
validated against the actual transport mass/volume operator as well; passing the
current field-integral tests is not sufficient for a calibrated dose claim.

This is an extension of the constrained operation specified in the replacement
contract. It does not authorize general boundary smoothing or deletion of an
oxide/silicon interface point. Extended fresh-flow results belong in the
manifest, separately from the captured activation milestone.

The old interior-only fresh-flow controls reached their 900-second wall-clock
limits during pad oxidation (1/16 checkpoints). The extended candidate passes
pad oxidation and reaches later stages. Fresh PMOS reaches LOCOS, then fails
during the initial Kooi-oxide step; it is not a successful fresh PMOS process.
The bounded replay from its `locos.str` reproduces the failure in
`sol_side_boundary_threshold_control`.

In that failure, triangle 7395 has positive exact orientation, but `back_pedal`
subtracts float coordinates before assigning its edge vectors to doubles. The
rounded vectors become parallel. Two stored vertices lie at x=5.5169353e-23
and 5.6843424e-21 cm with identical y; the third is at x=4.9980287e-8 cm.
The panic partition independently passes 7,946 triangles. Evidence is saved in
`sol_side_boundary_fresh/pmos/back_pedal_diagnosis.json`.

The repaired calculation in `dbase/grid_upd.c::back_pedal` converts each stored
coordinate to double *before* subtracting, including the displacement and 1D
paths. On triangle 7395 the double-vector determinant is 1.57378e-28 cm²,
where the former float-vector determinant was zero. A new frozen `pmos_kooi`
fixture in `tests/fixtures/rit_pmos_kooi/` fails with the preceding binary
`native64_n0paar46` and passes with `native64_6jztdvdl` and strict ASan/UBSan
`native64_32q9_n2u`. The ordinary saved result independently audits 7,900
triangles and 380 boundary edges. A longer ordinary continuation completes
Kooi and reaches LDD (11/16 accepted checkpoints), then fails at spacer etch
with `ambiguous clipped boundary connectivity`. This proves the arithmetic
crash was real; it does not establish a stable fresh process.

`sol_back_pedal_spacer_diag` records two retained clipping arcs whose end/start
positions map ambiguously at the left boundary. The initial failed run's panic
mesh passes an independent partition audit of 11,777 triangles; the failure is
in constructing the next boundary, not in that stored predecessor. The
`pmos_spacer_clip` fixture freezes the LDD input and reproduces the failure.
The source/clip intersection trace shows competing crossings at x approximately
5.45e-27 and 5.68e-21 cm near y=-1.0807e-5 cm. These near-zero coordinates
already occur after LOCOS; repeated left-side x drift is visible in its saved
structure. The physical side-boundary constraint must be represented through
motion and clipping, with canonical shared intersection identity and a
well-defined arrangement of retained arcs. Do not resolve this by selecting an
arbitrary successor, snapping by epsilon, or skipping the panic.

This exposes two distinct requirements: geometric calculations must consistently
promote operands before subtraction, and the upstream construction must preserve
authoritative boundary constraints and shared point identity. Establish why the
near-coincident boundary vertices are created before deciding the physical
operation. Do not merely remove the assertion, snap a tiny coordinate to zero,
or enlarge an epsilon. A CAD-style constraint representation must determine
which points truly belong to the same line/interface and when they are identical.
The present mesh alone does not encode all of that construction intent.

The final fresh NMOS run writes an annealed structure with a passing independent
partition audit (13,853 triangles, 26,535 candidate overlap pairs). It then
panics during ILD deposition. It also fails the existing body-polarity checks
at source/drain and annealed stages: the latter samples net donor concentration
8.9629e18 cm^-3 where an NMOS acceptor body is required. Consequently only
12/16 checkpoints are accepted by the runner, despite 14 stage files being
written. This is a failed full flow, not a functioning NMOS device.

## Next geometric failure: ILD offset boundary

The successfully annealed checkpoint reloads cleanly. The next command,
`deposit oxide thick=0.5 divisions=10 space=0.06`, fails before producing ILD.
Its generated skeleton contains 1,305 vertices. Independent exact orientation
tests identify a proper crossing between segments 746 and 710; an atan2 turning
sum is -4*pi. This is a genuinely self-intersecting constructed boundary, not
an overly strict angle tolerance or a misleading panic export.

The full before-state is captured in `tests/fixtures/rit_pmos_ild_boundary/`.
`run_rit_locos_strip.py --case pmos_ild_boundary` reproduces the failure. The
crossing endpoints and audit are in
`decks/subcmos2014/results/sol_stable_weights_ild_boundary/`. The attempted
post-activation finish correctly reports FAIL (14/16 checkpoints).

The next deposition design must treat the offset as a planar arrangement:
construct and node all segment intersections, classify resulting faces against
the deposited-material domain, preserve physical interfaces, and emit simple
boundary loops (or explicit holes for a capable backend). Prepare that complete
domain before modifying the live database. Verify the union and area against an
independent reference. Do not delete an inconvenient segment, enlarge the angle
tolerance, or let the triangulator guess which crossing loop is intended.

## Modern TCAD and FEA guidance for the next failure

[DEVSIM](https://devsim.net/CommandReference.html) exposes positive and
log-damped variable updates. This is a useful reference for treating a
concentration's domain explicitly. It does not justify clipping a negative
concentration and accepting an unconverged solve.

[PETSc SNES](https://petsc.org/release/manual/snes/) provides Newton line-search
and trust-region methods and distinguishes invalid function-domain evaluations.
[PETSc stage checking](https://petsc.org/release/manualpages/TS/TSAdaptCheckStage/)
separates successful stages from nonlinear failure or infeasible solutions.
Apply that separation to primary integration stages and the error estimator.

[MOOSE's adaptive timestepper](https://mooseframework.inl.gov/moose/source/timesteppers/IterationAdaptiveDT.html)
uses solve difficulty to adjust subsequent steps and cuts back failed steps.
For SUPREM, reducing a retry budget or increasing it is not a convergence
algorithm. A successful substep must inform the controller rather than force
completion of an increasingly deep recursive tree at an obsolete requested dt.

Required implementation sequence (first two diagnoses now have evidence):

1. Capture the first failed nonlinear stage with species, material, node,
   iterate, Newton correction, residual, dt and failure reason. Separate
   invalid geometry, inadmissible iterate, linear failure, stagnation and failed
   error-estimator solves. Preserve the unchanged input and model settings.
2. The failed-comparison acceptance regression and fix are implemented above.
   Retain this test: a failed estimator must not supply a timestep recommendation.
3. Specify one coupled trial state: coordinates/topology, accepted nodal
   solutions, active fields, old/mid/new buffers, velocity history, matrix/cache
   invalidation and simulation time. Prove reject/retry equivalence to a fresh
   smaller step. `mesh_trial` alone deliberately does not own all these fields.
4. If traces identify inadmissible Newton corrections, implement and test a
   residual-based damped update with explicit concentration bounds. Evaluate
   trial residuals without publishing fields; accept only a valid improving
   iterate. Follow the principle of
   [PETSc's Jacobian comparison](https://petsc.org/release/manualpages/SNES/SNESTestJacobian/)
   before compensating for a wrong Jacobian with damping. Signed fields need a
   separate policy. The captured infinity failure instead requires the geometry
   weight correction; do not replace it with a damped Newton workaround.
5. Replace recursive forced-interval subdivision with accepted-step control
   only after the coupled-state contract is tested. Recompute moving-interface
   quantities consistently on cutback, respect process-stage endpoints, and
   retain explicit work/minimum-step failure reporting.
6. Check constant-state preservation under prescribed mesh motion, closed-domain
   species balance, interface flux balance, and conservative transfer through
   every remeshing operation. Then run unchanged NMOS and PMOS flows, followed
   by mesh-density and mask-shift convergence studies.

No full-flow or launcher promotion is justified until these acceptance gates
pass. Geometry reliability and physical accuracy remain separate requirements.

## Remaining priorities after the current verification

1. The Kooi arithmetic regression is closed. Trace the left-side boundary
   drift from LOCOS and the new spacer clipping ambiguity; implement explicit
   boundary constraints and canonical intersection identity before reconnecting
   retained arcs. Require a passing `pmos_spacer_clip` control and fresh PMOS
   continuation, not only a Kooi replay.
2. Replace the ILD offset boundary construction with an arrangement/domain
   operation, using the captured crossing as a rejection and coverage oracle.
3. Trace the NMOS phosphorus inventory through deposition, transport and mesh
   changes; align conservative transfer with the actual transport volume
   operator. Do not tune implant doses to conceal the polarity error.
4. Only after unchanged fresh flows pass, run density/mask perturbations and
   physics convergence checks before promoting a launcher.

## Reproduce the latest activation milestone

From the repository root:

```sh
SUPREM4_MATERIAL_DOMAINS=oxide SUPREM4_MESH_PATCH_REMOVAL=1 \
python3 decks/subcmos2014/run.py \
  --binary suprem4/build/native64_n0paar46/src/suprem4 \
  --devices pmos --through annealed \
  --resume decks/subcmos2014/results/live_interface_fresh --after source_drain \
  --output /tmp/rit_domain_activation_new --timeout 900 --jobs 1
```

The output directory must be new. Expected result: PASS with `annealed.str`.
This resumes the archived source/drain checkpoint; it is not a fresh full flow.
The earlier `native64__faneiky` control stops with the retry failure described
above. Full fresh runs use the same candidate/options without `--resume` or
`--after`, and `--through metal1`; their currently failing results are recorded
in `sol_side_boundary_fresh` and in the machine-readable manifest.

To reproduce the newly closed Kooi arithmetic failure and the next open spacer
boundary failure, use fresh output directories:

```sh
SUPREM4_MATERIAL_DOMAINS=oxide SUPREM4_MESH_PATCH_REMOVAL=1 \
python3 suprem4/tests/run_rit_locos_strip.py --case pmos_kooi \
  --binary suprem4/build/native64_6jztdvdl/src/suprem4 \
  --output /tmp/rit_kooi_verified

SUPREM4_MATERIAL_DOMAINS=oxide SUPREM4_MESH_PATCH_REMOVAL=1 \
python3 suprem4/tests/run_rit_locos_strip.py --case pmos_spacer_clip \
  --binary suprem4/build/native64_6jztdvdl/src/suprem4 \
  --output /tmp/rit_spacer_open
```

Expected results: Kooi PASS with audited `after.str`; spacer FAIL with
`ambiguous clipped boundary connectivity`. The first result does not imply the
second is safe to bypass.
