# Open-source TCAD comparison and implementation direction

Reviewed 2026-11-15 against primary project documentation and source. These
projects inform the design; they do not independently validate this process.

2026-11-16 follow-up: [TCAD repair-method research](tcad_repair_methods_research_20261116.md)
adds thickness-aware pre-etch refinement, sub-grid oxidation, conservative
transfer and transient-history evidence. The earlier dyadic-grid backend
results below establish bounded topology successes; they do not establish
deposit/equal-etch construction identity, which subsequently failed.

## What transfers to this code

* [DEVSIM numerical methods](https://devsim.net/solver.html) provide BDF1,
  TRBDF and BDF2 integration and configurable precision in model evaluation,
  equation assembly and linear solves. Its [model documentation](https://devsim.net/models.html)
  separates edge flux integration, nodal volume integration and interface
  surface area. Apply that separation here: audit the measures before tuning
  Newton tolerances, and use a damped implicit method only with an error estimate.
* [Genius solver configuration source](https://github.com/cogenda/Genius-TCAD-Open/blob/master/src/solver/solver_specify.cc)
  explicitly distinguishes lower-order BDF startup, LTE rejection and timestep
  limits. It also exposes Voronoi truncation. This supports treating positivity,
  volume construction and temporal accuracy as separate acceptance conditions.
  Device-carrier formulas cannot simply replace process dopant equations.
* [ViennaPS](https://github.com/ViennaTools/ViennaPS) represents moving process
  surfaces using level sets and separates surface evolution from volume fields.
  Its [oxidation solver](https://github.com/ViennaTools/ViennaLS/blob/master/docs/OxidationSolver.md)
  couples oxidant transport and mechanical deformation on a Cartesian grid.
  The applicable principle is explicit geometry/field ownership and controlled
  motion. A level-set migration is an architectural alternative, not a small
  replacement for SUPREM's existing moving finite-volume equations.
* [Clipper2 offset documentation](https://www.angusj.com/clipper2/Docs/Units/Clipper/Functions/InflatePaths.htm)
  provides polygon dilation and joins. A deposition replacement should construct
  the deposited domain by dilation and Boolean difference, validate every loop,
  then triangulate. It must handle holes, side constraints and fields. Merely
  removing a crossing from the legacy offset path is insufficient.

## Findings from applying these principles

1. Reflecting normals: cancellation in normal subtraction creates false lateral
   geometry. Direct tangent projection repairs construction before `SplitNode`.
   See `reflecting_boundary_repair.md`, with failing controls and sanitizer passes.
2. Interface quadrature: `dbase/couple.c::sum_coeff` formerly included every
   interface face of a node's incident triangles, including the opposite face
   which does not contain that node. This attributes transport to unrelated
   nodes and overcounts interface measure. The same incidence omission exists
   in `diffuse/moving.c::moving_bound`. Both now filter for face incidence.
   `tests/interface_measure.c` checks a hand-computed interface and exposed-face
   measure; the old source fails. `tests/moving_flux.c` adds an opposite-face
   exclusion case to its analytic swept-flux and immutable-coordinate checks.
3. Legacy circumcentric mass: an obtuse Delaunay patch can have positive assembled
   edge couplings while its unsigned element volume weights exceed actual area.
   An experimental `SUPREM4_MASS_LUMPING=p1` option uses triangle area/3 at each
   vertex, including after the legacy constrained obtuse treatment. It matches
   linear-field integration and the coarsener's P1 field-transfer measure.
   Twenty-eight rotated/scaled analytic patches pass; the legacy source fails.
   Fresh P1 runs still fail and this option is not accepted as a CMOS fix.
4. Temporal positivity: the fresh PMOS activation failure after the side repair
   includes severe negative TR updates in thin oxide, and later BDF rejection.
   A restart can resolve some of these with much smaller timesteps. Increasing
   retry limits alone would not establish accuracy, dose balance or robustness.
   Rejection-stage diagnostics now distinguish motion, primary and comparison
   failures. The next temporal implementation must carry an LTE estimate and
   use the same conservative spatial measures as field transfer.

## Acceptance order

First verify face incidence and interface measure, then repeat fresh NMOS/PMOS
with the unchanged deck and recorded runtime options. Track total chemical dose
by material before/after each diffusion and remesh operation. Any purported
body-polarity repair must explain where phosphorus moved, not alter the implant.

For moving diffusion, test the discrete geometric conservation law with a
constant field under prescribed motion. Compare the actual nodal mass changes,
swept boundary flux and remapping inventory. A successful nonlinear solve does
not prove this identity. Implement a controlled lower-order implicit method
only after that audit; verify temporal convergence and positivity without
concentration clipping.

For ILD, replace the offset/loop-deletion domain builder and preserve all source
surface constraints. Validate the proposed domain before any GAS/material
mutation. The captured crossing fixture, planar/step/trench reference areas,
round-corner cases, side translations and holes are required regressions.

Finally require complete fresh process checks, ordinary/sanitized replays,
mesh scaling, mask shifts and dose/timestep convergence. Keep the launcher on
the validated binary until these gates pass.

## Results of the domain construction implementation

The opt-in Clipper2 backend is implemented, not just proposed. Build it with
`--domain-offset` and select `SUPREM4_DEPOSITION_GEOMETRY=domain`. It constructs
an offset occupied domain, validates its Boolean difference, and passes the
outer chain to the existing constrained boundary. Source material vertices
and fields are retained. Unsupported holes/multiple components fail explicitly.

A fine-grid prototype failed upon conversion to binary32. The corrected backend
constructs its output directly on a dyadic grid exactly representable in the
current point database, with geometric resolution derived from coordinate
magnitude. This avoids a second, topology-damaging conversion. Arc chord error
and grid resolution are reported independently.

The archived PMOS ILD fixture passes independent partition checks on 49,207
triangles, ordinary and ASan/UBSan (leak detection disabled). All four existing
region areas are unchanged. The added oxide area is 5.374806130246482 um2.
The captured fresh spacer case also passes, ordinary and sanitized.

The archived PMOS `annealed` continuation now completes contacts and metal1:
16/16 runner checks pass. Separate final partition audits pass 21,159 contact
triangles and 37,062 metal1 triangles. Evidence is under
`decks/subcmos2014/results/astra_domain_pmos_finish/pmos/`.

Full fresh CMOS remains open. Fresh PMOS fails during LOCOS, with a genuine
nitride boundary crossing (segments 232 and 234, 601 boundary vertices,
turning 4*pi). The stable `atan2` angle formula reproduces the same crossing;
it is not an angle-sum false alarm. The saved mesh independently fails overlap
checks on triangles 11739 and 11740. Evidence is
`decks/subcmos2014/results/astra_domain_stable_angles/pmos/locos_boundary.json`.
This must be traced to its first moving-boundary/remeshing mutation, with a
complete constrained replacement and conservative field transfer. A different
triangulator alone cannot repair an already self-intersecting material domain.

The old `dbase/grid_loop.c::DetectLoop` is not an acceptable recovery algorithm:
it assigns both `pcc` and `pc` from `ncc`, making its first edge vector zero;
its angle normalization also adds a dimensional constant 1e-10 to products of
edge lengths. Correcting those lines would activate heuristic material-swallowing
branches, some explicitly unfinished. Do not enable that as a claimed cure.
Use explicit contact/topology events and validated material-domain operations,
as motivated by ViennaPS's geometry/field separation.

The recursive velocity-history repair is separately verified: temporary
midpoint arrays were overwritten when child solvers reloaded nd[].sol. History
now lives in that source and follows stable node objects across reordering.
The old source fails a nested-halving/reordering regression; the correction
passes. It does not cure the captured nitride crossing.
