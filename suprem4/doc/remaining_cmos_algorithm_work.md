# Remaining CMOS algorithm work after the domain-backend repair

Current 2026-11-16 evidence: [connectivity transfer repair](conservative_connectivity_repair.md)
and [deposit/etch construction identity](deposition_etch_construction_identity.md).
The priorities below are historical. The opt-in
[construction repair](process_surface_construction_repair.md) passes bounded
identity tests; full fresh CMOS remains under validation.

2026-11-15 historical review. Read `open_source_tcad_repair_research.md` and the status
in `rit_spacer_boundary_repair.md`. Do not infer fresh acceptance from the
archived annealed-to-metal continuation, which passes.

## First priority: locate the first nitride embedding failure

Captured evidence:
`decks/subcmos2014/results/astra_domain_stable_angles/pmos/locos_boundary.json`.
The represented nitride boundary has 601 vertices and a proper crossing
between segments 232 and 234:

* 232: (0.0005582667072303593, 1.5819057352928212e-7) to
  (0.0005582962185144424, 1.581956894369796e-7) cm.
* 234: (0.0005582961020991206, 1.581956894369796e-7) to
  (0.0005583299207501113, 1.5819567522612488e-7) cm.

The turning number is 2 (turning angle 4*pi). The saved mesh also has overlapping
triangle interiors, 11739 and 11740. Stable atan2 angle evaluation reproduces
this failure. The input LOCOS-mask checkpoint is available in the same run.
The failure occurs during wet oxidation after the initial dry segment, around
781.526 s total before the next attempted macrostep.

Save the last accepted state before that macrostep and instrument these
operation boundaries: oxide velocity calculation, DetectLoop, grid_add/new_layer,
each tentative mesh move/flip, and each grid_sub removal. Use exact stored-coordinate
segment/triangle predicates in a small window around x=0.0005583 cm. Determine
which *completed mutation* first changes a valid embedding into overlap.
Do not begin with downstream polygon cleanup.

If motion introduces it, the accepted-step condition must include material
embedding/contact, not only individual positive triangle orientations. The
physical topology event must trigger a complete local constrained replacement
with positive volumes and conservative chemical-field transfer. Refusing an
invalid move is a safety gate, not completion of that event algorithm.

If node insertion/removal introduces it, construct its entire material patch in
private storage and validate the old/new boundary chain, all constraints,
nonlocal overlap and fields before commit. Follow `mesh_patch_replacement_contract.md`.

Do not simply correct `DetectLoop`'s `pc=pt_nd(ncc)` typo and activate it.
The zero first edge vector currently defeats its angle test; its normalization
also adds an absolute dimensional 1e-10. Several branches swallow/reassign
material and are explicitly incomplete. They require an algorithmic replacement.

## Second priority: NMOS spacer side topology

Frozen input: `tests/fixtures/rit_nmos_domain_spacer/`.
Use `tests/run_rit_locos_strip.py --case nmos_domain_spacer` with the domain
backend and the runtime options in fixture.json. The input partition passes
independent checks (9,484 triangles). The downstream failure is a reversing
side-boundary quad with stored points:

```
(0.0010000000474974513, -2.2899841496837325e-5)
(0.00099994894117116928, -2.5622528482927009e-5)
(0.0010000000474974513, -2.5622462999308482e-5)
(0.0010000000474974513, -2.5622492103138939e-5)
```

Its last side overlaps the prior collinear side. Evidence:
`decks/subcmos2014/results/astra_ordered_projection_nmos/nmos/suprem4.log`.
Exact quad admissibility correctly refuses this input. Ordering projected
rectangle refinements repairs a separate proven defect but reproduces the
same captured quad; it is not the cause of this case established so far.

Capture the complete region immediately after clipping, before each split and
before any projected-point insertion. Bracket whether the reversing chain
already exists at Boolean output or is created while triangulating. The
current fuzzy ON classification is a hypothesis to test, not a proven cause
of this second case. A previous fully exact ON experiment broke rounded shared
intersections; do not revive it without consistent construction provenance.

For a replacement Boolean pipeline, construct shared intersections once for
all affected materials, preserve physical point identity and constraint ownership,
regularize exactly lower-dimensional results, then triangulate validated domains.
Original boundaries and chemical inventory must survive the transfer. The
successful offset backend demonstrates domain construction for deposition; it
does not yet replace the shared multi-material clipping/transfer operation.

## Modern triangulator evaluation

[CDT](https://github.com/artem-ogre/CDT) is a candidate constrained Delaunay
backend with robust predicates, explicit constraint intersections and hole
removal. Use a mode which rejects intersecting constraints. It can replace
angle-based splitting of *valid* domains, but cannot make either captured
self-overlapping input physically correct. Require every boundary vertex and
constraint to survive, an independent partition audit, and sizing/convergence
checks. Do not reduce node counts silently and claim robustness.

## Physics and integration gates

NMOS body polarity has failed earlier candidates after poly activation. It
must be rerun on a successful fresh flow. Preserve the PDF implant doses and
schedule. Compare chemical inventory by material and flux before/after every
remesh, and bracket the first unexplained transfer.

The optional P1 mass model fixes analytic volume overcounting on obtuse patches,
and the optional BE controller passes analytic temporal/conservation tests.
Both remain experimental because their process replays fail. A moving-domain
acceptance test must include geometric conservation and restart/timestep
convergence; a converged Newton iteration alone is insufficient.

## Reproduction and status

Build ordinary candidates with `tools/build_native64.py --ordering amd --domain-offset`;
add `--sanitize` for ASan/UBSan. The optional deposition backend is selected by
`SUPREM4_DEPOSITION_GEOMETRY=domain`, alongside `SUPREM4_MATERIAL_DOMAINS=oxide`
and `SUPREM4_MESH_PATCH_REMOVAL=1`. Runtime options and build/source/data hashes
are recorded in result files. Full-engine sanitizer runs currently disable
leak detection; this is not a passed leak gate.

The launcher remains on `native64_2r55t88h`. Full fresh CMOS, mesh/mask
perturbations and physical convergence must pass before promotion.
# Latest construction and transfer evidence (2026-11-16)

Read [conservative_connectivity_repair.md](conservative_connectivity_repair.md)
for the reproduced automatic-flip inventory loss and bounded repair. Read
[deposition_etch_construction_identity.md](deposition_etch_construction_identity.md)
for the newly isolated planar deposit/equal-etch identity failure. Neither a
longer retry tree nor deletion of tiny oxide regions fixes that construction.
Full CMOS remains unvalidated; the launcher remains unchanged.

