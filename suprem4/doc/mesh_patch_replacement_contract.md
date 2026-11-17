# Replace the mesh operation, not its failure symptoms

Decision: 2026-11-15. This is the required design and acceptance contract for
further coarsening work. It is not a claim that the replacement exists.

Implementation update: [material domains and solver reliability](material_domain_reliability.md)
records the candidate interior-removal work and the upstream artificial oxide
constraint diagnosis. The per-region preservation contract below still applies
to removal. A separately selected, model-aware domain-unification operation may
identify equivalent regions before meshing; it must preserve protected physical
interfaces and fields. It is not an implicit exception inside point removal.

## Why the approach changes

The old path is `grid_sub -> rem_pt -> kill_node -> triangulation -> geom`.
`rem_pt` mutates material nodes sequentially. `kill_node` orders neighbors by
floating-point angle, edits global mesh state and uses the legacy triangulator.
The current implementation also contains duplicate-triangle cleanup inherited
from earlier work. These are coupled topology, geometry and bookkeeping changes;
fixing their downstream symptoms one at a time has not established correctness.

The withdrawn rewrite replaced angular ordering with directed topology and then
exact ear clipping. On the captured activation replay its panic mesh passed the
independent partition audit, but the engine still reported an exposed oxide
point without a GAS node. That is a failed operation, not partial acceptance.
The precise GAS-node mutation responsible is not yet established. Do not assert
that the triangulator alone explains the entire failure.

The experiment is preserved in `../experiments/rit_cavity_20261115/`, with hashes
and baseline-relative patches. Four active source files were restored from
`build/native64_b43fm60y/src/`, preserving prior workspace changes. The old
coarsener remains defective; restoring it withdraws an unvalidated replacement.
Existing binaries are historical artifacts, not rebuilt by this restoration.

## Evidence to explain before changing the engine

The [first invalid removal diagnosis](rit_first_invalid_removal.md) now traces
`grid_sub -> rem_edg(12550) -> rem_pt(3796) -> kill_node(5102)` and extracts
the crossing open fan. It also records GAS deletion before the attempted
oxide replacement. Use that reduced fixture as the initial control.

Use `tests/fixtures/rit_pmos_activation_seed/` and
`tests/run_rit_locos_strip.py --case pmos_activation_seed`. Its tiny dry-oxygen
step reproduces failure without the full furnace schedule. Keep its recipe,
input coordinates and model settings unchanged.

Archived stage snapshots show valid partitions after native oxide creation,
loop processing and grid addition; corruption appears during node removal.
The original failure is a self-crossing proposed boundary cavity. The earlier
claim of overlapping live oxide triangles came from serializing retired
triangles in the panic export and has been withdrawn. The withdrawn
candidate `native64_44j683ve` instead exits with an exposed-point invariant
failure. Evidence is in
`decks/subcmos2014/results/pmos_exact_ear_replay/` from the repository root.

First capture one removal's complete before-state and every mutation of its
points, material nodes, triangles, edges and exposure. Identify the first
invariant violation, using physical identity and coordinates as well as IDs:
compaction changes IDs. Extract that patch and its external neighbors as a
standalone failure fixture. Geometry audits alone cannot verify GAS-node or
field ownership. An audit of the panic file cannot prove earlier states valid.

## One operation and one authoritative patch

A removal proposal must cover the complete physical point, all of its material
nodes (including GAS bookkeeping), incident live leaf triangles, and external
neighbors. Region identity must be preserved separately from material identity.
Do not apply independent per-material removals to the global mesh.

Construct a read-only local representation from live triangles and physical
point identities. Reconstruct directed incidence and external attachment faces;
validate it against stored adjacency before proceeding. Contradictory input is
an explicit invalid-input result, not an invitation to guess connectivity.

Produce a replacement in private storage. No global allocation, deletion,
renumbering, flag change, field update or cache invalidation occurs while
building or checking the proposal. The candidate includes triangles, adjacency,
material ownership, exposure/GAS changes and transferred fields, not just a
polygon triangulation. Exact orientation refers to the stored binary32 input
coordinates; it does not make arbitrary floating-point constructions exact.

For the first implementation, preserve the patch domain and all material
interfaces exactly. Internal point removal can preserve the old domain.
Boundary removal is admissible only when it preserves the same segments and
region partition (for example an exactly collinear subdivision point).
Pin corners, noncollinear interface points and unsupported junctions with an
explicit unsupported reason. Approximating a curved moving boundary requires
a later, separate geometric-error and conservative-transfer contract. Do not
silently straighten it to obtain a simpler cavity.

## Required invariants

* Every candidate triangle has three distinct vertices and positive exact
  orientation. Directed incidence has one owner on exterior faces and two
  oppositely directed owners on interior faces. No nonmanifold edge, duplicate
  triangle, crossing or interior overlap is allowed.
* The candidate union equals the old patch domain, per region, with the same
  directed exterior attachment faces and constrained interfaces. Check coverage
  and area as well as intersections; absence of crossings does not exclude gaps.
* Triangle-to-node, node-to-triangle, edge and neighbor references agree. No
  live reference targets a deleted object. Physical points retain precisely the
  required material-node ownership; distinct regions are not merged by material.
* Exposure follows the resulting boundary conditions. An exposed point has the
  required GAS node and matching flags; covered material interfaces do not gain
  artificial exposure. Audit all affected boundary points, not only the removed
  point. Define back boundaries separately from exposed surfaces.
* Species inventories use the existing nodal finite-element representation.
  Removal transfers chemical inventory conservatively within each material
  unless an explicit physical flux applies. Positivity and bounds must be
  specified per field: chemical concentrations differ from signed net doping,
  velocities, defects and derived active concentrations. Copying/interpolating
  every field identically is not an acceptable transfer rule. Report transfer
  residuals separately from physical diffusion or surface loss.
* An unsuccessful proposal leaves the semantic mesh state unchanged, including
  list membership, ownership, flags, field values and dirty-cache state. Optional
  coarsening refusal is counted and reported; required geometry work may fail
  explicitly. Refusal does not count as a successful repair of this replay.

## Commit and integration

Allocate and populate the complete replacement before commit; prepare an
explicit attachment/retirement map. Validate it against the unchanged patch
version. Publish all material changes together, then retire the old entities.
The commit must have no fallible geometry or allocation step. If that cannot be
made true in the legacy database, use a complete owned-state rollback that also
covers entity creation/deletion and all field buffers. Prove failure restoration
with injected failures before relying on it.

`mesh_trial_save/restore` is NOT that rollback: its header forbids adding,
removing or renumbering entities and its implementation rejects changed counts.
It was built for moving-boundary diffusion trials and edge flips. Do not extend
its use to point removal without redesigning and testing the ownership model.

Initially rebuild derived connectivity from candidate triangles and compare it
to the committed database. Optimize local updates only after correctness is
established. Never delete duplicate triangles, invent boundary labels or create
missing GAS nodes after a failed audit to force an operation through.

## Implementation order and acceptance

1. Extract the first failing removal and demonstrate the exact pre/post
   invariant violation. Add an independent topology/ownership checker alongside
   the existing geometric partition checker. Keep this reproducer small enough
   to inspect without a complete process run.
2. Build a standalone proposal kernel for internal single-region removals.
   Test convex and concave stars, near collinearity, high valence, invalid input,
   permutations and allocation failure. Compare union coverage independently.
3. Extend the same proposal to straight boundary subdivisions and complete
   multi-material patches, including GAS ownership and constrained interfaces.
   Include same-material/different-region boundaries and point-only contacts.
   Refused operations must preserve semantic state exactly.
4. Implement field-specific conservative transfer and atomic database commit.
   Check inventory residuals, derived incidence and rejection restoration with
   failure injection. Then integrate behind an explicit candidate option; keep
   the validated launcher unchanged.
5. Require the captured activation replay to complete with engine self-tests,
   independent partition checks, ownership checks and strict ASan/UBSan. Require
   demonstrable valid coarsening, not success obtained by disabling removal.
6. Run fresh NMOS and PMOS through the unchanged process, then mesh-density,
   mask-shift and controlled perturbation matrices. Separate geometry reliability
   from NMOS phosphorus transport, dopant polarity and convergence. Timeouts,
   missing outputs and interrupted runs are failures/incomplete results.

Do not replace these gates with a picture, a single passing deck or sanitizer
silence. Full-flow acceptance and launcher promotion remain open.
