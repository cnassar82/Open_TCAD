# Sol: implement a correct mesh-removal operation

Prepared 2026-11-15. Read this before modifying the coarsener. This guide turns
[the replacement contract](mesh_patch_replacement_contract.md) into bounded
implementation tasks. The contract defines acceptance; this guide defines the
order of work. Full RIT CMOS acceptance is still open.

## First assignment completed, 2026-11-15

The [causal trace and reduced fixture](rit_first_invalid_removal.md) identify
edge 12550, point 3796, its GAS/silicon deletion order, and the crossing
proposal. Start implementation from this evidence. A standalone exact-rational Stage A
reference is in `tests/mesh_patch_oracle.py`, with five focused tests in
`tests/test_mesh_patch_oracle.py`. Read the
[reference review](mesh_proposal_reference_review.md) before porting it: the
actual RIT case needs a constrained edge operation, quality acceptance,
ownership, conservative field transfer and atomic commit. The opt-in boundary refusal
trial failed later on neighbor reciprocity and was withdrawn. Do not repeat
that narrow gate as a claimed repair.

## Original first assignment

Produce a small, deterministic reproducer of the **first invalid removal** in
the captured PMOS activation step, together with a before/after explanation.
Do not begin with another triangulator rewrite or another full process run.

The deliverable must identify:

1. Which operation and branch ran, including the caller and removal sequence.
2. Whether the input to that operation was already invalid.
3. The first mutation that made the result invalid, or the first completed
   operation with invalid output if intermediate mutations are intentionally
   inconsistent.
4. Which surviving points, material nodes, regions and boundary faces changed.
5. A failing automated test that exercises this operation without rerunning the
   whole CMOS process.

A trace that ends at `bd_connect` or `geom` is insufficient: these routines may
be detecting damage created earlier. Intermediate checks must distinguish an
intentionally incomplete update from an invalid committed state.

## Establish facts before choosing the repair

| Observation | What it establishes | What it does not establish |
| --- | --- | --- |
| The original captured activation produces overlapping oxide triangles during removal | The existing removal path can corrupt the partition | Angle sorting is the only cause |
| The withdrawn exact-ear candidate has an audited panic mesh but fails an exposed-point/GAS-node self-test | Geometry checks alone do not cover mesh ownership | The missing GAS node should simply be recreated |
| `rem_ele` can delete unused nodes; `rem_1nd` recursively deletes a remaining GAS node | Deletion has cascading ownership effects | This cascade caused the recorded missing GAS node |
| `rem_edg` can move nodes between physical points | Tracking numeric point IDs alone can miss ownership changes | Every failure goes through `kill_node` |
| NMOS body phosphorus remains wrong | Physical/numerical validation remains necessary | Fixing mesh geometry will necessarily fix doping |

Important correction to an overly simple diagnosis: the restored `kill_node`
currently calls `grid()` **before** deleting the old incident triangles. Do not
claim it always deletes old triangles before creating replacements. Trace the
actual ordering, including cleanup deletions and edge-collapse branches.

## Read the actual call graph

Relevant files are relative to `suprem4/src/`:

* `dbase/grid_upd.c::grid_sub`: first calls `pre_oxide_remesh` conditionally,
  then runs an edge-removal pass selected by `CknMcLen`, then a point-removal
  pass selected by `CknMcObt`. Both passes use the same self-test message.
* `dbase/remove.c::rem_edg`: chooses a point from endpoint velocities; either
  calls `rem_pt`, or transfers material nodes to the other physical point and
  calls `rem_edg_end`. Trace both branches.
* `dbase/remove.c::rem_edg_end`: substitutes vertices, deletes triangles and
  removes a node; it is a separate topology-changing path.
* `dbase/remove.c::rem_pt`: loops through a saved list of material-node IDs and
  mutates them sequentially. Earlier deletions can affect later ownership.
* `dbase/kill_node.c::kill_node`: angle-sorted skeleton, region selection,
  global triangulation, duplicate cleanup, then retirement of old triangles.
* `dbase/remove.c::rem_ele`, `rem_1nd`: cascading deletion and boundary effects.
* `dbase/check.c`: point/node reciprocity, boundary/GAS and triangle checks.
* `dbase/mesh_trial.c` and `include/mesh_trial.h`: rollback is restricted to
  unchanged entity counts; it cannot roll back general removal/replacement.

Instrument operation boundaries around all these paths. Never conclude which
pass failed from the text “after removing silicon nodes” alone.

## Reproduce the baseline

Run from `/home/chris/Open_TCAD`. This candidate was built after withdrawing the
experimental coarsener; it is expected to retain the known activation defect.
Use a new output directory for every run; the runner rejects populated ones.

```sh
python3 suprem4/tests/run_rit_locos_strip.py \
  --case pmos_activation_seed \
  --binary suprem4/build/native64_6r5i3c40/src/suprem4 \
  --output decks/subcmos2014/results/sol_activation_baseline_01
```

If that historical candidate is absent, build from current source, record the
printed directory, and substitute its `src/suprem4` path:

```sh
python3 suprem4/tools/build_native64.py --ordering amd --jobs 4
```

Fixture: `suprem4/tests/fixtures/rit_pmos_activation_seed/`. The decompressed
input hash is
`054e0c83f97c055d03712217b8f1a0bc9612042738873bf38d0b7fc279b5aaa0`.
The command is `diffuse time=0.000001 temp=825 dryo2`, with the settings recorded
in `fixture.json`. Preserve them. The tiny step isolates a defect; it is not
acceptance of the full anneal.

Historical failed experiment:
`decks/subcmos2014/results/pmos_exact_ear_replay/`, candidate
`native64_44j683ve`. Experimental source and baseline-relative patches are in
`suprem4/experiments/rit_cavity_20261115/`. Do not copy the experiment back into
the active source as an accepted fix.

## Capture without changing behavior

Add opt-in diagnostics, default off. Assign a monotonically increasing operation
number and record operation type, caller phase, selected edge/point, endpoints,
coordinates, region/material IDs, flags and node ownership. Coordinates must
retain the exact stored float values (hexadecimal or equivalent bit patterns).
Record ID remapping at compaction and node transfers between physical points.

At entry and exit, snapshot the affected triangles and their external neighbors,
all participating physical points and material nodes, boundary codes, incidence
lists, relevant fields, and dirty-state flags. Include old/dead/refinement-parent
status explicitly. A `.str` export alone is insufficient if it omits internal
state or renumbers entities. Diagnostic writing must not rebuild or compact the
mesh; doing so can hide the original failure.

Add a read-only checker that reconstructs incidence from live leaf triangles.
Compare it with stored connectivity and point/node ownership. For local checks,
include the external attachment faces so a cropped patch is not mistaken for
an exposed surface. Check the complete mesh at pass boundaries to detect effects
outside the expected patch. Run geometry and ownership checks independently.

Trace GAS creation, deletion and transfer, but also SURFACE/BACK flag changes:
a missing-GAS report can mean either lost ownership or an incorrectly set flag.
Record cascades through `rem_ele` and `rem_1nd`. The expected GAS requirements
must follow the engine's boundary contract, including back boundaries.

Reduce the failing patch while keeping its external boundary conditions and
region identities. Store both an accepted input and the invalid output. If the
input is already invalid, move the capture upstream until the first invalid
transition is found. Do not redesign a downstream operation to accommodate
corrupt input.

## Implement in bounded stages

### Stage A — standalone geometry proposal

Start with internal single-region removal. Build the cavity from directed
incident-triangle connectivity, preserving region identity and external faces.
Do not infer adjacency from coordinate proximity or angular order. Reject
nonmanifold or disconnected inputs explicitly. Produce candidate triangles in
private storage; do not call mutating legacy database constructors during
proposal generation.

Use exact orientation signs for the stored coordinates. A cavity triangulation
must cover exactly the old triangle union with positive triangles, no duplicates,
no crossings and no holes. Ear clipping is acceptable only if its assumptions
are proved for the input cavity and its output is independently checked. A
Delaunay implementation is not automatically a solution to constrained material
interfaces or ownership.

Test concave and convex stars, collinear vertices, narrow slivers, high valence,
permuted IDs/order, and invalid incidence. Translation/scale tests must account
for binary32 rounding: transformations that change the represented geometry
cannot be expected to produce identical connectivity. Compare valid partitions,
not a particular diagonal choice.

### Stage B — complete physical-point proposal

Extend the proposal to all material nodes at a physical point, and to the
edge-collapse path identified by the reproducer. Include surviving endpoint
ownership, GAS nodes, region membership and boundary labels. Keep one coherent
operation rather than calling single-material proposals sequentially on globals.

Initially permit boundary removal only when the old region domains and boundary
segments are exactly preserved. Preserve corners, curved interface vertices and
unsupported junctions with a reported reason. This is a limited supported domain,
not evidence of general robustness. If the captured failure requires a genuinely
moving/collapsing interface, specify and implement its permitted geometric change
and transfer law before claiming that case repaired. Merely declining every
removal cannot satisfy the activation acceptance gate.

### Stage C — field transfer and commit

Inventory every field stored on nodes. Define which are conserved chemical
quantities, interpolated state, and recomputed derived quantities. Use the
existing finite-element integration convention when measuring inventory. For
unchanged region domains, conservative projection must preserve the specified
species totals within a documented roundoff/error budget and retain required
positivity. Interface fluxes and physically disappearing material need explicit
rules; do not redistribute dopants between materials just to balance a total.

Build replacement incidence and ownership before retiring old entities. Reserve
all storage before commit, validate that the source patch has not changed, and
publish the complete replacement as one operation. Avoid generic removal helpers
whose cascading behavior violates the prepared ownership map. Any new retirement
primitive needs an explicit contract and tests. Do not silently weaken existing
helpers for all other callers.

Inject allocation/preparation failures and prove that rejected proposals leave
coordinates, fields, ownership, incidence, flags and cache state unchanged.
Use semantic comparisons rather than struct `memcmp` with padding/pointers.
The existing fixed-count `mesh_trial` snapshot is insufficient for this work.

### Stage D — integration and acceptance

Replace both relevant removal entry paths through the new operation, or clearly
identify any still-unsupported path. Revalidate each proposal against the current
mesh: selections made before earlier removals may become stale. Count committed,
rejected and unsupported operations and report reasons.

Require the local failing control to fail on the original code and pass on the
replacement. Then require the captured activation to complete with actual valid
removal activity, engine self-tests, independent partition and ownership checks,
field-inventory checks, and strict ASan/UBSan. No duplicate deletion or synthetic
boundary cleanup may be needed to pass.

Only after these gates run fresh NMOS/PMOS through unchanged process recipes.
Keep numerical/physical validation separate: correct polarity, inventory trends,
mesh/time-step convergence and furnace-budget completion are additional gates.
Repeat mesh-density and mask-perturbation cases after the baseline succeeds.

## Preserve established work

Keep the independently reproduced parser, mesh-checker bounds, finite-segment
intersection, live-interface ownership and diffusion-initial-step repairs.
Do not reset the dirty repository. The restored old coarsener still contains
previous duplicate-cleanup code; remove it as part of replacing the operation,
not by blindly reverting unrelated source changes.

Relevant focused checks:

```sh
python3 suprem4/tests/test_float_orientation.py
python3 suprem4/tests/test_partition_audit.py
sh suprem4/tests/run_key_reader.sh
sh suprem4/tests/run_etch_interface.sh
sh suprem4/tests/run_etch_corner.sh
```

Build sanitizer candidates with `build_native64.py --ordering amd --jobs 4
--sanitize`. Use `UBSAN_OPTIONS=halt_on_error=1`; record ASan options and whether
leak detection was enabled. A LeakSanitizer environment failure is not a passing
leak check. Full-engine leak checking remains a separate open item.

Run `run_rit_nmos_spacer_regression.py` with `--device nmos` and `--device pmos`,
plus the LOCOS replay (`run_rit_locos_strip.py --case nmos_locos_strip`), against
the actual new candidate. Supply `--binary` and a fresh `--output` to each.
Use `run_diffusion_restart.py --help` for its timing regression interface.
The final process runner is `decks/subcmos2014/run.py`; its default target is
`metal1`. A resumed checkpoint is useful diagnosis but cannot replace fresh-flow
acceptance. Do not promote the launcher until the release gates are satisfied.

## What to report at each checkpoint

Record source/build/input hashes, exact command, exit status, operation counts,
checker results, sanitizer options and field residuals. Preserve failing controls.
Distinguish PASS, FAIL, TIMEOUT, INTERRUPTED and NOT RUN. Do not label a run passing
because a picture looks plausible or because a panic mesh passes geometry checks.

Update the repair report with the proven cause, the changed operation, evidence,
and remaining limitations. Your first update should deliver the reduced failing
removal and causal trace. Your final update must state separately whether local
mesh correctness, captured activation, fresh process execution, and physical
validation have passed.
