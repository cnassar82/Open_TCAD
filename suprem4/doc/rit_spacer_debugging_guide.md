# Sol: fix the RIT spacer boundary at its source

This companion to `rit_cmos_reliability_plan.md` narrows the next assignment to
one defect. Complete this assignment before expanding into PMOS doping or the
full reliability matrix. Prepared 2026-11-14 from source and the captured replay.

## Repair status

The captured failure now has a tested repair. Read
[the repair report](rit_spacer_boundary_repair.md) and
`rit_spacer_fix_results.json` before repeating diagnosis. The guide below
preserves the original evidence and debugging procedure. Both shared junction
identity and exposure of the uncovered gate required repairs; accepting a
successful but partial resist deposit would have missed the second defect.

## Assignment

Find the first mutation that creates overlapping exposed edges during NMOS
spacer formation, fix its geometric/topological cause, and make the archived
replay pass through photoresist deposition with independently valid geometry.

The previous session successfully added detection. It did not fix the engine.
Do not spend another session rediscovering the downstream crash. Start from the
fixture and the coordinates below. A local fix is appropriate if it establishes
a general invariant; a special case for these coordinates is not.

## What the evidence establishes

The persistent input is:

`suprem4/tests/fixtures/rit_nmos_spacer/ldd.str.gz`

Its uncompressed SHA-256 is:

`06f0e640fbe652c274905f3b22690ee032342b0e84b17fd1655a94fbddc1c38e`

The three operations in `run_rit_nmos_spacer_regression.py` are:

```text
deposit oxide thick=0.5 divisions=10 space=0.04
etch oxide dry thick=0.5
deposit photoresist file=resist_top.dat divisions=8
```

The resist file defines a flat top from (0, -2) to (10, -2), in micrometers.
On the diagnosed build, the first operation produces a saved structure that
passes the current audit. The second exits successfully but produces overlapping
exposed edges. The third aborts with `Region crosses itself.`

**Passing the current audit is limited evidence.** It checks triangle orientation,
connectivity, fields, and exactly horizontal/vertical boundary overlap. It does
not establish absence of diagonal crossings, triangle interior overlap, missing
material, or incorrect exposure labels. The earlier deposited structure still
needs the local partition check described below.

The malformed spacer mesh contains these points, in micrometers:

| Label | x | y |
| --- | --- | --- |
| A | 4.570329911075532 | 0.07342216122196987 |
| B | 4.570341552607715 | 0.07342216122196987 |
| C | 4.570333403535187 | 0.07342216122196987 |
| D | 4.5706500532105565 | 0.07342216122196987 |

Spatial order is **A, C, B, D**. The exposed traversal instead contains
**A → B → C → D**. Thus interval C–B is traversed repeatedly. Its length is
about 0.000008149 um, or 0.008149 nm. Small size does not make the topological
contradiction valid.

In one saved replay, A/B/C/D were point IDs 3879/6972/4192/4725. Edge A–B belonged
to oxide region 36; B–C and C–D belonged to oxide region 21. IDs are unstable
across rebuilds, compaction and structure reloads. Match coordinates and record
operation-local identities rather than hard-coding these IDs.

Region 21 was a very thin remnant of deposited oxide. Region 36 was a larger
oxide region. This suggests a shared-interface or rounding problem, but does
not yet identify which function is wrong. No reliable interior-overlap area
has been established; avoid drawing conclusions from tiny areas computed by
subtracting large shoelace terms.

## Reproduce once, then work on the first bad operation

Run from the workspace root:

```sh
cd /home/chris/Open_TCAD
python3 suprem4/tools/build_native64.py --ordering amd --jobs 4
```

Use the exact executable printed by the builder:

```sh
python3 suprem4/tests/run_rit_nmos_spacer_regression.py \
  --binary /absolute/path/printed/by/builder/src/suprem4
```

The default creates a fresh output directory. Read its `result.json`,
`replay.log`, `spacer_deposit.str`, and `spacer_etch.str`. A failing replay is the
expected starting result. Keep this build as the control for your eventual fix.
For detailed boundary diagnostics, set `SUPREM4_MESH_DIAG=1` on the same command.

Do not reuse a populated output directory: the replay runner now rejects it so old `.str` files cannot survive a failed
run and be accepted as new output. Separate imported
checkpoints from newly produced checkpoints in reports.

Archive the build manifest and run manifest alongside the minimized case.
The fixture metadata's captured binary identifies a replay consumer; because
previous runs resumed across builds, do not assume it created every upstream
checkpoint. Fresh full flows remain necessary later.

## Trace the mutation, not the panic

Relevant call path in the current source:

```text
user_etch (refine/etch.c, dry branch)
  find_surf → construct translated etch polygon
  etch, once per exposed matching region
    ck_clock(subject and clip)
    sub_skel (refine/skel.c)
    for each retained component: grid → interp
    remove old triangles/region
    bd_connect_material (dbase/make_db.c)
      rebuild neighbors/edges and exposure-related state
```

The eventual resist deposition calls `find_surf`, `build_skel`, and `ck_clock`.
Those functions reveal the bad input surface; they are not yet implicated as
the creator of the overlapping spacer edges.

Add diagnostics gated by `SUPREM4_MESH_DIAG`. At each boundary below capture
coordinates with 17 significant digits, directed endpoints, region/material,
edge owner triangles, boundary codes and flags. Dump only the local neighborhood
around x=4.57034 um plus a map to the complete operation. Avoid giant unstructured
logs that hide the first transition.

1. **Before dry etch:** capture the deposited mesh's local interfaces, including
   both adjacent oxide regions. Check their triangle interiors and interface
   connectivity, not only exposed edges.
2. **After building the etch polygon:** capture the local polygon segments and
   their unrounded double intermediate values before storage as floats.
3. **Before and after each `sub_skel`:** capture the subject ring, clip ring, all
   crossings, retained arcs and output rings. Include already-processed neighbors
   of the current region: individually valid rings can assemble incorrectly.
4. **Before and after `grid`:** determine whether subdivision or triangulation
   moved a point or replaced one shared segment with inconsistent pieces.
5. **Before and after `bd_connect_material`:** determine whether geometry is
   already wrong or only interface/exposure classification changes.

Write down the first invariant violation as a before/after pair. Do not make a
production fix until you can name the operation that creates it. It is possible
that an earlier undetected invalid partition is merely exposed by the etch.

## Distinguish these hypotheses with evidence

| Hypothesis | Discriminating observation | Appropriate repair direction |
| --- | --- | --- |
| Independent clipping gives neighboring regions different subdivisions | One side has A–B while the other has A–C–B, before triangulation | Share intersection identities and propagate splits to all incident regions |
| Collinear overlap is handled as an arbitrary endpoint hit | Reordering equivalent clip edges changes retained arcs or intersection choice | Represent overlap intervals explicitly and node both endpoints before traversal |
| Rounding changes crossing order or collapses a thin face | Double intermediate coordinates have a different order/orientation from stored floats | Define representability and shared rounding rules; evaluate geometry storage precision |
| Refinement introduces a fold | Output rings are valid before `grid`, invalid afterward | Preserve constrained segments and update both sides of a shared interface |
| Exposure rebuild exposes an internal interface | Partition geometry is conforming before and after, but internal edges acquire gas/exposed status | Derive exposure from the assembled material partition and correct adjacency |
| Input already contains overlapping material | Independent local audit fails before dry etch | Move the minimal replay upstream to the operation that first violates the partition |

Do not choose a hypothesis solely because it matches a familiar earlier bug.

## Specific source hazards to inspect

* `node_multiple_crossings` in `refine/skel.c` handles `det != 0` intersections
  and retains groups with more than one crossing on a subject edge. That is not
  a complete noding algorithm for collinear overlaps or all endpoint contacts.
  Determine whether this case actually reaches the missing case before editing.
* `edg_crs` tries true intersections, then falls back to an endpoint accepted by
  `clip_on_segment`. A collinear overlap is an interval, not a unique crossing.
  Returning the first matching endpoint may depend on edge iteration order.
* `same_crossing`, `same_clip_point` and `clip_on_segment` have related but distinct
  tolerances. Increasing one can make topology decisions disagree with another.
* `sp_edge_keep` preserves endpoint positions, but preserving endpoints alone
  does not guarantee that all incident regions share the same subdivisions.
* `bd_connect_material` rebuilds exposure using triangle boundary codes. A
  geometrically overlapping or unnoded interface can look like two exposed
  edges to that code. Removing gas nodes cannot repair such a partition.
* `dbase/geom.c::intang` uses `acos` and a cross-product sign. The recorded
  backtrack contains 180-degree turns and a signed-zero ambiguity. Replacing
  this with `atan2` may improve angle numerics, but cannot make overlapping
  boundary segments valid. Keep the overlap regression independent of angle sum.
* The translated determinant edit in `geom/misc.c::area_tri` concerns a separate
  interpolation cancellation case. It is present in the build that still
  reproduces this spacer failure. Do not attribute a spacer fix to it.

## What a coherent repair must guarantee

A shared interface is represented by the same ordered geometric segments on
both sides. Every intersection and overlap endpoint is inserted consistently
before face classification. A segment has a well-defined material on each side;
a same-material internal interface is not an exposed gas surface.

For a collinear overlap, project endpoints onto a common segment parameter,
compute the overlap interval, split both participants at its endpoints, and
classify the resulting intervals using adjacent faces. Preserve distinct nodes
where material-specific fields require them while sharing geometric points.
Do not blindly merge nearby points or nodes with different physical meaning.

A robust local repair can satisfy this contract within the existing engine.
If independent per-region clipping cannot preserve it, isolate a planar
arrangement/material-partition operation behind a narrow interface. First
prototype it against the minimized polygons; require a clear mapping back to
materials, boundary conditions and transferred fields before integrating it.
An external Boolean library does not supply those semantics automatically.

Double arithmetic helps predicates; double coordinate storage may be necessary
for some features. Neither is proof of correct topology. If storage conversion
collapses a genuine feature, define an explicit resolution failure or a validated
geometry policy with error bounds. Do not silently delete slivers by area.

## Smallest useful tests and acceptance gates

After locating the first bad mutation, extract its two adjacent subject rings
and clip polygon into a small fixture. Keep the full LDD replay as integration
coverage. A useful reduced test must exercise the same failing code path.

Required cases for the repaired operation:

* The captured shared-interface/collinear-overlap case, with independently known
  expected material intervals and area.
* Reversed ring orientation, rotated starting-edge index and swapped processing
  order of adjacent regions. The assembled partition must be equivalent.
* Whole-geometry translations and consistent unit/scale variations within the
  supported resolution. Distinguish these from physical mask shifts.
* Endpoint touch, partial collinear overlap, complete shared edge, and a small
  positive gap that must remain open. Do not tune a tolerance to one fixture.
* A deliberately invalid partition that still triggers the independent audit.

Acceptance sequence:

1. Demonstrate the reduced regression fails without the fix and passes with it.
2. Verify output polygons, shared segmentation, no crossings/interior overlaps,
   and material-area accounting with an independent oracle. Use translated
   determinants; preserve the actual stored-coordinate geometry in the oracle.
3. Make the full spacer replay pass all three stages. Check the resist covers
   the intended full lateral domain, not just a valid-looking partial surface.
4. Run the focused sanitizer suite and existing LOCOS/STI checks appropriate to
   changed code; keep the current audit active.
5. Resume NMOS through heavy implants and report geometry separately from doping.
   Only after this gate advance to the next item in the reliability plan.

A successful simulator exit, a smoother viewer image, or an angle sum of 2*pi
is insufficient evidence by itself. Preserve the original recipe during these
checks, including equal 0.5 um spacer deposition and etch depths.

## Report back in this format

* First invalid operation and exact violated invariant.
* Before/after minimal geometry and why the repair fixes the general case.
* Failing-control result, repaired result, and independent audit results.
* Exact source/build/data identifiers and commands to reproduce.
* Remaining failures, including dopant polarity and full-flow validation.

The immediate deliverable is a proven spacer-boundary repair. The broader plan
remains the release standard; this document provides the next bounded task.
