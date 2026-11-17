# SUPREM-IV Mesh Robustness Notes

Date: 2026-05-12

This note records the recent attempt to make the oxidation/deposit/etch mesh
handling more robust, what worked, what failed, and the approach that should be
used next. The short version is that the failure is a topology corruption during
oxidation mesh update, not a simple crash in the later geometry or neighbor code.

## User-visible problem

The `rit_pmos.in` deck is sensitive to the starting surface mesh spacing:

- `line y loc=0.0 spac=0.0005 tag=surf` had been able to complete.
- `line y loc=0.0 spac=0.0006 tag=surf` crashed during a later oxidation step.

The goal was not just to avoid the crash. The generated `.str` structures also
need to remain correct when compared against the known-good structure files.

## Related work done before the mesh debugging

Several earlier improvements were made or explored:

- Added isotropic/wet etch support to the etch command. Wet and isotropic should
  be treated as the same etch mode.
- Added mesh diagnostics guarded by `SUPREM4_MESH_DIAG` in parts of the refine
  path, including etch/deposit-related triangulation diagnostics.
- Added or used regression tooling under `tools/`, including:
  - `tools/run_pmos_mesh_regression.py`
  - `tools/str_compare.py`
- Rebuilt with 32-bit flags. The installed binary was verified as a 32-bit ELF:
  `bin/suprem4`.

Those pieces are useful context, but the main robustness issue described below
is in oxidation grid update/node removal.

## What the failure actually was

The crash path was narrowed to oxidation grid update:

```text
grid_sub(dt)
  CknMcLen(...) / CknMcObt(...)
  rem_edg(...) or rem_pt(...)
  rem_pt(...) -> kill_node(...)
  bd_connect("after removing silicon nodes")
  geom(...) / nxtel(...) / tnabor(...)
```

Relevant files:

- `src/dbase/grid_upd.c`
- `src/dbase/remove.c`
- `src/dbase/kill_node.c`
- `src/dbase/make_db.c`
- `src/dbase/element.c`
- `src/geom/tnabor.c`

The strongest diagnostic finding was that the failing mesh contained two live
triangles with the same three physical points:

```text
tnabor overmatch ie=7245 in=7290 shared_points=3
regs=4,4 mats=1,1 dead=0,0
pts=(3833,3857,3832)/(3832,3833,3857)
nd=(658,758,607)/(607,658,758)
```

That is an invalid 2D triangulation. Two live triangles in the same material
cannot occupy the same three points. Once that happens, `tnabor()` is right to
fail: those triangles are not merely neighbors sharing one edge; they are
duplicates sharing all three vertices.

## What was tried and rejected

Several downstream recovery ideas were tried. They were useful diagnostically,
but they are not acceptable fixes.

### Guarding `flip()` and geometry consumers

I added checks around edge flips and edge geometry so stale or inconsistent
edge/triangle connectivity would be skipped instead of crashing. This moved the
failure but did not fix the cause. It also risked letting bad topology survive.

Conclusion: do not fix this by making `flip()` or geometry consumers tolerant of
corrupt topology.

### Changing `tnabor()`

I temporarily changed `tnabor()` to avoid stack overwrite when more than two
shared points were found and to print diagnostic details. That proved the
duplicate-triangle condition, but returning false or recovering there is wrong.

Conclusion: `tnabor()` should remain a detector, not a repair layer.

### Rebuilding neighbors more aggressively

I tried clearing/rebuilding neighbor state more broadly in `nxtel()` and
`bd_connect()`. This did not solve the underlying duplicate-triangle creation.

Conclusion: rebuilding connectivity cannot make an invalid triangle set valid.

### Duplicate cleanup after or during rebuild

I tried removing duplicate triangles either inside the database rebuild or just
before `bd_connect()`. One version let `spac=0.0006` finish, but the structure
comparison still failed by about 2.5 percent silicon area from `4000_ox.str`
onward, and `spac=0.0005` later segfaulted.

Conclusion: duplicate cleanup is a band-aid. It may hide the crash, but it does
not preserve the physical structure.

## Important lesson

The failure is upstream of the crash. The crash site is not the bug site.

The real problem is likely in the topology-changing operation that removes
points/nodes/edges during oxidation and retriangulates the local cavity. In
particular, `kill_node()` is suspicious because it builds a replacement skeleton
from neighboring nodes, calls `grid(sr, reg, FALSE, FALSE)`, and only then
removes old triangles containing the killed node:

```c
grid(sr, reg, FALSE, FALSE);

while (num_tri_nd(in) > 0)
    rem_ele(tri_nd(in, 0), FALSE);
```

That ordering may be valid only if the cavity boundary is perfectly identified
and the replacement triangles cannot overlap any existing live triangle. The
observed duplicate triangle suggests that assumption is sometimes false.

## Why the previous approach was wrong

The previous approach drifted into trial and error:

- Add a guard.
- Run the deck.
- See the crash move.
- Add another guard.
- Run again.

That can make the program finish while silently changing the structure. For TCAD,
that is worse than a crash. A finished but wrong structure is not a fix.

The acceptance standard must be:

- `spac=0.0005` completes.
- `spac=0.0006` completes.
- The generated `.str` files compare acceptably against known-good references.
- No topology invariant is violated.
- The solution is in the topology creation/removal path, not in later consumers.

## Recommended next approach

### 1. Keep the baseline trusted

Before making more mesh changes, keep production code close to the known-good
state. Diagnostics and tools are useful, but recovery changes in geometry,
neighbor, `tnabor()`, and flip code should not be left installed unless they are
purely diagnostic and guarded.

### 2. Build a focused reproducer around the failing oxidation update

The full `rit_pmos.in` run is too large for repeated blind experiments. The next
step should be to capture or construct the smallest deck/state that reaches the
bad oxidation node-removal step.

Useful targets from the observed failure:

- The bad physical points were around `3832`, `3833`, and `3857` in one run.
- The duplicate triangles were in region `4`, material `1`.
- The failing phase was `bd_connect("after removing silicon nodes")`.

The repro should stop just before the failing oxidation grid update if possible,
then run only that update.

### 3. Add provenance to triangle creation and deletion

The next diagnostic should answer this exact question:

> Which function created the second triangle with the same three points, and
> which old triangle should have been removed but survived?

To answer that, instrument triangle creation/deletion, not the later crash:

- In `mk_ele()` / `mk_ele_nd()`, record or print:
  - new triangle id
  - region/material
  - node ids
  - point ids
  - current high-level operation, e.g. `kill_node`, `grid`, `triang`
- In `rem_ele()`, record or print:
  - removed triangle id
  - node ids
  - point ids
  - whether it was part of the active node-removal cavity
- In `kill_node()`, print:
  - node being removed
  - point being removed
  - selected region
  - sorted neighbor loop
  - old incident triangles before `grid()`
  - new triangles created by `grid()`

The instrumentation should be gated by `SUPREM4_MESH_DIAG` and preferably by a
small coordinate/point/region filter to avoid enormous logs.

### 4. Check invariants at the creation boundary

The right invariant checks are:

- No two live triangles in the same material have the same three point ids.
- No live triangle has duplicate point ids.
- Every replacement triangle created by `kill_node()` lies inside the intended
local cavity.
- No triangle outside the cavity is duplicated by the replacement triangulation.

The checks should run immediately after one local topology mutation, before the
global rebuild, so the first bad operation is identified.

Avoid all-pairs checks after every node removal in a full deck; they are too
expensive and can make the run appear stalled. Use focused checks around the
active cavity or around a filtered point set.

### 5. Fix the cavity algorithm, not the symptoms

The likely real fix is to make node removal explicitly cavity based:

1. Identify the full set of triangles that form the cavity around the node or
   edge being removed.
2. Identify the cavity boundary loop.
3. Remove or mark the entire old cavity before inserting replacement triangles.
4. Retriangulate only inside that cavity.
5. Reject any replacement triangle that duplicates or overlaps a live triangle
   outside the cavity.
6. Rebuild connectivity once the triangle set is valid.

This is more invasive than a guard, but it matches the actual failure.

### 6. Validate with structure comparisons, not just completion

Every candidate fix must be run through:

```sh
tools/run_pmos_mesh_regression.py --spacings 0.0005 0.0006 --mesh-diag --keep
```

Success means both decks complete and the `.str` comparisons pass or have a
well-understood, physically justified difference. A crash-free run with a 2-3
percent material area mismatch is a failure.

## Current recommendation

Do not continue with downstream repairs. The next serious attempt should be a
short diagnostic branch focused on triangle provenance during oxidation node
removal. Once the first invalid triangle creation is known, make one local
algorithmic change in the cavity-remesh path and validate against both spacings
and the stored `.str` files.

