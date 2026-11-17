# First invalid PMOS activation removal

Status: causal diagnosis and reduced fixture, 2026-11-15. The simulator is not
repaired or accepted. This result supersedes the earlier speculation that a
later removal at point 4575 first corrupted the mesh; that inference came from
`panic.str`, which serializes triangles marked for deletion. An independent
audit of that raw panic export therefore gave a misleading overlap.

## Reproduction and evidence

The unchanged fixture `tests/fixtures/rit_pmos_activation_seed/` runs
`diffuse time=0.000001 temp=825 dryo2`. The captured original-code replay is
`decks/subcmos2014/results/sol_activation_baseline_01/`: exit -6, `Region
crosses itself.` The trace build `native64_j7tth8sm` adds opt-in diagnostics
without changing the removal decision. Its output is
`decks/subcmos2014/results/sol_activation_causal_final/`.

From that trace, the failing call path is:

```text
grid_sub: CknMcLen edge-removal pass
  rem_edg(edge 12550): remove point 3796, survivor point 3567
    rem_pt(point 3796): exposed, GAS node 5101 and oxide node 5102
      rem_1nd(5101): GAS node is removed first
      kill_node(5102): old star has three live triangles in region 14
        angle-ordered neighbor polygon: 3567, 3566, 3793, 3797
        new closing edge 3793-3797 crosses old outer edge 3567-3566
        grid -> triang -> ck_clock: Region crosses itself
```

Point 3796 is an open boundary fan, not an interior vertex. Its two boundary
rays end at points 3793 and 3797. Exact orientation of (3793, 3796, 3797) is
nonzero, so connecting the endpoints changes the represented region boundary.
The old three triangles each have positive exact orientation. The proposal is
self-crossing before the triangulator runs; the triangulator correctly detects
it. This is an algorithmic defect in the removal decision and cavity model,
not a floating tolerance or a missed cleanup step.

At `rem_pt` entry, point 3796 has live GAS node 5101 and live oxide node 5102
and is marked SURFACE. `rem_pt` stores both node IDs and processes GAS first.
After that deletion, the still-exposed point has no live GAS node while the
oxide cavity is being constructed. This is an intermediate ownership breach.
The whole physical-point removal never completes because the cavity panics.
Moving GAS deletion later could avoid this intermediate state, but would not
make the crossing cavity valid. Both geometry and ownership need one prepared
operation.

The input fixture is geometrically valid: 12,688 triangles and 31,254
candidate triangle pairs pass the independent partition audit. An opt-in
live-triangle export at `kill_node` entry also passes the triangle-only audit:
12,737 live triangles and 30,473 candidate pairs. Its boundary codes are stale
mid-operation, so the full boundary-code audit is not meaningful there. A
completed-operation ownership audit reported no missing GAS point before this
failing operation. These observations support the local cause; they do not
prove every other mesh invariant or field quantity at this intermediate stage.

## Reduced reproducer

`tests/fixtures/rit_pmos_first_bad_removal/patch.json` contains the five exact
binary32 point coordinates (hex notation), three old oriented triangles,
angle-ordered proposal, call path and ownership sequence.

```sh
python3 suprem4/tests/test_removal_patch_fixture.py
python3 suprem4/tests/test_removal_patch_fixture.py --legacy-control
```

The normal test passes after checking positive old triangles, open-fan boundary
endpoints, noncollinearity, ownership and the crossing proposed edge. The
`--legacy-control` command intentionally exits 1 because the legacy polygon
crosses itself. The reduced test is an independent geometric/ownership oracle;
it does not yet call a replacement production kernel. A future production test
must call the new proposal and require an explicit rejection or a supported
conservative interface move, with unchanged input state on rejection.

## Standalone proposal reference

`tests/mesh_patch_oracle.py` now provides a read-only exact rational reference
for **internal, single-region triangular stars**. It reconstructs the boundary
from directed incidence, validates old and candidate geometry, and returns an
explicit unsupported result for open boundary fans.
`tests/test_mesh_patch_oracle.py` passes convex, concave, high-valence, near
collinear, invalid-incidence and captured PMOS open-fan cases. This is a Stage A
reference for production implementation; it has no field transfer, GAS
ownership or database commit and cannot make the simulator pass activation.

## Rejected trial and implementation consequence

An opt-in exact-domain boundary preflight declined this removal, but the replay
later failed a neighbor reciprocity self-test (`12618 <== 12505`). Its source
and failed result are archived in `../experiments/rit_boundary_preflight_20261115/`.
The gate was withdrawn from active source. Rejecting this one cavity cannot
repair the broader removal algorithm.

Next implement a complete physical-point proposal that constructs geometry,
material/GAS ownership, adjacency, and field transfer in private storage. An
edge-collapse operation must explicitly define how the represented boundary
may move and how fields transfer. If the first supported version preserves the
exact region domain, it must refuse this noncollinear boundary removal before
any node deletion and report the refusal. Other accepted removals must be
independently validated; declining all removals is not success.

See [Sol’s implementation guide](sol_mesh_repair_implementation_guide.md) and
the [mesh replacement contract](mesh_patch_replacement_contract.md) for the
remaining stages and release gates.

## Subsequent candidate work

The interior proposal/atomic commit now passes this bounded replay, with real
accepted removals, under the explicit candidate option. Longer activation
revealed separate stored-coordinate motion failures caused by constrained
oxide slivers already present after spacer formation. See
[material-domain investigation](material_domain_reliability.md) for the
upstream cause, CAD-derived replacement approach, and remaining gates.
