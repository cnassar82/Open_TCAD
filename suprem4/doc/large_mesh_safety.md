# Large-mesh memory safety

2026-11-08

## Reproduced failure and changes

The previously installed `bin/suprem4` segfaults on
`tests/large_mesh.in`, a 160 by 160 micron silicon rectangle. GDB locates
the failure in `add_list`, called by `build_reg` during mesh creation,
with item 32760. Connectivity capacity grew in increments of ten and
wrapped its signed 16-bit counter on the next allocation.

- `list_str` now uses `int` counters, matching mesh indices. Capacity grows
  geometrically, with checks for integer/byte-size overflow and allocation
  failure. A failed reallocation preserves the previous pointer and counters.
- `kill_node` sizes its neighbor and angle scratch arrays to the actual
  neighbor count. Previously these buffers held 100 and 50 entries,
  respectively. Both allocations are freed on normal and low-valence exits.
- Node and point neighbor collectors use growing lists instead of static
  200-entry buffers. Point collection uses separate material-node and
  triangle loop counters, so it visits every material node at an interface.
  Empty lists stay empty during duplicate elimination.
- `tnabor` reserves space for all nine vertex-pair comparisons. Duplicate
  or degenerate triangles still stop the simulation, but cannot overwrite
  its two-entry stack array before the topology check. Its diagnostic
  includes triangle, region, and physical-point indices.

These changes do not repair or discard invalid topology. Existing local
duplicate-removal code in `kill_node` predates this change.

## Build and checks

Rebuild all consumers of `geom.h`; its in-memory structure layout changed.
The current compiler needs an explicit dialect for the legacy C sources:

```sh
make -C suprem4/src suprem4 COPT='-O1 -fcommon -std=gnu89'
sh suprem4/tests/run_mesh_memory.sh
python3 suprem4/tests/run_large_mesh.py --keep
```

The build produces `suprem4/src/suprem4`. The installed binaries and
`suprem4gs` launcher were not replaced.

Validation results:

- AddressSanitizer, LeakSanitizer, and UndefinedBehaviorSanitizer pass for
  100,000-entry lists, removal/reuse, forced reallocation failure, integer
  capacity exhaustion, node valences 0/2/3/50/51/100/101/512, and valid plus
  invalid triangle-neighbor inputs. Tests exercise the real neighbor
  collectors, including multi-material point connectivity and empty lists.
  The node tests mock skeleton operations;
  they check scratch storage and ordering, not cavity triangulation physics.
- The real executable creates and writes 25,921 points and 51,200 triangles.
  The output has no duplicate/zero-area triangles, all referenced coordinates
  exist and are finite, and the total area is 25,600 square microns.
- The existing installed executable fails that same large-mesh test with
  SIGSEGV; the debugger identifies connectivity-list growth as the crash site.

## Oxidation follow-up

The subsequent [oxidation cavity fix](oxidation_mesh_fix.md) resolves the
0.0006 crash described below. Both spacings now finish; the saved-reference
oxide-area discrepancy remains. The following records the state before
that fix.

### Original oxidation findings

The full RIT PMOS regression is **not passing**. Both the installed binary
and this candidate finish surface spacing 0.0005, but miss the stored oxide
area tolerance from `700_ox.str` onward (approximately 3%). The installed
binary aborts at `tnabor` for spacing 0.0006; the candidate still segfaults
during oxidation node removal for that spacing. Completion at 0.0005 is
therefore not a physics acceptance result. Tight geometry and solution
comparisons between the installed binary and candidate also differ; those
binaries are not a controlled same-source comparison.

Logs and structures from this investigation are in:

- `/tmp/suprem4_pmos_regression_edqg8iva` (installed executable)
- `/tmp/suprem4_pmos_regression_wwl3r6q0` (candidate executable)
- `/tmp/suprem4_pmos_regression_zvor8zq2` (final candidate, same failures)
- `/tmp/suprem4_large_mesh_dyl1o95o` (large-mesh failure)
- `/tmp/suprem4_large_mesh_ebqfmwvp` (large-mesh success)
- `/tmp/suprem4_large_mesh_h60joi94` (final candidate large-mesh success)

A full-run GDB trace also reproduces the candidate oxidation failure:
`ck_flip -> flip -> geom -> bd_connect -> grid_sub -> solve_time -> diffuse`.
The log is `/tmp/suprem4_oxidation_full_debug/gdb.log`. This identifies a
downstream crash site; it does not establish which earlier operation first
corrupted the topology. No recovery guard was added to `ck_flip`.

Restarting from the candidate's `4000_ox.str` and running the remaining
deck completed under GDB; it did not reproduce the full-run failure. A
restart changes the numerical/mesh trajectory and is not yet a reliable
minimal reproducer. The next topology investigation needs local cavity
provenance upstream of that backtrace, as described in
`mesh_robustness_notes.md`.

The executable remains 32-bit and region counts and other legacy limits
remain. This patch removes specific memory corruptions, not all scale limits.
