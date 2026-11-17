# Native 64-bit scaling and reliability

The native candidate is built with:

```sh
python3 suprem4/tools/build_native64.py
```

The tool prints the candidate's absolute path. Each build uses a fresh source
snapshot under `suprem4/build`, copies its runtime data, records source/data
SHA-256 hashes and compiler arguments in `build.json`, and retains `build.log`.
It uses `Makefile.proto` so generated dependencies on 32-bit system headers
cannot leak into the native build. Pointer/integer conversions and incompatible
pointer types are compiler errors. It checks that the result is a 64-bit ELF.
The builder leaves source-tree objects and installed launchers untouched.

## Installed launcher

The `suprem`, `suprem4`, and `suprem4gs` commands under `suprem4/bin` now link
to `tools/launch_suprem4.sh`. This launches `build/validated-native64/src/suprem4`
with the accepted build's model, key, and implant files. It resolves the
installation from its own location and preserves quoted command arguments.
Both `SUP4KEYFILE` and the older `KEYFILE` override point to the matched key file.

Original binaries and the previous launcher are backed up, with SHA-256 hashes,
under `suprem4/build/launcher_backup_20261109T121141Z`. To restore one entry point,
replace its symlink with the corresponding backed-up file.

All three entry points pass a launch/mesh smoke test from directories and deck
paths containing spaces, even with stale runtime-data environment variables.
The installed launcher also passes the 25,921-point / 51,200-triangle construction
and connectivity regression. Launch smoke artifacts are retained in
`/tmp/suprem4 launcher smoke iz_258wd`; the larger launch regression is in
`/tmp/suprem4_large_mesh_5bfqyu8v`.

## Build options

For the accelerated 2D ordering path, use:

```sh
python3 suprem4/tools/build_native64.py --ordering amd
```

This requires SuiteSparse AMD headers and the AMD library (available on this
machine). The default remains `--ordering legacy`, with no new library
dependency. AMD changes the matrix permutation, not the process equations or
physical model coefficients. It is used only for 2D ordering; 1D retains its
existing ordering. The adapter converts SUPREM row pointers, inverts AMD's
new-to-old permutation into SUPREM's old-to-new map, checks the permutation,
and supplies AMD's fill estimate to symbolic allocation.

On the captured initial fourfold-density graph (41,237 nodes and 283,909
adjacency entries), legacy ordering takes 28.395434 CPU seconds and AMD takes
0.016690 CPU seconds, about 1,700 times faster for ordering alone. Both outputs
are checked as complete permutations. This is not a full-process speedup claim.
The captured map is `/tmp/suprem4_ordering_4x.bin`; results are
`/tmp/suprem4_ordering_legacy.json` and `/tmp/suprem4_ordering_amd.json`.
`tests/benchmark_ordering.c` can replay either backend on that native-int map.

The AMD baseline completes with all ten saved meshes valid. Its nine process
checkpoints pass comparison against the final legacy native control at 0.1%
material-area and 0.005 µm bounding-box tolerances. The ordering regression also
passes sanitizers when compiled with AMD (`SUPREM_TEST_AMD=1`).

## Corrections

* Parser objects were passed to 59 command/helper parameters declared as `int`.
  On 64-bit systems their addresses were truncated. Those parameters now use
  `struct par_str *`; command callbacks and parameter accessors have typed
  interfaces.
* Shared allocation macros cast byte counts to `unsigned`, truncating requests
  above 4 GiB. They now preserve `size_t` and check multiplication overflow.
* Sparse-matrix neighbor IDs were signed 16-bit values. Fourfold-density PMOS
  initialization exceeded 32,767 nodes, and negative IDs crashed matrix ordering.
  Diffusion, finite-element, elastic-flow, and circuit neighbor buffers now use
  `int`, matching the mesh database. Their shared interfaces are declared in
  `sparse.h`.
* The diffusion adjacency callback ignored its caller's buffer capacity. It now
  reports insufficient capacity before writing. Matrix construction grows its
  map before writes, checks index limits, and releases neighbor scratch storage.
* The minimum-fill ordering pool used `realloc` even though live adjacency lists
  held pointers into it. Pool growth now allocates stable blocks; all blocks are
  released after ordering. This is a separate latent defect from the observed
  16-bit-ID failure.
* Rectangular triangulation could clip a right-angle corner whose triangle
  contained the remaining boundary vertices on its diagonal. That left a flat
  quadrilateral and dereferenced a missing split. Rectangle decomposition now
  requires positive-area halves. Ear clipping uses coordinate-based interior
  tests and preserves every boundary vertex. Invalid geometry is still rejected.
* The accelerated left-shifted mask case exposed a second path: angle/proximity
  heuristics accepted a cut that left a self-intersecting, nearly collinear
  quadrilateral. All heuristic ears now require a coordinate-based empty-triangle
  test. General and rectangular diagonal candidates must leave positive-area
  halves and avoid crossing nonincident boundary edges. The original failing
  accelerated mask run is retained under
  `/tmp/suprem4_mask_perturbations_w5ehc819/all_left`; its exact transient quad is
  in `trace/capture.log`. The updated accelerated left-shifted case completes.

The geometric regression uses real allocation, skeleton splitting, and triangle
creation. Its five-vertex right triangle with a subdivided hypotenuse crashes
the previous `rect_tri()` at the same location as the full-process failure.
The corrected code emits three positive triangles with the original area and
all five boundary vertices, tested from every cyclic starting vertex.

## Validation

Validated accelerated executable:
`/home/chris/Open_TCAD/suprem4/build/validated-native64/src/suprem4`
(points to `native64_2r55t88h`). The final legacy-ordering comparison build is
`suprem4/build/native64_s0_rs198/src/suprem4`.

Complete PMOS runs with the final accelerated candidate:

| Target initial density | Initial points | Final points | Final triangles | Saved meshes |
| --- | ---: | ---: | ---: | --- |
| 1× | 9,815 | 5,639 | 10,598 | 10/10 valid |
| ~2× | 19,809 | 9,984 | 19,078 | 10/10 valid |
| ~4× | 40,033 | 17,702 | 34,336 | 10/10 valid |

All three processes returned zero. Calibration targets are within 2% of the
requested initial point-count ratios; adaptive processing changes later counts.
The fourfold case exercises sparse ordering at 41,237 nodes and recovers from
small-time-step retries during post-implant diffusion. Final oxide area changes
by -0.8369% at doubled density and -0.6591% at fourfold density relative to the
current control; silicon area changes by less than 0.034%.

Observed process elapsed times were 88.2 s, 250.1 s, and 643.1 s respectively,
in a mixed parallel validation workload; these are not controlled full-process
speedup measurements. Scaling results and area measurements are retained in
`/tmp/suprem4_mesh_density__y54u82m/results.json` and `summary.json`.
Together with the final mask-shift suite, all 57 saved meshes pass the existing
finite-coordinate/field, positive-area, duplicate-triangle, edge-ownership,
and reciprocal-neighbor audit. Superseded longer runs were explicitly stopped
after the geometric acceptance fix and are marked `superseded.json`; they are
not included in these acceptance results.

All nine focused regressions pass AddressSanitizer and UndefinedBehaviorSanitizer:
allocation sizing, mesh memory, quad splitting, trial rollback, thin polygon
ears, retained etch corners, ordering-pool growth, sparse capacity, and
rectangular ears. Allocation testing verifies requests above 4 GiB with recording
allocators rather than physically allocating gigabytes. The sparse regression
uses 40,001 IDs, a high-degree node, undersized buffers, and repeated construction.

The final native candidate constructs, writes, and validates a 160 µm square
with 410,881 points and 819,200 triangles. Simulator elapsed time was 3.03 s and
peak RSS was 560,620 KiB (about 547 MiB), excluding the Python audit. All triangles
are positive, connectivity checks pass, and domain area remains 25,600 µm².
This is a mesh construction benchmark, not an oxidation benchmark.

Final-candidate large-mesh result:
`/tmp/suprem4_large_mesh_n_qc1cno/results.json`.

A clean 32-bit compatibility build also passes the 25,921-point,
51,200-triangle construction regression. The 64-bit executable is an additional
validated build path; the legacy build defaults are retained.

Two independent final-candidate PMOS control runs produce byte-identical
`pmos.str` files (one runner saves an additional initial-mesh checkpoint).
The final control also passes a comparison against the previous 32-bit control
with 0.1% material-area and 0.005 µm bounding-box tolerances. These are portability
checks, not a replacement for the separate stored-reference accuracy criteria.

The final candidate completes a control and translations of all mask groups by
0.1 µm left/right. All 27 saved meshes pass the topology and finite-field audit.
Final oxide area changes by -0.3501% and -0.2704% relative to this candidate's
control. Results: `/tmp/suprem4_mask_perturbations_gxsye8hh/results.json`.

Repeat large construction and complete PMOS scaling tests, substituting the
candidate path printed by the builder:

```sh
python3 suprem4/tests/run_large_mesh.py --binary /absolute/path/to/suprem4 \
  --cells 160 320 640 --keep
python3 suprem4/tests/run_mesh_density.py --binary /absolute/path/to/suprem4 \
  --ratios 1 2 4 --jobs 2
python3 suprem4/tests/run_mask_perturbations.py --binary /absolute/path/to/suprem4 \
  --groups all --jobs 2
```

Mesh IDs and matrix offsets remain `int`; 64-bit addressing does not make their
range unlimited. Sparse-ordering and solver costs also remain scaling limits.
These tests establish the audited mesh properties and tested process behavior,
not physical calibration. The older stored-reference oxide discrepancy remains
a separate issue.

## Eightfold initial mesh, 2026-11-10

The installed native AMD build completes the full PMOS flow with 79,900 initial
points (8.14 times the 9,815-point control, within the runner's 2% targeting
tolerance). The final structure has 33,374 points and 65,382 triangles. All ten
saved meshes, including the initial mesh, pass the topology and finite-field
audit. Simulation time was 2,496.3 seconds (41.6 minutes); timing is not a
controlled benchmark against earlier runs.

Final oxide area is 49.72513749 square microns, 1.333% below the control's
50.39676068. Completion and valid topology do not establish mesh convergence
or physical calibration. This result uses the installed build, before the
separate experimental LOCOS/STI clipping changes.

Full artifacts: `/tmp/suprem4_mesh_density_42mhzdux/points_8x`.
Machine-readable result: `doc/mesh_8x_results.json`.
