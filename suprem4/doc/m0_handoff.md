# M0 geometry core: handoff (state at 2026-11-22)

Read this first if you pick up the geometry-core work.  Companion documents:
- `doc/roadmap_2d_cmos_25nm.md`: why a new geometry core (goal: 2D Si CMOS down to 25 nm).
- `doc/plan_m0_geometry_core.md`: the M0 plan (WP0-WP8), acceptance criteria, and **section 4b,
  the revised WP6 plan; start there**.
- `doc/poly_oxidation/README.md` and `research_oxidation_meshing.md` (with the 2026-11-22 addendum):
  why SUPREM-IV's mesh fails in poly re-oxidation, and how Sentaurus/Synopsys handle it.
- `core/README.md`: the core's pipeline, checks and measured results per WP.

Branch `geometry-core` (pushed).  Commits: WP0-2 `5954ada`, WP3 `b40f2f2`, WP4 `8bedee0`,
WP5 `8ea880a`, WP6 work in progress `fc6b2ac`, `18fbe1f` and the commit that adds this file.

## 1. What exists and works

`suprem4/core/` is a standalone C++17 library and set of tools that read and write SUPREM `.str`
files:

| WP | Module | Status |
|---|---|---|
| 0 | CMake build; Mmg 5.8 (2D, LGPL-3) and CDT vendored; Clipper2 from `src/third_party` | done |
| 1 | `str_file`, `mesh`: exact .str I/O, semantic mesh, invariants | 108/108 checkpoints round-trip byte-identical |
| 2 | `brep`: region loops and chains (twin interfaces, pinch vertices, closed chains) | 108/108 |
| 3 | `remesh`: B-rep -> chains -> thin-film layers -> CDT -> Mmg2d -> Delaunay flips -> Ruppert splits | 108/108; worst box-method coupling >= -0.0087 (inputs down to -2.5e5) |
| 4 | `transfer`: per-material P1 transfer, P1-integral inventory restoration (min-change weighted) | inventory error <= 2e-14; SUPREM reloads and anneals remeshed structures; remeshed results mesh-converged |
| 5 | `evolve`: motion with step limits, events (thin, crumb, gap, cavity), Clipper2 repair, rebuild | shrink, merge and slit synthetic tests pass |
| 6 | SUPREM hook (`src/oxide/geometry_core.c`, `core/src/capi.cpp`) | runs, **no poly-line deck completes**; see 3. |

The 108 checkpoints are SubCmos2014 0.5 um (two runs, NMOS and PMOS, 16 stages each) and AdvCmos2014
0.25 um (18 stages x 2).

## 2. Build and test

```
# core library and tools
cmake -S suprem4/core -B suprem4/build/core_wp5 -DCMAKE_BUILD_TYPE=Release
cmake --build suprem4/build/core_wp5 -j8

# SUPREM with the core linked in (same options as the production build, plus the core)
python3 suprem4/tools/build_native64.py --jobs 8 --ordering amd --fe-solver umfpack \
        --domain-offset --geometry-core suprem4/build/core_wp5
#   -> prints the new build directory (suprem4/build/native64_XXXX)
#   A build without --geometry-core is unchanged: the hook is a weak symbol.

# core tests
C=suprem4/build/core_wp5
python3 suprem4/core/tests/roundtrip.py   $C/core_str    decks/subcmos2014/results/m0_reference_20261122
python3 suprem4/core/tests/remesh_suite.py $C/core_remesh decks/subcmos2014/results/m0_reference_20261122
$C/core_evolve_test decks/subcmos2014/results/m0_reference_20261122/nmos/initial.str

# SUPREM with the core on the six poly-line decks (about 5-10 min until they fail)
suprem4/core/tests/run_polyline_core.sh suprem4/build/native64_XXXX /tmp/polyline_run
```

Tools: `core_str IN.str [--write OUT.str]` (summary, invariants, B-rep, inventories);
`core_remesh IN.str OUT.str [options]` (remesh + transfer + all checks; `--worst-max N` lists the
worst triangles).  Plotting helper used during the work (not in the repo): a 30-line matplotlib
script that draws `c`/`t` records of two .str files in a window.  Write one; it pays off quickly.

## 3. WP6: what happens when the core runs inside SUPREM

### Mechanism
- `SUPREM4_GEOMETRY_CORE=1` turns the hook on (needs a `--geometry-core` build).  At the top of every
  `oxgrow` call SUPREM writes its structure, `suprem_core_oxidation_check` decides whether to rebuild,
  and the result is reloaded like `structure infile=` (`ig2_read` + `create_db`).  The velocity history
  slots `sol[MAXIMP-5/-6]` are not in the file and are restored from XVEL/YVEL.
- While active: SUPREM's `DetectLoop`, `grid_add` and `ChooseKillNodes` are skipped (the core owns the
  topology); the prototype `oxide_repair` is off.  SUPREM still moves the points inside its time step
  (`upd_diff` / `mesh_trial_move` in the TR-BDF stages; that motion is part of the diffusion solve).
- The core rebuilds on: first call, events (on a 1e-6 um Clipper grid), triangles compressed below
  30 % of their built area, triangles grown 2x and coarser than 2x the target size, and requests from
  SUPREM's step limiter.  A rebuild: union of triangles -> events -> `clean_shapes` (sub-delta slits and
  needles) -> `simple_shapes` (no holes) -> `mesh_from_shapes` -> remesh with points on SUPREM's
  float grid -> transfer -> float-orientation repair by flips -> checks.
- Step limits added in `oxgrow`: no point moves more than delta/4, no triangle loses more than half
  its area, no exposed vertex-edge pair closes more than half its clearance (`geometry_core_collision_dt`).

### Environment switches
| Variable | Default | Meaning |
|---|---|---|
| `SUPREM4_GEOMETRY_CORE` | off | turn the core on |
| `SUPREM4_CORE_DELTA` | 1e-3 um | event size; motion limit delta/4 |
| `SUPREM4_CORE_EVENT_DELTA` | = DELTA | event size only (tests) |
| `SUPREM4_CORE_HMAX`, `_SCALE` | 0.5 um, 1 | mesh size |
| `SUPREM4_CORE_SOLUTION_DECADES` | 0.25 | dopant-gradient refinement (0 = off) |
| `SUPREM4_CORE_COARSEN_TOL` | delta/10 | boundary simplification tolerance |
| `SUPREM4_CORE_MIN_COUPLING` | -0.5 | (unused since "first step" rebuild) |
| `SUPREM4_CORE_FORCE_FIRST`, `_EVERY=N` | off | force rebuilds (tests) |
| `SUPREM4_CORE_RELOAD_ONLY=1` | off | write + reload unchanged, no core (tests) |
| `SUPREM4_CORE_SPEED_DIAG=1` | off | print step limits (`[core-speed]`, `[core-area-limit]`, `[core-collision]`) |
| `SUPREM4_CORE_KEEP_BAD=1` | off | keep `core_embedding_fail.str` when a reload fails SUPREM's embedding check |
| `SUPREM4_CORE_TMP` | `.` | directory for the exchange files |
| `SUPREM4_CORE_DUMP=DIR` | off | save every mesh the hook receives as `DIR/core_in_NNNN.str` (test data) |
SUPREM's own `SUPREM4_MOTION_TRACE=1` prints why `mesh_trial_move` rejects a step (orientation or
embedding) and is the most useful diagnostic.

### Results (base deck: 0.4 um poly on 100 A oxide, dry O2 1000 C, 51 min = 3060 s)
| Configuration | Reached | End |
|---|---|---|
| plain SUPREM motion, no repair | 47 s | stalls (1000+ tiny steps) |
| prototype `oxide_repair` | 29 s | duplicate-triangle panic (known open bug) |
| core, first versions | 25-50 s | regeneration every step, points 16-42k, solver retries exhausted |
| core, last version | 126-145 s (other decks 58-400 s) | "Region has malformed boundary" on reload, or SUPREM rejects all motion (embedding) |

### Findings that hold regardless of the design
1. Reload mid-oxidation is exact: a viscous LOCOS run stopped at 30 of 65 min, saved, reloaded and
   finished matches the uninterrupted run to 1e-4; 201 reloads do not change a run.
2. SUPREM stores coordinates as **float cm**.  Any point closer than ~1e-7 um to an edge becomes a
   contact or a flat triangle; the core must create points on that grid and validate in floats.
3. SUPREM's validity rules for a reloaded mesh: exact float orientation; `mtest1` (a gas node exactly
   on exposed and backside points); `ck_clock` (each region one simple boundary loop: no holes, and
   probably no self-touching); `mesh_embedding_check` (no boundary vertex on any other element's
   edge, checked again for every trial motion).
4. In the poly re-ox decks every failure is at the poly/oxide/ambient corners: the gap between
   sidewall oxide and gate oxide closes as a wedge, and SUPREM's own in-step motion (compress model,
   singular stress at the corner) folds the surface there.
5. Size fields must come from a fixed reference (the grid at the start of oxidation), never from a
   mesh the core produced (thin-film layers and compressed elements make sizes ratchet down).
6. Mesh dependence is real: the original SUPREM meshes are under-resolved (B dose in Si after an anneal
   -3.6 % vs. converged remeshed results; P loss through an exposed oxide surface not converged even
   at 3x the points).

### Why WP6 does not work yet, and what to do (see plan section 4b)
The WP6 design rebuilt the *whole* mesh from a rounded polygon union every few steps, re-interpolating
every value and perturbing the whole boundary.  The production tools do the opposite: Sentaurus 2D moves
the mesh with TSUPREM-4-style local operations and remeshes globally only after the oxidation step; the
Synopsys 3D MovingMesh (US 9,477,798) keeps Lagrangian points and their values, repairs the boundary
only where a defect is detected, and interpolates only inserted points.  Revised order:
- **6a** port SUPREM's validity checks (orientation, embedding, single-loop regions, mtest1) into the
  core and run them on every output; cross-validate each by reloading into SUPREM.
- **6e** (cheap, parallel) give SUPREM the two missing local operations: interface short-edge removal
  and a per-step growth cap.  This alone may fix poly re-ox and makes the core a fallback.
- **6b** Lagrangian update instead of regeneration (keep points and values; remove/add points near
  moving interfaces; interpolate only new points).
- **6c** local event repair (only inside a window around each event).
- **6d** hook policy: 6b after each step, 6c on triggers, full regeneration only after a step or as a
  fallback.
Stop after two attempts on any single symptom and re-read section 4b; this work has a strong pull
towards chains of local patches.

## 3b. Progress on the revised WP6

**6a done.**  `core/src/suprem_check.cpp` (`check_suprem`) and the tool `core_check_suprem` hold
SUPREM's reload/motion contract as exact ports: no flat triangle in float cm (negative ones are
swapped on load by `sclk_tri`, so only flatness is fatal), `nxtel` topology, `mesh_embedding_check`,
one simple border loop per region with turn sum +-2 pi (`skel_reg` + `ck_clock`; an empty region
segfaults), and `mtest1` (gas nodes; `create_db` adds missing ones at exposed points, not at
backside points).  Strict-only rules (SUPREM loads the file anyway): counter-clockwise order as
written, one component and no holes per region, gas nodes written at every exposed point.
Cross-validation, `core/tests/suprem_contract.py` (SUPREM side: `SUPREM4_LOAD_CHECK=1` makes
`structure infile=` print the float-orientation count and the embedding verdict): core verdict =
SUPREM's verdict on all 1429 readable `.str` files under `decks/` (1382 OK, 10 WARN, 37 FAIL, among
them the saved panic states of earlier work) and on 5 mutations.  The WP6 hook refuses to write a mesh
that fails it (`capi.cpp`, rc -4); `core_remesh` reports it as the `suprem` section of the suite.
Suite on `m0_reference_20261122` + AdvCmos `polyreox_geo_20261122`: 68/72; NMOS spacers/source_drain
outputs have triangles flat in float cm (0.2 A slivers; SUPREM would panic on reload), PMOS
spacers/source_drain fail the coupling threshold (same with the WP5 build).  See section 4, slivers.

**6b done (core side).**  `update_mesh` / `update_fields` (`core/src/remesh.cpp`, `transfer.cpp`;
tool `core_update`, test `core/tests/update_suite.py`; details in `core/README.md`) keep SUPREM's
points and their values and change the mesh only where the spacing left the reference (per-material
size fields from the grid at the start of the oxidation): short interface segments and interior
points crowding a chain are removed, long segments / grown triangles get points, CDT with the chains
as constraints, encroachment repair limited to what the update introduced.  Nodes outside the
changed set keep their values bit for bit; the conserved inventory is restored on the changed points
only.  Test data: `SUPREM4_CORE_DUMP=DIR` makes the hook save every mesh it receives
(`experiments/wp6b_update_20261123/base_dump/dump`, 36 states of the base deck): 36/36 pass, 0.2 s
per call; 116/118 static structures (the 2 failures are inputs already flat in float cm).  The hook
does not call it yet: that is 6d.  Lesson: a single size field from the oxidation-start grid
stripped the film points of the growing oxide (it measured the oxide against the silicon it
replaced; SUPREM's box dose then moved 1 % per call); per-material fields fixed it.

**Finding for 6e.**  `ChooseKillNodes` (`src/dbase/grid_upd.c`) is an empty stub in the current
tree: production SUPREM removes **no** nodes during oxidation (the old McCoy passes `CknMcLen` /
`CknMcObt` are dead code; `rem_edg` / `mesh_patch_remove_point` are reachable only from other
callers).  So 6e is not "add interface points to an existing removal pass" but "design a removal pass
for oxidation" (candidate selection by approach rate, interface-aware edge collapse validated with the
6a predicates, local dose conservation).  That is more than the 2 days in the plan and touches the
fragile database; decide the scope before starting it.

## 4. Known issues and loose ends
- `core_evolve_test merge` produced rounding crumbs with the 1e-6 um grid until the wedge-tip fix in
  `significant()`; keep an eye on event counts in the tests (one gap + one cavity is correct).
- Spacer-stage structures contain 0.2 A oxide slivers from SUPREM itself; the remesher meshes them
  (Ruppert floor 1e-5 um), but they should be removed by events before meshing.
- Solution-adaptive refinement adds ~40 % points; SUPREM runs ~2.6x slower on a remeshed structure with
  20 % more points: profile before WP7.
- AdvCmos geometric re-ox run: PMOS punches through (poly 0.355 um after re-ox); rerun with the gate
  etched to ~0.445 um for a 0.40 um final gate.
- Experiment outputs of this session were in a temporary scratch directory and are gone; the numbers
  above and in `core/README.md` are the record.
