# SUPREM geometry/mesh core (M0)

Plan: `../doc/plan_m0_geometry_core.md`.  **Handoff / current state: `../doc/m0_handoff.md`.**  C++17, standalone of the SUPREM C code;
it reads and writes SUPREM `.str` structure files.

```
cmake -S suprem4/core -B suprem4/build/core -DCMAKE_BUILD_TYPE=Release
cmake --build suprem4/build/core -j8
```

Third party: Mmg v5.8.0 (2D library only, LGPL-3, unmodified sources in
`third_party/mmg`), CDT 1.4.5 (MPL-2.0, header-only), Clipper2 (Boost, from
`../src/third_party/clipper2`).

## Status

| WP | Content | Validation |
|---|---|---|
| 0 | CMake project, vendored libraries | builds without warnings |
| 1 | `str_file` exact reader/writer, `mesh` semantic model + invariants | `tests/roundtrip.py`: every SubCmos2014 checkpoint byte-identical |
| 2 | `brep` loops/chains, closure/simplicity/area/twin checks | 36/36 + negative tests (broken repair/merge meshes flagged) |
| 3 | `remesh` mesh generation from the B-rep | `tests/remesh_suite.py`: 72/72 SubCmos2014 (two runs) + 36/36 AdvCmos2014 checkpoints |
| 4 | `transfer` field transfer + inventory restoration; solution-adaptive size field | suite (transfer checks) + SUPREM reload/anneal and mesh-convergence test |
| 5 | `evolve` motion, events (thin / crumb / gap / cavity), Clipper2 repair, rebuild from polygons | `core_evolve_test`: shrink, merge, slit pass |
| 6a | `suprem_check` SUPREM validity contract (float orientation, nxtel topology, embedding, region loops, mtest1); every core output is gated by it | `tests/suprem_contract.py`: core verdict = SUPREM's own load verdict on 1429/1429 structures + 5 mutations |
| 6b | `update_mesh` / `update_fields`: Lagrangian update (keep points and values, local removal/insertion, CDT, inventory restored on the changed points only) | `tests/update_suite.py`: 36/36 oxidation states of the base poly deck, 116/118 static structures (2 = inputs already flat in float cm) |
| 6c | `local_event_repair` / `apply_events_locally` (`src/local_repair.cpp`): events applied to the mesh locally (relabel / remove / fill, gap mouths), sub-resolution needles, contacts (zipper), flat triangles, pinches | `update_suite.py --events 1e-3`: every dumped oxidation state of the poly decks passes (3000+) |
| 6d | hook policy (`src/capi.cpp`): update on every call, local repair, regeneration only as a fallback; old triangles kept; oxide sizes >= `grid.oxide` | **all six poly-line decks complete** (0 regenerations, 4-11 min each); LOCOS field oxidation (viscous and compress) completes and matches the core-off run |

## WP6a SUPREM validity contract (`src/suprem_check.cpp`, `tools/core_check_suprem`)

Exact ports of what SUPREM does to a structure it reloads (`ig2_read` +
`create_db`) and to every trial motion, on the float cm coordinates SUPREM
holds (`(float)(x_um*1e-4)`):

| Rule | SUPREM source | On reload |
|---|---|---|
| orientation: no flat triangle (exact `orient2d_float_sign`) | geom/misc.c, geometry.c | panic "triangle is flat" |
| (strict) ordering: counter-clockwise as written | clktri.c `sclk_tri` | swapped silently |
| topology: a point-edge in at most two triangles | element.c `nxtel` | panic in `tnabor` |
| embedding: no overlap, no outer-boundary contact | include/mesh_embedding.h | every trial motion rejected |
| loops: each region border simple, turn sum +-2 pi, region not empty | refine/skel.c `skel_reg`, triang.c `ck_clock` | panic / segfault |
| (strict) components: one loop per region (no holes, one component) | `skel_reg` keeps only the first loop | loads, skeleton wrong |
| nodes: mtest1 (gas node at backside points, none at unexposed points, one node per material and point, every material node in a triangle) | dbase/check.c | panic "Mesh failed self test" |
| (strict) gas_added: gas node at every exposed point | make_db.c `create_db` | added on load |

`core_check_suprem A.str ...` prints OK / WARN (reloads, a strict rule fails) /
FAIL (SUPREM panics, crashes or rejects the embedding).  The cross-validation
`tests/suprem_contract.py CORE_BUILD SUPREM_BUILD [--mutate] PATH...` loads
every file in SUPREM too (`SUPREM4_LOAD_CHECK=1` makes `structure infile=`
print the float-orientation count and the embedding verdict) and requires the
same verdict.  Result over every .str under `decks/` (1429 structures,
including the saved panic states of earlier work): 1382 OK/loads, 10
WARN/loads (two-component oxide regions, rit_pmos_remesh), 36 FAIL/panic or
crash, 1 FAIL/embedding rejected; no disagreement.  Three further files are debug dumps with a gap in their
triangle numbering, which SUPREM renumbers silently and the exact WP1 reader
refuses ("triangles out of order"); the test lists them separately.  The mutation set
(inverted order, missing exposed gas node, interior gas node, empty region,
coincident points) agrees as well.  `core_remesh` reports the contract as
the `suprem` section; the WP6 hook (`capi.cpp`) refuses to write a mesh that
fails it.

Remesh suite with the contract (`m0_reference_20261122` + AdvCmos
`polyreox_geo_20261122`, 72 structures): 68/72.  New finding: NMOS `spacers`
and `source_drain` outputs contain triangles that are **flat in float cm**
(the 0.2 A oxide slivers near x = 9.96, y = -0.24 um), i.e. SUPREM would panic
reloading them although the WP3 checks passed.  PMOS `spacers` and
`source_drain` fail the coupling threshold (worst -9.2 / -12.4); identical
with the WP5 build, so this predates WP6a (the earlier 72/72 was on other
SubCmos runs).  Both point to the same open item: remove sub-delta slivers by
events before meshing.

## WP6b Lagrangian update (`update_mesh` in `src/remesh.cpp`, `update_fields` in `src/transfer.cpp`, `tools/core_update`)

Remesh mode "update": instead of regenerating, the mesh SUPREM moved is kept
and changed only where its spacing no longer matches the reference (the grid
at the start of the oxidation, as in the hook).  The reference gives three
size fields **per material** (finest altitude, area-equivalent edge, longest
edge; a material is measured against its own reference triangles, so growing
oxide is compared with the reference oxide, not with the silicon it
replaced, and SUPREM's anisotropic elements are their own target):

1. chain vertex (not a junction) removed when a segment is shorter than
   1/3 finest and the chord stays within `tol` (1e-4 um); chains that would
   intersect or collapse are restored (`restore_invalid_chains`, shared
   with the WP3 coarsening);
2. chain segments longer than 1.5 x longest subdivided (exact geometry);
3. interior points closer than 1/3 finest to a chain or to a kept interior
   point removed (the encroaching interface of a consumed material);
4. CDT of kept + inserted points with the chains as constraints; triangles
   larger than 2.5 x the reference area get a point (longest-edge midpoint or
   centroid), up to 8 rounds;
5. encroachment repair only where the update introduced it (an edge whose
   apex already faced it keeps the old quality), no split creating a segment
   step 1 would remove, at most 5 % x points splits.

Fields (`update_fields`): nodes at kept points copy their values exactly,
inserted nodes are interpolated (WP4), and the conserved P1 inventory per
material and species is restored by one factor applied only to the nodes at
**changed points** (inserted points and vertices of triangles the old mesh did
not have).  Every other node keeps its old values bit for bit; `core_update`
checks this (`retained`) together with invariants, B-rep, geometry, gas nodes,
the WP6a contract and the inventory (1e-10).

Results (`tests/update_suite.py`):
- 36 oxidation states of the base poly deck (every mesh the hook received
  with `SUPREM4_CORE_DUMP`, reference = the first): 36/36 pass, 0.2 s each
  (7-11k points); obtuse boundary edges 0-4, SUPREM's box-volume dose
  changes ~1e-4 per call in oxide and silicon (with a single global size
  field it was 1 %, the growing oxide lost its film points);
- static structures (each its own reference; SubCmos m0_reference, AdvCmos
  polyreox_geo and the WP3 outputs of both): 116/118 pass; on core meshes the
  update changes nothing or a few hundred triangles of 10^5; the 2 failures
  are WP3 outputs that are already flat in float cm (the NMOS 0.2 A sliver,
  an event for 6c).  Large meshes (120k points) take 30 s, mostly the
  repeated CDT and flip passes; the hook meshes are far smaller.

Also fixed on the way: `SegGrid` overflowed its int cell count for zero-length
segments (cell 1e-9 um), degrading to one cell and O(n^2) queries.

## WP6c/6d hook policy and local repair

Every call of `suprem_core_oxidation_check` (one per SUPREM oxidation step):

1. `apply_events_locally` (repeated up to 3 rounds), on the mesh itself:
   - events of `repair_shapes` (delta): thin / crumb pieces relabel the
     triangles whose centroid lies in the piece (or remove them, to the
     ambient); gap / cavity pieces are filled from a local CDT of the
     boundary points around them, after splitting the faces at the gap mouths
     (the only points 6c creates);
   - below the shapes' resolution (Clipper grid 1e-6 um), on the mesh:
     needles (two or more edges on the material's border, altitude
     < delta/8: to the material across, monotone so it terminates),
     contacts (an exposed vertex within delta/2 of a facing exposed edge,
     also of its own folded surface: zipped into that edge; ambient wedges
     with ends closer than delta/2 closed), flat triangles (altitude
     < delta/8: removed, the triangle across split at the apex);
   - then pinches (a region touching itself at a vertex) and detached
     components smaller than delta^2 are resolved.
   No point moves; everything outside the pieces is unchanged.  Pieces
   smaller than the triangles they lie in are left until they become a
   needle, a flat triangle or a contact.
2. `update_mesh` on the result (6b) with the grid at the start of the
   oxidation as reference; old triangles are kept unless something changed
   there or one of their edges has a coupling below -0.01 (their edges are
   CDT constraints and are not flipped).  SUPREM's force request (predicted
   squeeze / collision) means a full Delaunay triangulation.
3. `update_fields` from the original mesh, restoring the inventories in
   SUPREM's box measure on the changed points.
4. Float repair, invariants, gas nodes, the WP6a contract; nothing is
   written when nothing changed.
5. Fallback: full regeneration (WP3) only when the B-rep is invalid or the
   local repair / update fails (`SUPREM4_CORE_REGENERATE=1`: old policy).

On the SUPREM side the global delta/4 speed cap is off by default
(`SUPREM4_CORE_SPEED_LIMIT=1`): it throttled LOCOS to ~1 s steps; the
collision prediction (relative velocities, pairs closer than delta/2 left
to the core) and the area-halving limit remain.

Findings that shaped this:
- re-triangulating the whole mesh every call (CDT of SUPREM's moved,
  non-Delaunay mesh: ~25 % new triangles) made the viscous stress solve
  fail after ~30 calls; keeping the old triangles fixed it;
- every poly-deck failure after 6d was at the poly corners (gate-oxide /
  sidewall wedge, poly top cusp): full regeneration there left couplings of
  -4 to -800 and meshes SUPREM could not move; the local rules above
  replaced every one of them (no regeneration in the runs since);
- SUPREM's embedding check rejects every motion once a surface vertex sits
  within float resolution of a facing edge (a closed gap); the zipper.

Results (2026-11-23, `experiments/wp6d_hook_20261123/run8`, `locos*/core_g6`):

| Deck | Reached | Wall | Hook calls | Local repairs | Regenerations |
|---|---|---|---|---|---|
| poly base (51 min dry 1000 C) | 3060 s, complete | 323 s | 68 | 29 | 0 |
| coarse | complete | 416 s | 70 | 40 | 0 |
| dry1000_120 (120 min) | 7200 s, complete | 679 s | 128 | 56 | 0 |
| narrow025 | complete | 260 s | 61 | 21 | 0 |
| wet900_20 | 1200 s, complete | 554 s | 101 | 44 | 0 |
| wide080 | complete | 365 s | 55 | 23 | 0 |

Every final structure passes the WP6a contract; the field oxide of the base
deck is 45.4 nm against 46.1 nm of plain SUPREM on a flat wafer.  Without
the core these decks stall at 47 s; with the 6a/6b-era hook they failed at
30-500 s.

LOCOS (SubCmos2014 NMOS field oxidation, 65 min wet 1100 C + ramps):
viscous: field oxide 564.1 nm vs 565.5 nm core off, oxide area -0.005 %,
628 s vs 705 s; compress: 615.2 vs 617.1 nm, area -0.12 %, 909 s vs 549 s.
Boron moves ~2 % from the oxide to the silicon (the interface
discretisation; the hook itself conserves it to 1e-15 per call).

Speed: SUPREM recomputes the oxide flow (7-node triox FE, direct solve)
periodically; its cost grows steeply with the oxide mesh.  The update kept
the oxide at the reference film's spacing (4.5 nm from a 10 nm film in 2
layers), 5x more elements than SUPREM's own `grid.oxide` (10 nm) would
add, and one flow solve took 180-280 s.  `grid.oxide` is now passed to the
core (`suprem_core_set_param("oxide_grid_um")`) as the lower bound of the
oxide target sizes: the flow solve takes ~50 s on the late base mesh and the
decks run 10-20x faster.

## WP3 remesh pipeline (`src/remesh.cpp`)

1. **Canonical chains.** Each interface chain is stored once and used by both
   regions (reversed), so the new mesh is conforming by construction.
2. **Chain discretisation.**  Vertices are removed only if the replacing chord
   stays within `coarsen_tol` (default 1e-5 um = 0.1 A) of them and is not
   longer than the target size; vertices shared by several chains stay;
   a simplification that would create an intersection is undone.  Segments
   longer than the target size are subdivided on the straight segment.
   `--coarsen-tol 0` keeps every vertex (exact geometry).
3. **Size field** (`SizeField`): area-weighted mean of the old triangles'
   area-equivalent sizes (degenerate triangles carry no weight), scaled and
   clamped.  (A tensor average was tried first: old sliver triangles dominated
   it and refined boundaries 50x.)
4. **Thin films** (thickness < `film_threshold`, 0.04 um): one projection pass
   along inward normals over a snapshot of the chain vertices; the hit is
   matched to an end vertex within a tolerance of min(0.2 t, 0.3 x spacing on
   both sides) or a vertex is inserted.  Points inserted by the pass are never
   projected again (an iterated version cascaded to 0.1 A segments in film
   wedges).  Pairs carry `film_layers - 1` interior points (layers at least
   `film_min_layer` thick), giving aligned near-rectangular elements, which
   are Delaunay with right angles at any aspect ratio.
5. **CDT** of chains + film points, regions labelled by flood fill from the
   region loops (a conflict is an error).
6. **Mmg2d** (2 passes) with all region-boundary edges/vertices and the film
   triangles required; isotropic metric from the size field.  Required edges
   must survive exactly (checked); coordinates of kept vertices are restored
   bit-exactly (Mmg rescales internally).
7. **Delaunay flips** inside regions, then **local Ruppert splitting** of
   region-boundary edges whose opposite angle exceeds 91 deg (midpoint inserted
   in the mesh and the chain, flips restore Delaunay), until none is left.
8. Assembly: neighbours, boundary codes from the chains, one node per
   (point, material) and a GAS node on every exposed and backside point
   (SUPREM's `mtest1` rule `is_surf(p) || is_back(p)`; missing backside gas
   nodes made SUPREM panic on load).
9. Size field, solution part (WP4): where a conserved dopant varies by more
   than `solution_decades` (0.25) across an old element, the target size
   shrinks in proportion, but not below the old mesh's finest local size.

## WP4 field transfer (`src/transfer.cpp`)

- Per material P1 interpolation from the old triangles of that material
  (interface values stay per material, segregation jumps are kept); points
  with identical coordinates copy the old node exactly; gas nodes copy the
  nearest old gas node (ambient values).  Velocities, defects, oxidant,
  stresses and potential are interpolated.
- Conserved species (As, P, Sb, B, Au, Cs, ... totals) are restored per
  material to the old **P1 integral** by the minimum change weighted by
  each node's interpolation uncertainty (spread of the old values it came
  from; exact copies unchanged, values kept non-negative); the active
  concentration follows its total.  A uniform scale (`--restore-scale`) is
  the alternative; the weighted one gives lower round-trip errors.
- Why P1 and not SUPREM's mixed-Voronoi box volumes: on SUPREM's sliver /
  obtuse triangles the box sum is a poor quadrature (pad oxide in `wells.str`:
  box inventory 25 % below the field integral, which the remesh preserves
  to 1e-4); restoring the box value would remove real dopant.  The box-volume
  change is reported.
- Identity transfer (same mesh) is bit-exact (checked by `core_remesh`).

Results (SubCmos2014 reference run, 36 structures): all pass; inventory error
after interpolation median 1.3 %, max 14.5 % (P at an exposed oxide surface,
2.7 decades across one old element); after restoration <= 2.2e-14; transfer
<= 0.6 s.  AdvCmos2014 (36 structures): all pass; 1.55M -> 2.15M points;
inventory error after interpolation median 1.8 %, max 20 %; after
restoration <= 1.7e-14; transfer <= 1.8 s.

SUPREM test (`pmos/gate.str`, then 950 C 20 min N2): SUPREM loads the
remeshed structure and anneals it (invariants ok).  Mesh convergence after
the anneal (original 30k points; remesh at scale 1 / 0.7 / 0.5):

| | original | remesh 36k | 60k | 103k |
|---|---|---|---|---|
| B in Si | 1.856e17 | 1.9259e17 | 1.9255e17 | 1.9248e17 |
| B in oxide | 9.88e16 | 9.31e16 | 9.35e16 | 9.36e16 |
| P in oxide | 6.7e16 | 2.26e17 | 2.91e17 | 3.11e17 |
| P in Si | 8.8546e17 | 8.8547e17 | 8.8547e17 | 8.8547e17 |

The remeshed results converge (B in Si within 0.06 %) and the original
SUPREM mesh is the outlier (-3.6 %).  P loss through the exposed oxide
surface is still mesh dependent at 103k points: exposed surfaces with
steep profiles need dedicated refinement (WP7 item).  Runtime: 102 s
on the 36k-point remesh vs 39 s on the original 30k mesh (to be looked at
in WP6/WP7).

## Checks (`core_remesh` prints all of them)

- invariants (orientation, neighbour symmetry, nodes per material, no duplicates)
- B-rep of the output (closure, simplicity, area, twin chains)
- geometry: two-sided Hausdorff distance old/new boundary <= `coarsen_tol`,
  region areas within `coarsen_tol` x perimeter
- gas nodes exactly on exposed points
- quality: angles, and box-method couplings (cot sum / 2 per edge and region).
  The suite requires no coupling below -0.01 (a boundary edge opposite angle
  of about 91 deg).

## Results (SubCmos2014 0.5 um, all 16 checkpoints x NMOS/PMOS, two runs)

- 72/72 pass; total points 977k -> 1.03M (+5 %) with the geometric size
  field alone, 1.34M (+38 %) with the WP4 solution adaptation; 0.2-80 s per structure.
- Worst coupling: input down to -2.5e5 (degenerate SUPREM triangles), output
  >= -0.0087 everywhere.  Max angle 180 -> <= 165 deg (interior obtuse
  triangles are Delaunay, hence harmless for the box method).
- Gate stages: 30k -> 15k points with max angle 131 deg.
- Spacer/source-drain stages: points double (aligned film layering takes the
  union of both film sides' vertex spacing); these inputs also contain oxide
  slivers about 0.2 A thick (e.g. `spacers.str` region 3), which cannot be
  meshed well while preserving geometry; removing them is WP5 event handling
  (thickness < delta).

AdvCmos2014 0.25 um (geometric poly re-ox run, 18 checkpoints x 2): 36/36 pass,
1.55M -> 1.70M points (+10 %), worst coupling input -9.9e3 -> output >= -0.0087;
up to 85 s for the 120k-point metal structures (Mmg dominates; a runtime
target belongs to WP7).

## WP5 motion, events, repair (`src/evolve.cpp`)

- `evolve(mesh, velocity(p, t), time)`: sub-steps with travel <= min(max_step,
  delta/4) (a gap closing from both sides is seen before the fronts cross), no
  triangle losing more than half its area per step, boundaries kept simple
  (checked; the step is halved otherwise).
- After each step the geometry (per material, Clipper2 integer grid 1e-7 um)
  is checked for events of size delta (1 nm):
  - thin: part of a region narrower than delta (opening) -> neighbour with
    the longest shared boundary (removed if that is the ambient);
  - crumb: connected component narrower than delta/8 everywhere -> as thin;
  - gap: ambient gap narrower than delta -> highest-priority adjacent
    material (SiO2 > OxNi > Si3N4 > poly > Si ...);
  - cavity: ambient pocket enclosed by material -> as gap.
  Each piece of space changes owner at most once per repair (no ping-pong
  of slivers between two materials).
- Rebuild (on events, or when a triangle has lost 70 % of its area since the
  last rebuild): planar graph of all polygon edges (T-junctions split exactly,
  crossings resolved by the CDT), point-in-polygon labels, regions = connected
  components per material (touching same-material regions merge), then the
  WP3 remesh with the target sizes of the mesh *as built at the last rebuild*
  (sizes from the deformed mesh ratchet: compressed material got finer at
  every rebuild, 2.8k -> 8k points before this was fixed), then the WP4
  transfer.

Synthetic tests (`core_evolve_test REFERENCE.str [OUTDIR]`, every step checked):

| test | what happens | events | result |
|---|---|---|---|
| shrink | poly square in oxide consumed from all sides | thin poly -> oxide | one oxide region, area 1, 3.4 s |
| merge | silicon slot closing from both sides, top first | gap + trapped cavity -> Si | one Si region, no ambient below the surface, 1.7 s |
| slit | slit through oxide on Si closing | gap + cavity -> oxide | one oxide + one Si region, 1.1 s |

Fixes in WP2/WP3 found by these tests: closed interface chains (an island
inside one region) compared independent of start point; chains split at
pinch vertices (boundary-graph degree > 2) so both sides split identically;
chain coarsening also rejected when it duplicates a segment or collapses /
flips a loop (regions thinner than the tolerance); Ruppert splitting stops at
edges of 2 x 1e-5 um (it does not terminate near small input angles). The
floor is also set by SUPREM, which stores coordinates as float cm (~1e-7 um
resolution): a 1e-9 um floor produced triangles that were flat after
rounding.  Sub-0.1 A slivers (spacer structures) keep some negative
couplings; on the oxidation path they are removed as thin/crumb events.

## Tools

- `core_str IN.str [--write OUT.str]`: summary, invariants, B-rep, round trip.
- `core_remesh IN.str OUT.str [--size area|min] [--scale S] [--hmin H] [--hmax H]
  [--hgrad G] [--coarsen-tol T] [--film-threshold T] [--film-layers N] [--passes N]
  [--encroach-rounds N] [--no-mmg] [--no-flips] [--worst N | --worst-max N]`.
- `core_update IN.str OUT.str [--reference REF.str] [--float] [--timing] [--near F] [--short F] [--tol T]
  [--long F] [--grow F] [--scale S] [--hmax H] [--solution-decades D] [--box] [--retriangulate]
  [--keep-coupling C] [--events DELTA]` (WP6b; `--events`: 6c local repair first, as the hook).
- `core_check_suprem A.str ...` (WP6a).
- `core_evolve_test REFERENCE.str [OUTDIR]` (env `CORE_EVOLVE_ONLY=name`, `CORE_EVOLVE_TRACE=1`).
- `tests/roundtrip.py CORE_STR DIR...`, `tests/remesh_suite.py CORE_REMESH DIR... [--out D]`,
  `tests/suprem_contract.py CORE_BUILD SUPREM_BUILD [--mutate] PATH...`,
  `tests/update_suite.py CORE_BUILD [--reference REF.str] [--events DELTA] PATH...`.
