# SUPREM geometry/mesh core (M0)

Plan: `../doc/plan_m0_geometry_core.md`.  C++17, standalone of the SUPREM C code;
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
   (point, material) and a GAS node on every exposed point.  Node values are
   zero: field transfer is WP4.

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

- 72/72 pass; total points 977k -> 1.03M (+5 %); 0.2-18 s per structure.
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

## Tools

- `core_str IN.str [--write OUT.str]`: summary, invariants, B-rep, round trip.
- `core_remesh IN.str OUT.str [--size area|min] [--scale S] [--hmin H] [--hmax H]
  [--hgrad G] [--coarsen-tol T] [--film-threshold T] [--film-layers N] [--passes N]
  [--encroach-rounds N] [--no-mmg] [--no-flips] [--worst N | --worst-max N]`.
- `tests/roundtrip.py CORE_STR DIR...`, `tests/remesh_suite.py CORE_REMESH DIR... [--out D]`.
