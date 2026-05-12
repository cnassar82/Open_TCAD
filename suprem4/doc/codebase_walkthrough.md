# SUPREM4 Codebase Walkthrough

This note records the current understanding of the reverted SUPREM-IV.GS
codebase. Its purpose is to give future mesh/refine work a factual reference
before changing behavior.

## Build And Runtime Entry Points

The top-level build is driven by `Makefile`, which delegates to
`src/Makefile`. `make depend` regenerates `src/Makefile` from
`src/Makefile.proto`, so durable build defaults belong in both the top-level
Makefile and the prototype.

The executable starts in `src/main.c`. Startup does the following:

1. Resolves help, key, model, and implant file locations from environment
   variables or compiled-in defaults.
2. Initializes the command parser from `suprem.uk`.
3. Initializes diffusion and, when enabled, device code.
4. Initializes plotting and rectangular mesh helpers.
5. Sources `modelrc`, user rc files, and any command-line input decks.
6. Enters the parser loop for interactive input if no deck was supplied.

Commands are registered in the `command` table in `src/main.c`. The process
geometry commands of interest are:

- `initialize` -> `initialize()`
- `structure` -> `structure()`
- `line`, `region`, `boundary` -> mesh construction input
- `deposit` -> `user_deposit()`
- `etch` -> `user_etch()`
- `diffuse`, impurity commands, and `oxide` -> diffusion/oxidation flow
- `implant` -> implant flow
- `regrid` -> explicit regrid command

The shell layer is in `src/shell`. `do_source()` and `do_string()` feed input
into the parser. `do_command()` maps parsed command names to command-table
entries. `do_exec()` validates command parameters with the `check` subsystem
before calling the command function.

## Major Source Directories

- `src/check`: command keyword/parameter checking. The generated key file is
  produced by `src/keyread`.
- `src/shell`: command input, parsing support, macros, source files, and shell
  fallback.
- `src/mesh`: mesh initialization, structure read/write, profile loading, and
  rectangular mesh setup.
- `src/dbase`: persistent mesh database, connectivity rebuild, geometry,
  neighbor links, edge/node/region management, consistency checks, and moving
  grid support.
- `src/refine`: process-geometry refinement: deposit, etch, skeleton regions,
  polygon subtraction, edge splitting, and skeleton triangulation.
- `src/diffuse`: impurity/defect diffusion setup and solve.
- `src/oxide`: oxidation, moving boundaries, oxidant model, viscous/stress
  finite element coupling.
- `src/implant`: implant profile and Pearson/Gaussian distribution handling.
- `src/device`: electrical device solve support.
- `src/finel`: finite element assembly/solve support used by oxidation/stress.
- `src/math`: sparse matrix and linear algebra routines.
- `src/plot`, `src/gpsup`, `src/xsupr4`: plotting backends and output helpers.
- `src/imagetool`: image/grid helper routines.
- `src/misc`: utility functions, panic handling, model/key file reading, CPU
  logging, and small portability helpers.

## Core Database Model

The central data structures are declared in `src/include/geom.h` and accessed
mostly through macros in `src/include/dbaccess.h`.

The main objects are:

- `pt`: a physical coordinate point. It stores coordinates, velocities, local
  spacing, flags, and a list of material-specific nodes at that location.
- `nd`: a material-specific node attached to one point. It stores solution
  values for impurities/defects, material id, lists of triangles and edges, and
  the owning point.
- `tri`: an element. In 2D this is a triangle. It stores node ids, neighbor
  element ids, edge ids, geometry coefficients, region id, tree refinement
  metadata, and flags.
- `edg`: an edge between two nodes. It stores connected elements, connected
  skeleton regions, length, coupling, and flags such as exposed surface.
- `reg`: a material region. It stores material id plus lists of triangles,
  edges, boundary faces, nodes, and an associated skeleton region.
- `sreg`: a temporary skeleton region, declared in `src/include/skel.h`, used
  as a polygonal boundary during deposit/etch/triangulation operations.

A physical interface is represented by multiple nodes at the same point, one
per material. For example, an oxide/silicon/gas surface point can have a
silicon node, an oxide node, and a gas node all sharing one coordinate. This
is fundamental to the code: changing point/node relationships casually can
break material interfaces.

Important flags include:

- Point flags: `SURFACE`, `BACKSID`, `SKELP`
- Triangle flags: `TRIPTS`, `GEOMDN`, `CLKWS`, `NEIGHB`
- Edge flags: `REGS`, `ESURF`, `EBACK`, `MARKED`
- Region flags: `EXPOS`, `SKEL`, `ETCHED`

The dirty flags in `dbaccess.h` tell the database rebuild logic which derived
connectivity is stale.

## Database Rebuild Path

`bd_connect()` in `src/dbase/make_db.c` is the main consistency rebuild. It is
called after major geometry edits such as deposit and etch.

The rebuild sequence is:

1. Free existing skeleton regions.
2. Remove dead points/nodes/elements/edges/regions if `need_waste` is set.
3. Ensure triangle orientation with `clock_tri()` if needed.
4. Rebuild node-to-triangle lists with `node_to_tri()`.
5. Rebuild triangle neighbors with `nxtel()` if needed.
6. Rebuild edges with `build_edg()`.
7. Recompute geometry and coupling data with `geom()`.
8. Rebuild each region with `build_reg()`.
9. Run `mtest1()` and `mtest2()`.

`build_reg()` reconstructs each region's triangle list, edge list, boundary
face list, node list, exposure flag, and skeleton boundary.

`mtest1()` and `mtest2()` in `src/dbase/check.c` verify database consistency:
allocated objects, point/node back-links, gas node vs surface consistency,
triangle materials, loose nodes, node/triangle/edge reverse links, neighbor
symmetry, and exposed point/edge consistency. These tests catch corruption,
but they do not prove the physical shape is correct.

## Deposit Flow

`user_deposit()` is in `src/refine/deposit.c`.

The current deposit path is:

1. Reject invalid current meshes via `InvalidMeshCheck()`.
2. Increment `process_step` and save grid state with `GridSave()`.
3. Parse material, thickness, divisions, optional spacing, optional file,
   optional square behavior, and optional dopant concentration.
4. Create a new material region with `mk_reg()`.
5. Extract the current exposed surface using `find_surf()`.
6. Generate an offset surface using `gen_offset()`, unless a file supplies it.
7. Build a skeleton polygon between the existing surface and the offset line
   with `build_skel()`.
8. Triangulate that skeleton into the new region with `grid()`.
9. Free the temporary skeleton.
10. Rebuild connectivity with `bd_connect("after deposit")`.
11. Initialize dopant values in new nodes when requested.
12. Run mesh self-tests.

Deposit does not directly deform old triangles. It builds a new strip-shaped
polygon between the old surface and the offset surface, triangulates that
polygon, then reconnects the database.

`find_surf()` in `src/refine/surface.c` walks exposed `ESURF` edges and returns
an ordered surface polyline. In 2D it starts from a leftmost exposed edge and
walks edge-to-edge through shared nodes. Correct `ESURF` flags and edge
connectivity are prerequisites.

`gen_offset()` in `src/refine/offset.c` computes a parallel curve at the
requested thickness. It adds extra points around corners using local spacing,
then attempts to remove loops. This is a high-risk stage for self-intersection
or tiny-edge creation.

## Etch Flow

`user_etch()` is in `src/refine/etch.c`.

Supported modes include:

- `dry`: finds the exposed surface and shifts it by the requested thickness.
- `left` / `right`: constructs a side-removal polygon.
- `start` / `continue` / `done`: constructs a user-defined polygon.
- `all`: removes all exposed regions of the selected material.
- `physical`: uses rates and time, implemented in `src/refine/rate.c`.
- `file`: etches from a string-defined polygon.

The main `etch()` function uses polygon subtraction:

1. Clear `ETCHED` flags.
2. For each exposed region matching the selected material, build a temporary
   skeleton from the etch polygon.
3. Take the target material region's current skeleton.
4. Call `sub_skel(target_skeleton, etch_skeleton, result_skeletons)`.
5. If the target is fully removed, remove the region.
6. If one or more pieces remain, create new regions, triangulate each piece
   with `grid()`, interpolate solution data from the old region, remove old
   triangles, and delete the old region.
7. Rebuild connectivity with `bd_connect()`.
8. Clear temporary skeleton flags.

Etch is therefore not a triangle-by-triangle deletion algorithm. It is a
boundary skeleton subtraction followed by retriangulation.

## Skeleton Subtraction

Skeleton logic is in `src/refine/skel.c`.

`skel_reg()` builds a skeleton boundary from an existing filled material
region. `sub_skel()` subtracts one skeleton from another.

The subtraction path:

1. Classify target boundary nodes as `IN`, `OUT`, or `ON` relative to the etch
   skeleton using `check_in()`.
2. Return early if the target is entirely inside or outside the etch polygon.
3. In 2D, call `sub_2dskel()`.
4. Walk boundary edges looking for transitions between inside/on and outside.
5. Find edge crossings with `edg_crs()`.
6. Split crossing edges with `sp_edge()`.
7. Follow surviving target boundary segments.
8. Add newly exposed segments along the etch skeleton.
9. Return one or more output skeletons.

This code depends on tight geometric classification. `on_bound()` uses a
tolerance around `1e-8` cm, and `pt_in_skel()` uses a ray-crossing
point-in-polygon test. Near-coincident edges, tiny edges, duplicate points, and
self-crossing polygons are likely failure sources.

## Skeleton Triangulation

`grid()` in `src/refine/grid.c` dispatches to:

- `triang()` for 2D
- `lineseg()` for 1D

`triang()` in `src/refine/triang.c` takes a skeleton polygon and produces
triangles for a material region.

The algorithm:

1. Check and normalize boundary orientation with `ck_clock()`.
2. Duplicate the skeleton so the working copy can be split destructively.
3. Prefer rectangular decomposition with `rect_div()` for long regions.
4. Optionally subdivide exposed boundary edges.
5. Repeatedly process the top skeleton region:
   - if it has 3 edges, create a triangle with `cr_tri()`;
   - if it has 4 edges, split as a quadrilateral;
   - otherwise try to chop off a good acute triangle;
   - otherwise divide the region;
   - otherwise chop the least bad triangle.

This is an ear-clipping/region-splitting style triangulator. It is not a
general constrained Delaunay remesher. It assumes a valid, simple, continuous,
counter-clockwise boundary.

## Edge Splitting

`sp_edge()` in `src/refine/sp_edge.c` performs low-level edge splitting. It
may move an existing endpoint when the split point is very near an endpoint,
or create a new point/node/edge when the split is interior. It updates skeleton
references and connected triangle/edge relationships.

Because many routines use existing edge identity and skeleton edge lists,
edge splitting is dangerous to modify without proving all reverse links remain
valid.

## Diffusion And Moving Grid

The diffusion command flow is primarily in `src/diffuse`.

Broadly:

- Impurity-specific files configure active concentration, boundary conditions,
  coupling, and model constants.
- `diffuse.c`, `setup.c`, `prepare.c`, `solve.c`, `solve_time.c`, and `time.c`
  manage assembly and time stepping.
- `Interst.c`, `Vacancy.c`, and `defect.c` manage point-defect species.
- `moving.c` and related code interact with moving boundaries.

The database includes `grid_upd.c` and `grid_loop.c` for moving-grid updates,
triangle inversion checks, loop handling, and connectivity rebuilds. This is
separate from the deposit/etch skeleton pipeline but shares the same mesh
database and invariants.

## Oxidation And Finite Elements

Oxidation code lives in `src/oxide`.

The important responsibilities are:

- Oxidation coefficients and model setup: `coeffox.c`, `oxrate.c`,
  `Oxidant.c`
- Oxide growth and velocity calculation: `oxgrow.c`, `oxide_vel.c`
- Material and boundary-condition handling: `mater.c`, `FEbc.c`,
  `FEconvert.c`
- Viscous/stress support: `viscous.c`, `elast.c`, `triox.c`, `vert.c`

Finite element support is in `src/finel`, with assembly, solve, and triangle
element routines. The FE code interacts with the same mesh database, so
triangle quality and boundary consistency matter for oxidation/stress as well
as pure geometry operations.

## Implant, Device, Plot, And Structure I/O

`src/implant` handles implant profiles and distributions. It uses the current
surface/material geometry to place dopants.

`src/device` handles electrical solves after process simulation. It relies on
the same impurity solution arrays stored on nodes.

`src/mesh/structure.c`, `ig2_meshio.c`, `pi_meshio.c`, and `save_simpl.c`
handle structure import/export. Structure output should be used as a
regression checkpoint for future mesh work.

`src/plot` and plotting backends consume mesh, material, and solution data for
visualization.

## Critical Invariants For Mesh Robustness Work

These invariants should be treated as contracts:

1. A physical point may have multiple material nodes; do not merge nodes just
   because coordinates are equal.
2. Gas nodes should exist only at exposed/backside points and should not belong
   to triangles.
3. Material nodes should belong to at least one triangle.
4. Triangle nodes must all match the triangle's region material.
5. Edge/node/triangle/region reverse links must be symmetric after edits.
6. Exposed `ESURF` edges must correspond to point `SURFACE` flags and gas
   nodes.
7. Region skeletons must be simple, continuous, and correctly oriented before
   triangulation.
8. Skeleton subtraction assumes stable `IN`/`OUT`/`ON` classification near
   boundaries.
9. `bd_connect()` is the authoritative rebuild after geometry edits.
10. `mtest1()` and `mtest2()` catch database corruption, but not necessarily
    wrong physical shape.

## Where Robustness Should Be Added First

Before changing algorithms, add diagnostics around the existing pipeline:

1. Dump the surface returned by `find_surf()`.
2. Dump the offset line from `gen_offset()`.
3. Dump deposit skeletons from `build_skel()`.
4. Dump etch polygon skeletons before subtraction.
5. Dump result skeletons from `sub_skel()`.
6. Record triangulation inputs and outputs from `triang()`.
7. Record region skeletons rebuilt by `bd_connect()`.
8. Add geometry metrics: minimum edge length, minimum area, minimum angle,
   aspect ratio, duplicate consecutive points, self-crossing skeletons, and
   material area.

The immediate goal should be to identify the first operation where a known
deck diverges geometrically from expected output. Only after that should we
change meshing behavior.

## Risk Areas

The highest-risk code for incorrect structures is:

- `find_surf()`: can start from or follow the wrong exposed boundary.
- `gen_offset()`: can create self-crossing or poorly mapped offset curves.
- `sub_skel()` / `sub_2dskel()`: can misclassify nodes or construct incorrect
  remaining polygons.
- `edg_crs()` and `on_bound()`: sensitive to tolerance and near-coincident
  geometry.
- `sp_edge()`: updates many coupled data structures.
- `triang()`: assumes the input skeleton is already valid.
- `bd_connect()` / `build_reg()`: rebuilds the persistent interpretation of
  regions and exposed surfaces.

Future robustness work should preserve behavior until a diagnostic proves the
exact failure stage.
