# Oxidation retry and thin-polygon triangulation

The small-perturbation sweep exposed two independent failures in the full
PMOS process. These fixes build on `oxidation_mesh_fix.md`.

## Rejected time steps retained tentative edge flips

At surface spacing 0.000594, the first gate-oxidation step was rejected by
the diffusion solver. `timed_trbdf()` moved coordinates backward but left
the edge flips performed by `geom()` in place. One new diagonal lay outside
its cavity at the restored coordinates. The next attempt aborted in
`tri_geom()` before it could solve the smaller time step.

The diagnostic captured triangle 7058 changing from nodes
`(1612,2335,451)` to `(2335,451,449)`. After rollback, its doubled signed
area was `-3.0452528519686606e-17` cm², even though its point coordinates
had returned to their original values. This establishes connectivity
rollback as the cause, rather than simply an excessively large step.

Each moving-mesh diffusion trial now saves point state, triangles, edges,
and ordered node/edge connectivity after initial node reordering. A failed
TR or BDF stage restores that state, invalidates matrix assembly through
`UpdateSymbolic`, and restores the simulation clock. Successful trials
retain their new mesh. Snapshots are freed on both paths. This adds linear
memory and copying cost per moving-mesh trial. It assumes trials do not
add, remove, or renumber mesh entities after the checkpoint; a count
change remains an error.

`tests/mesh_trial.c` exercises the production edge flipper and geometry
code on a moving concave cavity. Coordinate-only rollback produces an
inverted triangle; the snapshot restores valid geometry and reciprocal
connectivity. Repeated rollback, growing connectivity buffers, and
successful-step commit are covered under ASan and UBSan.

## Thin valid polygons were rejected during contact processing

With correct rollback, both spacing perturbations reached later process
steps but exposed `error in chop` during triangulation. A captured
five-sided contact-etch polygon was simple and had positive area, but its
cached angles included values rounded to π. The existing heuristic
excluded the vertices at which valid ears could be removed.

When the existing mandatory-cut heuristic finds no candidate, `chop()`
now checks the actual stored coordinates in double precision. It rejects
zero-length edges and intersecting boundaries, then chooses a positive-area
ear containing no other boundary vertex. It preserves every boundary
vertex and retains the existing quality heuristic on its normal path.

`tests/ear_clip.c` includes the captured polygon and cached angles, all
boundary starting positions, a concave polygon, and crossed, flat, and
duplicate-point boundaries. It checks positive ear area and conservation
of polygon area. The captured fixture fails with the previous splitter
and passes with the fallback.

## Repeating validation

```sh
make -C suprem4/src suprem4 COPT='-O1 -fcommon -std=gnu89'
sh suprem4/tests/run_mesh_trial.sh
sh suprem4/tests/run_ear_clip.sh
sh suprem4/tests/run_quad_split.sh
sh suprem4/tests/run_mesh_memory.sh
python3 suprem4/tests/run_oxidation_perturbations.py
python3 suprem4/tools/run_pmos_mesh_regression.py \
  --binary "$PWD/suprem4/src/suprem4" --spacings 0.0005 0.0006 --keep
```

Saved-mesh integrity and process completion do not establish physical
accuracy. The earlier oxide-area discrepancy against stored reference
structures remains a separate acceptance issue; reference tolerances
have not been relaxed. Installed binaries and the viewer launcher are
unchanged.

## Final results

- All four perturbations complete: surface spacing 0.000594 and 0.000606,
  and gate-oxide opening shifts of ±0.1 micrometers at spacing 0.0006.
- Both original spacings, 0.0005 and 0.0006, also complete.
- All 54 checkpoint meshes pass positive-area, finite-coordinate,
  duplicate-triangle, edge-ownership, and reciprocal-neighbor checks.
  All saved solution values are finite.
- The rollback, ear-clipping, quadrilateral-splitting, and mesh-memory
  sanitizer tests pass.
- Across the four perturbations, final silicon area differs from the
  current 0.0006 baseline by less than 0.003%, and oxide area by less
  than 0.14%. These are geometry sensitivity measurements, not proof of
  physical accuracy.
- The stored-reference comparison still fails its 2% oxide-area tolerance:
  final differences are approximately 3.1% at spacing 0.0005 and 2.6% at
  spacing 0.0006. No reference files or tolerances were changed.

Final retained results:

- `/tmp/suprem4_oxidation_perturbations_flgx69f2/results.json`
- `/tmp/suprem4_pmos_regression_6efstlbh` (spacing 0.0005)
- `/tmp/suprem4_pmos_regression_jr_18cgf` (spacing 0.0006)

Debugging artifacts are under `/tmp/suprem4_perturb_debug`; the captured
thin polygon is `chop/polygon.json` there. The focused source tests embed
the reproducing geometry and do not depend on these temporary files.
