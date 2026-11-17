# Retained-corner fix for malformed etch boundaries

The doubled-density PMOS run (19,809 initial points, versus 9,815 at the
baseline) aborted with `Region has malformed boundary` while preparing
the gate-oxide region after `4000_ox.str`. The same failure reproduces
from that checkpoint without repeating oxidation.

## Cause and correction

`sub_2dskel()` classified a material-boundary vertex as OUT: it belonged
to material that should survive the etch. It then passed an intersecting
edge to `sp_edge()`, whose mesh-quality policy could move an endpoint
to a nearby intersection instead of inserting a point.

In the captured failure, that policy moved retained oxide node 9256
about 0.069 micrometers to the first cut. This changed the adjacent
edge before its intersection was computed. The next cut selected the
same node, so `first_node == last_node == 9256`. Closing the component
created a one-edge loop. `ck_clock()` correctly rejected it.

The correction is in clipping and edge splitting. The two 2D clipping
calls now identify their retained endpoint to `sp_edge_keep()`. When
snapping would move that endpoint, the routine inserts the intersection
through the existing material-aware edge-splitting code instead. An
intersection exactly at the retained endpoint still reuses it. This
also applies below the old minimum-grid-spacing threshold. Other
callers retain the existing snapping policy through `sp_edge()`.

Thus the original corner and both cut locations survive. Shared material
interfaces use the same new physical point, and nodal fields are
interpolated by the existing insertion routine. No region, edge, or
triangle is discarded to bypass an error, and boundary validation is
unchanged. Some meshes can contain additional small edges because a
minimum-spacing preference can no longer erase a retained corner.

## Regression tests

```sh
make -C suprem4/src suprem4 COPT='-O1 -fcommon -std=gnu89'
sh suprem4/tests/run_etch_corner.sh
python3 suprem4/tests/run_mesh_density.py
python3 suprem4/tests/run_oxidation_perturbations.py
python3 suprem4/tools/run_pmos_mesh_regression.py \
  --binary "$PWD/suprem4/src/suprem4" --spacings 0.0005 0.0006 --keep
```

The focused C regression uses local edges captured from the failure,
with a straight cut along the first etch segment's supporting line.
It demonstrates the legacy corner displacement and checks the corrected
operation with both edge orientations, two materials sharing the boundary,
field interpolation, retained area, and real skeleton validation. It
also checks exact endpoint reuse and cuts below the minimum spacing.
Point allocation, edge insertion, and topology routines are not mocked.
AddressSanitizer and UndefinedBehaviorSanitizer pass.

An additional analytic triangle/mask case exercises complete
`sub_skel()` subtraction, verifies retained area, and keeps the original
corner in place. With the previous clipping calls it aborts with
`Region has malformed boundary`; with the fix it passes sanitizers.

## Full-process validation results

| Case | Initial points | Final points | Result |
| --- | ---: | ---: | --- |
| Half density | 4,982 | 3,807 | Complete; 10/10 meshes valid |
| Double density | 19,809 | 10,240 | Complete; 10/10 meshes valid |

Both original surface spacings (0.0005 and 0.0006) and all four earlier
spacing/gate-position perturbations also complete. All 74 saved meshes
pass finite-coordinate/solution, positive-area, duplicate-triangle,
edge-ownership, and reciprocal-neighbor checks. The fast restart's five
subsequent checkpoints pass as well. All five focused mesh sanitizer
suites pass.

Relative to the current 0.0006 baseline, final oxide area changes by
+0.48% at half density and -0.76% at double density; silicon area changes
by less than 0.03%. These quantify mesh sensitivity, not physical accuracy.
The original-spacing comparisons still fail the stored-reference 2%
oxide-area tolerance, by approximately 3.1% and 2.6% respectively.

Retained full-run results:

- `/tmp/suprem4_mesh_density_v1flj8vk/results.json`
- `/tmp/suprem4_oxidation_perturbations__euozzqj/results.json`
- `/tmp/suprem4_pmos_regression_f4sib8n2`

Debugging artifacts and the fast restart deck are under
`/tmp/suprem4_boundary_debug`. The original failing run is under
`/tmp/suprem4_mesh_density_tm9aomqn/double_points`.

Mesh integrity and completion do not establish physical accuracy. The
stored-reference oxide-area discrepancy remains a separate acceptance
issue. Installed binaries and launchers are unchanged.
