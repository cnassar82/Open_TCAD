# Oxidation cavity triangulation fix

The full RIT PMOS deck at surface spacing 0.0006 crashed during oxidation.
The final backtrace was `ck_flip -> flip -> geom -> bd_connect -> grid_sub`:
an edge referred to a triangle that did not contain that edge.

Tracing replacement cavities upstream identified `quadsplit()` in
`src/refine/triheur.c`. Its quality score multiplied absolute distances to
each diagonal. For a concave quadrilateral, this could favor the exterior
diagonal and create overlapping or inverted triangles. Subsequent node
removal then corrupted connectivity.

The splitter now requires the other two vertices to lie on opposite sides
of a candidate diagonal. It preserves the existing quality ranking when
both diagonals are valid. If neither is valid it returns NULL, using the
caller's existing failure path. No downstream crash suppression or new
triangle deletion was added.

## Validation

```sh
make -C suprem4/src suprem4 COPT='-O1 -fcommon -std=gnu89'
sh suprem4/tests/run_quad_split.sh
sh suprem4/tests/run_mesh_memory.sh
python3 -u suprem4/tools/run_pmos_mesh_regression.py \
  --binary "$PWD/suprem4/src/suprem4" --spacings 0.0005 0.0006 --keep
```

The focused test calls the production splitter with the observed oxidation
cavity and convex, concave, collinear-corner, and crossed quadrilaterals.
It checks every boundary start and both windings, requiring interior
diagonals and area conservation. The old splitter fails this test; the
fix passes AddressSanitizer and UndefinedBehaviorSanitizer.

Both full PMOS simulations now complete and write all nine checkpoints.
Temporary cavity diagnostics also passed throughout both runs; they were
removed from production code. All 18 output meshes were checked for
positive triangle area, finite coordinates, valid references, duplicate
triangles, nonmanifold edges, and reciprocal face neighbors. All passed.
The final 0.0006 rerun without diagnostics produced byte-identical
checkpoints (`/tmp/suprem4_pmos_regression_2m6287z3`). The mesh memory
sanitizer suite also passes.

The saved-reference comparison still fails its 2% oxide-area tolerance
from `700_ox.str` onward: approximately 3.0% at spacing 0.0005 and 2.6%
at spacing 0.0006. The 0.0005 discrepancy existed before this fix; its
final oxide area changes by only about 0.035% compared with the preceding
same-source candidate. This fixes the reproduced crash, but does not
establish physics acceptance or explain the reference discrepancy.
The tolerance and stored references were not changed.

Investigation artifacts (temporary files, not required to run the tests):

- `/tmp/suprem4_oxidation_trace/gdb.log`: downstream crash trace.
- `/tmp/suprem4_pmos_regression_jt5cr1el`: first invalid cavity diagnostic.
- `/tmp/suprem4_pmos_regression_ftxf_ilw`: fixed 0.0006 run.
- `/tmp/suprem4_pmos_regression_xe4d0vqq`: fixed 0.0005 run.

The tested executable is `suprem4/src/suprem4`; the installed executable
and launcher have not been replaced. This is a targeted triangulation
fix, not a guarantee that every legacy mesh limit has been removed.

## Small-perturbation sweep

The failures recorded in this section are resolved by the subsequent
[retry and thin-polygon fixes](oxidation_retry_fix.md). All four
perturbations now complete and pass saved-mesh audits. The table below
records the findings before those fixes.

Run `python3 suprem4/tests/run_oxidation_perturbations.py` to repeat four
full-process cases around surface spacing 0.0006. The gate cases translate
only the gate-oxide opening, preserving its 40 micrometer width. The script
retains decks, logs, checkpoint structures, and a JSON mesh-audit summary.

| Change | Completion and saved-mesh audit |
| --- | --- |
| Surface spacing -1%, to 0.000594 | Abort (SIGABRT); first 4 checkpoints valid |
| Surface spacing +1%, to 0.000606 | Completes; all 9 checkpoints valid |
| Gate-oxide opening shifted -0.1 micrometers | Completes; all 9 checkpoints valid |
| Gate-oxide opening shifted +0.1 micrometers | Completes; all 9 checkpoints valid |

Results are in `/tmp/suprem4_oxidation_perturbations_j8600xzk/results.json`.
The failing case reaches the 700-Angstrom oxidation step and reports
`triangles are not clock wise, data base corrupted` before writing
`700_ox.str`. This is an additional reproducible robustness failure; its
upstream cause has not yet been established. Passing saved-mesh audits
does not establish validity at every internal solver step or physics
accuracy. The simulator is therefore still sensitive to small mesh changes.
