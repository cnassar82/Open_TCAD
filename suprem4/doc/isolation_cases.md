# LOCOS and STI reliability exercises

These decks exercise process completion and geometry. They are not calibrated
TSUPREM-4 process replacements. STI uses polygonal trench etching, viscous liner
oxidation, conformal oxide deposition, and a geometric planar cut at y=0;
there is no CMP or deposition-void calibration and no STI hard-mask sequence.

## Decks

- `decks/locos/locos.in`: based on bundled example 9, a 20 nm pad oxide,
  150 nm nitride mask, and 90 minutes of stress-dependent viscous wet oxidation
  at 1000 C, followed by nitride stripping.
- `decks/sti/sti.in`: 600 nm sidewall spacing and 350 nm trench depth,
  10 minutes of dry liner oxidation at 1000 C, oxide filling and planar cut.
- `decks/sti/sti_rounded.in`: the same process with 30 nm upper-corner and
  50 nm lower-corner radii (12 polygon segments per quarter circle). Initial
  substrate spacing is 5 nm around the trench edges and grows toward 100 nm
  in the deeper substrate. Oxide spacing is 10 nm and deposition has 24 layers.

Refinement lines initialize the substrate mesh; they are not persistent interior
constraints during etch/retriangulation. Boundary refinement therefore does not
guarantee a uniformly refined band beneath the trench floor. That would require
interior sizing constraints during remeshing.

## Geometry corrections

The original clipping algorithm independently closed each retained boundary arc.
A planar cut through a U-shaped oxide liner consequently closed its inner arc
across the trench void. It then tried to interpolate oxide solution values in
space that had never belonged to that oxide region and aborted. The correction
collects retained arcs, joins them in clipping-perimeter order with material on
the same side of every boundary, and emits the actual connected components.
Point contacts between retained pieces are not treated as nonzero interior
connections. Clipping corners are preserved, except for exact coordinate coincidences with
intersection endpoints, which would introduce zero-length edges.

Edge insertion can append multiple material nodes; the classification array now
grows with the node list. Intersection search checks all clipping edges for a
true double-precision segment intersection before considering legacy ON
tolerances. This avoids selecting an earlier near-corner match instead of the
actual intersected segment. All-in/all-out paths also free their scratch array.

Viscous oxidation also exposed an invalid edge -1 access in layer insertion.
When all neighboring edges point along or tangential to interface motion, none
qualifies for splitting. Layer insertion now leaves that node alone in this case.

## Validation and limitations

Candidate: `suprem4/build/native64_e0r_18r5/src/suprem4` (AMD ordering).
The installed launcher remains on its previously validated build.

This candidate also completes the complete PMOS flow at 9,815 and 19,809 initial
points: all 20 saved meshes pass. The final normal-density geometry passes a
comparison against the installed-build control at 0.1% material-area and 5 nm
bounding-box tolerances. Results:
`/tmp/suprem4_mesh_density__p6nf0r5/results.json`.

Both STI decks pass the baseline, finer initial spacing (factor 1/sqrt(2)), and
translations of the entire trench polygon by +/-25 nm. All 40 saved meshes pass
positive area, finite fields, duplicate-triangle, edge-ownership and reciprocal
neighbor checks. Trench area is independently calculated from the etch polygon;
the final filled cross-section has area 2 square microns.

Results: `/tmp/suprem4_isolation_6x8tbdbq/results.json`.
The rounded baseline finishes with 3,035 points and 5,782 triangles; the finer
rounded case has 4,274 points and 8,174 triangles.
Final baseline structures and provenance are retained in `decks/locos/output`,
`decks/sti/output`, and `decks/sti/rounded_verified_output`. The earlier rounded
structure shown in the desktop viewer is retained in `decks/sti/rounded_output`;
its `result.json` identifies the exact build and tests for that artifact.

LOCOS baseline and +/-25 nm mask shifts complete with valid saved meshes.
The baseline silicon surface is about 0.238 um deep in the field, 0.000865 um
under the distant mask, and 0.140 um at a probe 50 nm inside the mask edge.
These probes demonstrate masked growth and lateral encroachment, not calibration.
The substrate-only finer LOCOS case still encounters an oxidation mesh inversion
and remains a known failing regression. A later study demonstrates passing 2x
and 4x configurations when oxide insertion spacing and deposited pad/nitride
layers are refined along with the substrate; see `doc/locos_refinement.md`.
Experimental motion-rejection and coordinate-normal
changes did not resolve it and are not included in this candidate.

Run the supported STI exercises:

```sh
python3 suprem4/tests/run_isolation.py \
  --binary suprem4/build/native64_e0r_18r5/src/suprem4 \
  --cases sti sti_rounded --jobs 2
```

The runner's default also includes LOCOS, including the known failing fine case.
Focused sanitizer checks include `tests/run_etch_corner.sh` (concave multi-arc
cuts, rotations, disconnected components and exact/near corner intersections)
and `tests/run_new_layer_corner.sh` (ineligible/tangential/stationary motion and
an eligible edge split). Existing mesh/ordering sanitizer checks also pass.

The retained clipping algorithm still assumes crossings are resolved by the
subject boundary sampling. A narrow etch can cross an edge twice while both
endpoints remain outside; the older all-out shortcut does not detect that case.
This is a separate remaining clipping limitation, not covered by these passes.
