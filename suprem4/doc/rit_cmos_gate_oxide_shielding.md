# Directional contact-etch shielding repair

**Subsequent presentation audit:** retained gate oxide is verified, while
profile and electrical agreement remain unqualified. See the
[dated audit](rit_cmos_presentation_audit_20261117.md) for the demonstrated
2D poly-implant resolution error, oxide/well deviations and 3.3 V results.

2026-11-17. **Gate-oxide repair verified:** both fresh nominal processes finish
with 16/16 audited checkpoints. All five gate-oxide thicknesses are exactly
unchanged between activation and first metal. Center values: NMOS 16.67454 nm,
PMOS 16.97010 nm. Every downward-facing gate boundary has an oxide neighbor;
all four final contacts connect to the intended conductors.

This is scoped process/geometry verification. General launcher promotion and
PDF-calibrated electrical behavior are not claimed. Machine-readable evidence:
`gate_oxide_validation_20261117.json`.

## Root cause and removal geometry

The old `etch <material> dry thick=...` translated the free surface downward
and subtracted that cutter from every exposed-connected target region.
Connectivity to the ambient was mistaken for directional accessibility.
The RIT contact step (`etch oxide dry thick=1.5`) consequently removed the
buried gate dielectric under polysilicon, even though the polysilicon itself
was not etched. Changing the requested depth would hide this defect.

The new 2D implementation snapshots the whole material partition before
removing any region. Its vertical front is the piecewise-linear lower envelope
of the depth-limited free surface and the shadows of all non-target solid
regions. Every solid material other than the selected target is opaque; no
special rule names a gate, oxide thickness, region area, or CMOS deck.

For each opaque region, split its x projection at boundary-vertex abscissas.
On each open slab, sort boundary crossings by depth and pair successive
crossings to obtain every occupied vertical interval. Entry and exit are
linear functions of x. Terminate the cutter at the midpoint of each chord.
Take the lower envelope over
all these segments, splitting at segment endpoints and actual line crossings.
There are no sampled rays or mesh-size displacement parameters.

This interior termination has a set-theoretic justification. For each ray,
the interval between the first opaque entry and an interior termination lies
inside non-target material. Therefore its intersection with the target is
empty. Moving the cutter endpoint inside that interval cannot change which
target material is removed. The union of region shadows retains this property
for nonconvex solids and is independent of internal triangulation. In particular,
a mask protects target material in a disconnected region below it too.
The cutter may geometrically enter silicon or poly, but the selected-material
subtraction never removes those materials.

Stopping directly on either entry or exit interfaces instead introduced
coincident clipping edges and unnecessary subdivision of material interfaces
in the legacy skeleton Boolean implementation. Independent audits rejected
those preliminary candidates. A triangle-based intermediate implementation
also let internal mesh vertices imprint isolated contacts with the entry
interface. Region-boundary decomposition eliminates that dependency.
Interior chord termination avoids that
representation problem while preserving exactly the same target-removal set;
it is not an offset of the physical etch surface or a post-etch fill operation.

The final cutter is projected once to the engine's stored float coordinate
type. Existing exact duplicate canonicalization and skeleton subtraction then
apply. This change applies to the **2D vertical dry etch command**. It does not
add angular plasma flux, finite mask erosion, or repair the separate 1D model.

## Polygon-subtraction repair

A second root cause appeared in an analytic stack: an existing subject edge
had both endpoints ON two different vertical cutter walls, with its entire
interior outside the cutter. No individual subject edge had multiple proper
crossings. `node_multiple_crossings` returned before inserting an interior
classification sample. The retained arc walker saw only ON endpoints, lost
the real subject arc, and followed the recessed cutter boundary through the
opaque gate, creating material outside the original target. Interpolation
correctly rejected the manufactured material as outside the source mesh.

Interior classification of ON-to-ON edges now runs independently of whether
proper crossings were discovered. Exactly collinear boundary intervals need
no extra interior point. `check_on_on_retained_arc` in `tests/etch_corner.c`
checks the independently known retained rectangle area in four orientations.
The former source snapshot fails that regression; the repaired source passes
the complete existing clipping suite under sanitizers.

## Independent acceptance checks

* `run_directional_etch.sh`: analytic finite/sloped/overlapping shadows,
  endpoint reversal, narrow masks, and opaque-cell triangulation invariance;
  address and undefined-behavior sanitizers.
* `run_directional_stack.py`: oxide/poly/oxide stack at 5/10/20 divisions;
  all exposed oxide clears, the complete buried rectangle retains its input
  represented area, silicon/poly areas are unchanged, and the output partition
  has positive nonoverlapping triangles and noded boundaries. The nominal
  buried area is 0.04 square microns; its reference uses the input's actual
  float-projected rectangle rather than assuming exact decimal deposition.
* `run_rit_locos_strip.py --case nmos_contact_shielding`: frozen real failing
  contact input, unchanged 1.5-micron etch, independent partition audit.
* `run_contact_shielding.py`: real contact input, +/-2 nm whole-structure
  origin translations and +/-1 percent etch-depth changes; requires preserved
  dielectric AND opened body/source/gate/drain windows after resist removal.
  These are origin translations, not independently shifted masks.
* `run_process_surface_regression.py`: planar cancellation, retained real
  films, restart, construction invalidation and sloped-film control.

`decks/subcmos2014/validate.py` now requires actual oxide material in the gate
gap and an oxide neighbor along **every downward-facing polysilicon boundary
edge**, in addition to the five thickness probes. A regression explicitly
places a void between probes to prevent a return to sampled-only acceptance.
Final metal checks also require four islands and contact to their conductors.

## Candidate and records

Current ordinary candidate: `suprem4/build/native64_rxq2pnrb/src/suprem4`.
Experiment root: `suprem4/experiments/directional_etch_20261117`.
Fresh full-process run: `decks/subcmos2014/results/directional_final_fresh_20261117`.
Earlier intermediate restart comparison from previously saved annealed structures:
`decks/subcmos2014/results/directional_cell_resume_20261117`.
Current full electrical rerun starts only after fresh process audits pass:
`suprem4/experiments/directional_etch_20261117/electrical`.

Earlier `native64_ec3ihno0` (entry-surface caps) passed the contact replay and
resumed final structures but failed a fresh upstream strip. The subsequent
`native64_wpdiiq88` (exit-surface caps) failed independent contact and analytic
stack partition checks. Neither is accepted. Their logs are retained.

The previous electrical comparison in `rit_cmos_electrical_pdf_comparison.md`
remains a valid diagnosis of the old invalid inputs, not a prediction for
these corrected structures. Fresh electrical results must be identified by
input hashes and model assumptions before making electrical claims.

## Additional strict inventory probe

`final_non_target_inventory.json` records a deliberately strict 1e-12 relative
before/after comparison in nonetched silicon/poly. All three analytic stacks
pass (largest error 6.99e-15). The captured real contact etch has chemical
inventory differences up to 2.12e-9 relative and **does not meet 1e-12**.
That result is retained, not corrected or accepted by widening the criterion.
The operation refines interfaces stored in native float coordinates; this is
not the fixed-domain, double-coordinate conservative-rezoning case. A general
release still needs a separate coordinate/transfer error budget and full-flow
mesh-convergence qualification. Gate-film occupancy, opened contacts and
partition validity remain independently checked.

## Reproduce

```sh
SUPREM4_PROCESS_SURFACE=1 SUPREM4_DEPOSITION_GEOMETRY=domain \
SUPREM4_MATERIAL_DOMAINS=oxide SUPREM4_MESH_PATCH_REMOVAL=1 \
python3 decks/subcmos2014/run.py \
  --binary suprem4/build/native64_rxq2pnrb/src/suprem4 --through metal1
```

The recorded full fresh runtimes and every checkpoint audit are in the fresh
run's per-device `result.json` files. Bounded final checks all pass: three stack
mesh densities (ordinary and sanitizers), five real contact origin/depth cases,
29 construction/restart controls, seven captured failures, the analytic
clipping suite and directional-region tests. No recipe depth was reduced and
no post-etch oxide was inserted.

## Electrical follow-up on the corrected fresh structures

The unchanged baseline model from the earlier diagnostic was rerun: Genius
DDML1, 300 K, 1-micron device width, active B/P, Fermi statistics, Lucent mobility,
ideal n+ gate work function 4.17 eV, explicitly zero fixed interface charge.
Both devices finish 67 converged gate points, 0 to +/-3.3 V, with drain held at
+/-0.05 V. No ERROR or DIVERGED markers appear. Export manifests identify exact
source meshes and preserve every vertex/triangle (poly plus gate metal are
one ideal conductor, so poly depletion is omitted).

| Device | Drain-current magnitude at Vg=0 | At 3.3 V gate drive | On/off ratio | Maximum terminal KCL residual |
| --- | --- | --- | --- | --- |
| NMOS | 0.181 nA | 11.3749 uA | 62,926 | 2.20e-15 A |
| PMOS | 4.878 nA | 5.2331 uA | 1,073 | 4.77e-16 A |

These are unfitted model results. They demonstrate restored gate control,
not agreement with the PDF's nominal +/-0.75 V thresholds. No threshold is
extracted in this check. Numerical current-conservation residuals do not prove
mesh convergence or process-model calibration.

Records: `suprem4/experiments/directional_etch_20261117/electrical/summary.json`,
per-device transfer inputs/logs/data/CSV and `export_manifest.json`.
Plots: `electrical/corrected_transfer.png` and `final_gate_oxide.png` under the
same experiment root. Corrected final `nmos.str` and `pmos.str` are open in the
viewer. The general launcher remains unchanged pending broader qualification.
