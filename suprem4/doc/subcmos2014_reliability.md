# Sub-CMOS 2014 development and reliability coverage

**2026-11-17 update:** the gate-oxide shielding and missing retained-arc repairs now pass both fresh 16-checkpoint processes. See [current repair and limits](rit_cmos_gate_oxide_shielding.md). Earlier missing-oxide structures and withdrawn acceptance below are historical evidence.

**Electrical follow-up: final acceptance withdrawn.** The contact etch removes
buried gate oxide in both devices. Process completion and the earlier checks
missed this material void. See [actual electrical/PDF comparison](rit_cmos_electrical_pdf_comparison.md).
The validated RIT pin is removed and the acceptance manifest is marked false;
completion claims below describe the earlier, insufficient validation.

2026-11-16 update: both fresh decks complete all 16 checkpoints on the
construction candidate. See [implementation and validation](process_surface_construction_repair.md)
for current evidence and remaining scope. The investigations below are historical.

This work adds `decks/subcmos2014/` from the local RIT process PDF. It is a
process-deck implementation and an engine development candidate, not an
assertion of calibrated transistor performance or general engine robustness.
The installed launcher has not been promoted.

The current source candidate is `suprem4/build/native64_4vwld4yj/src/suprem4`.
Both fresh-prefix/resumed device runs passed all eight checkpoints through
`gate_oxide` in `decks/subcmos2014/results/gate_oxide_check/`. This is partial
process validation, not end-to-end success. The complete 16-stage processes
are **not accepted**. A reference PMOS etch replay still fails with `error in
chop`, and a later full-flow run also found a region-discontinuity failure.
An NMOS source/drain probe found the assumed body tie had the wrong net doping
polarity. These findings prevent promotion of the development build.

The latest candidate was also checked in
`decks/subcmos2014/results/dopant_check/`: NMOS completed through source/drain
implantation but failed the body-tie polarity check (net donor concentration
about 1.03e20 cm⁻³, where a p-type tie is required); PMOS still failed during
spacer formation. The numerical exit status alone is therefore not an
acceptance criterion. `subcmos2014_validation.json` records the exact candidate,
passed focused checks, 27 local perturbations, partial successes and failures.

## Engine changes motivated by these decks

- **Matrix assembly and directed boundaries:** `assem_map` orders an edge's
  endpoints for matrix indexing. It now reverses the direction flag in every
  referring skeleton when swapping those endpoints. Otherwise a subsequent
  native-oxide operation traverses the wrong vertices and creates repeated-
  vertex triangles. A two-region shared-edge test verifies both traversals and
  repeated assembly.
- **Implant tables:** the shipped table body uses SUPREM-IV material IDs, as
  recorded by its 1992 modification note; its older introductory header still
  lists SUPREM-III IDs. Oxide must look up material **1**, not material 5. The
  BF₂ ion must look up **6**, not elemental boron 5. Native P-in-oxide data
  support 125 keV. The lookup handles exact first rows and skips malformed
  records rather than using stale/uninitialized moments. Seven-column BF₂
  records use a documented isotropic lateral-straggle estimate.
- **Failed implants:** `do_implant` returns failure before transforming the
  mesh. The caller now respects that return instead of restoring coordinates
  from an uninitialized/stale transform buffer. A missing-table test compares
  geometry, connectivity and the original boron values before/after failure.
- **Clipping predicates:** boundary classification and endpoint intersection
  use a common finite-segment predicate at float storage precision, replacing
  the inconsistent fixed one-Angstrom ON band. Collinear endpoints have a
  defined intersection. Near-coincident clipping corners do not introduce a
  reversing microscopic edge.
- **Multiple crossings:** an individual source edge can cross a curved etch
  perimeter more than once. Such edges are split at all crossings, and both
  perimeters receive the same stored intersection coordinates before arc
  classification. ON-to-ON intervals receive a midpoint sample. The regression
  includes a notch that crosses an edge twice while every original subject
  vertex lies outside, checking the retained corners and analytical area.
- **Refinement status:** a negative `sp_edge` return identifies an existing
  endpoint. Boundary subdivision no longer treats it as another edge ID.
  A nonmoving split request does not move the retained endpoint when minimum
  spacing prevents subdivision.
- **Structure restart:** coordinates and solution values are written with
  sufficient precision; coordinate units are converted in double precision
  before rounding to stored float coordinates. Six-digit coordinate output
  flattened/inverted narrow triangles in the spacer case. A narrow-triangle
  save/reload test checks geometry and solution preservation.
- **Native oxide joins:** new seed oxide tapers to existing oxide at exposed-
  run endpoints rather than displacing and inverting a residual oxide wedge.
  Gas nodes on neighboring boundary points are not mistaken for covering solid
  material. A one-edge gap with oxide at both ends receives an interior seed
  vertex; connectivity is rebuilt before computing its normal. Tests cover a
  bare domain, an existing oxide wedge, and a gap between two oxide wedges.
- **Junction collapse:** `rem_edg` checks the actual edge endpoint nodes before
  its material search. The previous check used a search result that could be
  `-1`, dereferencing `nd[-1]` at a three-material junction.
- **File-defined deposition:** initialize the side-squaring option for an
  explicit deposition string; that path previously passed an uninitialized
  option to skeleton construction.

These changes preserve the requested implant energies and species. No
synthetic initial oxide, mirrored silicon moments for oxide, or numerical
change to spacer thickness is needed to bypass the failures above.

## Reproduce focused checks

From the workspace root, with a freshly built candidate:

```sh
sh suprem4/tests/run_assembly_boundary.sh
sh suprem4/tests/run_implant_tables.sh
sh suprem4/tests/run_etch_corner.sh
python3 suprem4/tests/run_implant_failure.py --binary /absolute/path/to/suprem4
python3 suprem4/tests/run_structure_roundtrip.py --binary /absolute/path/to/suprem4
python3 suprem4/tests/run_native_oxide_join.py --binary /absolute/path/to/suprem4
```

The C checks use AddressSanitizer and UndefinedBehaviorSanitizer. The clipping
suite includes shared interfaces, retained corners, concave components,
tangencies, coincident and near-coincident corners, parallel and collinear
edges, multiple crossings, minimum-spacing statuses, and junction collapse.

`run_etch_snapshots.py --binary ... --input ... --depth ...` exercises a saved
strip/spacer case at thickness factors 0.99, 1 and 1.01 and origin shifts
−0.01, 0 and +0.01 µm. These are **local geometry perturbations**, not full-flow
mask perturbations or mesh-convergence evidence. The captured NMOS field strip,
PMOS spacer and later NMOS spacer each passed all nine combinations during
investigation. Logs and exact inputs are retained in the printed run folders.
Existing LOCOS and rounded-STI baseline cases also passed their geometry audits
on the accumulated candidate.

## Limits and open issues

The full viscous stress model has still reported nonlinear convergence failures
on later thin-oxide steps. Its legacy continuation code can announce "Giving
up" and continue with a partial state; this has **not** been repaired here. The
new runner rejects such messages instead of reporting success. A passing local
geometry regression does not establish convergence of that model.

The default deck retains stress-dependent viscous LOCOS, then uses vertical
thin-oxide growth. This omits lateral thin-oxide growth and later mechanical
stress feedback. `build_decks.py --post-locos-model viscous` preserves the fuller
model for further investigation; it is not a validated variant. Later implant
and contact masks use an explicit flat-top resist approximation. The conformal
offset routine can still create a malformed coating over spacer topography.
Using a flat resist surface is a process-model assumption appropriate to a
spin coat; it does not repair that general conformal-offset limitation.

The older RIT PMOS reference comparison has shown about a 2.6–2.7% oxide-area
change after gate oxidation with corrected implant lookup. This is outside its
stored geometry tolerance and must not be presented as a passing golden-file
regression. Geometry auditing and physical calibration are distinct checks.
With subsequent changes, the older PMOS flow has also failed during its gate
oxide strip. Its captured replay still reproduces the triangulation failure
on the current candidate. This regression remains open; the passing focused
tests do not supersede it.
