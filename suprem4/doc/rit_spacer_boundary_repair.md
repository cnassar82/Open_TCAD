# RIT CMOS mesh and process repair report

**2026-11-17 update:** the gate-oxide shielding and missing retained-arc repairs now pass both fresh 16-checkpoint processes. See [current repair and limits](rit_cmos_gate_oxide_shielding.md). Earlier missing-oxide structures and withdrawn acceptance below are historical evidence.

**Electrical follow-up: final acceptance withdrawn.** The contact etch removes
buried gate oxide in both devices. Process completion and the earlier checks
missed this material void. See [actual electrical/PDF comparison](rit_cmos_electrical_pdf_comparison.md).
The validated RIT pin is removed and the acceptance manifest is marked false;
completion claims below describe the earlier, insufficient validation.

2026-11-16 update: both fresh decks complete all 16 checkpoints on the
construction candidate. See [implementation and validation](process_surface_construction_repair.md)
for current evidence and remaining scope. The investigations below are historical.

Implementation handoff: [Sol’s mesh repair guide](sol_mesh_repair_implementation_guide.md).
The later numerical and material-domain investigation, with build hashes and
fresh-flow results, is in [material domains and coupled-solver reliability](material_domain_reliability.md)
and its [result manifest](material_domain_reliability_results.json).

## Current status, 2026-11-15 — open-source TCAD follow-up

The [research and implementation report](open_source_tcad_repair_research.md)
compares DEVSIM, Genius and ViennaPS and records the remaining algorithmic gaps.
The [reflecting-boundary report](reflecting_boundary_repair.md) gives the new
construction-level diagnosis and failing controls.

The Clipper2 deposition-domain backend now passes the captured PMOS ILD and
spacer cases, ordinary and strict ASan/UBSan with leak detection disabled.
It replaces raw offset-path loop deletion with occupied-domain dilation and
Boolean difference. Its output uses a dyadic grid exactly representable in the
float database; original material vertices remain unchanged. Unsupported holes
and multiple coating components are explicitly refused.

The final ordinary candidate is `native64_l3ozbabw`; its ASan/UBSan counterpart
is `native64_92cvoyp7`. Both source manifests match the working engine. Thirteen
focused numerical/geometry checks pass; the domain backend also passes its
analytic-area, translation and storage-representability sanitizer tests.

The archived PMOS continuation from `annealed` through contacts and `metal1`
now completes: **16/16 runner checks pass**. Separate partition audits pass
21,159 contact triangles and 37,062 final metal triangles. This is a resumed
success, not fresh full-process acceptance.

| Acceptance scope | Latest completed evidence | Status |
| --- | --- | --- |
| Captured reflecting-side construction | Old source leaves its prescribed plane; corrected sanitizer build passes | PASS, bounded |
| Captured finite spacer clip and exactly collinear residual | Ordinary and sanitized partition checks pass | PASS, bounded |
| Captured ILD deposition with domain backend | 49,207 positive, nonoverlapping triangles; existing region areas unchanged | PASS, bounded |
| Archived PMOS `annealed` → `metal1` | `astra_final_pmos_finish/pmos`, 16/16 checks; final outputs match the independently audited preceding continuation | PASS, resumed |
| Fresh PMOS with domain backend | Fails during LOCOS; actual nitride boundary crossing and overlapping triangles | FAIL |
| Fresh NMOS with domain backend | Reaches LDD, then a spacer quad has overlapping collinear side edges | FAIL |
| Recursive velocity history | Original fails nested-halving/reordering regression; corrected controller passes | PASS, bounded |
| Experimental P1 mass and BE integrator | Analytic conservation/convergence tests pass; process replays still fail | EXPERIMENTAL, not accepted as CMOS fixes |
| Full fresh CMOS, physics and scaling convergence | Not achieved | OPEN |

Evidence is under `decks/subcmos2014/results/astra_*`; build/source hashes and
runtime options are retained in the individual result files. The validated
launcher still resolves to `native64_2r55t88h`. The experimental domain backend
requires a build with `--domain-offset` and `SUPREM4_DEPOSITION_GEOMETRY=domain`,
plus the existing oxide-domain and point-removal candidate options.

The remaining LOCOS boundary is genuinely self-intersecting: 601 vertices,
turning 4*pi, crossing segments 232 and 234. The stable angle formula reproduces
it, and an independent mesh audit finds overlapping triangle interiors.
The NMOS spacer failure also retains a side boundary that doubles back. Do not
repair either by loosening winding tolerances or discarding a finite-area
fragment. Trace the first invalid motion/construction mutation, then implement
a complete constrained domain replacement with conservative field transfer.
The old `DetectLoop` code has a zero edge-vector bug and unfinished material
swallowing branches; simply enabling that heuristic is not a valid repair.

The sections below retain historical implementation and control evidence.
Earlier build/status claims apply to their recorded candidates.

The recent coarsening rewrite described below is **withdrawn and not accepted**.
Its history is retained as a negative control; the newer opt-in constrained
point-removal operation is documented in the linked material-domain report.
Independent geometry checks passed on the withdrawn rewrite's captured panic
mesh, but the engine failed an exposed-point/GAS-node invariant. Four active
files were restored to the
pre-experiment `native64_b43fm60y` source snapshot; experimental sources, tests,
patches and hashes are preserved in `suprem4/experiments/rit_cavity_20261115/`.
The old coarsener is still defective.

The [first invalid removal diagnosis](rit_first_invalid_removal.md) now identifies
a noncollinear open boundary fan selected by edge collapse. The legacy cavity
closes across an existing edge; it also deletes GAS before constructing the
silicon replacement. A narrow opt-in refusal avoided that crossing but failed
later on neighbor reciprocity, so it was withdrawn.

Further coarsening work must follow the
[mesh replacement contract](mesh_patch_replacement_contract.md): extract the
first invalid mutation, construct a complete multi-material replacement in
private storage, validate geometry, ownership and fields, and commit atomically.
Do not add downstream cleanup to make a failed mesh appear acceptable.

The sections below retain implementation history. Earlier remaining-defect
statements, gate results and probe values describe their recorded candidates,
not the latest result table above.
Use the [reliability plan](rit_cmos_reliability_plan.md) for release scope.

## Follow-up implementation, 2026-11-14

The parser startup UB and a subsequent mesh-checker bounds read are repaired.
`misc/read.c` now passes the valid `cards` base and a logical last index of -1;
`misc/reader.c` advances integer indices within the array and retains a NULL
terminator. A real binary-key regression checks first/root cards, nested
parameters and siblings. Reintroducing the before-array pointer fails strict
UBSan. The focused test also passes with leak detection enabled.
`dbase/check.c::mtest2` now begins its reverse vertex scan at `num_vert(t)-1`,
rather than reading index 3 of a three-entry triangle array.

The full spacer replay passes ASan and strict UBSan on `native64_qxu4kodi`,
with leak detection disabled. See [sanitizer results](rit_parser_sanitizer_results.json)
and `decks/subcmos2014/results/meshcheck_sanitized_spacer/` for persistent
provenance, logs and all three outputs. This supersedes the startup blocker
recorded in the original verification below, without changing its historical
results. Full-engine leak checking remains open. The builder now supports
`--sanitize` directly.

Per-command chemical-species traces are available through
`decks/subcmos2014/trace_poly.py`; its linear-field inventory/probe test passes.
The resumed NMOS trace first develops substantial body phosphorus during poly
activation: 1e5 cm^-3 immediately after implantation, 6.67e15 after the second
800 C segment, and 1.31e20 after the full activation sequence at the 50 nm
probe. Its body oxide gap is 4.64 nm. These measurements locate the change;
they do not establish conservation, numerical convergence or correct transport.

The PMOS heavy-implant trace reveals a separate **remaining coating defect**:
after spacer etching, the unpatterned resist ends at x=5.54049 um. The drain at
x=7.75 is uncovered even before resist patterning. Its phosphorus rises from
1.60e17 to 1.80e20 cm^-3 during the body-tie phosphorus implant. Thus the NMOS
fixture's full-width success does not generalize to PMOS. Extend the exposure
repair and add a PMOS coverage regression before claiming the spacer work is
complete for both devices. Persistent traces are under
`decks/subcmos2014/results/{poly_trace_nmos,poly_trace_pmos,heavy_trace_nmos_complete,heavy_trace_pmos_complete}/`.
The first heavy-trace attempts omitted `resist_top.dat` and failed before
coating; those diagnostic harness failures are retained separately. The runner
now copies that input.

### Extended PMOS exposure repair

The PMOS failure is now repaired. Refinement left neighbor arrays referring to
obsolete unsplit parent faces. The exposure snapshot now finds matching live
faces through shared physical point identities, including distinct material
nodes at the same point. It requires owners on opposite sides of the face;
old and replacement triangles coexist during remeshing and their coincident
outer edges must not be mistaken for interfaces. A broad owner match without
this side test was rejected because it exposed artificial domain sides.
Exact stored-float orientation signs now support that side test.

The captured PMOS replay fails on `native64_88fw7xz8` at x=5.54049 um and passes
on the repaired candidate. Both NMOS and PMOS full-width replays pass strict
ASan/UBSan on `native64_3ur90zja` (leak detection disabled). The focused interface
ownership test passes with leak detection enabled. The resumed PMOS heavy
implants now pass all three polarity probes: drain P=1.603e17 and B=2.523e19
cm^-3, body net donor=3.234e20 cm^-3. This closes the captured PMOS mask failure,
not fresh full-process acceptance.

### Fresh-flow diagnosis and exact intersection repair

Fresh NMOS execution on `native64_mad7zzna` failed at the 50 nm dry oxide strip
after LOCOS. A single-command fixture preserves its predecessor. Double
intersection parameters accepted the extension of a mask segment beyond its
endpoint near the coordinate origin, inserting an outside corner and a
reversing boundary. A trial tightening coordinate range checks did not fix it
and was removed.

The engine now checks exact orientation signs before accepting finite segment
intersections, including preliminary multiple-crossing detection. Products of
two stored binary32 coordinates fit exactly in binary64; error-free expansion
sums preserve their determinant sign across large exponent differences. This
adds no geometric snapping tolerance. A rational oracle checks 10,003 cases.
The captured LOCOS strip fails on the control and passes on `native64_sca2gh2a`,
including independent partition checks on 9,642 triangles; strict ASan/UBSan
also passes on `native64_1ikkuqm2`. The focused clipping suite still passes.
Use `suprem4/tests/run_rit_locos_strip.py` to reproduce it.

Fresh PMOS on the preceding candidate reaches `source_drain` with passing
polarity checks, then aborts during activation with an inverted triangle.
The saved `pmos_activation_trace` isolates the onset to the first dry-oxygen
segment after two nitrogen segments. A separate restart regression proves
that `diffuse ... continue` erased the configured initial step: two independent
anneals requested at 0.125 s instead started at 0.125 and 0.6 s. The proposed
repair retains the method setting and passes zero only for the individual
continuation call. Its verification and fresh-flow results remain pending.
This timing defect is established independently; whether it explains the
activation inversion is still being tested.

The planar reference tools predict P around 2.0–2.6e17 cm^-3 at 50 nm beneath a
4.64 nm oxide, versus 1.31e20 in the archived RIT body. Halving the fixed step
from 30 to 15 s changes the refined 1D probe by about 0.084%; the 2D comparison
and full RIT profiles are not yet converged. The archived mesh has coincident
Si/poly points without shared Si/poly edges. Their role in transport is an
open hypothesis, not established direct contact. These same-engine references
are diagnostic comparisons, not independent physical calibration.

## Causes and repairs

Two defects were involved.

1. **A shared junction was recomputed after rounding.** Clipping the first oxide
   region inserted a point into its interface with the next region. That point
   classified ON the etch boundary. The next region nevertheless solved the
   nearly parallel intersection again, using the rounded, subdivided edge. It
   moved the junction by 28 float ULPs in x, producing exposed traversal A–B–C–D
   instead of spatial order A–C–B–D. This caused an overlapping collinear edge.

   `refine/skel.c::classified_crossing` now preserves the classified endpoint
   when it belongs to the clipping segment actually intersected. It uses the
   existing finite-segment predicate; it does not enlarge a tolerance or move
   existing vertices. When an ON endpoint belongs to another mask side, the
   actual crossing is retained. Existing coincident-corner tests exercise this
   distinction.

2. **The uncovered poly gate was not exposed.** Removing the old covering oxide
   did not expose interface points classified ON instead of IN. The subsequent
   surface walk stopped at x=5.5 um. Fixing only the junction let deposition
   return successfully while covering only half the structure.

   `refine/etch.c` now snapshots interfaces to other regions after refinement
   and before deleting the old triangles. After the assembled mesh rebuild, a
   previously covered face is exposed only if it has no solid neighbor. Retained
   cover stays internal. Snapshot coordinates survive database compaction; IDs
   are not assumed stable. The additional logic is restricted to 2D etching.
   Marking all deleted triangle vertices as exposed was rejected because it
   incorrectly exposed artificial domain sides.

## Verification

Candidate:

```text
suprem4/build/native64_mnxbw7bf/src/suprem4
SHA-256 4348b0446e0c073374d5eac6aa6ea013b183feb034900355f2edede28d1cdb02
```

The candidate's C/header hashes match the working engine source at report time.
Full machine-readable provenance, stage statistics, and remaining failures are
in `rit_spacer_fix_results.json`.

| Check | Result |
| --- | --- |
| Captured LDD → spacer deposition → spacer etch → resist deposition | PASS |
| Independent triangle-interior, boundary crossing, overlap and unnoded-contact checks on all three outputs | PASS |
| Resist extent | x=0 to 10.000000475 um, top=-1.999999949 um |
| Coated material area | 99.999924779 um²; nominal rectangle 100 um², preselected replay tolerance 0.0001 um² |
| Reduced shared-junction test with preservation disabled | FAIL as expected: loses C and inserts B, 28 ULPs away |
| Junction-only candidate without exposure repair | FAIL as expected: resist ends at x=5.5 um |
| Focused AddressSanitizer, UndefinedBehaviorSanitizer and LeakSanitizer suite | PASS |
| Shared-cut orientation, start-edge, scale, origin and 90-degree variants | PASS, 70 cases |
| Independent-audit and boundary-audit unit tests | PASS, 8 tests |
| LOCOS, STI and rounded STI baseline simulations | PASS |
| Resumed NMOS through heavy implants | Exit 0; geometry PASS; body polarity FAIL |
| Resumed PMOS through heavy implants | Exit 0; geometry PASS; drain polarity FAIL |

The independent audit uses filtered orientation signs with exact rational
fallback on stored coordinates. It detects overlapping triangle interiors,
including containment, plus crossing/overlapping boundary edges and unnoded
boundary contacts. It does not prove the intended process geometry, correct
material labels, absence of every possible hole, or dopant conservation. The
flat-coat replay has a separate reference extent/area check.

The reported 100 um² is the sum of **all material areas after coating**, not
the photoresist area alone. Its deficit is 0.0000752214 um², about 75% of the
preselected absolute tolerance. It passes that bounded check; the discrepancy
is not zero and has not been assigned a demonstrated cause. Retain the measured
value and investigate its convergence rather than increasing the tolerance.

The 70 cases exercise the reduced shared-cut geometry. Coordinate scaling by
0.5x/2x is not mesh-density refinement, and these cases do not establish full
RIT stability at 0.5x, 2x or 8x point counts. The disabled-junction failing
control establishes the captured baseline defect, not a recorded failing-control
run for all 70 variants. The isolation results are baseline cases only.

A full-engine sanitizer attempt also found an existing startup bounds violation
at `misc/read.c:58`, `read_list(&(cards[-1]), fd)`. With strict UBSan it stops
before geometry execution. A clean full-engine sanitizer result is therefore
still open; the focused kernel sanitizer suite above passes. The replay runner
now rejects reported runtime sanitizer errors even if the process exits 0.

## Reproduction

From `/home/chris/Open_TCAD`:

```sh
python3 suprem4/tools/build_native64.py --ordering amd --jobs 4
python3 suprem4/tests/run_rit_nmos_spacer_regression.py \
  --binary /absolute/path/printed/by/builder/src/suprem4
sh suprem4/tests/run_etch_corner.sh
PYTHONPATH=suprem4/tests python3 -m unittest \
  suprem4/tests/test_partition_audit.py \
  suprem4/tests/test_boundary_axis_overlap.py
```

The replay rejects populated output folders to prevent stale structures from
being mistaken for new output. Its compressed LDD fixture and checksum remain
in `suprem4/tests/fixtures/rit_nmos_spacer/`.

For a local etch trace, set `SUPREM4_ETCH_TRACE_X` to two coordinates in centimeters.
For example, `0.000457032 0.000457035` captures the affected shared interface.
The trace records subject/clip rings and retained rings before/after meshing.

## Implementation order and success gates

Keep the successful spacer fixture as a mandatory regression for subsequent
changes. Work through the following gates in order, recording each as PASS,
FAIL, or NOT RUN. An unresolved gate stays open even when a later diagnostic
run completes. Do not change doses, energies, furnace budgets, or acceptance
tolerances to turn a failure into a pass.

### 1. Restore full-engine sanitizer coverage

Repair the startup parser's before-array pointer construction and verify both
first-card and nested parameter-list parsing. Add a regression that detects
the original undefined behavior. Run the complete spacer replay under ASan and
strict UBSan after the parser fix, with timeouts and explicit model/data paths.
Report leak checking separately: the recorded full-engine attempt disabled it.
Passing the focused clipping suite does not close this gate.

Exit evidence: parser regression, sanitizer build manifest, clean replay log,
and all three audited outputs. Current status: **PASS for the captured spacer
replay** on `native64_qxu4kodi` and later bounded sanitized replays. Full-engine
leak detection remains **NOT RUN** as an accepted gate. The startup failure in
the original verification section is a historical negative control.

### 2. Locate and explain the first wrong phosphorus transfer

Instrument both device flows immediately before/after the poly implant and
each activation segment, then gate patterning, LDD and heavy implants. Keep
chemical and electrically active concentrations separate. For each checkpoint,
record species inventory by material/region, body and active-area profiles,
oxide thickness and poly/oxide/silicon connectivity. Use consistent units for
2D inventory per unit out-of-plane depth; include implant inputs, boundary
losses and interface transfers in the balance.

Find the first unexpected increment before changing transport code. Check
implant windows and field transfer if it appears during implantation or a
topology change; check interface flux and thermal transport if it appears
during annealing. Reproduce that segment in a planar poly/oxide/silicon stack
with the same thicknesses, coefficients and thermal history. Independently
refine the oxide mesh and halve time steps. Define comparison tolerances before
evaluating the fix. A model-predicted penetration problem needs a documented
model/recipe explanation, not an arbitrary concentration clamp.

Exit evidence: identified source/path, inventory balance, converging profiles,
and a minimized failing-control regression for any engine correction. Require
correct body/source/drain polarity across declared sampling regions and junction
checks, not just the current isolated probes. Current status: **FAIL** for fresh
NMOS body polarity. The archived PMOS heavy-implant polarity regression passes,
but it is not a fresh-flow or transport-conservation result. The first
unexpected NMOS phosphorus increment is bracketed during poly activation;
its source, material/interface balance and mesh/time convergence remain **OPEN**.

### 3. Complete fresh flows, then test robustness

Run fresh NMOS and PMOS from initialization through all 16 checkpoints to
`metal1` on one recorded source/build/data snapshot. Audit each checkpoint for
geometry and fields; verify intended material coverage, dielectric separation,
contact connectivity and dopant profiles in addition to exit status. Extend
the audit where needed: the current independent partition checker alone does
not establish all those properties.

Then execute the reliability plan's work-package-6 matrix: actual 0.5x/1x/2x/8x
point counts, local interface refinement, individual mask shifts of ±10/25 nm,
origin translations, ±1% film changes, time-step refinement and restart
comparisons. Record actual counts, memory and elapsed time. Coordinate-scale
unit tests cannot substitute for these runs. Include the remaining geometry,
implant accounting, rollback/failure-injection, conformal-offset and full
viscous-oxidation gates in work packages 1–5 before release.

Exit evidence: fresh nominal flows and required variants pass the declared
geometry, conservation, polarity and convergence criteria; unsupported
resolution fails explicitly without accepting corrupt output. Current status:
**FAIL** for fresh nominal NMOS and PMOS on `native64_n0paar46`. The newer build
passes a captured Kooi replay but its PMOS continuation fails at spacer etch;
a fresh run from initialization has **NOT RUN** on that build. The robustness
matrix is **NOT RUN** on an accepted fresh-flow candidate. PMOS spacer clipping,
ILD boundary construction and NMOS polarity are current blockers.

### 4. Preserve evidence and promote only an accepted build

Archive logs, decks, relevant input/output structures and manifests under a
persistent results location before temporary folders disappear. The persistent
fixture and JSON summary are useful but do not replace all referenced `/tmp`
artifacts. Every new fix needs a failing control and passing regression tied
to its source/build identifiers. Keep historical failures in the report and
add new outcomes with provenance rather than overwriting them as successes.

Only after the release gates pass, switch the launcher to the accepted build,
verify its binary/data resolution and run a launcher smoke case. Preserve the
previous build for rollback. Report engineering reliability separately from
physical-model calibration; process completion does not validate transistor I–V.
Current status: **NOT RUN**; launcher promotion remains withheld.

## Current blocking evidence and next repair units

1. **Fresh PMOS side boundary and spacer clipping:** the Kooi arithmetic
   failure is fixed and captured in `tests/fixtures/rit_pmos_kooi/`. Start from
   `tests/fixtures/rit_pmos_spacer_clip/`, the saved LDD predecessor, and
   `sol_back_pedal_spacer_diag/pmos/suprem4.log`. Trace left-side coordinate
   drift from LOCOS and the repeated near-zero intersections at spacer etch.
   Encode the side boundary as a persistent geometric constraint and preserve
   canonical point/intersection identity through motion and clipping. Require
   one unambiguous directed boundary arrangement, independent partition and
   coverage checks, and a fresh PMOS replay beyond spacers. Do not pick a
   successor by index or snap coordinates by an arbitrary tolerance.
2. **PMOS/NMOS ILD:** use
   `suprem4/tests/fixtures/rit_pmos_ild_boundary/` as a failing control. The
   generated 1,305-vertex offset has a proper crossing between segments 746 and
   710. Construct an intersection-noded planar arrangement, classify deposited
   faces, and emit valid simple loops or explicit holes before committing any
   database mutation. Require independent area/coverage and interface checks,
   the captured ILD replay, and fresh-flow replays for both devices. A larger
   tolerance or segment deletion is not an accepted fix.
3. **Fresh NMOS doping:** trace chemical phosphorus and boron inventories and
   boundary fluxes through the first poly activation increase, with the same
   transport-volume operator used by diffusion. The latest fresh NMOS sample
   has net donor 9.9083e19 cm^-3 at source/drain and 8.9629e18 cm^-3 after
   anneal. Determine whether the increment comes from a model, interface flux,
   or remap defect before changing code or recipe. Require a failing control,
   balance closure, mesh/time convergence and correct declared body polarity.

The bounded activation and Kooi passes and geometry audits establish useful local progress,
not full RIT CMOS reliability. After these repair units pass, rerun both
unchanged fresh recipes to `metal1` and the density, mask, film and restart
matrix above. The launcher remains on its prior validated build.
