# Construction geometry and RIT CMOS validation, 2026-11-16

**2026-11-17 update:** the gate-oxide shielding and missing retained-arc repairs now pass both fresh 16-checkpoint processes. See [current repair and limits](rit_cmos_gate_oxide_shielding.md). Earlier missing-oxide structures and withdrawn acceptance below are historical evidence.

**Electrical follow-up: final acceptance withdrawn.** The contact etch removes
buried gate oxide in both devices. Process completion and the earlier checks
missed this material void. See [actual electrical/PDF comparison](rit_cmos_electrical_pdf_comparison.md).
The validated RIT pin is removed and the acceptance manifest is marked false;
completion claims below describe the earlier, insufficient validation.

This is the implementation follow-up to
[the TCAD research](tcad_repair_methods_research_20261116.md) and
[the planar cancellation defect](deposition_etch_construction_identity.md).
Both fresh CMOS cuts now complete all 16 checkpoints through first metal
on the final construction/traversal/overlap-noding build. All 32 saved stages
pass independent partition, chemistry, topology and applicable device checks.
The accepted build/options and evidence are recorded in
[the validation manifest](process_surface_validation_20261116.json).

## Geometry repair

`SUPREM4_PROCESS_SURFACE=1`, with `SUPREM4_DEPOSITION_GEOMETRY=domain`, retains
a construction-precision representation of the last deposited free surface,
separately from the binary32 simulation mesh. The domain backend evaluates
this curve on a finer dyadic grid. Horizontal facets retain their analytic
source-plane offset through integer generating-feature identity. The existing
coarser mesh backend and its chord-error certificate are unchanged.

Dry etching translates the retained construction, then projects the cutter
once into stored mesh coordinates. Deposition and etching now use the same
double-precision conversion of recipe thickness to cm. A complete coordinate
snapshot and exposed material ownership guard the descriptor: any altered
surface refuses its use. Hexadecimal binary64 construction records in saved
structures preserve it through restart without a cm/micron round-trip loss.

This is deliberately limited to the last free-surface construction. It is not
a complete CAD representation of buried interfaces, and later depositions
currently derive from the current mesh surface. A stale descriptor safely
falls back to the mesh path, which does not guarantee analytic cancellation.
The option remains explicit; acceptance is scoped to the documented process cuts and controls.

No finite region is removed by area, no dry-etch depth is enlarged, and no
species concentration is reset. A sloped conformal film must retain its real
residual after an equal nominal vertical etch.

## Etch representation and ownership repairs

Fine construction curves can have adjacent vertices that become identical
when projected to binary32. Creating an edge between them manufactured a
zero-length cutter edge in the current NMOS spacer case. Etching now
canonicalizes exact adjacent/closing coordinate identities before creating
its skeleton. Every nonzero edge is retained, however small. Vertex spacing
constraints are combined using the smallest positive request.

The same larger cutter exposed existing memory ownership errors in
`user_etch`: polygon growth reallocated the y-coordinate pointer as the
spacing buffer, `done` failed to grow spacing, and completed polygon/line
etches leaked their spacing array into later dry etches. Each completed
operation now frees and clears its buffers; all three arrays grow together;
dry etches have no spacing array. Invalid start/continue/done sequences are
rejected before touching buffers.

## Bounded evidence

Candidate `native64_jbar3ce0`, sanitized counterpart `native64_fiwjr_56`:

* `experiments/process_surface_20261116/canonical_extended` and
  `extended_sanitized`: all 29 construction cases pass, including translations,
  12/100 micron widths, unequal/equal depths, separate restarts, stale-cache
  rejection, retained sloped film, independent partition, and substrate boron
  inventory.
* `experiments/process_surface_20261116/etch_buffers` and
  `etch_buffers_sanitized`: 501- and 1003-vertex masks followed by a dry etch
  pass. `etch_buffers_old_control` fails on the pre-repair engine.
* `decks/subcmos2014/results/astra_construction_spacer_canonical` and
  `astra_construction_spacer_sanitized`: the newly captured NMOS spacer
  reproducer passes, including independent partition audit.
* `tests/run_domain_offset.sh`: existing mesh construction tests plus 81
  analytic nominal-plane cancellation cases pass under sanitizers.

## Native-oxide boundary traversal repair

A fresh PMOS run exposed another root cause after stripping the well-mask
oxide. `NativeOxide` began its scan at the skeleton list head even if that
head was in the middle of an exposed run. The skipped leading edge retained
an exposed silicon point without a supporting exposed face after the next
vertex was coated. The scan now finds the actual beginning of the maximal
exposed run and counts its edges. A completely exposed loop receives one
vertex operation per edge; open runs receive the additional endpoint.
This also removes the old unbounded end search on a fully exposed loop.

`astra_native_seed_runfix` passes the captured fresh PMOS reproducer on
`native64_uu2_2yx9`. `native_closed_tiny` and `native_closed_sanitized` pass
an independently checked, fully exposed silicon island. The old engine
hangs until the 20-second test timeout (`native_closed_old_control`).
The first attempted closed-loop test used a forbidden zero process time;
its failed result is retained in `native_closed`. The accepted test uses
1e-9 minutes, checks command errors, and requires every silicon edge coated.

The focused P1 mass test originally depended on an external environment
setting. It now selects `SUPREM4_MASS_LUMPING=p1` itself. Its original
failed invocation and explicit-mode rerun are both retained in
`experiments/process_surface_20261116/focused`; the engine default has not
changed. All 16 focused checks have successful evidence after that correction.

## Noding an endpoint contact and an interior crossing

The +2 nm interior perturbation of the captured NMOS spacer exposed a
collinear-overlap case. A source edge starts at a clip corner, overlaps a
straight clip segment, then exits at a different corner in its interior.
The old prepass counted only proper interior crossings, so this appeared to
have one crossing. The arc walker selected the starting endpoint, retained
the overlapping edge, then closed back through it. The first invalid geometry
was already present at Boolean output; triangulation correctly refused it.

The prepass now recognizes an exact existing endpoint contact plus a proper
interior intersection as multiple contacts. It splits both perimeters at
the shared interior intersection before classifying arcs. Endpoint contact
uses the robust stored-coordinate orientation predicate and exact interval
bounds, with no new proximity band. No edge is deleted after clipping.

Two new analytic overlap tests check component count, signed area, and
absence of collinear reversals, including rotation. The pre-repair source
fails their component-count assertion; the correction passes them and all
70 existing rounded shared-cut variants. Control source/test hashes and
logs are in `experiments/process_surface_20261116/overlap_noding_control.json`.

Final candidates: ordinary `native64_ducjxrq7`, ASan/UBSan `native64_rgomhyjr`.

* `overlap_perturbations_inventory`: all four ±2 nm interior-geometry cases
  pass, with reflecting sides fixed. Both NMOS spacer cases retain silicon
  boron and phosphorus inventories exactly at saved-file precision. This
  perturbs the represented geometry, not the complete mask recipe.
* `overlap_focused`: all seven captured RIT cases pass under sanitizers.
* `overlap_extended_sanitized`: all 29 construction/restart/retained-film
  controls pass.
* `overlap_buffers_sanitized` and `overlap_closed_sanitized` pass.

The RIT runner now performs independent partition and nonnegative finite B/P
checks at every checkpoint, in addition to its topology, finite-field,
polarity, dielectric and contact checks. Chemical inventories are reported;
they are not asserted constant across physical implantation, etching or
transport steps. Rezoning-only tests separately assert conservation.

## Historical process-completion validation (acceptance withdrawn)

`decks/subcmos2014/results/astra_surface_overlap_fresh` runs both unchanged
generated decks from initialization on `native64_ducjxrq7`, using:

```
SUPREM4_PROCESS_SURFACE=1
SUPREM4_DEPOSITION_GEOMETRY=domain
SUPREM4_MATERIAL_DOMAINS=oxide
SUPREM4_MESH_PATCH_REMOVAL=1
```

No alternate time integrator, mass matrix, positivity predictor or vertical
constraint option was enabled. Both runs exit zero with no runner errors:

| Device | Checkpoints | Wall time | Gate dielectric probe | Metal islands |
| --- | --- | --- | --- | --- |
| NMOS | 16/16 | 588.56 s | 0.0166745 microns | 4 |
| PMOS | 16/16 | 552.65 s | 0.0169701 microns | 4 |

Both pass well, poly, source/drain and body-tie polarity, gate width, remaining
gate dielectric, and all four contact-center gap checks. Every saved stage
passes independent partition and finite/nonnegative B/P checks. Inventories
are recorded as process results, not assumed conserved across physical steps.
`results.json` contains the original complete process results;
`complete_checkpoint_audit.json` applies the current strengthened checkpoint
checker to all 32 outputs. Earlier successful continuations and the initial
fresh pair remain archived separately; acceptance uses this final fresh pair.

The former workspace pin `suprem4/build/validated-rit-cmos` selected this build;
it has now been removed following the electrical material-occupancy failure.
Use the four options above; [the deck README](../../decks/subcmos2014/README.md)
gives a directly runnable command. The general launcher remains on
`native64_2r55t88h` and its defaults are unchanged.

This validates these process cuts and bounded geometry perturbations. Electrical
calibration, a mesh-convergence study, the full-viscous variant, and arbitrary
large decks remain separate work. ASan/UBSan tests use leak detection disabled
for the legacy engine; the complete fresh flows use the ordinary candidate.
