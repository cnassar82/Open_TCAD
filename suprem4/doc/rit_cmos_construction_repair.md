# RIT CMOS construction repair, 2026-11-15

**2026-11-17 update:** the gate-oxide shielding and missing retained-arc repairs now pass both fresh 16-checkpoint processes. See [current repair and limits](rit_cmos_gate_oxide_shielding.md). Earlier missing-oxide structures and withdrawn acceptance below are historical evidence.

**Electrical follow-up: final acceptance withdrawn.** The contact etch removes
buried gate oxide in both devices. Process completion and the earlier checks
missed this material void. See [actual electrical/PDF comparison](rit_cmos_electrical_pdf_comparison.md).
The validated RIT pin is removed and the acceptance manifest is marked false;
completion claims below describe the earlier, insufficient validation.

2026-11-16 update: both fresh decks complete all 16 checkpoints on the
construction candidate. See [implementation and validation](process_surface_construction_repair.md)
for current evidence and remaining scope. The investigations below are historical.

Full fresh acceptance is still being tested. Do not equate the completed LOCOS
replay or captured spacer success with a completed CMOS process.

## NMOS spacer: a constructed clipping corner lost its identity

`astra_embedding_spacer` brackets the first invalid skeleton at `triang input`,
immediately after `sub_2dskel`. The clipping result already reverses along the
right reflecting side; triangulation does not introduce the defect.

`classified_crossing` recomputes an intersection at an actual clipping vertex,
then substitutes a nearby subject endpoint because that endpoint is classified
ON by the supporting-line roundoff allowance. The clip walker subsequently
visits the actual corner before returning to the substituted endpoint. That
creates a reversed collinear boundary edge.

The repair preserves the constructed intersection when it is an existing clip
vertex. Edge-interior rounded shared intersections retain their established
junction policy. No fuzzy classification band, winding check or triangle
admissibility test was relaxed.

`tests/etch_corner.c::check_constructed_corner` checks the captured geometry,
half/double scale and a 90-degree rotation. It verifies no reversing collinear
triple and the analytic retained area. The old source fails this assertion;
the repair also passes all 70 existing rounded shared-cut cases.

The frozen full spacer replay passes with ordinary `native64_v1173sxh` and
ASan/UBSan `native64_bru9u5y4`: 17,769 triangles pass independent partition
checks. A later replay with the simplified deposition construction and final
sanitizer candidate `native64_j5pmyekc` also passes. Leak detection was disabled
for these sanitizer runs.

## PMOS LOCOS: redundant offset edges precede the crossing

The mutation trace in `astra_first_crossing_locos/pmos` locates the first
proper boundary crossing after `mesh_trial_move`, before geometry/edge flips.
All local triangle orientation tests pass. The corresponding fresh trace,
`astra_first_crossing_pmos`, reproduces the same event.

`astra_unrounded_motion_locos/pmos/embedding_motion.dat` records the old
coordinates, point velocities, unrounded candidates and dt=27.675354722035426.
The crossing exists in the unrounded candidate too. It is not solely a float
storage artifact. The LOCOS input already contains a one-grid backward edge
on the deposited nitride surface near x=5.58284 microns. Subsequent smooth
deformation turns that redundant join into a loop.

The [Clipper2 offset documentation](https://www.angusj.com/clipper2/Docs/Units/Clipper.Offset/Classes/ClipperOffset/_Body.htm)
explicitly describes redundant short segments as a source of offset blemishes
and recommends simplification between offset operations. The deposition bridge
had omitted that construction stage.

`domain_offset.cpp` now allocates half the existing chord-error budget to round
arc construction and uses a quarter-budget simplification tolerance on the
newly grown domain. It independently measures every removed polygonal
subchain against its retained finite chord and refuses simplification exceeding
the remaining half-budget. The original occupied domain, original interface
coordinates and input fields are not simplified.

The captured 143-point pad-oxide surface is frozen in
`tests/fixtures/rit_pmos_offset_surface.dat`. Its dilation must remain a
single-valued surface. The old backend produces a backward segment and fails;
the repaired backend passes. Measured simplification error is
3.8718992315913946e-6 microns against a total chord budget of approximately
1.5e-5 microns. Analytic planar/step areas, translations, side constraints and
binary32 representability checks also pass under ASan/UBSan.

`astra_offset_budget_locos/pmos` completes from the unchanged wells checkpoint
through LOCOS: 6/6 runner checkpoint checks pass. A separate partition audit
passes 7,646 triangles and 380 boundary edges. This is a resumed result.

## Current integration evidence

Ordinary certified candidate: `native64_jjf7h6aw`. Sanitizer counterpart:
`native64_j5pmyekc`. Runtime geometry options:

```
SUPREM4_DEPOSITION_GEOMETRY=domain
SUPREM4_MATERIAL_DOMAINS=oxide
SUPREM4_MESH_PATCH_REMOVAL=1
```

The preceding fresh ordinary run `astra_offset_budget_fresh/nmos` reaches
source/drain implantation and passes all three source/drain/body polarity
checks. It then exhausts the unchanged retry budget at the first activation
nitrogen step. `astra_offset_activation_diagnosis` identifies a negative
phosphorus update in a tiny oxide control volume. An opt-in P1 mass replay
passes that startup but fails later oxidation and is not accepted.

Fresh certified legacy-mass runs are in `astra_certified_offset_fresh`; fresh
P1-mass comparisons are in `astra_p1_offset_fresh`. The legacy NMOS fails
activation and PMOS times out. P1 PMOS fails LOCOS; P1 NMOS times out and
fails the body-polarity check. None is a passed full-process result. The
launcher has not been promoted.

Subsequent work is recorded in [conservative_connectivity_repair.md](conservative_connectivity_repair.md).
It isolates a chemical-inventory loss in automatic diagonal swaps after a
shared-interface insertion. That bounded transfer now passes ordinary and
sanitized replay checks, but both resumed devices still fail activation.

The 13 focused numerical and geometry regressions pass. New clipping and
offset tests and their failing controls are recorded separately under
`experiments/corner_identity_20261115`.
