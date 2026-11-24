# M0 geometry core: WP8 go/no-go review (2026-11-24)

Scope: the six acceptance criteria of `doc/plan_m0_geometry_core.md` section 5, measured in WP7 on
branch `geometry-core` (commits `b78ba7f` .. `1096761`).  Decision asked for: proceed with roadmap
option C (a new geometry/mesh core under SUPREM-IV, `doc/roadmap_2d_cmos_25nm.md`) or not.

**Recommendation: GO for option C, with one condition.**  Before M1 starts, close the one failing
case (rand7_35) with the sub-resolution clean-up described below.  It is a known gap in a
specific mechanism (handoff section 4 already lists it), not a design limit.  Every other criterion
is met, or the shortfall comes from SUPREM itself and the core does not make it worse.

| # | criterion | result | verdict |
|---|---|---|---|
| 1 | six decks + 50 random variants complete | 6/6 decks, 49/50 variants, LOCOS 2/2 | **not met: 1 failure** (diagnosed) |
| 2 | physics vs compress reference, mesh convergence | ratio 0.4545 on all planar faces; gate edge within 3 %; refinement within 0.6 % | met |
| 3 | dose drift < 1e-4; areas vs expansion ratio | 41/55 below 1e-4, the other 14 below plain SUPREM's own flat drift; core changes areas < 0.02 % | met as far as SUPREM allows |
| 4 | invariants after every step | every mesh validated before hand-over; 0 refused in all calls | met |
| 5 | core off bit-identical; core on RIT flows pass | 72/72 checkpoint files identical; core on: 0.25 um 17/17 + 17/17, 0.5 um 16/16 + 16/16 | met |
| 6 | poly re-ox step <= 3x SUPREM vertical | 2.6x per step; 12.8x per run | met as worded; per-run cost is the lever |

## 1. Completion (criterion 1)

Harness: `tools/run_geometry_suite.py --matrix --random 50 --seed 7 --locos` (58 cases, 14 h CPU).
Random variants span line width 0.2-0.8 um, poly 0.25-0.5 um, gate oxide 50-200 A, dry
900-1050 C for 10-120 min, and wet 850-925 C for 5-25 min.  Each poly case also runs a flat wafer
of the same recipe with plain SUPREM (core off) as the field and dose reference.  Report:
`experiments/wp7_suite_20261123/report.md` (re-evaluate with `--report`).

57/58 pass.  Over about 4,250 core calls the update path handled all but two calls; those two fell
back to a full regeneration.

**rand7_35 fails** (0.636 um line, 0.276 um poly, 162 A oxide, dry 1050 C, 118.5 min).  SUPREM's
time step collapses at t = 2949 s of 7110 s (retry budget exhausted at dt = 0.017 s).
Diagnosis (deterministic, replayed offline from `SUPREM4_CORE_DUMP` inputs,
`experiments/wp7_physics_20261123/rand7_35_repro`):
- At the right poly bottom corner the gate-edge bird's beak lifts the poly, and a small gas gap
  under the overhang pinches off.  The core fills it (two `gap` events, 4e-7 and 4e-8 um2).
- The fill piece comes from a morphological opening on the integer grid and does not match the old
  oxide surface exactly.  A **hairline gas slit 0.01-0.07 A wide** remains between them, below the
  significance filter, so no event removes it.
- The local repair rejects the resulting B-rep and the call falls back to regeneration.  The
  regeneration meshes the slit faithfully, with 0.13-0.16 A edges and couplings down to -5.9, and
  SUPREM cannot advance on it.
- Ruled out: Mmg (same defect with Mmg off) and the film layer points (a clearance rule on them
  changed nothing).

The fix is the item already in the handoff: geometry below the event resolution must be removed
before meshing.  Concretely, close the material union morphologically at the clean-up tolerance
(delta/10 = 1 A) after events, or fill every gas piece under that width, not only the significant
ones.  Expected effort: small.  Risk: it also touches the cases that pass, so the suite must be
rerun.

## 2. Physics (criterion 2)

Base deck (0.4 um poly on 100 A oxide, dry 1000 C, 51 min) with the core, against the validated
compress reference (`doc/poly_oxidation/README.md`):

| quantity | core | reference / target |
|---|---|---|
| consumption ratio, Si field / poly top / poly sidewalls | 0.4545 / 0.4545 / 0.4545 | 0.44 +- 5 % (SUPREM: 1/2.2) |
| poly height | 0.400 -> 0.355 um | 0.355 um |
| sidewall oxide | 502 A | 500 A |
| gate-edge dielectric at 0.02 / 0.05 / 0.10 um | 280 / 162 / 108 A | 271 / 162 / 108 A (within 15 %) |
| field oxide | 454 A | flat plain SUPREM 454 A |

Mesh convergence: `SUPREM4_CORE_SCALE=0.5` (10,459 vs 5,948 points) changes every metric by less
than 0.6 %, and the `coarse` deck agrees within 1 %.  Over the suite, every planar face of every
case shows a ratio of 0.4545-0.461 (the native oxide is subtracted), and the field oxide is within
2 % of the flat wafer.  Faces shorter than 4x their oxide thickness, such as the top of a narrow line
under 0.11 um of oxide, are 2-D and are reported, not checked (rand7_05: 0.468).

## 3. Conservation (criterion 3)

Dose: boron drift (box measure, as SUPREM's `geometry.c`) over the full oxidation.
- Kept-triangle CDT constraints left non-Delaunay edges at the boundary between kept and new
  triangles.  These negative couplings cost SUPREM dose (5e-4 on a flat wafer).  Fixed with Delaunay
  flips at those edges (`b78ba7f`); the flat wafer then drifts -9.4e-5.
- Suite: 41 of the 55 completed poly cases drift less than 1e-4, and the largest drift with the core is 1.8e-4.  In
  every case above 1e-4, plain SUPREM on a flat wafer of the same recipe drifts more (up to 2.7e-4).
  The floor is SUPREM's own oxidation, not the core.
- LOCOS: viscous +5.4e-5 with the core vs -5.8e-4 without; compress -7.6e-5 vs -8.7e-4.
- Remaining: 1-6 obtuse boundary edges per call at the poly corners.  They are already in SUPREM's
  input and are kept on purpose.

Areas: oxide gained / (Si + poly consumed) ranges 2.16-2.19 for dry cases and 2.08-2.12 for wet
cases (native oxide subtracted).  To attribute this, every core call of two cases was replayed
offline: rand7_33 (wet, 2.10) and rand7_10 (thin dry, 2.37 before the native-oxide correction).
The core's calls change the oxide area by 1.9e-5 um2 of 9.47e-2 grown, and by 3e-7 of 2.1e-2
(< 0.02 %).  The deviation from 2.2 is entirely in SUPREM's motion: the compress model in fast wet
2-D growth, plus the 10 A native oxide on bare poly.  The suite therefore reports the ratio and does
not check it.

## 4. Invariants (criterion 4)

Before any mesh goes back to SUPREM, the hook checks it for orientation, conformity, embedding,
region consistency, surface flags, gas nodes and the SUPREM reload contract (`finish()` in
`core/src/capi.cpp`).  A failing mesh is refused and SUPREM continues with its own.  Across the
suite (about 4,250 calls), the RIT flows with the core on (403 calls at 0.25 um, 363 at 0.5 um),
**no call was refused and no call raised an error**.  The final structures pass `core_str`
(invariants, B-rep) and `core_check_suprem`.

## 5. Production identity and RIT flows (criterion 5)

Core off.  A build with the core linked in must not change production results.  At first it did:
the 0.5 um PMOS aborted at LOCOS and the 0.25 um NMOS crashed.  Cause: the linker resolved
SUPREM's own Clipper2 and predicate symbols to the core's copies, which are compiled differently.
Fix (`0182753`): the core and Mmg are pre-linked into one object that exports only the C API
(`ld -r --force-group-allocation` + `objcopy --keep-global-symbol`), and the core is built with
`-fno-gnu-unique`.  Result: all checkpoints of both flows, both devices, are **byte-identical to
`main`**:

| flow | NMOS | PMOS |
|---|---|---|
| 0.5 um (subcmos2014), 16 checkpoints | 18/18 files identical | 18/18 identical |
| 0.25 um (advcmos2014), 17 checkpoints | 18/18 identical | 18/18 identical |

(`main` built from a clean worktree; data in `experiments/wp7_identity_20261123`.)

Core on (`SUPREM4_GEOMETRY_CORE=1`):
- 0.25 um flow: PASS 17/17 NMOS and 17/17 PMOS, all checkpoint audits.
- 0.5 um flow: PASS 16/16 NMOS and 16/16 PMOS (363 core calls, 0 refused, 0 regenerations).

## 6. Runtime (criterion 6)

Same deck, same machine load: SUPREM's vertical model (core off) takes 93 s in 20 steps; compress
with the core takes 1,194 s in 100 steps.  That is **2.6x per step**, which meets the criterion as
worded, but **12.8x per run**.  The core costs about 0.7 s per call.  Most of the time is SUPREM's
compress solve and the five times larger step count that compress needs.  Plain compress without the
core cannot finish this deck (it stalls at 47 s), so there is no like-for-like plain number.
LOCOS compress wall times are inconclusive: 1,620 s with the core vs 4,241 s without in the suite
(4 parallel jobs), but 909 s vs 549 s when measured alone in WP6.

Levers, in order: the time-step control after mesh changes, an in-memory exchange instead of the
file write/reload per call, and the per-call update cost.

## 7. Decision

Option C is technically viable: an exact-geometry core with Lagrangian updates, local events and
validated hand-over runs real poly re-oxidation and the full RIT flows inside SUPREM-IV.  It
reproduces the validated physics, does not worsen SUPREM's dose conservation, and leaves production
bit-identical when it is off.  The criterion that is not met (one of 50 random variants) traces to
one identified mechanism with a contained fix.

Conditions for M1:
1. Sub-resolution clean-up after events (rand7_35); rerun the suite (expect 58/58).
2. Per-run runtime: treat the step count after mesh changes as an M1 work item; the 3x criterion
   should be restated per run for M1.
