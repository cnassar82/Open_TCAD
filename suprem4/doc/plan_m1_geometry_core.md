# Implementation plan: M1 geometry core in production (draft, 2026-11-24)

Follows M0 (`plan_m0_geometry_core.md`, review `m0_review_wp8.md`: GO for option C, with
conditions).  M0 showed that the core can carry SUPREM-IV through poly re-oxidation and through both
RIT flows.  M1 makes it **the geometry layer SUPREM-IV runs on by default**, for every step that
changes geometry, at a cost that does not make flows impractical.  That is the foundation the
roadmap's Phase 2 (materials, mechanics, physics port) and Phase 3 (topography) build on.

Not in M1: porting SUPREM's physics to a new solver (Phase 2 proper), new materials, mechanics.
Junction physics (roadmap Phase 1) proceeds in SUPREM-IV in parallel and is unaffected: it runs on
whatever mesh the core hands over.

## 1. Starting point (measured in M0)
| Item | State |
|---|---|
| Oxidation (compress, viscous) | Core runs on every step (Lagrangian update, local events, regeneration as fallback), opt-in via `SUPREM4_GEOMETRY_CORE=1` |
| Suite | 58/58 after the M1.1 fix (2026-11-24) |
| RIT flows | Core off bit-identical to `main`.  Core restricted to compress/viscous (M1.0a): 0.25 um flow byte-identical to production; 0.5 um differs from LOCOS on, Vt 13-20 mV, Ion < 3 %, SS < 5 %, DIBL 5.9 % (`m0_review_wp8.md` 5b) |
| Runtime | 2.6x SUPREM's vertical model per step, 12.8x per run (5x more steps in compress) |
| Etch / deposit | SUPREM's own grid code (Clipper2 backend for deposition); source of sub-angstrom slivers at spacer stages |
| Exchange | File write + reload per call (`.str` in `SUPREM4_CORE_TMP`) |
| Mesh quality | 1-6 obtuse boundary edges per call at poly corners, inherited from SUPREM and kept |

## 2. Principles (unchanged from M0, plus two)
1. Geometry is the source of truth (B-rep); the mesh is derived.
2. Every mesh handed to SUPREM passes the validity contract; failures stop with a dump.
3. Core off stays bit-identical to production until the default is switched (WP6), and remains
   available afterwards as the reference.
4. **New:** geometry below the event resolution (delta = 1 nm, tolerance delta/10) never reaches a
   mesh.  Every input (events, etch, deposit, SUPREM's own motion) goes through the same clean-up.
5. **New:** cost is measured per flow, not per step.  A change that is cheaper per step but needs
   more steps is not an improvement.

## 3. Work packages

Status 2026-11-25: M1.0(a), M1.0(b) and M1.1 done.  M1.0 result (`doc/m1_0_core_on_off.md`): the remaining core-on vs core-off differences are discretization errors of SUPREM's vertical-model steps (gate-oxide time step, near-surface P pile-up mesh), not the core; they vanish when those steps are converged.  The production RIT decks are themselves unconverged (gate oxide 8 A too thick); that is a Phase 1 item.  Next: M1.2.

| WP | Content | Deliverable | Effort |
|---|---|---|---|
| M1.0 | **Core on vs core off in the RIT flows.**  The core currently runs on vertical-model oxidations too (never validated).  (a) Restrict the hook to the models it was validated for, and measure what remains.  (b) For the vertical-model steps, establish the converged answer: mesh-refinement study of the core-off reference itself (pad, field, gate and re-oxidation; dose after each anneal), then decide whether the core should run there.  (c) Trace each remaining difference to its step. | A table of every process step with core-on vs core-off difference and the converged value; the electrical differences explained | 4 d |
| M1.1 | **M0 closure.**  Sub-resolution clean-up after events and before every mesh build: morphological close/open of each material at delta/10 on the integer grid, applied to the union so that slits and spikes below 1 A vanish without moving resolved geometry.  Regression case from the rand7_35 dump (replay of call 96). | rand7_35 completes; suite 58/58; no change > 0.1 % in any metric of the 57 passing cases | 2 d |
| M1.2 | **Runtime: measure first.**  Per-flow profile of where wall time goes (SUPREM solves, core update, file I/O, reload/`create_db`, symbolic factorisation) and *why compress takes 5x the steps of vertical*: SUPREM's step control (area loss, collision prediction, retry after mesh change) vs the physics.  Comparison with plain compress where it completes (LOCOS). | Profile note with the budget per component, base deck and both RIT flows | 2 d |
| M1.3 | **Runtime: act.**  Expected levers, confirmed or replaced by M1.2: (a) in-memory exchange through the C API instead of write/reload (keeps SUPREM's factorisation when the topology is unchanged); (b) step control that does not shrink after a core update that changed nothing SUPREM's limits depend on; (c) incremental update cost (only changed region of the mesh). | Per-run cost <= 3x SUPREM's vertical model on the poly matrix, and not slower than core off on the RIT flows | 5 d |
| M1.4 | **Boundary quality.**  Remove the obtuse boundary edges at corners that the core currently keeps (split or relocate within the B-rep tolerance), so that no negative box coupling reaches SUPREM.  Dose drift then has no geometric contribution. | 0 negative couplings after every call on the suite; dose drift <= plain SUPREM flat drift in every case | 2 d |
| M1.5 | **Etch and deposit through the core.**  After every `etch` / `deposit` the core rebuilds the affected region from the B-rep (exact Booleans on Clipper2's grid, M1.1 clean-up), instead of SUPREM's grid code splicing triangles.  Removes the spacer-stage slivers at the source; same contract and transfer as oxidation. | Both RIT flows with the core on for all geometry steps pass their audits; no structure in either flow has an edge < 1 A | 5 d |
| M1.6 | **Default on.**  Decision gate, with data: process metrics (junction depths, oxide thicknesses, poly CDs) and electricals (Vt, SS, DIBL, Ion, Ioff) core on vs core off for both RIT flows, both devices, with tolerances agreed up front (section 4).  If within tolerance, switch the default; re-baseline the RIT reference results; keep `SUPREM4_GEOMETRY_CORE=0` as the legacy path. | Default switched, or a list of what blocks it | 2 d |
| M1.7 | **Suite and regression.**  Extend `tools/run_geometry_suite.py`: etch/deposit cases (spacer, recess, narrow trench), LOCOS and recessed oxide at several pitches, the failure catalogue of `doc/poly_oxidation/README.md`, and the RIT flows with electricals as a nightly job. | One command for the whole M1 acceptance run, with report | 3 d |

Total about 25 working days; plan 5-6 weeks.  Order: M1.0 -> M1.1 -> M1.2 -> (M1.3 || M1.4) -> M1.5 ->
M1.7 -> M1.6.  M1.0 comes first: until it is resolved, the core cannot become the default (M1.6), and its outcome may change the scope of M1.3-M1.5.  M1.1 is next: every later change reruns the suite.

## 4. Acceptance criteria
1. Suite (M0 matrix, 50 random variants, LOCOS, new etch/deposit cases): all complete, physics checks
   as in M0.
2. RIT flows with the core on for oxidation, etch and deposit: all checkpoint audits pass.
3. Core on vs core off (tolerances to be confirmed before M1.6): oxide thicknesses and poly CD within
   2 %; junction depths within 5 nm; Vt within 20 mV; SS within 5 %; Ion and DIBL within 5 %;
   Ioff within a factor of 2.
4. Runtime per flow <= 3x SUPREM's vertical model on the poly matrix; RIT flows not slower with the
   core than without.
5. No mesh handed to SUPREM has a negative box coupling or an edge below 1 A.
6. Core off remains bit-identical to `main` until the default switch; afterwards the re-baselined
   references are reproduced bit-identically by the default build.

## 5. Risks
| Risk | Mitigation |
|---|---|
| The 5x step count is SUPREM physics (compress stress), not step control | M1.2 decides before M1.3 spends effort; then the criterion is restated against plain compress |
| In-memory exchange touches SUPREM's database code (the M0 file reload was chosen to avoid it) | Keep the file path as fallback and verify the in-memory path bit-identical to it |
| Etch/deposit through the core changes results that the RIT references encode | Expected; quantified in M1.6 with the same tolerances, not hidden |
| Clean-up at delta/10 removes real features (thin films) | Films are resolved far above 1 A; the regression requires < 0.1 % change on passing cases |
| Electrical differences core on vs off turn out larger than the tolerances | M1.6 is a gate: the default switch waits; the difference is traced to the process step that causes it |

## 6. Decisions for the reviewer
1. **Scope of M1**: geometry in production (this plan), or start the Phase 2 physics port (new
   FVM solver) now and leave etch/deposit to later.  Recommendation: this plan.  The physics port
   needs a stable geometry layer under every step, and the RIT flows are the regression for both.
2. **Tolerances** for the default switch (section 4, criterion 3).
3. **Default switch** inside M1, or keep the core opt-in until Phase 2.
