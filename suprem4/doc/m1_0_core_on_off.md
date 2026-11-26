# M1.0: core on vs core off in the RIT flows (2026-11-24/25)

Question (M1 plan, M1.0): why do the RIT devices differ with the geometry core on, and which answer
is right?  Flow: 0.5 um subcmos2014 (the 0.25 um flow has no compress/viscous step and is
byte-identical with the core on since M1.0a).  Build `native64_se0fhkjm`
(`build_native64.py --domain-offset --ordering amd --fe-solver umfpack --geometry-core build/core_wp6a`).
Data: `experiments/m10b_mesh_20261124`, `experiments/wp7_electrical_20261124/{sub_off,sub_m10,gdt_off,gdt_on}`.

## Result

The core is not the source of the remaining differences.  They are discretization errors in
SUPREM's own vertical-model steps after LOCOS.  The core's LOCOS hands those steps a mesh that is
slightly different.  Production's mesh is equally unconverged, so it lands on a different wrong
answer.  When those steps are converged, core on and core off agree.

| step | quantity | core off | core on | converged |
|---|---|---|---|---|
| LOCOS (viscous, core active) | field oxide at x = 2.5 um | 6657 A | 6659 A | - |
| LOCOS | active edges | 3.145 / 8.85 um | 3.145 / 8.855 um | - |
| LOCOS | n-well P at the channel surface | 6.389e16 | 6.388e16 | - |
| gate oxidation (vertical) | thickness, production settings | 149.1 A | 146.4 A | **141 A** |
| gate oxidation | with `oxide.gdt=0.002` | 141.3 A | 141.3 A | 141 A |
| sacrificial (Kooi) oxidation | PMOS P at 10 nm after gate ox, mesh 1 / 0.7 / 0.5 | 1.114 / 1.055 / 0.961e17 | 1.063 / 1.085 / 0.960e17 | not yet (still -9 % from 0.7 to 0.5) |

### Gate oxide: time-step error
- A flat 1-D wafer with the same recipe gives 155.7 A (0.04 um surface mesh, defaults) down to
  141.1 A (0.5 nm mesh, rel.err 1e-7).  The analytic Deal-Grove + thin-regime value with SUPREM's
  coefficients and the 10 A native oxide is 141.0 A.
- The error is in the step size, not the mesh.  SUPREM sizes oxidation steps by the dopant LTE
  (only `method boron rel.err` moves the answer; the `oxida` / `veloc` tolerances do not).  The
  interface limit `oxide.gdt` allows grid.oxide x 0.25 = 100 A of growth per step.  With
  `oxide.gdt=0.002` (0.8 A/step) the flat wafer gives 141.25 A on both a 0.5 nm and a 10 nm
  surface mesh, and both 2-D flows give 141.3 A.
- So the production RIT decks grow the 0.5 um gate oxide **8 A (6 %) too thick** for numerical
  reasons.

### Near-surface phosphorus: mesh error
- n-well phosphorus piles up at the silicon surface during the 1000 A sacrificial oxidation.  The
  production mesh resolves the pile-up with one or two silicon nodes.
- Core on vs off: -5.7 % at 10 nm right after that step, -2.2 % after the final anneal.  The
  mesh-refined runs (all mesh lengths x 0.7, x 0.5) move the value by 14 %, and the on/off gap
  shrinks from 4.6 % to 2.8 % to 0.03 %.

### Electricals (Genius, same protocol), with `oxide.gdt=0.002` after LOCOS

| device | quantity | core off | core on |
|---|---|---|---|
| NMOS | Vt max-gm / cc | 0.383 / 0.372 V | 0.383 / 0.373 V |
| NMOS | Ion, SS, DIBL | 257 uA, 85.5, 6.5 | 258 uA, 85.2, 6.2 |
| PMOS | Vt max-gm / cc | -0.827 / -0.851 V | -0.811 / -0.829 V |
| PMOS | Ion, SS, DIBL | 120 uA, 104.9, 52.6 | 121 uA, 107.6, 53.8 |

NMOS agrees.  The PMOS Vt still differs by 16-22 mV, because the PMOS is a buried-channel device
and responds to the under-resolved phosphorus pile-up.  The production-mesh devices also shift
by -36 mV (NMOS) and +34 mV (PMOS) between production settings and `oxide.gdt` (thinner gate
oxide).  That shift is larger than any core effect.

## Consequences
1. **M1.6 gate:** compare core on vs core off only with converged vertical-model steps (step
   limit on the oxide interface; near-surface silicon mesh fine enough for the pile-up).
   Otherwise the comparison measures SUPREM's discretization noise.  Core on vs off at the 0.5x
   mesh already agrees within 0.1 % on near-surface P and within 0.1 A on the gate oxide.
2. **RIT decks (Phase 1):** the production decks are not converged.  The gate oxide is 8 A too
   thick, and the near-surface well doping after the sacrificial oxide is still moving at 0.5x mesh.
   Fix with a principled step limit: interface motion per step as a fraction of the current
   oxide thickness, not of grid.oxide.  A fixed `oxide.gdt` cannot serve both the 5000 A LOCOS
   and a 150 A gate oxide.  Also refine the silicon mesh at the oxidizing surface.  Then
   re-baseline the electricals.
3. Still open: the converged near-surface P (a finer run beyond 0.5x, or a mesh that follows the
   interface), and the core-on vs off PMOS Vt at that resolution.  Addressed below and in
   `m1_interface_layer.md`.

## Converged reference decks (2026-11-25)
The 0.5 um flow with the vertical-model steps converged (`oxide.gdt` step limit after LOCOS,
near-surface silicon mesh refined to resolve the Kooi pile-up).  Data:
`experiments/m10b_mesh_20261124/conv_{off,on}`, electricals
`experiments/wp7_electrical_20261124/conv_{off,on}` (same Genius protocol).

| device | quantity | core off | core on | production mesh + `oxide.gdt` (off / on) |
|---|---|---|---|---|
| NMOS | Vt max-gm / cc | 0.412 / 0.404 V | failed (see below) | 0.383 / 0.383 V |
| NMOS | Ion, SS, DIBL | 252.9 uA, 86.2, 10.1 | - | 257 / 258 uA |
| PMOS | Vt max-gm / cc | -0.722 / -0.722 V | -0.721 / -0.722 V | -0.827 / -0.811 V |
| PMOS | Ion, SS, DIBL | 118.6 uA, 100.7, 60.6 | 119.3 uA, 100.1, 58.9 | 120 / 121 uA |

- **Core on vs off on converged decks: PMOS Vt within 1 mV, Ion 0.6 %, SS 0.6 %, DIBL 2.8 %.**
  That closes M1.0: the earlier 16-22 mV PMOS gap was the unresolved pile-up.
- Convergence itself moves the devices far more than the core: PMOS Vt by +105 mV, NMOS by
  +29 mV.  The production RIT decks are not converged (Phase 1 re-baseline).
- Caveat: the core-off reference has a one-element grown oxide in every vertical-model step
  (SUPREM never adds oxide nodes), which overstates oxide-to-silicon back-flux in inert anneals
  (`m1_interface_layer.md`: +2.3 % P at 10 nm in the 1-D Kooi test).  Both columns share that bias
  here, since the core only ran on compress/viscous.
- **Core-on NMOS failed** in LOCOS (viscous, build `native64_se0fhkjm`): SUPREM's finite-element
  continuation did not converge, worst coupling -0.379.  Not rerun yet with the current build.
