# RIT Advanced CMOS 2014 (0.25 um lambda): first full flow (2026-11-21)

Decks: `decks/advcmos2014` (README lists PDF pages and assumptions). Build `native64_u85__kd9`.
Flow `decks/advcmos2014/results/full_20261121`: NMOS and PMOS PASS 17/17.
Genius `suprem4/experiments/electrical_advcmos_20261121` (6/6 sweeps; n+ poly 4.17 eV, p+ poly 5.29 eV, QF = 0).

## Oxide calibration
| oxide | PDF target | recipe as given | calibrated |
|---|---|---|---|
| pad (recipe 250) | 500 A | 414 A (40 min) | 476 A (51 min) |
| recessed, each | 3650 A | (no recipe) | 3689 / 3554 A (1100 C wet, 22 min) |
| gate (recipe 213, N2O as N2) | 100 A | 144 A (30 min O2) | ~100 A as grown (12 min); 110 A after reox/anneal |

## Electrical
| | NMOS | PMOS | RIT Athena/Atlas (p.47-48) |
|---|---|---|---|
| Vt 100 nA, |Vd| = 0.05 V | 0.859 V | -0.457 V | 0.73 / -0.74 V |
| Vt max-gm | 0.906 V | -0.456 V | |
| SS | 103 mV/dec | 86 mV/dec | |
| DIBL | 10.9 mV/V | - | |
| Id(Vg=0, |Vd|=3.3) | 2.2e-14 A | 7.5e-6 A | |
| Ion | 253 uA | 260 uA | |

## Profiles versus RIT's Athena cuts (p.47-48, read from the plots)
- NMOS channel boron 4.2e17 (50-100 nm), 3.4e17 at 0.4 um; RIT ~4-5e17 flat. Agrees.
- PMOS: surface 2e17 after the P Vt adjust (RIT ~2e17), but the n-well below is 5e16 at 0.5 um and 3.4e16
  at 0.8 um; RIT ~1e17 at 0.5 um and ~7e16 at 1 um. The simulated n-well is ~2x lighter.
- Cause: the PDF implants both wells through the 1500 A nitride + 480 A pad oxide (steps 10-16), and the
  nitride is then stripped with the dopant it stopped. Gibbons tables: P 170 keV in Si3N4 Rp 0.140 um,
  dRp 0.041 um, so the peak sits inside the 0.150 um nitride. Retained in Si: P 6.1e12 of 3e13 (20 %),
  B 4.0e13 of 8e13 (50 %). Retention is a knife-edge on the nitride range table.
- PMOS consequence: S/D (B xj 0.79 um) approach to 0.155 um at 0.2 um depth under the shallow Vt-adjust
  layer with only 4e13-7e16 n doping: subsurface punchthrough, 7.5 uA off current at Vd = -3.3 V, and
  |Vt| reduced by short-channel roll-off (long-channel hand estimate -0.63 V).
- NMOS Vt 0.86 V vs RIT 0.73 V with matching doping: Leff 0.31 um, xj 0.49 um, 110 A oxide; extraction
  definitions differ (RIT method not stated).

## Layout effect on the n-well (2026-11-21)
The gate is the PDF minimum (Lpoly 0.5 um); it is formed after the well drive and cannot affect the wells.
The 3 um device active between recessed-oxide fields does: the field oxide blocks most of the 170 keV P,
and the 6 h 1100 C drive (2*sqrt(Dt) ~1.2 um) spreads the well laterally. P at x = 3.0 um after the drive:

| depth (um) | 0.02 | 0.2 | 0.5 | 1.0 |
|---|---|---|---|---|
| layout (3 um active) | 6.3e16 | 6.0e16 | 4.9e16 | 2.7e16 |
| no field oxide (1D) | 8.3e16 | 8.0e16 | 6.7e16 | 3.7e16 |
| RIT Athena (p.48, read) | ~1e17 | ~1e17 | ~9e16 | ~7e16 |

1D retained P in Si: 8.3e12 of 3e13. The layout costs ~25 %; the nitride costs the rest.

## Wide device active (7 um, 1.5-8.5 um; gate 4.75-5.25 um) — results/full_wide_20261121
Both decks PASS 17/17; Genius 6/6 (suprem4/experiments/electrical_advcmos_wide_20261121).

| | NMOS 3 um | NMOS 7 um | PMOS 3 um | PMOS 7 um | RIT |
|---|---|---|---|---|---|
| Vt_cc 100 nA | 0.859 | 0.984 | -0.457 | -0.518 | 0.73 / -0.74 |
| Vt max-gm | 0.906 | 1.034 | -0.456 | -0.517 | |
| SS (mV/dec) | 103 | 108 | 86 | 88 | |
| DIBL (mV/V) | 10.9 | 9.8 | - | 86.8 | |
| Id(Vg=0, abs Vd=3.3) | 2.2e-14 | 6.2e-14 | 7.5e-6 | 1.5e-8 | |
| Ion (A/um) | 253e-6 | 213e-6 | 260e-6 | 228e-6 | |

Well at the channel after the drive (abs net, 0.02/0.2/0.5/1.0 um): PMOS 8.1/7.7/6.4/3.5e16 (= 1D case);
NMOS 3.9/4.9/4.2/2.5e17. Metallurgical channel: NMOS 0.33 um at 20 nm, 0.56 um at 0.2 um; PMOS 0.25 um
at 20 nm, 0.17 um at 0.2 um (subsurface approach remains).
NMOS hand Vt with RIT's own ~4.5e17, n+ poly, 110 A, no Qf: ~1.0 V, so RIT's 0.73 V implies positive
fixed charge (~5e11 cm^-2 would move both devices about -0.27 V). The 0.5 um PDF Vt sheet uses Nss 3e11.

### Fixed oxide charge (Genius only, same wide structures; QF on all Si/SiO2 interfaces)
suprem4/experiments/electrical_advcmos_wide_qf{3e11,5e11}_20261121, 12/12 runs exit 0.

| QF (cm^-2) | NMOS Vt_cc | NMOS max-gm | PMOS Vt_cc | PMOS max-gm | PMOS DIBL | PMOS Id(Vg=0,Vd=-3.3) |
|---|---|---|---|---|---|---|
| 0 | 0.984 | 1.034 | -0.518 | -0.517 | 86.8 | 1.5e-8 |
| 3e11 | 0.817 | 0.867 | -0.675 | -0.671 | 83.8 | 8.9e-9 |
| 5e11 | 0.710 | 0.759 | -0.780 | -0.775 | 83.4 | 6.3e-9 |
| RIT | 0.73 | | -0.74 | | | |

QF = 5e11 reproduces RIT within 0.02-0.04 V for both devices (slope ~-0.055 V per 1e11, both devices).
Qf does not fix PMOS DIBL (subsurface, set by the n-well dose).

## N-well phosphorus dose balance (2026-11-22)
suprem4/experiments/nwell_dose_balance_20261122: restart from the wide PMOS recessed_oxide.str, P 3e13 170 keV,
nitride strip, recipe-11 drive; 1D integrals (cm^-2) at the channel x = 5.0 um.

| stage | nitride | pad oxide | silicon |
|---|---|---|---|
| after implant | 1.80e13 (60 %) | 7.1e12 (24 %) | 5.1e12 (17 %) |
| after nitride strip | removed | 7.1e12 | 5.1e12 |
| after 6 h drive | - | 0.17e12 | 8.07e12 (27 %) |

During the drive 3.0e12 moves from the pad oxide into Si and ~3.9e12 leaves through the oxide surface
(modelrc phosphorus oxide/gas Trn.0 = 1.66e-7 cm/s, i.e. a perfect sink over 6 h; P oxide sqrt(Dt) ~28 nm
vs 48 nm pad oxide). Silicon result equals the full-flow wells.str. Control (wells implanted after the strip):
3.0e13 in Si, fully conserved through the drive. RIT's ~1e17 n-well corresponds to ~1.1-1.2e13 retained,
about what results if the pad-oxide P were kept (12.2e12): consistent with a sealed oxide surface and/or a
slightly deeper nitride range in Athena.

### Oxide/ambient P transport versus TSUPREM-4
TSUPREM-4 default (App. A, interface transport table): P oxide/ambient TRANS.0 = 9.0e5 um/min = 1.5 cm/s,
TRANS.E = 1.99 eV (7.4e-8 cm/s at 1100 C); ours is 1.66e-7 cm/s constant. The loss is diffusion-limited
(D_ox/h ~0.05 nm << 48 nm), so the coefficient does not matter. Rerun of the through-nitride case (x = 5.0):

| oxide/gas transport | P in Si | P in pad oxide | P at 0.02/0.2/0.5/1.0 um |
|---|---|---|---|
| ours 1.66e-7 | 8.07e12 | 0.17e12 | 8.2/7.9/6.5/3.6e16 |
| TSUPREM-4 1.5*exp(-1.99/kT) | 8.07e12 | 0.17e12 | identical |
| sealed (0) | 9.61e12 | 2.46e12 (lost at pad etch) | 11.4/10.2/7.8/3.9e16 |
| RIT Athena (read) | | | ~1e17/~1e17/~9e16/~7e16 |

No model change: the setting already behaves like TSUPREM-4. Even sealed, the deep tail stays ~2x below
RIT's at 1 um, pointing at implant depth/dose (range in nitride or Athena structure), not the surface.

## SRIM-2013 check of the nitride range table (2026-11-22)
SRIM-2013 SR Module (srim.org, run under Wine in a no-network sandbox); outputs copied to
suprem4/experiments/nwell_dose_balance_20261122/srim/.

| | Gibbons table Rp/dRp (um) | SRIM Si3N4 3.44 | SRIM 3.1 | SRIM 2.9 | SRIM Si | SRIM SiO2 2.27 |
|---|---|---|---|---|---|---|
| P 170 keV in nitride | 0.1404/0.0411 | 0.1469/0.0394 | 0.1630/0.0437 | 0.1743/0.0468 | 0.2224 (Gibbons 0.2329) | 0.2338 (0.2185) |
| B 80 keV in nitride | 0.1972/0.0504 | 0.1769 | 0.1963/0.0532 | 0.2098 | 0.2609 (0.2473) | 0.2927 (0.2681) |

Gibbons' P-in-nitride row behaves like ~3.4 g/cm3 nitride; its B row like ~3.1. For LPCVD nitride (2.9-3.1)
SRIM puts the 170 keV P peak past the 0.150 um nitride. B is unchanged at 3.1, so the p-well is unaffected.

Rerun of the through-nitride n-well with only the P-in-nitride rows replaced by SRIM (test tables
sup4gs_srimP_nitride_{3p1,2p9}.imp), x = 5.0 um:

| P-in-nitride data | P in nitride at implant | P in Si after drive | P at 0.02/0.2/0.5/1.0 um |
|---|---|---|---|
| Gibbons | 1.80e13 | 8.07e12 | 8.2/7.9/6.5/3.6e16 |
| SRIM 3.1 g/cm3 | 1.18e13 | 1.38e13 | 13.8/13.2/11.1/6.3e16 |
| SRIM 2.9 g/cm3 | 0.94e13 | 1.63e13 | 16.1/15.5/13.1/7.5e16 |
| RIT Athena (read) | | | ~10/~10/~9/~7e16 |

SRIM at 3.1 g/cm3 brings the n-well to RIT's level including the deep tail; the PDF process is unchanged.

## Full flow with the SRIM P-in-nitride table (adopted 2026-11-22), QF 5e11
results/full_srim_20261122 (17/17 both), suprem4/experiments/electrical_advcmos_srim_20261122 (Genius 6/6).

| | NMOS before | NMOS SRIM | PMOS before | PMOS SRIM | RIT |
|---|---|---|---|---|---|
| Vt_cc 100 nA (V) | 0.710 | 0.710 | -0.780 | -0.919 | 0.73 / -0.74 |
| Vt max-gm (V) | 0.759 | 0.759 | -0.775 | -0.911 | |
| SS (mV/dec) | 107.6 | 107.6 | 88.7 | 92.5 | |
| DIBL (mV/V) | 9.6 | 9.6 | 83.4 | 37.7 | |
| Id(Vg=0, abs Vd=3.3) (A) | 1.2e-13 | 1.2e-13 | 6.3e-9 | 2.4e-14 | |
| Ion (A/um) | 258e-6 | 258e-6 | 198e-6 | 158e-6 | |

("before" = wide layout with QF 5e11, Gibbons table.) PMOS channel net donors at x = 5.0 after the D/S
anneal, 0.02/0.1/0.2/0.5/1.0 um: before 2.1e17/1.6e17/9.1e16/6.4e16/3.5e16; SRIM 2.7e17/2.2e17/1.4e17/
1.1e17/6.1e16. Punchthrough is gone; the n-well is now ~30 % above RIT's plotted level, and PMOS Vt is
0.18 V more negative than RIT's.

## Poly re-oxidation (2026-11-22)
RIT's process (PDF p.4): Lpoly 0.50 drawn, 0.40 after etch, 0.35 after the 500 A re-oxidation.
SUPREM-IV.GS never oxidizes poly (every oxidation model reacts only at Si/SiO2; the poly hooks in
oxgrow.c/FEconvert.c were never finished).  TSUPREM-4 oxidizes poly only in COMPRESS/VISCOUS/VISCOELA,
with the <110> linear rate and alpha 0.44.

Physics (implemented in src, compress model only): poly/oxide interface as a reacting sink at the <110>
rate (Oxidant.c, oxide_vel.c), poly consumed and held rigid in elast.c, oxide grid on poly (grid_upd.c).
Test (0.4 um poly line on 100 A oxide, dry 1000 C 51 min): 499 A on poly top, 222 A consumed (0.44 x),
sidewall oxide 500 A, poly 0.400 -> 0.355 um, gate-edge dielectric 270 A tapering over ~0.1 um.

Mesh: the Lagrangian moving mesh cannot follow the topology changes (convex poly corners collapse, the
sidewall and substrate oxides close a notch).  Local fixes reached 5/6 test variants but are not robust;
a Boolean-repair-and-remesh module (refine/oxide_repair.cpp with Clipper2 + vendored CDT 1.4.5,
oxide/oxide_repair.c) works on short runs but is not finished.  Neither is used by the decks. Everything (research, failure catalogue, patches, test decks, plan) is saved in suprem4/doc/poly_oxidation/.

Adopted: geometric re-oxidation in the deck (`build_decks.py --poly-reox geometric`, oxide on poly
default 0.05 um): poly deposited 0.4 - 0.44*t thick and etched to drawn - 2*0.44*t; sidewall oxide by
conformal deposition of t and vertical etch-back of t (the nitride-spacer mechanism); then the normal
recipe-250 oxidation with the vertical model for the substrate and the thermal budget.  Omitted: the
500 A top oxide on poly (only screens the S/D implant into the gate) and the gate-edge bird's beak.
Runs entirely on the unchanged production build (native64_u85__kd9).
