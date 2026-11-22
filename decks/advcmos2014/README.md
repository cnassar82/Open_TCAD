# RIT Advanced CMOS 2014 NMOS and PMOS (0.25 µm lambda)

Separate 2D transistor process cuts based on `rit_process/AdvCmos2014_0p25um.pdf`
(Lynn Fuller, slides dated August 17, 2014). λ = 0.25 µm and Lpoly = 0.5 µm
(p.4), so both decks use a 0.5 µm drawn poly gate. Targets (p.6): Vdd 3.3 V,
Vt = ±0.75 V, 100 Å gate oxide, n-well and p-well 1e17 cm⁻³.

Generate the decks with `python3 build_decks.py`. Run and audit them with
`python3 run.py --binary <suprem4>`, which uses the same runner and audits as
`decks/subcmos2014`. Export to Genius with `prepare_electrical.py`.

## Process (PDF step list p.20, with the slide that gives each value)

| Step | PDF | Deck |
|---|---|---|
| Substrate | p.6: p-type 10 Ω·cm | B 1.3e15 cm⁻³, (100) |
| Pad oxide 500 Å | p.22–23, recipe 250 (1000 °C dry) | recipe 250, soak 40 → **51 min** (calibrated: 40 min grows 414 Å; 51 min grows 476 Å) |
| Nitride 1500 Å | p.24 | `deposit nitride 0.15` |
| Recessed oxide | p.36–39: 3650 Å thermal oxide, strip, regrow 3650 Å | 1100 °C wet, **22 min** soak each (calibrated: 3689 / 3554 Å); vertical strip of the first oxide |
| N-well | p.20 step 11: P 3e13 170 keV | PMOS deck only, through nitride/pad (active) and recessed oxide (field) |
| P-well | p.20 step 14: B 8e13 80 keV | NMOS deck only |
| Nitride strip | p.40 hot phosphoric | `etch nitride all` |
| Well drive | p.42–43, recipe 11: 1100 °C N₂, 6 h | as given |
| NMOS Vt | p.51: B 3e12 30 keV through pad oxide | as given |
| PMOS Vt | p.53 and p.48: **P 3.5e12 60 keV** | as given; see discrepancy below |
| Pad oxide etch | p.55: 10:1 BOE 45 s | vertical etch 600 Å |
| Gate oxide 100 Å | p.56–58, recipe 213 (900 °C: ramp O₂ 30, N₂O 30, O₂ 30) | N₂O interval as N₂; final O₂ soak 30 → **12 min** (calibrated: 30/10/0 min grow 144/94/63 Å) |
| Poly 4000 Å | p.63, undoped | `deposit poly 0.4` |
| Poly gate etch | p.65–66 | gate centred at x = 5.0 µm; drawn length `ADVCMOS_LPOLY` (default 0.5 µm; RIT as-etched 0.40 µm) |
| Poly re-ox 500 Å | p.74, recipe 250 | same calibrated recipe 250 (vertical model: poly itself is not oxidized); `--poly-reox geometric` adds poly consumption and sidewall oxide, see below |
| p-LDD / n-LDD | p.76: B 4e13 50 keV; p.80: P 4e13 60 keV | as given |
| Nitride spacers | p.83–85: LPCVD 3500 Å, ~0.3 µm wide | deposit 0.35, vertical etch 0.37 |
| N+ / P+ D/S | p.88–92: P 4e15 60 keV; B 4e15 50 keV | gate open in its own D/S window: n+ poly NMOS, **p+ poly PMOS**; opposite implant on the body-contact island |
| D/S anneal | p.95, recipe 284: 1000 °C N₂ 20 min | as given |
| TiSi₂ salicide | p.96–103 | **not modeled**: no silicide material in this simulator build |
| ILD | p.105: PECVD TEOS 5000 Å | `deposit oxide 0.5` |
| Contacts, W plugs, Al | p.106–117 | contact cuts filled with aluminum, four metal islands |

Shared models: TSUPREM-4 P diffusivity in SiO₂ (see `decks/subcmos2014`),
`method fermi vertical`.

## PDF discrepancies and assumptions

- **PMOS Vt implant.** The step list (p.20 step 22) says 1.75e12 B11 60 keV.
  The implant slide (p.53) and RIT's Athena result (p.48) say 3.5e12 P31
  60 keV. Phosphorus is the right species for a p+ poly PMOS; the deck uses p.53.
- **Recessed-oxide recipe.** "OX04 Tube 1" has no recipe slide. The deck
  assumes the SubCmos2014 field-oxide form (1100 °C wet, recipe 406 ramps) and
  calibrates the soak to the 3650 Å Athena result (p.39).
- **N₂O.** SUPREM-IV has no N₂O ambient. The interval is run in N₂ and the
  final O₂ soak is calibrated to the 100 Å target. Nitrogen in the oxide is not modeled.
- **Channel stop.** The final cross-section (p.134) shows a channel stop that
  the step list does not include. It is not modeled.
- **Mask geometry** is assumed (`build_decks.py`: `FIELD`, `GATE`,
  `CONTACTS`, `METAL_GAPS`), not extracted from artwork.

## RIT reference results (Athena/Atlas, pp.47–48)

| | Without Vt adjust | With Vt adjust |
|---|---|---|
| NMOS | 0.61 V | 0.73 V (B 3e12 30 keV) |
| PMOS | −0.31 V | −0.74 V (P 3.5e12 60 keV) |

## Electrical defaults

`prepare_electrical.py` places a fixed positive charge, **QF = 5e11 cm⁻²**, on every Si/SiO₂ interface
(override with `--qf`). Without it, the simulated Vt are 0.98 / −0.52 V on the wide layout, and a hand
calculation with RIT's own plotted doping also gives about 1.0 V for NMOS. QF = 5e11 gives 0.71 / −0.78 V,
against RIT's 0.73 / −0.74 V. Each 1e11 cm⁻² shifts both devices by about −0.055 V. See
`suprem4/doc/rit_advcmos_first_flow_20261121.md`.

## Layout (assumed, `build_decks.py`)
Body-contact island 0.3–1.5 µm (1.2 µm: a 0.8 µm window is nearly filled by the 0.35 µm spacer
nitride), field 1.5–1.9 µm, device active 1.9–8.5 µm, gate centred at 5.0 µm, domain 9.5 µm.
A 3 µm device active was found to dilute the n-well by ~25 % (field oxide blocks the 170 keV P and the
6 h drive spreads it laterally), so the device active is wide.

## Poly re-oxidation (`--poly-reox geometric`)
SUPREM-IV.GS never oxidizes poly.  The geometric option deposits the poly 0.44·t thinner and etches it
0.88·t shorter (t = `--poly-reox-ox`, 0.05 µm), forms the sidewall oxide by conformal deposition and
vertical etch-back of t, then runs the normal recipe-250 oxidation for the substrate.  With
`ADVCMOS_LPOLY=0.40` the gate is 0.355 µm after re-ox (RIT: 0.35).  The vertical model still forms the
gate-edge "smile" (≈300 Å at the edge, tapering to 100 Å over ≈0.15 µm).  Limitations: no 500 Å top
oxide on the poly; small oxide slivers on field-oxide steps.  Background and the parked 2D poly
oxidation work: `suprem4/doc/poly_oxidation/`.

## Calibration notes
- PMOS n-well: phosphorus-in-nitride ranges from SRIM-2013 (3.1 g/cm³) instead of Gibbons; the well
  implant goes through the 1500 Å nitride (`suprem4/data/sup4gs.imp`).
- Electrical: fixed oxide charge 5e11 cm⁻² by default (`prepare_electrical.py --qf`), which reproduces
  RIT's 0.73 / −0.74 V.
- Results and analysis: `suprem4/doc/rit_advcmos_first_flow_20261121.md`.
