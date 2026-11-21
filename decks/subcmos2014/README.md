# RIT Sub-CMOS 2014 NMOS and PMOS

**Presentation audit (2026-11-17): agreement is not yet achieved.** The
[full audit](../../suprem4/doc/rit_cmos_presentation_audit_20261117.md) includes
224 stage/location profiles and operating-bias electrical sweeps. It identifies
a severe 2D poly-implant resolution failure, oxide/well-depth deviations,
threshold estimates of +0.246/−0.555 V versus ±0.75 V, and PMOS off-current
0.722 µA/µm at 3.3 V drain bias. A completed process run does not qualify these
profiles or electrical predictions.

**Gate-oxide shielding repair verified on fresh NMOS and PMOS runs (2026-11-17).**
Both complete all 16 checkpoints through first metal. The gate dielectric
survives the contact etch unchanged: 16.67 nm NMOS and 16.97 nm PMOS at the
center probe, with oxide adjacency along every downward-facing gate edge and
four connected conductor contacts. See `suprem4/doc/rit_cmos_gate_oxide_shielding.md`
and `suprem4/doc/gate_oxide_validation_20261117.json` for the algorithm,
regressions, exact builds and remaining qualification limits.

The old missing-oxide outputs and electrical curves remain historical failure
evidence in `suprem4/doc/rit_cmos_electrical_pdf_comparison.md`. General release
qualification and PDF electrical calibration are separate from this repair.

These are separate 2D transistor process cuts based on
`rit_process/SubCmos2014_0p5um.pdf` (Lynn Fuller, slides dated September 10, 2015).
The PDF's **0.5 µm is the layout lambda**; page 11 specifies a **1 µm minimum
channel length**. Both generated decks therefore use a 1 µm drawn gate.

The decks include wells, LOCOS, threshold adjustment, n+ polysilicon gates for
both transistor types, LDD implants, oxide spacers, source/drain and body-contact
implants, activation, contact openings and patterned first aluminum. They model
fabrication geometry and dopants; they do not calculate or validate transistor
I–V curves or reproduce the PDF's threshold-voltage targets.

## Run

**Fresh full flow (2026-11-19):** both devices complete all 16 audited checkpoints
through metal1 in one invocation. The installed launcher now uses this build with
the validated options as built-in defaults, so no environment variables are needed:

```sh
python3 decks/subcmos2014/run.py --binary suprem4/build/validated-native64/src/suprem4 --through metal1 --timeout 7200
```

**2026-11-20:** thin gate oxides now get at least 8 element layers before each
non-oxidizing diffusion. The NMOS channel is now p-type: before, poly phosphorus
crossed a one-element gate oxide. The remaining 8 nm surface skin from the S/D
anneal matches an independent 1D solution of the same model. Saved run:
`results/fresh_full_flow_thinfilm_20261120`. See
`suprem4/doc/rit_cmos_full_flow_repair_20261119.md`.

The older 2026-11-17 instructions follow.

To reproduce the corrected candidate through first metal:

```sh
SUPREM4_PROCESS_SURFACE=1 SUPREM4_DEPOSITION_GEOMETRY=domain \
SUPREM4_MATERIAL_DOMAINS=oxide SUPREM4_MESH_PATCH_REMOVAL=1 \
python3 decks/subcmos2014/run.py --binary suprem4/build/native64_rxq2pnrb/src/suprem4 --through metal1
```

The fresh audited outputs are in `results/directional_final_fresh_20261117`.
The frozen contact-shielding regression now passes without reducing its etch depth.

To rebuild, run `python3 suprem4/tools/build_native64.py --ordering amd
--domain-offset --jobs 4`, then substitute the printed executable in the command
above with the same runtime options. `build_decks.py` regenerates the recipes
if you deliberately change their translation settings.

Run from the workspace root. The runner prints an output directory and writes
`results.json`, per-device logs, audits and 16 saved checkpoints. `nmos.str` and
`pmos.str` are the first-metal results. Add `--output /absolute/path` to retain
results in a chosen location, `--devices nmos` or `--devices pmos` for one device,
or `--through locos` to stop early. Each device's output folder must be empty.

To resume an unchanged generated process from a saved checkpoint, use the same
runtime options as its original run:

```sh
python3 decks/subcmos2014/run.py --binary /absolute/path/to/suprem4 \
  --resume /absolute/path/to/previous/run --after locos
```

The runner checks the original process prefix, audits imported checkpoints, and
restores method/stress settings before reading the checkpoint. A resumed run is
identified in its result JSON. For verification of the entire process, use a
fresh run. Generated decks and `process.json` record the recipe segments.

The withdrawn `validated-rit-cmos` pin has been removed. The general installed
launcher has not been switched to this RIT candidate.
Older builds do not contain the boundary, implant and restart fixes needed by
these decks.

## Process translation

Page numbers below refer to PDF pages, counting the cover as page 1. Doses are
cm⁻², energies keV, lengths µm. Detailed recipe slides take precedence over the
short process-summary slide where they disagree.

| Operation | PDF pages | Deck translation |
| --- | --- | --- |
| Substrate | 11, 18 | (100) p-Si, B background 1e15 cm⁻³ |
| Pad oxide | 23–24 | Recipe 250: dry oxidation, 1000 °C / 40 min soak; nominal 0.05 |
| Well mask | 38, 41, 43 | P 9.5e12 / 125 for N-well; recipe 330 masking oxide, wet 1000 °C / 50 min; nominal 0.30 |
| Prestrip BOE | 48 | 0.025 oxide removal |
| P-well | 50, 54 | Default improved BF₂ 8e12 / 110; optional legacy B 2e13 / 35 |
| Well drive | 55–56 | N₂ 1100 °C / 600 min soak, then strip oxide |
| LOCOS mask | 58–63 | Second recipe-250 pad oxide; 0.15 nitride, patterned active and body islands |
| Channel stop | 66 | NMOS field windows only: B 8e13 / 100 |
| LOCOS | 67–69 | Recipe 406, wet 1100 °C / 65 min; nominal 0.65, followed by nitride strip and 0.05 oxide loss |
| Sacrificial oxide | 71–73 | Recipe 311, wet 900 °C / 40 min; nominal 0.10 |
| Threshold adjust | 75, 77–78 | Blanket B 7e11 / 60; additional PMOS B 2.6e12 / 60; strip 0.10 oxide |
| Gate oxide | 80–81 | Recipe 215, dry 900 °C / 45 min; nominal 0.015 |
| Poly gate | 85, 87–91 | Deposit 0.50 poly; P 2e16 / 50; recipe 120 N₂ 1000 °C / 15 min; pattern 1 µm gate |
| LDD | 99, 102 | NMOS P 2.5e13 / 60; PMOS B 4e13 / 50 |
| Spacers | 104–107 | Deposit 0.50 conformal oxide; directional oxide etch 0.50 |
| Heavy implants | 109, 112 | N+ P 2e15 / 60; P+ B 2e15 / 50; opposite-type body tie in separate implant window |
| Activation | 116 | Recipe 280, dry O₂ 1000 °C / 20 min soak |
| Contacts and metal 1 | 117–128 | 0.50 oxide ILD, four selective contact openings; 0.75 aluminum and pattern |

The furnace recipes include the listed loading, stabilization, flood, anneal,
cooldown and pull times, not just their headline soak. Temperature ramps are
approximated by constant-temperature midpoint segments of at most 50 °C.
Oxides are **grown by the simulator**, not forced to their nominal slide
thicknesses. Actual thicknesses and dopant profiles depend on its models and
require calibration against measurements.

`build_decks.py --pwell boron` selects the older p-well implant. Regenerate both
decks before running that variant; the default saved results use BF₂.

## Layout and model assumptions

The PDF contains process drawings, but no mask-coordinate file. Each isolated
cut is 10 µm wide and 8 µm deep. The main active island spans x=3–9; a separate
body-tie island spans x=0.5–2; the gate spans x=5.5–6.5. Contact centers are
x=1.25 (body), 4.25 (source), 6 (gate), and 7.75 (drain). These are assumed masks,
not a reconstructed layout or a complete CMOS circuit. The N-well mask is open
across the PMOS cut and closed across the NMOS cut. The subsequent P-well implant
is screened by the thicker N-well masking oxide in the PMOS cut.

Depositions are ideal conformal geometry operations. Their CVD thermal budgets,
plasma damage and deposition kinetics are not modeled. Etches are ideal material
selective operations. The spacer etch uses the nominal deposited thickness;
the PDF does not specify an exact overetch depth. The contact etch uses an oxide
removal depth of 1.5 µm with ideal stopping on silicon/poly, under patterned
photoresist. It is not a calibrated contact-etch rate model.

LOCOS uses SUPREM's viscous/stress-dependent model and coefficients
`Vc=300 Vr=30 Vd=25` from SUPREM example 9; those coefficients do not come from
this PDF. The legacy BF₂ range table has no lateral-straggle column; the lookup
uses its vertical straggle as an explicit isotropic lateral estimate. The
correct BF₂ ion table is selected, rather than silently substituting elemental
boron. Crystal channeling, implantation damage/TED calibration and interface
charge calibration are outside this first implementation.

After LOCOS the default uses vertical growth for the thinner oxides. This
omits lateral thin-oxide growth and subsequent mechanical stress feedback.
`build_decks.py --post-locos-model viscous` retains the full viscous model for
investigation; it has encountered nonlinear convergence failures and is not a
validated variant. The runner rejects those failures even if later commands
continue. This model limitation is separate from the geometry fixes.

Later implant and contact masks use a flat resist top at y=−2 µm through the
existing file-defined deposition option (`resist_top.dat`), approximating a
spin coat over the spacer topography. Earlier masks use a conformal resist
approximation. These resist geometries are assumptions; they are not derived
from mask artwork or a resist-flow model. The general conformal-offset routine
still has a known failure over this spacer topography.

The process stops at patterned first metal. The forming-gas passivation and
second-metal interconnect sequence are not represented. Electrical contacts
for a device solver are not assigned by these process decks.

## Acceptance checks

Every saved stage is audited for finite coordinates and fields, positive signed
triangle area, distinct vertices, duplicate triangles, reciprocal neighbors and
edge ownership. Independent geometric predicates reject interior overlaps,
crossing boundaries and unnoded contacts; B/P chemistry must be finite and
nonnegative, and integrated inventories are recorded. The runner treats
reported simulator errors as failures even
if the executable exits with status zero. It also checks well polarity, n+ poly,
1 µm drawn gate width, a remaining gate dielectric, source/drain and body-tie
polarity, four separate first-metal islands, and the absence of oxide gaps at
the four contact-center probes. These are structural checks, not electrical
calibration or a proof of mesh convergence.

See `suprem4/doc/subcmos2014_reliability.md` for engine fixes and regression
coverage. Generated simulation outputs under `results/` are local artifacts.

## Corrected-structure electrical smoke check (2026-11-17)

Both fresh final cuts complete 67-point Genius DD transfer sweeps at 300 K and
|Vds|=50 mV. For 1-micron width, the gate-driven endpoint currents are 11.37 uA
NMOS and 5.23 uA PMOS. These use an unfitted model; threshold/PDF agreement and
full mesh convergence are not established. See the current shielding repair
report and `suprem4/experiments/directional_etch_20261117/electrical/summary.json`.
