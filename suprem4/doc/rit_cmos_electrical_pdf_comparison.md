# Electrical follow-up: the final RIT structures are not accepted devices

**2026-11-17 update:** the gate-oxide shielding and missing retained-arc repairs now pass both fresh 16-checkpoint processes. See [current repair and limits](rit_cmos_gate_oxide_shielding.md). Earlier missing-oxide structures and withdrawn acceptance below are historical evidence.

2026-11-16. Actual Genius drift-diffusion transfer sweeps complete, but they
expose a missing gate dielectric in both final process structures. The earlier
claim of validated final CMOS structures was too strong and is withdrawn.
Process-command completion and nonoverlapping positive mesh elements did not
establish that every required material volume existed.

## PDF comparison

Reference: `rit_process/SubCmos2014_0p5um.pdf`, page 11. Gate-oxide recipe and
measurement targets are repeated on pages 81-82. The PDF specifies nominal
thresholds, rather than a measured transfer curve or a threshold extraction
criterion for this particular 2D cut.

| Quantity | PDF | NMOS | PMOS |
| --- | --- | --- | --- |
| Minimum drawn gate length | 1 micron | 1 micron | 1 micron |
| Nominal threshold | +0.75 / -0.75 V | Cannot validly extract | Cannot validly extract |
| Gate oxide after gate-oxide recipe | 15 nm | 16.6741 nm (+11.16%) | 16.9699 nm (+13.13%) |
| Gate oxide in final gate stack | Must remain between gate and channel | Missing in contact opening | Missing in contact opening |
| Field oxide, x=2.5 microns, after activation/before ILD | 600 nm nominal | 629.28 nm | 625.90 nm |
| Well net doping 0.5 microns below well surface, before threshold implants | P-well 1e16 / N-well 3e16 cm^-3 nominal | P-type 5.207e16 cm^-3 | N-type 5.909e16 cm^-3 |

The well probes describe a specified point in a nonuniform profile. The PDF's
nominal well concentrations do not define that same sampling location, so this
is a calibration warning, not an exact error percentage. Subsequent threshold
and source/drain implants are excluded from these well probes.

## What was actually run

The checked-in Genius executable could not start against the installed PETSc
and other library versions. An isolated snapshot was configured, built and
installed against the available PETSc 3.24.4; the existing installation was not
overwritten. The resulting executable SHA-256 is
`deefa548eb0d88b05f11bd0bbd0f17c50f5867b947da0b024a4ea9d00dd24c71`.
Build/configuration logs and the runtime are in
`suprem4/experiments/rit_electrical_20261116/`.

Both final `metal1.str` structures were exported without deleting any vertex
or triangle. Active B/P concentrations are used. Body/source/gate/drain are
assigned by the four existing metal islands. Physically contacting n+ poly
and gate metal are combined into one ideal conductor for electrical analysis.
This approximation omits poly depletion; the 4.17 eV work function is consistent
with this Genius Si-affinity model and its n+ gate example, and is used for
both devices. It was not fitted to obtain the PDF thresholds.

Baseline: 300 K, one-micron width, Lucent mobility, Fermi statistics, fully
ionized supplied active dopants, explicitly zero Si/oxide fixed charge and no
fitted interface traps. An initial equilibrium check using the solver's default
1e10 cm^-2 interface charge is archived separately; the reported transfer
sweeps explicitly override it to zero.

Each device completes 67 converged gate-bias points from zero to +/-3.3 V,
after a drain ramp to +/-0.05 V. The recorded terminal voltages confirm that
the drain bias is held and every gate boundary receives the requested gate
voltage. No `ERROR` or `DIVERGED` result occurs. Terminal KCL discrepancies are
at most 2.03e-15 A for NMOS and 1.18e-16 A for PMOS.

* NMOS drain current remains at numerical noise, maximum magnitude 9.03e-16 A.
  This must not be interpreted as excellent transistor off-current.
* PMOS current magnitude is 0.85713 microamps at Vg=0 and 0.85887 microamps
  at Vg=-3.3 V, with essentially no gate modulation.

These are diagnostic curves of invalid structures. They do not supply
meaningful threshold voltages or validate either +/-0.75 V PDF target.
The missing material separates the gate and silicon from the intended oxide
PDE domain; the converged algebraic solution is insufficient for a transistor
performance claim.

Data, assumptions and hashes:

* `experiments/rit_electrical_20261116/qf0/export_manifest.json`
* `experiments/rit_electrical_20261116/qf0/{nmos,pmos}/transfer_execution.json`
* `experiments/rit_electrical_20261116/comparison/comparison.json`
* `experiments/rit_electrical_20261116/comparison/{nmos,pmos}_transfer.csv`
* `experiments/rit_electrical_20261116/comparison/diagnostic_transfer.png`

## First failing process operation

Both gate stacks contain oxide at `gate`, `spacers`, and `annealed`, and lack
oxide in the gate opening at `contacts` and `metal1`. A new NMOS replay brackets
this more tightly:

| Saved operation | Required gate oxide present |
| --- | --- |
| Read annealed state | Yes |
| Deposit 0.5-micron interlayer oxide | Yes |
| Deposit/pattern contact resist | Yes |
| `etch oxide dry thick=1.5` | No |
| Strip resist | No |

Evidence: `experiments/rit_electrical_20261116/contact_bracket/occupancy_bracket.json`.
The input is frozen in `tests/fixtures/rit_nmos_contact_shielding`. Its replay
must require gate-oxide occupancy, not merely successful simulator exit.

The dry-etch implementation builds a globally translated surface cutter and
subtracts it from eligible oxide regions. Being part of an exposed connected
oxide region does not make every buried part visible to the incident etch.
The cutter reaches oxide beneath the nonetched polysilicon instead of stopping
there. The gate contact opening is x=5.55..6.45 microns, directly over the gate.

The algorithmic repair must account for directional accessibility and shielding
by nonetched materials. A selective vertical etch front may remove exposed
oxide until it reaches poly or silicon; it cannot then continue into buried
oxide beneath that blocking solid. Connectivity alone is not an accessibility
criterion. Required gates include intended contact openings plus unchanged
buried gate oxide, valid assembled domains, and conservative field transfer.
Reducing this particular etch depth until the example passes, filling the void
afterward, or applying an area cutoff would not establish the required behavior.

## Validation corrected now

`validate.py::gate_dielectric` checks actual oxide coverage between poly and
silicon at five positions across the gate. `gate`, `annealed`, `contacts`, and
`metal1` validation now require it. A regression distinguishes identical
positive poly/Si distances filled by oxide, partially filled, or empty. The
current saved final structures fail this corrected check in both devices.

The generic exporter's legacy dopant map also mislabeled this SUPREM fork's
holes/velocity and grain-size/tracer codes as dopants. That map is corrected;
`tests/test_genius_dopant_export.py` checks active-dopant precedence and excludes
these unrelated fields. `prepare_electrical.py` now rejects missing-dielectric
structures unless the explicit diagnostic override is supplied.

The RIT acceptance manifest is marked `accepted: false`, with the original
acceptance archived. The `validated-rit-cmos` pin has been removed. All source
builds, structures and logs remain available. The general launcher is unchanged.
The physical contact-etch repair remains open; no threshold fit or repaired
process result is claimed here.
