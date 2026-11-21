# RIT CMOS numerical repair continuation

**Full CMOS/electrical acceptance: false.** See the [repair report](../../doc/rit_cmos_fem_repair_20261118.md).

The launcher is unchanged. Numerical changes in the main source are tested through immutable builds, with explicit runtime options. Experimental work under ignored `work/` is archived as `experimental_solver.patch`; the overlay includes unqualified algorithms and must not be applied as a qualified bundle.

## Evidence to inspect first

Latest main candidate `native64_5njvfhlt` completes both saved-LOCOS continuations through gate oxide (8/8 each). The native-join cap fixes real dielectric pinch-outs; optimal cavity triangulation fixes the separate PMOS coarsening deadlock. Both remain subject to the full-process acceptance boundary below.

- `main_native_optimal_gate/results.json`: combined main source, NMOS/PMOS 8/8 through gate oxide, explicit checkpoint restarts and provenance.
- `native_join_main_sanitized/result.json`: positive-width native join, no Si/Poly point contacts, independent geometry and chemistry audits under ASan/UBSan.
- `main_sanitized_native_join_controls/`, `main_sanitized_seed_loop/result.json`: bare/wedge/gap and closed-loop controls. Conservative patch removal enabled for the join controls.
- `main_sanitized_join_shifts.json`: four join shifts, +/-5 and +/-25 nm, sanitizer checks and physical dielectric separation.
- `nmos_native_chain_activation/independent_audit.json`: complete unchanged 117-minute activation with positive-width joins, body phosphorus reduced about 814 times and p-type body restored. Gate-channel transport is still unqualified.
- `pmos_cavity_quality.json`, `../../tests/fixtures/rit_pmos_quality_deadlock/cavity.json`: exact captured coarsening deadlock; greedy minimum quality 5.07e-7 versus optimal 0.8397, with unchanged quality guard.
- `native_join_comparison.png`: old zero-width pinch versus added gas-side cap.

- `production_full_locos/{execution,independent_audit}.json`: completed full captured PMOS LOCOS sequence, main-source build, independent partition/chemistry checks; final structure byte-identical to `trace_original_full`.
- `fresh_masks_embedding/results.json`: fresh NMOS/PMOS through their five LOCOS-mask checkpoints, all audited.
- `cmos_mixed_trace/results.json`: incomplete full-process continuation. NMOS activation failure; PMOS post-LOCOS vertical-growth stagnation. Historical failures; superseded for the gate-oxide suffix by `main_native_optimal_gate`. Neither historical full-process attempt is accepted.
- `chemical_capture/{linear_analysis,tail_row}.json`: an independently accurate linear solve reproduces the negative dopant update; the edge flux loses monotonicity.
- `activation_sg_production_60s/{execution,independent_audit}.json`: main-source SG flux completes a bounded 60-second activation replay; topology, partition and dopants pass. This is not the original 12-minute hold.
- `activation_sg_production_12min/execution.json`: completed original first 12-minute hold (888.10 wall seconds), with independent geometry/dopant audit.
- `oxide_transport_resolution/convergence.json`, `oxide_transport_sg/comparison.json`: gate barrier resolution changes transmitted P dose by about 29 times. The finest 64/128 pair differs by about 0.4%.
- `normal_growth_*/result.json`: achieved planar-growth spacing controls. Simulator exit zero does not imply the requested normal metric was achieved. Sub-1-nm requests are reset to 100 nm by the current parser.
- `production_sanitized_spacing/results.json`: six achieved deposition-spacing/implant controls under sanitizers.

## Reproduce the activation regression

From the repository root, build an immutable candidate:

```sh
python suprem4/tools/build_native64.py --domain-offset --ordering amd --fe-solver umfpack --jobs 4
```

Use the emitted binary path in a **new** output folder:

```sh
SUPREM4_CHARGED_FLUX=sg python suprem4/experiments/repair_continuation_20261118/run_activation.py \
  --binary /absolute/path/to/build/src/suprem4 \
  --output /absolute/path/to/new/output --first-minutes 1 --timeout 600
```

Omit `--first-minutes 1` for the original first 12-minute hold. `--phase full` uses the complete unchanged activation recipe. The initial structure defaults to the saved, corrected NMOS poly implant; the runner records its SHA-256. Runtime limits are failure/incomplete bounds, not convergence tolerances. Do not relabel a timeout as an algorithmic crash or a process pass.

Run the charged-edge source regression with:

```sh
bash suprem4/tests/run_charged_flux.sh
```

## Mesh work still required

The isolated normal-layer prototype is rejected as a complete remesher: it only targets the interface-adjacent layer, and coarsening does not enforce the requested anisotropic metric. Forecast motion helps one planar control but the finer control still exceeds the target. The proper repair requires a shared, achieved normal metric across refinement, motion and coarsening, plus conforming material constraints and conservative field transfer. Do not fit process parameters to compensate.

To replay the remaining activation budget from its completed 12-minute checkpoint, use `run_activation.py --phase remaining --initial /absolute/path/to/12-minute/activated.str` with a new output folder. This explicitly restarts the integrator for the first remaining hold and records `checkpoint_restart: true`.

The remaining activation budget also completed (750.22 wall seconds), with independent topology/partition/P/B audits. Its output and the completed first hold form an explicit checkpoint chain in `nmos_activation_checkpoint/nmos/activation_provenance.json`. Downstream numerical qualification is recorded separately in `nmos_sg_downstream`; the gate-barrier mesh and full CMOS electrical acceptance remain false.

`probe_kooi.py` reproduces the stressed PMOS failure geometry controls with all inherited oxide coefficients restored. `--jcheck` is a diagnostic that deliberately stops after finite differences, never a process pass. The latest derivative probes remain experimental and need a perturbation-step study before attributing discrepancies to an analytic-Jacobian defect.
