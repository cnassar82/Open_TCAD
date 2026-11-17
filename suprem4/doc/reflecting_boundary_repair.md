# Reflecting boundary construction repair

2026-11-15. Full CMOS acceptance remains open; see the process result files.

## Root cause

`dbase/new_layer.c::local_normal` previously enforced a reflecting boundary by
subtracting the normal projection from the interface normal. With float
intermediates, cancellation left a nonzero constrained component. `SplitNode`
then used that component as a geometric displacement. This is a construction
error, not a physical small feature and not a reason to relax solver tolerances.

A controlled pre-LOCOS PMOS replay starts from the same wells checkpoint on
`native64_0o563h8b` and `native64_yuj0z1w7`. The old formula first introduces
x = 5.684342356405358e-17 um during the 825 C dry oxidation following the
oxide strip. The corrected formula retains x = 0 through the mask sequence.
The old drift subsequently propagates through LOCOS refinement. A downstream
source/drain checkpoint contains oxide triangles with nodal control volume
1.2924685484455974e-33 cm2. Activation's first trapezoidal stage produces
updates approximately -2 times the boron/phosphorus concentrations there,
and halving the timestep to the retry limit does not resolve it.

The repair constructs the projection directly from the reflecting face's
tangent, with double intermediates, rather than subtracting two nearly equal
vectors. An axis-aligned face therefore produces an exactly zero constrained
component before any new point is constructed. It does not snap coordinates,
remove geometry, change dopant concentrations, or change the time integrator.
It does not repair an already corrupted checkpoint: a continuation from the
old `locos_mask` still fails with nonfinite oxidation velocities.

## Verification

* `tests/run_reflecting_normal.sh`: 5,600 cases spanning quarter-turns,
  coordinate scales, slopes and both normal styles. ASan/UBSan passes with
  leak detection disabled. The original source fails the exact constrained
  component assertion.
* `tests/fixtures/rit_pmos_reflecting_normal`: frozen pre-825 C checkpoint.
  The old build exits normally but fails the reflecting-plane audit. The new
  sanitizer build `native64_70li8k7h` passes the plane and independent partition
  audits. This distinguishes silent construction corruption from a crash.
* Persistent control and sanitizer results:
  `decks/subcmos2014/results/astra_reflecting_{control,sanitized}/`.
* Command-by-command traces:
  `decks/subcmos2014/results/astra_locos_mask_{control,trace}/`.
* Longer clean-history continuation:
  `decks/subcmos2014/results/astra_tangent_from_wells/pmos/`.
  Its newly generated LOCOS mask and LOCOS checkpoints contain no nonzero
  near-origin side coordinates (the historical pair contained 1 and 13).
  Consult `result.json` for final process acceptance, not output presence.

## Separate finite-segment classification correction

`refine/skel.c::clip_on_segment` must reject a point outside the finite segment's
coordinate bounds before applying its supporting-line roundoff allowance.
Previously that allowance extended the clipping segment past an endpoint,
causing two retained spacer arcs to select the same successor.

The regression uses the captured endpoint and false-ON point under rotations
and scaling. The original source fails; the correction and the existing rounded
shared-cut tests pass. A fully exact collinearity experiment was rejected:
rounded, previously constructed intersections still require their existing
supporting-line allowance. Only finite segment membership was corrected.

The captured spacer replay passes independent geometry checks under strict
ASan/UBSan on `native64_70li8k7h`, with leak detection disabled:
`decks/subcmos2014/results/astra_finite_spacer_sanitized/`.
This bounded repair does not replace the upstream reflecting-boundary repair.

Negative-control and focused-test logs are in
`suprem4/experiments/reflecting_boundary_20261115/`.
