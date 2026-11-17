# LOCOS refinement study, 2026-11-10

The stress-dependent 90-minute, 1000 C wet-oxidation process, geometry, and
material parameters are unchanged. Candidate executable:
`suprem4/build/native64_e0r_18r5/src/suprem4`.

Two kinds of refinement are distinguished:

- Substrate-only: scale all initial line spacings by 1/sqrt(N).
- Full mesh: use the same initial refinement, also scale oxide insertion
  spacing from 30 nm by 1/sqrt(N), and increase deposited pad/nitride divisions
  by sqrt(N), rounded upward. Physical deposited thicknesses stay unchanged.

The control starts with 560 points and completes with 737 points and 1,363
triangles. Substrate-only refinement to 1,173 points (2.095 times the control)
fails with an inverted triangle around 52.18 simulated minutes. Refinement to
2,263 points (4.041 times the control) fails around 86.84 minutes. Both failure
structures also contain negative-area triangles in the independent mesh audit.

Full 2x refinement uses 3 pad-oxide divisions, 5 nitride divisions, and
21.213 nm oxide insertion spacing. It completes with 1,474 points and 2,785
triangles. All four saved meshes pass the topology and finite-field audit.
A second independent run produces byte-identical initial, masked,
field-oxidation and final structures. Both reach 5,400 simulated seconds.

Final oxide area changes from 0.61199116 to 0.60967436 square microns (-0.379%).
The silicon surface probe 50 nm under the mask edge changes from 0.140218 to
0.136018 microns. These differences are discretization sensitivity, not proof
of calibrated accuracy or full mesh convergence.

The passing 2x deck is `decks/locos/locos_refined.in`; all four saved structures,
the exact input, log and result are in `decks/locos/refined_output`.

Full 4x refinement also completes, with 4 pad divisions, 6 nitride divisions and
15 nm oxide insertion spacing. It starts with 2,263 points and finishes with
2,791 points and 5,353 triangles. All four saved meshes pass the audit and the
oxidation reaches 5,400 seconds. The first attempt reached its 180-second wall
time limit; a longer run completes in 565.5 seconds. A timeout is not counted
as a mesh inversion. Its deck and outputs are `decks/locos/locos_refined_4x.in`
and `decks/locos/refined_4x_output`.

| Refinement | Initial points | Final triangles | Result |
|---|---:|---:|---|
| Control | 560 | 1,363 | Pass |
| Initial substrate only, approximately 2x | 1,173 | — | Inversion at ~52.18 min |
| Initial substrate only, approximately 4x | 2,263 | — | Inversion at ~86.84 min |
| Substrate and layers, approximately 2x | 1,173 | 2,785 | Pass twice; byte-identical |
| Substrate and layers, approximately 4x | 2,263 | 5,353 | Pass |

The 4x oxide area is 0.61614164 square microns (+0.678% from control). Its
silicon surface probe 50 nm inside the mask is 0.144451 microns. Area and probe
changes are not monotonic across these meshes; these completion tests do not
establish numerical convergence. There was no physics or mesh-engine change
in this study, so the substrate-only failures remain unresolved.

Machine-readable results, including the timed-out attempt, the longer 4x run,
binary SHA-256, and repeated-2x checkpoint comparisons:
`suprem4/doc/locos_refinement_results.json`.

Run the study or a selected case:

```sh
python3 suprem4/tests/run_locos_refinement.py \
  --binary suprem4/build/native64_e0r_18r5/src/suprem4 --timeout 600
python3 suprem4/tests/run_locos_refinement.py \
  --binary suprem4/build/native64_e0r_18r5/src/suprem4 \
  --cases full_2x --timeout 300
```

This establishes reproducibility for the tested 2x configuration and successful
completion at 4x with coordinated layer refinement. The installed launcher is
unchanged. The viewer was opened on the repeated/verified 2x result.
