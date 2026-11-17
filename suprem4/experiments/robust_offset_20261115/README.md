# Robust deposition-domain experiment

Clipper2 1.5.4, Boost Software License. The unmodified dependency and pinned
archive SHA-256 are recorded in `suprem4/src/third_party/clipper2/PROVENANCE.md`.
The original standalone prototype is `offset_domain.cpp`; it uses a finer
decimal grid and does not model the simulator storage conversion.

The integrated implementation is `suprem4/src/refine/domain_offset.cpp`.
Build with `python3 suprem4/tools/build_native64.py --ordering amd --domain-offset --jobs 4`.
Enable it with `SUPREM4_DEPOSITION_GEOMETRY=domain`. It constructs an occupied
closed domain, dilates it, clips it to the existing side planes, and uses
Boolean difference to validate the coating. It refuses holes and multiple
components that the existing skeleton mesher cannot represent.

The integrated construction uses a dyadic grid derived from the maximum
coordinate magnitude and binary32 storage precision. Every new outer vertex
is exactly representable in the database; source vertices are unchanged.
The early fine-grid implementation failed because conversion to float merged
vertices and reversed short edges. Those failed controls are retained.

Analytic planar and step coating areas, translations, exact side identities,
float representability and the failure contract pass ASan/UBSan, with leak
detection disabled. Logs: `domain_offset_float_grid.log`.

Ordinary and sanitized captured ILD and spacer cases pass:
`decks/subcmos2014/results/astra_domain_{ild,spacer}_{float_grid,sanitized}/`.
The ordinary ILD mesh has 49,207 triangles; all four pre-existing region areas
remain identical. Added coating area is 5.374806130246482 um2.

`decks/subcmos2014/results/astra_domain_pmos_finish/pmos/` completes the archived
annealed checkpoint through contacts and metal1, with 16/16 runner checks.
This is a resumed success, not complete fresh CMOS acceptance.

Fresh PMOS on the initial domain candidate fails during LOCOS; fresh NMOS
reaches spacer deposition and fails in quadsplit. The stable-angle candidate
is being evaluated separately. Consult each result.json for final outcomes.
