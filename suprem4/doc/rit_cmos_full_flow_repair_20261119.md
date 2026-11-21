# RIT CMOS: fresh full flow through metal1, 2026-11-19

**Promoted (2026-11-20).** `bin/suprem4` (via `build/validated-native64`) now runs `native64_8som14gn` (`--ordering amd --domain-offset --fe-solver umfpack`; sources identical to the working tree). It adds fix 4 below, thin-dielectric through-thickness resolution, to the 11-19 build `native64_7hxx2b_8`. Validated runtime options are built-in defaults (`src/include/suprem_options.h`); set any one to `0`/`off`/`legacy` for the legacy algorithm. A fresh run with no SUPREM4_* variables passes 16/16 audited checkpoints on both devices: `decks/subcmos2014/results/fresh_full_flow_thinfilm_20261120`, re-audited in `revalidated.json`.

**Doping.** On the 11-19 build the NMOS channel was n-type to 79 nm. Poly phosphorus crossed a gate oxide meshed one element thick. This was a discretization defect, not physics: with identical parameters, P reaching Si fell 1.7e10 → 2.0e8 → 6.8e7 cm^-2 for 1, 8 and 16 element layers. With fix 4, the channel at 50 nm is p-type (net −4.7e16 cm^-3, previously +3.2e16), and the poly-activation P at the channel surface is 1.7e14 cm^-3. An 8 nm n-type skin (surface P 4e16, about 3e10 cm^-2, ΔVt ≈ 25 mV) forms during the S/D anneal. This is the model's physics, not a mesh error: an independent 1D implicit finite-volume solution with the same modelrc parameters gives a Si surface P of 7.4e16 for poly activation plus S/D anneal (√ΣD_ox·t rises from 4.0 to 6.2 nm). Whether D_ox(P) is calibrated for this process is a model question. Poly retention, junctions (0.64 / 0.92 µm), wells (+26% N-well depth) and gate oxides (17.1 / 18.0 nm, target 15) are unchanged.

**Audit correction.** The PMOS body probe (50 nm, before S/D) required n-type. But the blanket boron threshold adjust makes the n-well p-type there before poly exists (`gate_oxide`: B 1.25e17, P 9.6e16; a buried-channel PMOS). The 11-19 run passed only because leaked poly P added about 6e16. `validate.py` now requires the poly steps to preserve the polarity found in `gate_oxide.str`. That check passes the new run and **fails** the 11-19 run, which correctly flags its leak.

Continues [the 11-18 FEM repair](rit_cmos_fem_repair_20261118.md). That work reached gate oxide only by resuming from saved checkpoints. Today's first fresh run with the 11-18 source reached poly activation on both devices, then failed:

| Device | Last good checkpoint | Failure |
|---|---|---|
| NMOS | `spacers` | `triangles are not clock wise, data base corrupted` (Delaunay flip during the photoresist deposit) |
| PMOS | `source_drain` | `native oxide cap requires a global union` (first oxidizing step of S/D activation) |

Three independent defects were fixed at their root, in general form. No process parameter, deck or tolerance was tuned for this deck.

## Result

Saved run: `decks/subcmos2014/results/fresh_full_flow_20261119`. Build `native64_kkj_q2ee` (`build_native64.py --ordering amd --domain-offset --fe-solver umfpack`), fresh single invocation, run options:

```sh
SUPREM4_PROCESS_SURFACE=1 SUPREM4_DEPOSITION_GEOMETRY=domain SUPREM4_MATERIAL_DOMAINS=oxide \
SUPREM4_MESH_PATCH_REMOVAL=1 SUPREM4_FE_OXIDANT_P1=1 SUPREM4_FE_BBAR_P1=1 \
SUPREM4_FE_LINEAR_TRACE=1 SUPREM4_NATIVE_JOIN_CAP=1 SUPREM4_CHARGED_FLUX=sg \
python3 decks/subcmos2014/run.py --binary <build>/src/suprem4 --through metal1 --timeout 7200
```

| Device | Exit | Checkpoints passing all audits | Wall time |
|---|---|---|---|
| NMOS | 0 | 16/16 | 3150 s |
| PMOS | 0 | 16/16 | 2955 s |

Metal1 mesh, compared with a run that had only fixes 1 and 2:

| | Triangles with quality < 1e-3 | Minimum edge | Triangles |
|---|---|---|---|
| NMOS, without fix 3 | 42 | 5.8e-7 um | 97,097 |
| NMOS, final | 7 | 1.2e-6 um | 80,244 |
| PMOS, without fix 3 | 46 | 3.6e-7 um | 96,409 |
| PMOS, final | 8 | 1.2e-6 um | 72,267 |

This is process completion with audited topology, partition, chemistry and device checks. It is **not** electrical qualification; the 11-17/11-18 caveats on profiles and threshold still apply.

## 1. Native-oxide join with a flat or falling neighbour (`dbase/new_layer.c`)

The positive-width cap followed the old oxide free boundary until it crossed the new film's offset line. That assumes the oxide rises toward the film, as at a LOCOS bird's beak. A spacer top, or any oxide flush with or below the bare surface, never reaches the offset. On the PMOS gate the walk followed the spacer and the field oxide to the domain edge, then panicked.

Now the walk continues only while the oxide's offset-normal projection *strictly* increases. Otherwise the film ends with a positive-width end face at the last rising vertex M: X = M + (t − proj M)·n, a new gas-exposed vertex with a GAS node. The old edge leaving M is recorded. `native_cap_fill` certifies with exact predicates that this edge lies outside the cap's interior angle and crosses no cap edge. The existing global embedding check still runs. A flat plateau is not walked, because that would coat the old oxide's top with native oxide.

Test: `tests/run_native_oxide_end_face.py` covers a Si mesa between oxide fills that are flush, falling, rising, mixed and barely recessed. It checks full mesa coverage, that old fills are unchanged and uncoated, film area, and audits. The 11-18 build panics on the first case with the PMOS message.

## 2. Orientation decided inconsistently (`dbase/clktri.c`, `dbase/geometry.c`)

`sclk_tri` ordered vertices with a **float** determinant; `tri_geom` checks orientation in double. A Delaunay flip that `ck_flip` correctly accepted with exact predicates produced a valid sliver: two vertices one ULP apart, area 9.1e-16 cm². Its float determinant had the wrong sign, so the vertex order was reversed and `tri_geom` then panicked. Both now use `orient2d_float_sign`. `tri_geom` keeps double magnitudes and panics separately if the double area disagrees with the exact sign.

Test: `tests/run_sclk_orientation.sh` (ASan/UBSan) covers the captured sliver in all permutations, scalings and translations, plus 20,000 random near-collinear triangles. The original sources fail it. `run_geometry_weights.sh` and `run_p1_mass.sh` now link `geom/misc.c`.

## 3. Etch cutter coincident with an existing surface (`refine/etch.c`, `refine/skel.c`)

A spacer etch removes exactly the deposited thickness, so on flat faces its front coincides with the pre-deposit surface. A mask edge placed on a mesh line does the same. After storage in binary32 they agree only to rounding. The exact Boolean then found crossings one ULP beside existing vertices, leaving one-ULP edges in the poly top and field oxide. Those created the NMOS sliver above and many poor triangles downstream.

- **Snap rounding of the cutter only** (`snap_cutter`, after Hobby 1999 and Guibas–Marimont 1998). The mesh is never moved:
  - a cutter vertex that agrees with a subject vertex to storage rounding takes that vertex's coordinates;
  - a subject vertex ON a cutter edge is inserted into it;
  - cutter detours ON a subject edge are dropped;
  - consecutive cutter vertices closer than the deposition coat grid 2^(exponent(extent)−24) are merged.

  The snapped cutter is used only if exact predicates prove it simple.
- **One tolerance.** "Agrees to rounding" is the Boolean's own ON predicate, now shared through `skel.h` as `rounding_same_point` and `rounding_on_segment`. It is per axis and relative to that axis's coordinate magnitudes. Exactly the contacts the Boolean would call ON become exact. A real film is kept unless it is thinner than storage rounding of its own coordinates, and the real t(1−cos θ) offset of a vertical etch on a sloped face is never snapped. An absolute tolerance was tried first and **rejected**: it deleted the deliberate 1e-5 µm film in `run_process_surface_regression.py`.
- **`rounding_on_segment`** used an exact bounding box before its lateral band, which gave axis-aligned segments (mask edges) no lateral allowance at all. It now uses a perpendicular slab, 0 ≤ (p−a)·(b−a) ≤ |b−a|², which is exact in double, plus the unchanged lateral band. Points beyond an endpoint are still not ON (`etch_corner` `check_finite_clip_endpoints`).
- **Coincident boundary edges in `sub_skel`.** A subject edge running along a clip edge cannot be classified by a midpoint sample: the sample is ON, or rounding puts it on an arbitrary side and fragments the boundary into hundreds of arcs. Both skeletons share the `ck_clock` orientation, so the edge is removed if it runs the same way as the clip edge and retained if it runs the opposite way. An ON vertex whose incident edges share one class is not a junction and takes that class. An isolated retained coincident edge gets a midpoint with that exact class.

Test: `tests/run_etch_snap_rounding.py` builds a poly gate on gate oxide, deposits a spacer of thickness T and dry-etches exactly T, with T ∈ {0.1, 0.25, 0.5} µm, two origins, two y offsets and two mesh densities. It requires no vertex pairs closer than 1e-5 µm, no residue on flat faces, gate oxide and both spacers retained, and silicon and poly unchanged. Final build 24/24; 11-18 build 4/24. On the captured RIT spacer step, near-duplicate vertices fall from 28 to 1; the remaining one is a genuine 3-ULP offset on a reflecting side.

## 4. Thin dielectrics one element thick (`dbase/thin_film.c`, `diffuse/diffuse.c`)

A film meshed one element thick couples its two interfaces directly, so barrier transport through it is grossly wrong. Before every non-oxidizing diffusion, `thin_film_conform(L)` (default L = 8, `SUPREM4_THIN_FILM_LAYERS`) enforces at least L element layers across every SiO2/oxynitride film that separates two **distinct semiconductor regions**: gate oxide, interpoly, buried oxide. Gas carries no impurities and its surface is no-flux. Nitride, metal and resist are not dopant reservoirs. So a film facing only one reservoir cannot leak and is not refined. (An earlier version also refined against gas: every native/pad oxide went from 4.7k to 27k nodes and LOCOS stalled.)

For each neighbour region c, a film coordinate s_c = d_c/(d_c + d_other) is computed from Euclidean distances to interface segments. It depends only on film geometry, so the construction is idempotent. Edges with |Δs_c| > 1/L are split at the levels j/M, M = ceil(3L/2), by `sp_gredg`. That splits conformingly and linearly interpolates every field, so the P1 dopant representation is unchanged exactly. Splits that would leave a child triangle non-positive are refused (exact predicate). Up to three rounds follow Delaunay reconnection.

Test: `tests/run_thin_film_barrier.py` (RIT gate stack). Results, P in Si (cm^-2):

| Case | P in Si (cm^-2) | Edges across oxide | Surface type |
|---|---|---|---|
| One-layer oxide | 2.7e8 | 12 | p-type |
| Grown oxide | 6.6e7 | 12 | p-type |
| 16-layer reference | 7.9e7 | – | – |
| One-layer oxide, pass disabled (control) | 1.9e10 | 1 | inverted (n-type) |

**These numbers are superseded by fix 5.** They were computed with the legacy control volumes. Those inflated the layered oxide ~3×, so the "converged" layered values were ~100× too low. The control's agreement with the exact reference is a coincidence of two errors.

## 5. Control volumes of obtuse triangles (`dbase/geometry.c`)

`tri_geom` gave each vertex the circumcentric pieces `edge² · |cot(opposite angle)| / 8`. The original SUPREM-IV code is equivalent: `sqrt(fabs(r² − d²/4))`. For a non-obtuse triangle these are the Voronoi cells and sum to its area. For an obtuse triangle the circumcentre lies outside the triangle, and the |cot| pieces sum to *more* than the area: unboundedly more as the triangle flattens. `repair_obtuse`/CHP corrects only triangles whose *edge coupling* sum goes negative. Obtuse triangles whose couplings stay positive kept the inflated volumes. That inflates the time-derivative (mass) term and breaks dose conservation as the solver sees it.

The effect is concentrated where flat elements fan out:
- A refined interface over a coarser bulk. In the barrier deck, 65 interface nodes over 3 silicon nodes have volumes inflated 13–26×, and the phosphorus there was diluted to 1e13 between the original grid lines.
- Layered thin dielectrics. Gate oxide volumes were 2.9× too large after poly deposition and 3.2× in the barrier deck, which slowed and stored dopant inside the barrier.

Fix: keep the Voronoi pieces for non-obtuse triangles, unchanged. For obtuse triangles use Meyer's mixed area (Meyer, Desbrun, Schröder, Barr 2003): the obtuse vertex gets s/2 and the others s/4. Volumes are positive, each triangle contributes exactly its own area, and the formula is continuous at 90° (a right-angle vertex gets s/2 either way). Couplings (`ehed`) and the CHP repair are unchanged. Test `tests/geometry_weights.c` now checks area conservation *before* any CHP repair, including a flat fan element; it fails on the legacy formula.

**Exact reference.** The gate stack built from `line`/`region` statements has only right triangles. Its control volumes are exact under either formula, so legacy and fixed builds give identical results: Si P 1.95e10 cm^-2 and a flat interface at 9.96e15 cm^-3, above the surface B of 4.4e15. Against it (P in Si, cm^-2):

| Case | Legacy volumes | Mixed area |
|---|---|---|
| Right-triangle reference (thin-film pass off) | 1.95e10 | 1.95e10 |
| One-layer oxide | 2.7e8 | 1.23e10 |
| Grown oxide | 6.6e7 | 1.65e10 |
| 16-layer oxide | 7.9e7 | 1.17e10 |
| One-layer oxide, pass disabled | 1.9e10 | 2.4e11 (12× high: fix 4 is still needed) |

Lateral convergence in the laterally uniform one-layer deck (interface P max/min):

| Lateral spacing | Legacy | Mixed area |
|---|---|---|
| 0.5 µm | 126× | 38× |
| 0.1 µm | 18× | 2.5× |
| 0.05 µm | 5.2× | 1.4× |

The median rises from 2.0e15 to 6.0e15 with refinement, toward 1e16. The remaining spread at 0.5 µm is genuine silicon under-resolution: 3 silicon nodes span 1 µm at 5 nm depth. In the device, channel silicon volumes are already exact (ratio 1.00), so there the error was in the gate oxide.

`tests/run_thin_film_barrier.py` now uses the right-triangle reference. It requires the one-layer, grown and 16-layer cases within 2× of it and with matching surface polarity, and the disabled control above 5× of it. The old "p-type surface" criterion encoded the legacy undercount: at this test's 1e16 B, the exact answer is a weakly n-type surface.

Other observations from this audit, not bugs:
- Poly loses ~55% of its phosphorus to the ambient during the 117 min N2 activation. This is the stock `phosphorus poly /gas Trn.0=1.5 Trn.E=1.99` evaporation model, identical under every discretization.
- The interface segregation coupling (half-lengths of the faces incident at each node, `couple.c`) is correct.
- Neither the field term nor the background doping causes the lateral pattern. Legacy charged flux and B 1e13 give identical results.

## Regression battery (final build vs 11-18 build `5njvfhlt`)

| Test | Final | 11-18 |
|---|---|---|
| deposit_etch_identity, contact_shielding, etch_buffer, construction_perturbations, native_seed_loops, RIT NMOS/PMOS spacer replays | pass | pass |
| process_surface_regression | only the six 100 µm-wide cases exceed the test's 30 s limit | identical |
| directional_stack | three cases exceed the test's 120 s limit | identical |
| native_oxide_join | wedge vertex check fails | identical (pre-existing) |
| native_oxide_end_face (new) | pass | panic |
| etch_snap_rounding (new) | 24/24 | 4/24 |
| C unit tests (36 `run_*.sh`) | pass | – |
| thin_film_barrier (new) | pass | – |

directional_stack case 5 now finishes. It was one float ULP short of exact buried-oxide area (2.7e-9 of 0.04 µm^2). Snap rounding routes the poly-edge cut through an existing Si-interface vertex 2.9e-7 µm away, instead of creating a sliver edge there. The test now bounds the area by the Boolean's ON band times the cut height (5.7e-8 µm^2) instead of 1e-12. The 11-18 build was exact here but created the 3e-7 µm interface edge.

`run_mesh_memory.sh` fails to link at HEAD (`kill_node.c` needs `ig2_write`, `orient2d_float_sign`). That predates this work.

Diagnostics: `SUPREM4_SNAP_DIAG=1` prints cutter snap counts and the boundary class runs per subject. `SUPREM4_MESH_DIAG=1` also traces the native-cap walk.

## Full flow with fix 5 (2026-11-20)

Build `native64_u85__kd9` (fix 5 included). Resumed from the saved `gate_oxide` checkpoint (`run.py --resume results/fresh_full_flow_mixedarea_20261120 --after gate_oxide --timeout 10800`). Output: `decks/subcmos2014/results/fresh_full_flow_mixedarea_resume_20261120`. NMOS and PMOS both PASS, 16/16 checkpoints. Elapsed 3894 s (NMOS) and 3571 s (PMOS). Poly activation alone is ~184 steps on ~95k nodes, ~42 min, at the same per-step cost as the previous build. The `run.py` default `--timeout 1800` is therefore too short for the full flow.

Electrical results (Genius DD, `experiments/electrical_mixedarea_20261120`) against the previous build (`electrical_goxcal_20261120`):

| Metric | NMOS before | NMOS fix 5 | PMOS before | PMOS fix 5 |
|---|---|---|---|---|
| Vt max-gm (V), target ±0.75 | 0.373 | 0.126 | −1.057 | −0.658 |
| Vt at 100 nA, Vd=0.05 V (V) | 0.363 | 0.122 | −0.736 | −0.342 |
| SS (mV/dec) | 86.4 | n/a | 253.6 | 225.8 |
| Id at Vg=0, Vd=0.05 V (A) | 8.4e-12 | 4.6e-9 | 6.5e-10 | 4.6e-9 |
| Ion at Vg=Vd=3.3 V (A) | 2.48e-4 | 2.74e-4 | 1.17e-4 | 1.27e-4 |

NMOS body P 50 nm below the gate after poly activation: 4.6e13 → 1.8e15 cm^-3 (40× more). This matches the barrier test, where fix 5 raised P penetration through the gate oxide by ~45× to agree with the exact right-triangle reference. The NMOS Vt drop is phosphorus from the poly penetrating the gate oxide: legacy control volumes had been hiding it. The numerics now agree with the exact reference, so the remaining Vt gap is set by the oxide P transport model (stock diffusivity and segregation) and the recipe, not by the mesh.

## Phosphorus diffusivity in oxide set to TSUPREM-4 default (2026-11-20)

Stock SUPREM-IV uses D(P in SiO2) = 7.6e-3 exp(-3.5 eV/kT) cm^2/s. The TSUPREM-4 2002.4 default (Appendix A, Table A-3, from Ghezzo & Brown, JECS 120, 146, 1973) is 4.338e10 um^2/min with 4.44 eV, i.e. 7.23 cm^2/s with 4.44 eV. That is 5.5x lower at 1000 C and 27x lower at 800 C. Shimakura et al. (SSE 18, 991, 1975) measured poly-source P through SiO2 in N2 at ~3e-18 cm^2/s at 1100 C, about 100x below even TSUPREM-4. Segregation (m=30 at Si/oxide and poly/oxide) and transport (0.1 um/min = 1.66e-7 cm/s) already equal TSUPREM-4's, which notes reported n-type values span 10 to 1000.

`build_decks.py` now emits `phosphorus oxide Dix.0=7.23 Dix.E=4.44` after the `method` line. `run.py --resume` now keeps impurity-parameter lines (boron/phosphorus/arsenic/antimony/interstitial/vacancy) in the resumed header.

Gate-stack barrier deck with the override: P in Si falls from 1.95e10 to 4.6e5 cm^-2 for the right-triangle reference, and to 1e4–2e4 for the one-layer, grown and 16-layer cases. All keep a p-type surface. The one-element oxide without thin_film_conform still leaks 1.0e10 cm^-2, so fix 4 remains necessary.

Full flow (`results/fresh_full_flow_tsupox_20261120`, build `native64_u85__kd9`): both PASS 16/16, NMOS 4916 s, PMOS 4592 s. NMOS P 50 nm under the gate after poly activation: 1.8e15 -> 5.5e9 cm^-3. NMOS channel surface after metal1: P 1.3e17 -> 2.7e14 cm^-3, surface net B 4.3e16.

| Metric | NMOS fix 5 | NMOS + TSUPREM-4 oxide | PMOS fix 5 | PMOS + TSUPREM-4 oxide |
|---|---|---|---|---|
| Vt max-gm (V), target ±0.75 | 0.126 | 0.419 | −0.658 | −0.234 |
| Vt at 100 nA, Vd=0.05 V (V) | 0.122 | 0.410 | −0.342 | not reached (always above 100 nA) |
| SS (mV/dec) | n/a | 86.5 | 226 | n/a |
| DIBL (mV/V) | 8.4 | 9.9 | – | – |
| Id at Vg=0, Vd=0.05 V (A) | 4.6e-9 | 2.3e-12 | 4.6e-9 | 3.6e-7 |
| Id at Vg=0, Vd=3.3 V (A) | 8.5e-9 | 4.0e-12 | 1.0e-6 | 5.4e-6 |
| Ion at Vg=Vd=3.3 V (A) | 2.74e-4 | 2.43e-4 | 1.27e-4 | 1.48e-4 |

PMOS no longer turns off. The n+ poly PMOS is buried-channel: the blanket plus PMOS boron Vt adjust (0.7e12 + 2.6e12 at 60 keV, PDF p.75–78) forms a p-layer from 35 to 215 nm below the surface already at the threshold and gate-oxide checkpoints, under an n-type surface held by P pile-up. After poly activation and the S/D anneal it now merges into one p-layer 0–230 nm deep (surface net −2.3e16, B ~1.1e17 over n-well P ~7–8e16), which the gate cannot pinch off at Vg=0. In every earlier run, phosphorus leaking from the gate kept the PMOS surface n-type (+1.8e17 cm^-3 in the legacy build), and that leak supplied the off-state. The open question is the n-well surface doping and depth (and the boron retained through the Kooi oxide), not the gate stack.

### PMOS Vt-adjust dose versus the PDF's design assumptions (2026-11-21)

PDF p.74 (MOSFETVT.XLS) sizes the Vt adjust with Nd = 1e17 cm^-3, Nss = 3e11 cm^-2,
Xox = 150 A, n+ poly, and **"assume 1/2 dose in Si"** (x2 factor). For a desired
-1 V it gives 2.65e12, so the PDF's 0.7e12 + 2.6e12 = 3.3e12 B is intended to leave
~1.65e12 in silicon.

Simulated (tsupox run, PMOS, x = 6 um):

| stage | B in Si (cm^-2) | B in oxide | net donor at 2/50/100/200/300 nm |
|---|---|---|---|
| locos | 7.7e11 (background) | - | +6.4e16 / +6.3e16 / +6.3e16 / +6.2e16 / +6.1e16 |
| threshold | 3.97e12 | 1.5e11 | +1.3e17 / -4.4e16 / -1.2e17 / -4.4e16 / +6.0e16 |
| gate_oxide | 3.77e12 | 1.9e11 | +9.7e16 / -4.4e16 / -1.1e17 / -3.0e16 / +5.6e16 |
| annealed | 3.72e12 | (gate ox) | -2.3e16 / -4.0e16 / -4.5e16 / -1.6e16 / +2.7e16 |

- 3.2e12 of the 3.3e12 implant (97 %) enters silicon (60 keV B, Rp ~0.19 um, through
  ~104 nm Kooi oxide); 2.95e12 remains after anneals. The "1/2 in Si" assumption does not hold.
- The n-well surface (6.3e16) is lighter than the spreadsheet's 1e17 but heavier than the
  p.10 summary's 3e16; it is not the cause.
- Genius uses QF = 0, whereas the spreadsheet uses Nss = 3e11 (-0.21 V).
- First-order estimate: Vt0 (Nd 6.3e16, QF 0) ~ -1.55 V; +q*2.95e12/Cox = +2.05 V gives
  Vt ~ +0.5 V (depletion mode), consistent with the PMOS that does not turn off. With
  the spreadsheet's assumptions (1.65e12 in Si, Nss 3e11) the same estimate gives about -0.6 V.
- A -0.75 V PMOS with full retention needs ~1.15e12 total B in Si (QF = 0),
  i.e. PMOS-only implant ~0.45e12, or ~0.8e12 if QF = 3e11 is added.

### PMOS Vt adjust 0.45e12 variant (2026-11-21)

`build_decks.py --pmos-vt-dose 0.45e12 --output variants/pmos_vt_0p45e12` (default decks unchanged,
PDF 2.6e12). PMOS resumed after `locos` from the tsupox run:
`results/pmos_vt_0p45e12_20261121` (PASS 16/16). Genius: `suprem4/experiments/electrical_pmosvt045_20261121`.

| PMOS | 2.6e12 (PDF) | 0.45e12 |
|---|---|---|
| Vt max-gm (V) | -0.23 | -1.31 |
| Vt 100 nA, Vd=-0.05 (V) | +0.35 | -1.34 |
| Vt 100 nA, Vd=-3.3 (V) | - | -1.26 |
| SS (mV/dec) | - | 81.4 |
| DIBL (mV/V) | - | 22.5 |
| Id(Vg=0, Vd=-0.05) (A) | 3.6e-7 | 6.6e-16 |
| Ion (Vg=Vd=-3.3) (A) | 1.48e-4 | 7.9e-5 |

Surface-channel enhancement PMOS, but overshoots the -0.75 V target. Vt moves 0.78 V per 1e12 cm^-2
of PMOS dose (q/Cox = 0.70 V per 1e12), so -0.75 V needs about 1.2e12. The hand estimate of Vt0 was too
optimistic because P pile-up leaves the n-well surface near 1e17.

### PMOS Vt adjust 1.2e12 variant (2026-11-21)

`results/pmos_vt_1p2e12_20261121` (PASS 16/16); Genius `suprem4/experiments/electrical_pmosvt12_20261121`.

| PMOS dose (60 keV) | 0.45e12 | 1.2e12 | 2.6e12 (PDF) |
|---|---|---|---|
| Vt 100 nA, Vd=-0.05 (V) | -1.34 | -1.10 | +0.35 |
| Vt max-gm (V) | -1.31 | -1.07 | -0.23 |
| SS (mV/dec) | 81.4 | 77.7 | - |
| DIBL (mV/V) | 22.5 | 28.9 | - |
| Id(Vg=0) (A) | 6.6e-16 | 9.0e-16 | 3.6e-7 |
| Ion (A) | 7.9e-5 | 1.01e-4 | 1.48e-4 |

Vt vs dose is strongly nonlinear: 0.32 V per 1e12 between 0.45e12 and 1.2e12, then 1.0 V per 1e12 up to 2.6e12.
At 60 keV through the 105 nm Kooi oxide the boron peak sits ~130 nm deep, beyond the surface depletion
width (~0.1 um at 1e17), so it barely shifts Vt until it forms a buried p-layer. A shallower implant
(30-35 keV: 52-67 % retained, peak at 5-15 nm) acts as a surface charge sheet and is the regime the
PDF p.74 spreadsheet assumes.

### PMOS Vt adjust energy variants, 2.6e12 at 30 / 35 keV (2026-11-21)

`build_decks.py --pmos-vt-energy` (default 60). Results `results/pmos_vt_2p6e12_{30,35}kev_20261121` (both PASS 16/16);
Genius `suprem4/experiments/electrical_pmos_vt_2p6e12_{30,35}kev_20261121`.

| PMOS Vt implant | 0.45e12/60 | 1.2e12/60 | 2.6e12/30 | 2.6e12/35 | 2.6e12/60 (PDF) |
|---|---|---|---|---|---|
| Vt 100 nA (V) | -1.34 | -1.10 | -1.11 | -0.885 | +0.35 |
| Vt max-gm (V) | -1.31 | -1.07 | -1.08 | -0.86 | -0.23 |
| SS (mV/dec) | 81 | 78 | 80.5 | 97 | - |
| DIBL (mV/V) | 22.5 | 28.9 | 26.5 | 43 | - |
| Id(Vg=0) (A) | 6.6e-16 | 9.0e-16 | 6.9e-16 | 1.2e-15 | 3.6e-7 |
| Ion (A) | 79e-6 | 101e-6 | 97e-6 | 117e-6 | 148e-6 |

The 1000 C poly and S/D anneals spread the Vt boron over ~100-150 nm whatever the implant depth, so a shallow
implant gives no special Vt efficiency. At 35 keV a weak p layer (-2.3e16, 0-150 nm) forms, but the n+ gate
still depletes it at Vg=0: the device turns off, but SS and DIBL degrade (onset of buried-channel behaviour).
Reaching -0.75 V with n+ poly needs a deeper counter-doped layer, trading SS/DIBL for Vt.

**Decision (2026-11-21):** default `pmos.in` now uses the PDF dose 2.6e12 at 35 keV (`--pmos-vt-energy 60` restores
the PDF energy). `nmos.in` unchanged.

### NMOS Vt diagnosis (2026-11-21)

NMOS channel boron (tsupox run, x = 6 um): after anneals 4.7e16 at 2 nm, ~7e16 over 20-100 nm, 5.4e16 at 300 nm
(p-well alone 4.3-4.9e16 after LOCOS). Hand Vt with n+ poly, Na = 7e16, tox 15 nm, QF = 0:
phi_ms = -0.96 V, 2phi_F = 0.80 V, Qd/Cox = 0.59 V -> Vt = 0.43 V, matching Genius 0.42 V. +0.32 V to the
0.75 V target needs roughly +0.5-0.6e12 more boron in Si. The only NMOS Vt implant is the shared blanket
7e11 at 60 keV, so raising it also shifts the PMOS, which then needs its own dose reduced.

### Blanket Vt dose, both devices (2026-11-21)

Shared blanket B at 60 keV varied in both decks (PMOS-only 2.6e12 at 35 keV). `variants/blanket_{1p3e12,1p6e12}`,
`results/blanket_*_20261121` (all four PASS 16/16), Genius `suprem4/experiments/electrical_blanket_*_20261121` (6/6 each).
The 7e11 column uses the tsupox NMOS and the 35 keV PMOS run, which have the same process.

| blanket (60 keV) | 7e11 (PDF) | 1.3e12 | 1.6e12 |
|---|---|---|---|
| NMOS Vt 100 nA (V) | 0.41 | 0.51 | 0.56 |
| NMOS SS / DIBL | 86.5 / 9.9 | 90.8 / 9.1 | 92.9 / 8.9 |
| NMOS Id(Vg=0) (A) | 2.3e-12 | 3.5e-13 | 1.5e-13 |
| PMOS Vt 100 nA (V) | -0.885 | -0.414 | -0.31 |
| PMOS SS / DIBL | 97 / 43 | 107 / 81 | 113 / - |
| PMOS Id(Vg=0) (A) | 1.2e-15 | 3.4e-11 | 3.8e-10 |

NMOS gains only ~0.16 V per 1e12: the 60 keV boron spreads past the depletion width and Vt grows as sqrt(Na).
PMOS loses ~0.8 V per 1e12 because the same boron counter-dopes the n-well. NMOS cannot reach 0.75 V by the
blanket alone without overshooting PMOS; the PMOS-only dose must come down to compensate. Estimate: PMOS at
-0.75 V needs ~2.6e12 B retained in Si (blanket x0.95 + PMOS-only x0.67 at 35 keV), so PMOS-only reaches 0
near blanket ~2.7e12, where NMOS is only ~0.65-0.7 V.
