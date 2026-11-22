# Roadmap: a 2D TCAD flow for silicon CMOS down to 25 nm (draft, 2026-11-22)

## 1. Goal and scope
Simulate complete 2D silicon CMOS process flows and the resulting devices, from the RIT 0.5 um and
0.25 um processes down to 25 nm planar bulk and SOI technology (the 32/28 nm generation: high-k/metal
gate, strain engineering, millisecond anneals), with results calibrated against published data.

In scope: planar bulk, PD/FD-SOI, 2D cross-sections (W assumed), process + device + calibration.
Out of scope: FinFET, gate-all-around, random dopant fluctuation and any other inherently 3D effect;
these need a 3D code base and are a separate decision.

## 2. Where we are
| Asset | State |
|---|---|
| SUPREM-IV.GS (modified, native 64-bit) | 2D process: analytic Pearson implant, Fermi/point-defect diffusion, vertical/compress/viscous oxidation, etch/deposition (Clipper2 deposition backend). Full RIT 0.5 um flows run; AdvCmos 0.25 um flow running. Fixes to date: mesh/FE repairs, TSUPREM-4 P-in-oxide, SRIM P-in-nitride ranges |
| Genius (open) | Device: drift-diffusion (ddm1/ddm2), density gradient (dg), energy balance (ebm3), AC, mixed-mode, stress solver; Si mobility Lombardi/Philips/HP/Darwish/Lucent/Sdevice; Schenk band model, generation models, gate boundary tunnelling hooks; materials incl. PolySi, TiSi2, W, Nitride (no HfO2, SiGe, NiSi) |
| Workflow | Deck generators, runners with checkpoint audits, SUPREM->Genius export, Vt/SS/DIBL extraction, regression folders |
| Parked | 2D poly oxidation physics + repair/remesh experiments: suprem4/doc/poly_oxidation/ |

## 3. The central lesson so far
Adding physics to SUPREM-IV is fast (poly oxidation physics, SRIM tables, P-in-oxide were hours each).
The **geometry/mesh core** is the bottleneck: its moving mesh, region bookkeeping and triangulator were
built for flat, vertically oxidized structures and fail on non-trivial topology (poly re-oxidation,
narrow windows, slivers).  Production tools use a Lagrangian moving mesh with disciplined local
operations plus whole-mesh regeneration from a boundary representation as fallback, and even Synopsys
documents unsolved cases (suprem4/doc/poly_oxidation/research_oxidation_meshing.md).  Nearly every
25 nm feature (epi, silicide, STI, eSiGe, spacers, CMP, gate-last) changes geometry, so this core
must be solved first.

## 4. Base-framework decision (needed before Phase 1 of the geometry work)
| Option | Description | For | Against |
|---|---|---|---|
| A. Extend SUPREM-IV | Replace its mesh/geometry core in place; add physics to its solver | Keeps working flows; physics reuse is direct | 1990s data structures (fixed MAXMAT=10, region/skeleton DB, K&R C); every core change risks regressions; mechanics and topography are new code anyway |
| B. Adopt an open framework | ViennaPS/ViennaLS (topography, level sets, GPL-3), FLOOXS (academic licence), DEVSIM (device, Apache) | Modern level-set topography exists (ViennaPS) | No complete open 2D CMOS process simulator; licences constrain reuse; their diffusion/defect physics is thin |
| C. New core, ported physics | Modern C++ core: B-rep geometry (Clipper2/CGAL), mesh regeneration (Triangle/Gmsh/CDT), FVM diffusion and FEM mechanics on PETSc; port SUPREM/TSUPREM-4 models as specifications; keep SUPREM-IV as the regression reference and Genius for devices | Designed for the target; geometry robust by construction; no material limit | Largest up-front cost; months before it reproduces today's flows |

**Recommendation:** C, staged so value arrives early: build the geometry/mesh core first and drive it
with SUPREM-IV's existing physics through an interface (or file exchange) so the RIT flows keep running,
then move physics modules over one at a time, each validated against SUPREM-IV and literature.
Continue adding junction physics (Phase 1) to SUPREM-IV meanwhile, because those models port directly.
Revisit this decision after a 2-4 week geometry-core prototype (milestone M0).

## 5. Feature inventory
Status: have / partial / missing.  Effort: S (days), M (1-3 weeks), L (1-2 months), XL (several
months) of focused work; rough, for planning only.

### 5.0 Foundation
| Feature | Status | Effort | Notes |
|---|---|---|---|
| Robust geometry/mesh core (B-rep, local ops + regeneration, adaptive refinement) | missing | XL | Prerequisite for 5.3-5.7 |
| Mechanics solver with stress history (intrinsic, thermal mismatch, lattice mismatch, relaxation) | partial (oxidation flow only) | L | Prerequisite for strain |
| General materials (SiGe(x), Si:C, NiSi, TiN, HfO2, SiON, nitrides) | missing (MAXMAT 10) | M | Part of core |
| Time integration for spike/flash/laser anneals (us-ms) | partial | M | |

### 5.1 Implantation
| Feature | Status | Effort |
|---|---|---|
| Analytic Pearson tables (SRIM-based for missing pairs) | partial (P-in-nitride done) | S |
| Tilt/rotation, halo/pocket (dual-Pearson or MC) | missing | M |
| Monte Carlo / BCA implant with channeling | missing | L |
| Damage, amorphization, pre-amorphization (Ge/Si PAI) | partial (damage option) | M |
| Species: In, F, N, C, BF2/molecular | partial | S-M |

### 5.2 Diffusion and activation
| Feature | Status | Effort |
|---|---|---|
| OED/ORD (oxidation injection, surface recombination) | in code, untested | S |
| Basic TED (+n damage into interstitial field) | partial | M |
| {311} clusters, dislocation loops, boron-interstitial clusters | missing | M-L |
| Solid-phase epitaxial regrowth with activation, end-of-range defects | missing | M |
| Transient activation/deactivation (As/P/B, metastable after ms anneals) | partial (As clustering) | M |
| Co-dopants: C traps interstitials, F, N | missing | M |
| Diffusion in SiGe (Ge-dependent B/As) | missing | M |
| Interface trap/segregation dose loss | missing | S-M |

### 5.3 Oxidation and dielectrics
| Feature | Status | Effort |
|---|---|---|
| Thin-oxide regime (Massoud), N2O/NO nitridation | missing | S-M |
| Stress-dependent viscoelastic oxidation (STI liner, corners) | partial (viscous) | M |
| Poly oxidation | physics done, blocked on 5.0 | S after 5.0 |
| Deposited films with intrinsic stress (SiN liners, high-k stack) | missing | S after mechanics |

### 5.4 Epitaxy
Selective Si/SiGe/Si:C, raised source/drain, in-situ doping, facets by orientation, loading.
Status missing; effort M-L; needs 5.0.

### 5.5 Silicide
Ni/Co/Ti deposition and reaction, Si consumption, phases, dopant segregation and snow-plow, stress.
Status missing (TiSi2 skipped in AdvCmos); effort M-L; needs 5.0 and mechanics.

### 5.6 Strain engineering
eSiGe source/drain, stressed nitride cap (CESL), stress memorization, STI stress.
Status missing; effort M once mechanics + epi exist.

### 5.7 Topography
Anisotropic/isotropic/recess/facet etch, conformal and directional deposition, CMP, STI, gate-last
replacement gate.  Status partial (vertical etch, conformal deposition); effort L; level-set approach
recommended, building on the Clipper2 backend.

### 5.8 Device (Genius)
| Feature | Status |
|---|---|
| Quantum confinement (density gradient) | have (dg solver); calibrate for ~1 nm EOT |
| Stress-dependent mobility (piezoresistance / deformation potential) | partial (stress solver, some strain terms); verify |
| Inversion-layer mobility, high-k remote phonon | partial (Lombardi/Darwish); high-k missing |
| Velocity saturation, energy balance | have (ebm3) |
| Gate tunnelling, BTBT/GIDL | partial (hooks/band models); verify |
| Metal-gate work function, interface traps, fixed charge | have (QF used), traps partial |
| Materials HfO2, SiGe, NiSi, TiN | missing | 
| Process-to-device export with stress and new materials | partial |

### 5.9 Calibration and validation
Reference data set (literature SIMS for OED/TED/activation, Vt roll-off, Ion/Ioff, C-V for 0.25 um,
90 nm, 45 nm, 32 nm class devices), automated regression of every flow, parameter-extraction tools, a
versioned calibration database.  Status partial (RIT regression flows); effort M, ongoing.

## 6. Phased plan and milestones
| Phase | Content | Milestone (demonstration) |
|---|---|---|
| M0 (2-4 wk) | Geometry-core prototype: B-rep + regeneration + field transfer; drive one SUPREM-IV oxidation case; decide A vs C | Poly re-ox test matrix (suprem4/doc/poly_oxidation/test_decks) completes robustly |
| 1 | Junction physics in SUPREM-IV: OED, TED with clusters, BIC, activation, tilt/halo, PAI, co-dopants | RIT 0.25 um junctions vs literature-type SIMS; RSCE trend |
| 2 | Core + materials + mechanics; port diffusion/oxidation to the new core | RIT flows reproduced on the new core (bit-level vs SUPREM-IV where models agree) |
| 3 | Topography (STI, etch/dep/CMP), epi, silicide, stressed films | 90 nm class planar CMOS cross-section with STI, spacers, silicide |
| 4 | Strain (eSiGe, CESL, SMT), spike/flash anneal kinetics | 45 nm class device with strain; Ion gain vs literature |
| 5 | Device: quantum, stress mobility, high-k/metal gate, tunnelling, GIDL | 32/28 nm planar bulk and FD-SOI demo; Vt roll-off, Ion/Ioff, C-V vs published |
| ongoing | Calibration database, regression, documentation | Every milestone re-runs all previous flows |

Rough total: a multi-year effort at full scope; Phases 0-1 deliver usable value within a few months.

## 7. Risks
- Geometry core proves harder than the prototype suggests (mitigation: M0 go/no-go).
- Calibration data for 25-45 nm processes is scattered and partly proprietary (mitigation: build the
  reference set from published papers early; accept trend-level accuracy where data is thin).
- Numerics: coupled point-defect/cluster models and millisecond anneals stress the solvers (mitigation:
  1D benchmarks before 2D flows).
- Licences of third-party components (GPL ViennaPS, MPL CDT, academic FLOOXS) constrain reuse.
- Scope creep into 3D; keep it an explicit separate decision.

## 8. Immediate next steps
1. Finish the AdvCmos 0.25 um flow (geometric re-ox run) and commit the production changes.
2. OED benchmark in SUPREM-IV (point-defect model on, buried boron marker under dry oxidation).
3. M0 geometry-core prototype scoping note (data model, mesher choice, field-transfer method) for
   review before any code.
