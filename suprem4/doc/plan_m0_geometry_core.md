# Implementation plan: M0 geometry/mesh core (draft, 2026-11-22)

Milestone M0 of `roadmap_2d_cmos_25nm.md`.  Goal: a robust geometry/mesh core that keeps the mesh valid
through moving-boundary oxidation, demonstrated on the poly re-oxidation cases that defeated the
SUPREM-IV core, and a go/no-go decision on building the rest of the roadmap on it.

## 1. Design principles
1. **Geometry is the source of truth.**  A boundary representation (B-rep) holds the material regions
   and interfaces; the mesh is derived and disposable.  Topology changes are exact Boolean operations
   on the B-rep (Clipper2, integer coordinates), never surgery on the mesh database.
2. **Mature libraries do the meshing.**  Constrained Delaunay generation (CDT 1.4.5, already vendored)
   plus Mmg2d (LGPL) for quality remeshing and adaptation with required interface edges and region
   references.  No hand-written triangulation or edge-collapse code.
3. **SUPREM-IV stays the physics engine for M0.**  The core replaces only the geometry step.  The mesh is
   handed back to SUPREM through its own structure reader (whole-mesh reload), not by editing its
   database in place.
4. **Invariants are checked every step**; failures stop with a dump, never silently continue.
5. **Off by default.**  `SUPREM4_GEOMETRY_CORE=1` (and only for oxidation in compress/viscous) until
   accepted; every existing deck must stay bit-identical with it off.

## 2. Architecture
```
SUPREM-IV time step (oxgrow)                     suprem4/core (C++17, C ABI)
  velocities, fields  --- export (in memory) --> 1. mesh + fields + velocities
                                                 2. B-rep from mesh (region loops, interface chains)
                                                 3. advance B-rep by velocities (step-limited)
                                                 4. detect events (gap/thickness < delta, inversion)
                                                 5. repair B-rep (Clipper2 open/close/difference)
                                                 6. mesh: CDT + Mmg2d (boundary layers, quality)
                                                 7. transfer fields (conservative per material)
                                                 8. invariant checks
  reload mesh via structure reader <--- import --- new mesh + fields (+ velocity history)
  re-initialise diffusion state, continue
```
Steps 3-7 run only when needed: normal steps keep SUPREM's own Lagrangian motion; the core is invoked
when the step-limited motion would degrade quality below a threshold or an event is detected
(the Sentaurus "Plan A local / Plan B regenerate" split, with Mmg2d providing the local part).

## 3. Data model (core)
- `Point` (int64 x,y on a 1e-11 cm grid = Clipper2 grid; double view).
- `Region` {material, outer loop, holes} - holes allowed in the B-rep; split only at mesh time if the
  consumer (SUPREM) needs simply connected regions.
- `Chain` - polyline shared by two regions or by a region and the ambient/backside/side wall, with a
  boundary type (interface, exposed, backside, reflecting).  Shared chains guarantee conformity.
- `Mesh` {points, triangles with region id, per-(point,material) node values}; multi-valued nodes at
  interfaces as in SUPREM.
- `Fields` - solution slots by name (dopants, oxidant, point defects, velocity, velocity history),
  conservation class per slot (conserved / interpolated / ambient-valued).

## 4. Work packages
| WP | Content | Deliverable | Effort |
|---|---|---|---|
| WP0 | Build setup: `suprem4/core/` CMake project; vendor Mmg (pinned release >= 2 weeks old, LGPL) and CDT; hook into `tools/build_native64.py` as `--geometry-core`; unit-test framework (Catch2 or plain asserts) | Core library builds and links into suprem4 | 2 d |
| WP1 | Mesh I/O: read/write SUPREM .str in the core (C++), plus in-memory export/import from SUPREM's database; round-trip test (str -> core -> str identical) | Round-trip passes on all RIT checkpoints | 3 d |
| WP2 | B-rep extraction and validation: region loops, shared chains, boundary types, hole handling; invariant checks (closed, simple, conforming, no overlap) | B-rep of every RIT checkpoint validates | 3 d |
| WP3 | Mesh generation: chain discretization from a spacing field (existing vertices kept where requested), anisotropic boundary layers for films thinner than a threshold, CDT + Mmg2d quality remesh with required edges and region references; minimum angle / quality targets | Remesh of every RIT checkpoint: quality >= target, area/interface conformity exact | 5 d |
| WP4 | Field transfer: P1 interpolation per material from the old mesh, per-material inventory restoration for conserved species, ambient values for new gas nodes, velocity-history transfer; transfer error metrics | Remesh changes dose by < 1e-6 relative, profiles within interpolation error | 3 d |
| WP5 | Motion + events: apply velocities with step limit (<= 0.3 x local element height, <= 1 nm growth), event detection on the B-rep (gap / thickness < delta = 1 nm, predicted inversion), Clipper2 repair (opening of consumed materials, closing of oxide, priority differences), port and clean up the repair logic from `doc/poly_oxidation/patches` | Synthetic moving-boundary tests (shrinking square, merging fronts, slit closure) pass | 5 d |
| WP6 | SUPREM integration: call site in `oxgrow` (behind `SUPREM4_GEOMETRY_CORE`), in-memory export, reload through the structure reader, re-initialisation of diffusion/oxidation state (`init_diffuse`, symbolic factorisation, velocity history), compress then viscous | Poly-line base deck completes with the core active | 5 d |
| WP7 | Test suite: poly-line matrix (6 decks), LOCOS and recessed-oxide corners, narrow windows, thin films, 50+ randomized geometries (gate width, film thickness, time, temperature), RIT flows with core off (bit-identical) and on; invariants after every step; runtime budget | `tools/run_geometry_suite.py` with pass/fail report | 4 d |
| WP8 | Review: results vs acceptance criteria, measured effort, go/no-go recommendation for roadmap option C | Review note | 1 d |

Total about 31 working days; plan 4-6 weeks of focused work including integration surprises.
Order: WP0 -> WP1 -> WP2 -> WP3 -> WP4 (WP7 harness grows alongside) -> WP5 -> WP6 -> WP7 -> WP8.
WP1-WP4 are independent of oxidation and are useful on their own (a robust remesh command for any
SUPREM structure, e.g. before device export).

## 5. Acceptance criteria (go/no-go)
1. All six poly-line decks (`doc/poly_oxidation/test_decks`) complete, plus 50 randomized variants.
2. Physics agrees with the validated compress results: consumption = 0.44 x oxide on every face
   (+-5 %), poly 0.400 -> 0.355 um, sidewall oxide 500 A, gate-edge dielectric within 15 % of the
   reference profile; results converge under mesh refinement (two mesh densities within 5 %).
3. Conservation: dopant dose drift < 1e-4 relative over a full oxidation; material areas consistent with
   the volume expansion ratio.
4. Invariants hold after every step (orientation, conformity, embedding, region consistency, surface
   flags and gas nodes).
5. With the core off, all RIT 0.5 um and 0.25 um decks are bit-identical to the current production
   build; with it on, they run and pass their checkpoint audits.
6. Runtime of a poly re-ox step at most 3x SUPREM-IV's own (vertical-model) time.

## 6. Risks and mitigations
| Risk | Mitigation |
|---|---|
| Reloading the mesh mid-diffusion leaves stale solver state in SUPREM | WP6 starts with a "reload and continue" test with no geometry change; compare to an uninterrupted run |
| Mmg2d does not preserve thin layers or interface vertices as needed | Required edges/vertices and per-region references; boundary layers generated before Mmg; fall back to CDT + own quality refinement |
| Viscous model uses quadratic elements (FEconvert) | M0 targets compress first; viscous in WP6 only after compress passes |
| Field transfer at interfaces (segregation jumps) | Per-material transfer only; interface nodes keep per-material values |
| Effort overrun | WP8 review with measured progress; WP1-WP4 still deliver a standalone remesher |

## 7. Licences
Clipper2 (Boost), CDT (MPL-2.0, file-level), Mmg (LGPL-3, link dynamically or keep sources unmodified
and separate).  Avoid Triangle (non-commercial) and GPL meshers in the core unless the project
licence is decided to allow them.

## 8. What is reused
- `doc/poly_oxidation/patches/*repair_module*`: Clipper2 opening/closing, sliver filtering,
  hole cutting, global vertex thinning, gas-node ambient values, inventory restoration.
- `doc/poly_oxidation/README.md` failure catalogue (items 1-15) becomes WP7 regression cases.
- The poly oxidation physics patch (compress model) is applied in WP6.
