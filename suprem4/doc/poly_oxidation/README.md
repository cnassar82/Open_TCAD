# 2D polysilicon oxidation in SUPREM-IV.GS: what was learned (2026-11-22)

Status: **parked**.  The AdvCmos2014 flow uses the geometric re-oxidation in
`decks/advcmos2014/build_decks.py --poly-reox geometric`, which runs on the unchanged production code.
This folder keeps everything needed to resume a real 2D implementation.

| File | Content |
|---|---|
| `research_oxidation_meshing.md` | How Sentaurus/TSUPREM-4/MovingMesh/diffuse-interface tools keep the mesh valid (sources, quotes) |
| `patches/repair_module_tracked.patch` | Poly oxidation physics + Boolean repair/remesh hooks (diff against commit 8d48f8e) |
| `patches/repair_module_new_oxide_repair_{c,cpp}.patch` | The repair/remesh module (new files) |
| `patches/third_party_cdt_1.4.5.tar.gz` | Vendored CDT 1.4.5 (MPL-2.0) used by the module |
| `patches/local_fixes_v1.patch` | First attempt: local mesh fixes (relabel, merges, notch/contact zippers) |
| `test_decks/poly_line_*.in` | 0.4 um poly line on 100 A oxide; base + 5 variations (the robustness matrix) |

Apply with `git apply patches/repair_module_tracked.patch`, the two new-file patches, and untar CDT into
`suprem4/src/third_party/`; build with `tools/build_native64.py --domain-offset --ordering amd --fe-solver umfpack`.
Run test decks with `SUP4KEYFILE/SUP4MODELRC/SUP4IMPDATA` pointing at the build's `data/`, the
`compress` model is selected in the deck.

## 1. Background
- SUPREM-IV.GS never oxidizes poly: every oxidation model reacts only at Si/SiO2.  Hooks for poly exist
  (`oxgrow.c` moves points with a Poly node, `FEconvert.c` reads the poly Fermi level) but nothing was
  finished; the Stanford/Cogenda tree (github.com/cogenda/Suprem4) is the same.
- TSUPREM-4 oxidizes poly only in COMPRESS, VISCOUS, VISCOELA (not VERTICAL).  Defaults (App. A): poly
  linear rate = <110> values (L.LIN.0 8.65e4 um/min O2, 2.87e4 H2O; E 2.0/1.6 eV), parabolic as Si,
  ALPHA poly/oxide 0.44, grain-boundary enhancement GB.FAC.0 10 (<10 % for 50 nm grains).
- RIT AdvCmos2014 (PDF p.4): Lpoly 0.40 um after etch -> 0.35 um after the 500 A re-oxidation (recipe 250).

## 2. Physics implemented (compress model only) and validated
`Oxidant.c` (poly/oxide boundary = reaction sink, <110> rate via an `Orr_poly` flag), `oxide_vel.c`
(interface velocity from the poly node), `elast.c` (Si and Poly both rigid and consumed:
`OXSEMI(M)`), `grid_upd.c` (oxide grid layers on poly).  All gated by `oxide_model == O_ELASTIC`, so
vertical/viscous decks are unchanged.

0.4 um poly line, dry O2 1000 C 51 min (patched build, completed run):

| quantity | result | check |
|---|---|---|
| oxide on Si | 454 A total (354 grown), 161 A Si consumed | 0.44 x 354 = 156 |
| oxide on poly top | 499 A, 222 A poly consumed | 222/0.44 = 505; <110> faster than <100> |
| sidewall oxide | 500 A per side | |
| poly width at mid-height | 0.400 -> 0.355 um | RIT 0.35 |
| gate-edge dielectric | 271 A at 0.02 um, 162 A at 0.05, 108 A at 0.1 um in | vertical model: 406/330/260/183 A (longer smile) |

Robustness matrix with the local-fix build (`native64_40c23ifv`): base, 120 min dry, wet 900 C 20 min,
0.25 um line, coarse mesh completed; 0.8 um wide line timed out at 1159 s of 3060 s.

## 3. Failure modes seen (the catalogue a new implementation must handle)
1. Poly sliver at the bottom gate corner whose three vertices are all on the oxide interface: default
   mesh-patch removal pins every multi-material point, so it can never be removed (stall at 47 s).
2. Convex consumed corners: interface points from two faces converge and must merge.
3. Concave corner beside the gate: sidewall oxide and substrate oxide grow into a narrow gas notch; its
   walls collide (embedding check), and with mismatched wall discretization the notch is a slit whose tip
   angle is not small.
4. Interior consumed points overrun by the interface inside a *split* time step (the first half commits).
5. Lagging oxide surface points invert surface triangles; obtuse-triangle rule skips boundary points.
6. Legacy `kill_node` retriangulation can loop allocating memory (OOM) on large fans.
7. Point indices change across `bd_connect()` (`waste()` renumbers): anything queued must be stored by
   coordinates.
8. Rejected candidate meshes leave dead points that crash `point_vel` unless collected (+UpdateSymbolic).

Repair-module specific (Boolean + rebuild):
9. Oxide encloses the poly (gate oxide + native oxide): SUPREM regions cannot have holes -> cut pieces.
10. SUPREM's skeleton triangulator (`triang`) clips ears at near-collinear boundary vertices -> zero-area
    slivers and T-junctions; a constrained Delaunay triangulation (CDT) is required.
11. Clipper offsetting drops collinear vertices, and offset rounding leaves ~1e-11 cm slivers that pinch
    loops; filter pieces by a small opening and add only "what changed" to the original polygons.
12. Deleting oxide triangles deletes surface points together with their gas nodes; new gas nodes need
    the ambient values or the oxide/gas boundary becomes an oxidant sink (growth stops).
13. Thinning vertices at a scale comparable to the native oxide (~1 nm) cuts through it and exposes poly.
14. Exposure must be derived from the assembled triangle boundary (`bd_connect_material`).
15. Last open bug: duplicated triangles between a rebuilt region and an existing one ("tnabor ... 3 vertex
    matches"), i.e. a region was not fully deleted before rebuilding.

## 4. What the production tools do (see research_oxidation_meshing.md)
- Primary: Lagrangian moving mesh with *local* operations each step - grid time step so no element
  collapses; growth per step capped (Sentaurus `dThickness` 1 nm; interface crosses at most one element);
  split growing-oxide edges perpendicular to the interface (`perp.add.dist` ~2 nm); remove short edges
  in the shrinking material (`Remove.Dist` 2-3 A, incl. interface points); flips.
- Fallback ("Plan B"): when the grid-limited step collapses, build a boundary representation, clean it
  to 0.5-1 nm, protect the original surface points, and regenerate the **whole** mesh by Delaunay from
  the B-rep plus the good bulk points.
- Even Synopsys: front collisions with points inside Si/poly/gas are "not completely solved"; gate
  corners in poly re-oxidation are mesh-sensitive (mechanics singular at corners).
- Alternatives without moving boundaries: diffuse interface (Hollauer, TU Wien 2007) and level sets
  (Victory, ViennaPS - GPL-3); Synopsys considers level-set regeneration non-conservative for dopants.

## 5. Recommended plan if resumed
1. Add the two missing standard local operations to the compress/viscous path:
   (a) growth cap per step (~1 nm, plus "no more than one element") in `oxgrow`/`back_pedal`;
   (b) short-edge removal for consumed-material interface points by validated edge collapse (positive
       orientation of every surviving triangle; reject otherwise), applied in `grid_sub`.
   Test on `test_decks/` before anything else.
2. Fallback: whole-mesh regeneration from a cleaned B-rep (Clipper2 is already linked) with a real
   constrained Delaunay mesher (CDT or Triangle/Gmsh) and region assignment from the B-rep, instead of
   deleting and rebuilding regions inside the database.
3. Field transfer: interpolate from a snapshot, restore each material's dopant inventory, give new gas
   nodes the ambient values, reset velocity history.
4. Consider the viscous model for poly (poly deforms, gate edge lifts) once 1-3 work; TSUPREM-4 supports
   poly oxidation there.
5. Validation: the six test decks; consumption = 0.44 x oxide on every face; area/dose conservation;
   AdvCmos/SubCmos decks bit-identical (they do not use compress).
