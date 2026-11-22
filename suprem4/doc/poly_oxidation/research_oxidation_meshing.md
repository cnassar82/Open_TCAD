# How production simulators keep the mesh valid during oxidation (research, 2026-11-22)

Question: how do robust TCAD tools handle the moving Si/SiO2, poly/SiO2 and oxide/gas boundaries
(corner collapse, surface collisions, thin layers) during oxidation, and what does that imply for
poly re-oxidation in SUPREM-IV.GS?

## Sources read
- Sentaurus Process User Guide, Y-2006.06, ch. 8 "Oxidation and silicidation" (2D/3D mesh handling,
  3DOxSetup parameter file) and ch. 1 (poly re-oxidation example).
- D. Guoy, A. Gencer et al. (Synopsys), "3-D Simulation of Silicon Oxidation: Challenges, Progress and
  Results", SISPAD 2013 (MovingMesh).
- C. Hollauer, "Modeling of Thermal Oxidation and Stress Effects", dissertation, TU Wien 2007.
- TSUPREM-4 2002.4 manual, "Grid Structure: Oxidation and Silicidation" (read earlier, local copy).

## Findings
1. **Primary method everywhere: Lagrangian moving mesh with local operations each step.**
   Sentaurus 2D (ch. 8): a grid time step is chosen so no element collapses; edges in the growing
   oxide that are nearly perpendicular to the interface are split above `perp.add.dist` (e.g. 2 nm);
   short edges in the shrinking region (Si/poly) are removed below `Remove.Dist` (2-3 A, "should be
   kept above 2e-8 cm"); growth per step is capped (`dThickness`, 1 nm default); "the interface [may
   not] traverse more than one element thickness at a time"; new nodes are interpolated after the step
   and the time step growth is capped (`IncreaseRatio`).  TSUPREM-4 describes the same three operations
   (remove nodes in consumed layers, add nodes in growing layers, flip/remove/relabel in congested
   regions).
2. **Fallback: global regeneration from a cleaned boundary representation.**  When local cleanup fails
   and the grid-limited step collapses, Sentaurus builds a boundary representation (B-rep), cleans it to
   `BrepAccuracy` (0.5-1 nm), protects the original surface points (decimating only short edges), and
   calls a Delaunay mesher with the B-rep and the "good" bulk points ("Plan B"); a level set is used
   for native-layer construction.  This is the same idea as the Boolean-repair module tried here, but
   used only as a fallback, with whole-mesh regeneration and region assignment from the B-rep.
3. **Even Synopsys calls it unsolved in hard cases.**  "handling the frequent collision of the oxidation
   front with points inside the silicon, polysilicon, or gas regions can cause intractable problems for
   the local mesh operations ... These problems have not been completely solved ... various global
   remeshing strategies have been implemented to try to recover."  The 3D MovingMesh needed "many
   years of research".
4. **Gate corners in poly re-oxidation are specifically hard.**  "the quasistatic mechanics equations
   ... are prone to high levels of shape dependence. This is most frequently seen at gate corners during
   polysilicon reoxidation ... [and] a high sensitivity to the mesh modification algorithms at these
   corners."  At sharp corners the mechanics are singular.
5. **Alternatives avoid moving boundaries altogether.**  Diffuse-interface (Hollauer: normalized
   silicon concentration on a static mesh) and level-set methods (Victory, ViennaPS) have no topology
   problem; Synopsys rejected level-set regeneration for oxidation as "non-conservative" for dopants.
6. Sentaurus' own poly re-ox example is a thin stress-relief oxide (900 C, 10 min, 0.5 atm), not 500 A.

## What it means here
- The SUPREM-IV local operations miss two of the standard ones: interface short-edge removal (the
  mesh-patch removal pins every multi-material point) and a per-step growth cap.  Those, not exotic
  topology handling, are the first things a robust implementation needs.
- The repair-and-remesh module is the right fallback idea, but it should regenerate the whole mesh from
  a B-rep (as Sentaurus Plan B does) rather than delete and rebuild regions inside the database, which
  is where its bookkeeping failures came from.
- Expect mesh sensitivity of the gate-edge stress/shape even when it runs.
