# How other TCAD packages address the current repair issues

Researched 2026-11-16. Scope: the construction-identity and conservative-connectivity
repair reports, including the RIT oxidation/activation failures. This is a research
and implementation-direction document. No external package was benchmarked on
our RIT deck, and no simulator or launcher changes were made for this research.

## Recommendation

Repair the material representation and its conversion into a numerical mesh
before adding more recovery logic to the time integrator. Keep three separate
contracts: geometric accuracy and material ordering; conservative field transfer;
and nonlinear/time-integration acceptance.

The evidence supports separating evolving process geometry from the PDE mesh.
It does **not** establish that every successful package uses an exact CAD history,
that level sets automatically preserve thin layers, or that changing coordinates
to double alone will solve this problem. Construction provenance is one useful
implementation choice for our analytic operations, not the only valid architecture.

## Findings from other packages

### FLOOXS: compute geometry, then generate its mesh

The University of Florida process tutorial describes level-set structure updates
followed by generation of a new mesh. Deposition exposes separate level-set
spacing and geometric tolerance controls; diffusion then operates on the finite
element grid. It also declares whether a solution variable may be negative.
This is a close process-simulation reference for SUPREM's replacement architecture.
The tutorial does not document the conservative transfer implementation or prove
positivity of its transient algorithm. A negative-variable declaration is not
that proof. [FLOOXS process tutorial](https://www.flooxs.ece.ufl.edu/index.php/Process_Tutorial)

**Application here — recommendation:** stop making the current triangle boundary
our sole record of what deposition constructed. Give geometric operations their
own accuracy controls and produce a validated solver mesh from the resulting
material domains. Keep the original recipe unchanged.

### Victory Process: implicit layers, hierarchical geometry and hybrid oxidation

Silvaco documents process-mode geometry as material layers bounded by implicit
surfaces on hierarchical Cartesian grids. Its oxidation modes include empirical,
physical and hybrid treatments, with switching related to oxide thickness;
empirical treatment can resolve planar layers below the mesh size. The physical
mode couples oxidant transport, reaction, flow, deformation and stress.
These are documented capabilities, not disclosure of every internal algorithm.
[Victory Process brief, especially page 6](https://silvaco.com/wp-content/uploads/product/pdf/VictoryProcess_brief.pdf)

**Application here — recommendation:** a real very thin layer need not always
be represented by vanishingly small moving triangles. A sub-grid or thin-layer
model is a possible architectural option, but would require its own transport,
storage, coupling and convergence validation. Switching oxidation models merely
to complete RIT would not validate the original process physics.

### Published Victory Process research: refine thin layers before etching

Lenz and coauthors identify precisely the limitation relevant to our proposal:
level-set methods can represent flat thin films on coarse grids, yet Boolean
etching exposes inadequate resolution and produces artifacts. Their work develops
local refinement based on affected layer thickness and demonstrates it on an LED
process. [2023 paper](https://www.iue.tuwien.ac.at/publist/pdf/20.500.12708_153610.pdf)

The accompanying thesis describes detecting interfaces affected by the cutter,
finding their nearest outward interface separation, and setting a target spacing
`h_target = d_closest / N_min`. Refinement precedes the Boolean operation and is
localized to affected regions. Its benchmark uses Victory Process. `N_min` is an
accuracy parameter, not a universal constant to copy. Wrapping layers prevents
another class of voids but does not cure inadequate resolution.
[Lenz thesis, Chapter 7](https://www.iue.tuwien.ac.at/phd/lenz/Thin-Material-Layer-Refinement-Etching-Simulations.html)

**Application here — recommendation:** estimate interface separation and required
resolution before cutting or deforming a physical thin layer. Do not wait for a
crossing or collapsed triangle. First distinguish a real layer from the artificial
film produced by our deposit/etch inconsistency; endlessly refining the latter
would preserve the wrong construction at increasing cost.

### ViennaPS/ViennaLS: material ordering and sub-grid boundary positions

ViennaPS uses a wrapped multi-material level-set representation: successive
functions represent cumulative unions of material regions. Automatic wrapping
maintains that relationship when inserting a layer.
[ViennaPS surface/material representation](https://viennatools.github.io/ViennaPS/domain/surface.html)

ViennaLS's oxidation documentation describes a Cartesian oxide solve grid and
interface distances reconstructed from level-set values. These distances enter
the boundary coefficients instead of forcing every boundary onto a grid node.
It also documents floating-point classification tolerances, minimum boundary
distance guards and persistence of interior fields across advection. These
mechanisms have assumptions and error implications; their constants should not
be copied into SUPREM. Persisting or interpolating a field is not itself proof
that dopant inventory is conserved.
[ViennaLS oxidation design](https://github.com/ViennaTools/ViennaLS/blob/master/docs/OxidationSolver.md)

**Application here — recommendation:** store material ordering and shared
interfaces explicitly; use physical boundary positions in interface quadrature.
Do not generate the same contact independently for two material owners. If
introducing sub-grid geometry, account for its actual volumes and fluxes.
The current upstream documentation is a design reference, not a validated RIT
replacement or a promise about every released version.

### Sentaurus: distinguish interpolation accuracy from conservation

The historical Synopsys-authored Sentaurus Process Y-2006.06 manual describes
MGOALS remeshing and two transfer choices. Conservative interpolation integrates
control-volume dose and redistributes it by overlap with new control volumes;
linear interpolation is faster and has different profile-preservation behavior.
The manual explicitly discusses smoothing of extrema by conservative transfer.
This is evidence of an established design tradeoff, **not** a claim about current
Sentaurus defaults or implementation.
[Historical manual, Mesh Generation chapter](https://dunham.ece.uw.edu/protected/PDFManual/data/sprocess_ug.pdf)

**Application here — recommendation:** retain our inventory test for every split,
collapse and flip. For larger remeshes, evaluate overlap-based conservative
transfer with a bounded reconstruction. Preserve material ownership and measure
peak/profile error as well as total dose. A local conservative projection is a
useful primitive but is not automatically a globally accurate remapping scheme.

### DEVSIM: integration and precision are independent controls

DEVSIM documents BDF1, TRBDF and BDF2, plus independently selectable extended
precision for model evaluation, equation assembly and linear/Newton solves.
It is a device simulator; this does not provide an etch or moving-boundary
replacement. [DEVSIM numerical methods](https://devsim.net/solver.html)

**Application here — recommendation:** diagnose geometric rounding, cancellation
in volume/flux assembly and linear-system conditioning separately. Higher
precision in the linear solve cannot restore an interface already rounded away.
A higher precision geometry path likewise cannot make a non-positive time
integration scheme unconditionally positive.

### Genius: error control, history and positivity require explicit policies

I inspected the local upstream-derived implementation, rather than inferring
behavior from a method name:

- `src/solver/ddm_common/ddm_solver.cc:1428` distinguishes nonlinear divergence
  recovery from LTE rejection; it restores/re-predicts state and controls the
  next step separately.
- `src/solver/ddm1/ddm1.cc:496` checks whether BDF2 history would violate a density
  condition and requests lower order when that condition fails.
- `src/solver/ddm1/ddm1.cc:457` also replaces negative **carrier predictors** with
  a fraction of previous values. Do not transplant that adjustment as a repair
  for conserved process dopants.

[Upstream transient controller](https://github.com/cogenda/Genius-TCAD-Open/blob/master/src/solver/ddm_common/ddm_solver.cc),
[upstream density/history handling](https://github.com/cogenda/Genius-TCAD-Open/blob/master/src/solver/ddm1/ddm1.cc).
The conclusions above come from the local files; upstream web fetches were
incomplete. Local source SHA-256 values:

```
ddm_solver.cc  da38e5cb3b3802d05fe488334643f2e711261ed8470a50a147bbc7be7ff55863
ddm1.cc        ae5dfa54c95b1d6b749b9e3c4d443bb2d4534602556e8413221a26dc3deb1f27
```

**Application here — recommendation:** distinguish a bad mesh, negative stage,
Newton failure and LTE failure. After remeshing, explicitly rebuild or restart
history on the new supports. Use a documented lower-order startup or positivity
strategy with an error estimate and rollback, rather than a longer retry tree.

## Mapping these findings to our repair report

| Open issue | Implementation direction | Required evidence |
| --- | --- | --- |
| Offset quantizes the source differently from the retained interface | One authoritative material-domain construction and a declared conversion error | Translated planar deposit/etch tests; vary lateral extent independently of local feature size |
| Mixed float/double thickness and binary32 boundary storage | Unified units/types at the geometry API, with audited conversions and shared construction identity | Show where rounding occurs; certify material ordering after export/restart |
| Very thin real layers become unmeshable | Pre-operation interface-aware refinement; evaluate sub-grid treatment only as a separate model | Thickness and surface-position convergence; preserved layer inventory and flux |
| Neighboring boundaries follow incompatible sampled motion | A common deformation/geometry description, with a globally valid resulting partition | Prescribed-motion tests and shared-interface identity, not only positive triangles |
| Flips/remeshing alter chemical inventory | Conservative transfer coupled to every connectivity operation | Per-material dose, positivity, bounds, affine-field and profile-error checks |
| TR/BDF stages fail after a mesh change | Correct mass/flux measures and controlled history restart, followed by a positivity-aware integrator | Temporal convergence, conservation, restart agreement and explicit failure reasons |

These rows are our proposed engineering response. No external source proves that
our existing implementation satisfies them.

## Changes to the earlier repair recommendations

**Keep exact cancellation as the acceptance criterion for an analytic planar
construction mode.** It exposes our inconsistent deposition/etch path without
involving oxidation. But do not assert that every discretized level-set package
must return bitwise-zero material after paired operations. A numerical geometry
mode needs a stated positional/volume error bound and refinement convergence.
That distinction should be made before comparing packages.

**Do not equate every tolerance with a band-aid.** Other implementations use
classification guards, distance limits and resolution policies. A defensible
policy states its units, scale, accuracy assumptions, affected geometry and
conservation consequences. An unexplained epsilon, enlarged etch depth or
small-region deletion that merely makes RIT finish is still unacceptable.

**Our observed artificial planar film does not explain every RIT failure.**
Real sloped surfaces can retain material after directional etching. The current
reports show separate transfer, motion and transient-stage failures; trace each
against its own invariant rather than assigning all of them to the planar test.

**Steady oxidant transport is not an unimplemented cure.** Our default O2/H2O
fields already carry STEADY flags in `diffuse/diff_init.c`. The reported negative
phosphorus stages involve a transported dopant. Making that dopant algebraic or
clipping its concentration would change the problem.

## Recommended implementation order

1. **Specify the geometry contract.** Decide which operations are analytic and
   which are approximate, select a physical geometric accuracy target, and
   separate it from PDE mesh spacing. Preserve source/interface identity across
   paired constructions. Account for holes and multiple material owners.
2. **Fix the planar construction path first.** Use the existing three translated
   tests and add lateral-domain scaling, nonzero starting heights, unequal etch
   depths and a true retained-layer control. Test the geometry before meshing,
   after conversion to the legacy database, and after a structure-file restart.
   A float-to-double local edit is not sufficient evidence of completion.
3. **Add a thin-layer preflight.** Determine actual interface separations in the
   operation's affected region. Refine physical layers before Boolean/motion
   operations, or refuse an unsupported embedding with a diagnostic. Do not
   refine a known construction artifact into an enormous mesh.
4. **Complete geometric conservation tests.** Under prescribed motion, compare
   actual nodal volume changes, integrated swept flux and chemical inventory.
   Audit the legacy obtuse-element mass treatment against physical area before
   claiming that a successful transient solve has conserved physical dose.
5. **Then repair time integration.** Rebuild histories after topology changes;
   combine LTE control with a justified positivity strategy. Verify a static
   stiff-diffusion reference and a moving-domain reference before another full
   process run. A lower-order method must converge with time-step refinement.
6. **Run unchanged fresh NMOS and PMOS.** Require every process checkpoint,
   independent partition/dose checks, device polarity checks, mesh/mask
   perturbations and temporal convergence before launcher promotion.

This order targets the construction and discrete conservation problems. It
avoids making completion depend on hand-tuned retry counts or deleting the very
material whose evolution the simulation is meant to predict.

## Supporting numerical references

Exact predicates alone do not make rounded geometric constructions exact;
CGAL documents how inexact intersection constructions can still create
self-intersections. This supports auditing the construction-to-mesh conversion,
not only adding more orientation tests.
[CGAL corefinement documentation](https://doc.cgal.org/6.1.2/Polygon_mesh_processing/index.html)

Published TR-BDF2 analysis establishes conditional positivity/monotonicity and
studies hybridization with implicit Euler. Its assumptions still have to be
checked for our spatial operator and moving mass terms; it does not guarantee
that a method-name change will repair SUPREM.
[Bonaventura and Della Rocca](https://arxiv.org/abs/1510.04303)
