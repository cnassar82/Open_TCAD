// WP3: mesh generation from the B-rep.
//
// Pipeline (remesh):
//  1. B-rep -> canonical chains (each interface stored once, shared by both
//     regions, so the regenerated mesh is conforming by construction).
//  2. Chain discretisation: chain vertices are dropped only where the chord
//     stays within coarsen_tol of them and is not longer than the target
//     size (coarsen_tol = 0 keeps every vertex: exact geometry); segments
//     longer than the target size are subdivided on the straight segment.
//  3. Thin films (local thickness < film_threshold): boundary vertices are
//     projected across the film onto the opposite side; the matched pairs
//     carry film_layers - 1 interior layer points, giving aligned,
//     near-rectangular film elements (good for box-method discretisation).
//  4. Constrained Delaunay triangulation (CDT) of all chains + film points,
//     triangles labelled by flood fill from the region loops.
//  5. Mmg2d remesh of the non-film part with all chain edges and film
//     triangles required; isotropic size from the old mesh's size field.
//  6. Delaunay edge flips inside regions, then local encroachment repair
//     (Ruppert): boundary/interface edges whose opposite angle exceeds
//     encroach_angle are split at their midpoint (in the mesh and the chain)
//     and flips restore the Delaunay property; repeated until none is left.
// Node values are allocated (zero) but not transferred: that is WP4.
#pragma once
#include <map>
#include <string>
#include <vector>

#include "suprem_core/mesh.hpp"

namespace suprem_core {

struct RemeshOptions {
    std::string size_mode = "area";   // "area" or "min" (see SizeField)
    double scale = 1.0;               // multiplies the old mesh size
    double hmin = 2e-4, hmax = 0.5;   // um
    double hgrad = 1.3;
    double solution_decades = 0.25;   // max dopant variation (log10) per old element before refining (0: off)
    double coarsen_tol = 1e-5;        // um; chain vertices may be dropped within this (0: keep all)
    double film_threshold = 0.04;     // um
    int film_layers = 2;              // layers across a film ...
    double film_min_layer = 0.002;    // ... each at least this thick (um)
    double film_max_angle = 100;      // deg: film triangles above this are left to Mmg
    double encroach_angle = 91;       // deg
    int encroach_rounds = 50;
    double min_split_length = 1e-5;   // um: Ruppert termination guard; also ~100 float ulps (SUPREM stores float cm)
    int mmg_passes = 2;
    bool mmg = true;
    bool float_cm_grid = false;       // create points on SUPREM's coordinate grid (float cm)
    bool flips = true;
    int verbose = 0;
};

struct Quality {
    int points = 0, triangles = 0;
    double min_angle = 180, max_angle = 0;
    int below_20 = 0, below_10 = 0, above_120 = 0;
    int negative_couplings = 0;   // box-method edge couplings (cot sums per region) < 0
    int obtuse_boundary = 0;      // region-boundary edges whose opposite angle > 90.5 deg
    int significant_negative = 0; // couplings below -0.01 (e.g. a boundary edge opposite > 91.1 deg)
    double worst_coupling = 0;
};

struct RemeshReport {
    std::vector<std::string> notes;
    int coarsened = 0, chain_points = 0, film_pairs = 0, film_points = 0, film_triangles = 0;
    int encroach_splits = 0, flips = 0, mmg_failures = 0;
};

class SizeField;

// WP6b: Lagrangian update (remesh mode "update").  The points of m and their
// node values are kept; the mesh is changed only where the spacing no longer
// matches the reference mesh (the grid at the start of the oxidation).  The
// reference supplies three size fields, so its anisotropic elements are
// their own target (an update of the reference itself changes little):
// finest (smallest altitude), area (area-equivalent edge) and longest edge.
//  - boundary/interface: a chain vertex (not a junction) is removed when one
//    of its segments is shorter than short_edge x finest and the chord
//    replacing it stays within tol of it (chains that would then intersect
//    or collapse are restored); segments longer than long_edge x longest are
//    subdivided (exact geometry);
//  - interior: a point closer than near x finest to a boundary/interface
//    chain or to an already kept interior point is removed (the encroaching
//    interface of a consumed material); a triangle whose area exceeds
//    grow x the reference area there gets a point (midpoint of its longest
//    edge, or its centroid when that edge is a boundary);
//  - constrained Delaunay triangulation of the kept and inserted points with
//    the chains as constraints, then the WP3 encroachment repair, limited to
//    encroachment the update introduced and to splits not finer than the
//    reference's finest spacing.
// The B-rep is kept exactly apart from the tol-bounded vertex removals.
struct UpdateOptions {
    double near = 1.0 / 3;         // x finest: interior points this close to a chain / each other are removed
    double short_edge = 1.0 / 3;   // x finest: chain vertices on shorter segments may be removed ...
    double tol = 1e-4;             // um: ... when the chord stays this close (0: never)
    double long_edge = 1.5;        // x longest: chain segments longer than this are subdivided
    double grow = 2.5;             // x reference area: larger triangles get a point
    int rounds = 8;                // interior insertion rounds
    double scale = 1.0, hmin = 2e-4, hmax = 0.5;   // size fields (um)
    double solution_decades = 0;   // area field: dopant adaptation as in WP3 (0: off)
    double encroach_angle = 91;    // deg: only encroachment the update introduced is repaired ...
    int encroach_rounds = 50;
    double encroach_budget = 0.05; // ... with at most this fraction x points splits (min 100)
    double min_split_length = 1e-5;   // um; encroachment splits also stop at hmin (features below it are events)
    bool float_cm_grid = false;    // inserted points on SUPREM's float cm grid
    bool timing = false;           // stage times as notes
    std::map<int, double> min_size; // material -> lower bound of its target sizes (um), e.g. SiO2: SUPREM's grid.oxide
    bool keep_triangles = true;    // old triangles stay unless something changed there ...
    double keep_coupling = -0.01;  // ... or one of their edges has a coupling (cot sum / 2) below this
    // M1 interface boundary layer: on the consumed side (silicon, poly) of every
    // interface with oxide, rows of points at normal spacing interface_h up to
    // depth interface_depth, then graded by interface_growth per row until the
    // spacing reaches the material's own size there.  Along the interface the
    // rows keep the interface's vertex spacing (anisotropic elements, so the
    // point count grows with the number of rows, not with 1/h^2).  Lagrangian:
    // rows the interface has not reached are kept with their values; rows the
    // interface eats (closer than near x interface_h) are removed; new rows
    // are added at the deep end, where the profile is smooth.  0: off.
    double interface_h = 0;        // um
    double interface_depth = 0.05; // um
    double interface_growth = 1.5;
    // M1.2 look-ahead: the step SUPREM wants next (s), with node velocities
    // from the mesh's XVEL/YVEL solutions.  Interior points and interface
    // layer rows that a boundary of their material would approach within
    // twice that step's motion are removed now, so SUPREM's half-area step
    // limit does not cut the step (and again at every step while the gap
    // halves).  0: off.
    double step_dt = 0;
    // Flow quality: oxide/nitride triangles with an angle below this (deg)
    // are refined by longest-edge bisection (graded to the boundary), since
    // the flow's Newton solve stalls on needles.  0: off.
    double flow_min_angle = 15;
};

struct UpdateReport {
    int interface_points = 0, interface_inserted = 0;
    int kept = 0, removed_interior = 0, removed_boundary = 0, inserted_boundary = 0, inserted_interior = 0;
    int encroach_splits = 0, flips = 0, new_triangles = 0;
    std::vector<int> origin;       // per new point: old point index, -1 for inserted points
    std::vector<char> changed;     // per new point: inserted, or a vertex of a triangle the old mesh did not have
    std::vector<std::string> notes;
};

// Topology/geometry only: node values of the result are zero; use
// update_fields (transfer.hpp) with rep.changed.  Throws if m's B-rep is invalid.
Mesh update_mesh(const Mesh &m, const Mesh &reference, const UpdateOptions &o, UpdateReport &rep);

Quality measure_quality(const Mesh &m);
std::string format_quality(const Quality &q);

class SizeField;
// size: target size field (default: derived from m itself, with solution
// adaptation).  An external field is used when m is a coarse mesh built
// from repaired geometry (WP5) and the sizes come from the previous mesh.
Mesh remesh(const Mesh &m, const RemeshOptions &o, RemeshReport &rep, const SizeField *size = nullptr);

// Geometry preservation of a remesh: two-sided Hausdorff distance between
// the old and new boundaries <= tol, region areas equal within tol x
// perimeter (exact, up to rounding, for tol = 0).  Returns violations.
std::vector<std::string> compare_geometry(const Mesh &before, const Mesh &after, double tol = 0,
                                          size_t max_reports = 20);

// Gas nodes: a GAS node at exactly the points on EXPOSED or BACKSIDE faces
// (SUPREM mtest1: is_surf(p) || is_back(p)).
std::vector<std::string> check_gas_nodes(const Mesh &m, size_t max_reports = 20);

}  // namespace suprem_core
