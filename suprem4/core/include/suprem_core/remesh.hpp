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
