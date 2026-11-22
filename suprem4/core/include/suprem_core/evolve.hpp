// WP5: boundary motion, geometric events and repair.
//
// Geometry is represented per material as polygons on an exact integer grid
// (Clipper2, 1e-7 um).  Events are found by morphology with the event size
// delta (default 1 nm):
//   thin    part of a consumed-material region narrower than delta (opening;
//           protected growing materials such as oxide are exempt)
//           -> given to the neighbour sharing most of its boundary
//              (removed if that neighbour is the ambient)
//   crumb   connected component narrower than delta/8 everywhere (e.g. a
//           remnant of a removed piece) -> as thin
//   gap     ambient gap narrower than delta (opening of the ambient)
//           -> filled by the highest-priority adjacent material
//   cavity  ambient pocket enclosed by material (cannot be represented)
//           -> filled by the highest-priority adjacent material
// The repaired polygons are turned back into a mesh: planar straight-line
// graph of all polygon edges (T-junctions split, crossings resolved by the
// CDT), triangles labelled by point in polygon, regions = connected
// components per material, then the WP3 remesh with the previous mesh's size
// field and the WP4 field transfer.
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "suprem_core/mesh.hpp"
#include "suprem_core/remesh.hpp"

namespace suprem_core {

class SizeField;

struct Shape {
    int material = 0;
    std::vector<std::vector<Point>> loops;   // outer boundary and holes (um)
};

// Boundary code of a domain-boundary segment (EXPOSED, BACKSIDE, BC_OFFSET).
using BoundaryCode = std::function<int(const Point &, const Point &)>;

std::vector<Shape> mesh_shapes(const Mesh &m);   // one shape per region (needs a valid B-rep)
// Per material, the union of its triangles: robust to overlapping or
// crossing boundaries (a collision produced by the motion).
std::vector<Shape> mesh_shapes_union(const Mesh &m);
// Codes taken from the coded (non-exposed) boundary of m; EXPOSED elsewhere.
BoundaryCode boundary_codes_of(const Mesh &m);

struct RepairOptions {
    double delta = 1e-3;   // um
    std::vector<int> fill_priority = {SIO2, OXNI, SINI, POLY, SI, AL, PHRS, GAAS};
    // Growing materials: a film thinner than delta is real (native / early
    // oxide), never thinned away; only their crumbs smaller than delta^2 go.
    std::vector<int> protect = {SIO2};
};

struct RepairEvent {
    std::string kind;      // thin | crumb | gap | cavity
    int material = 0;      // material the piece came from (GAS for gaps/cavities)
    int to = 0;            // material it went to (GAS = removed)
    double area = 0;       // um^2
    Point at;              // a point of the piece
};

std::vector<RepairEvent> repair_shapes(std::vector<Shape> &shapes, const RepairOptions &o);

// Mesh of the shapes: metadata and species from `reference`, node values zero.
// size: target sizes (nullptr: from the coarse mesh itself).
Mesh mesh_from_shapes(const std::vector<Shape> &shapes, const Mesh &reference, const BoundaryCode &code,
                      const SizeField *size, const RemeshOptions &o, RemeshReport &rep,
                      std::vector<std::string> *notes = nullptr);

using Velocity = std::function<Point(const Point &, double t)>;   // um per time unit at time t

struct EvolveOptions {
    RepairOptions repair;
    RemeshOptions remesh;
    double max_step = 1e-3;          // um per step (also limited to delta / 4 so that gaps are seen before closing)
    double min_area_ratio = 0.5;     // no triangle loses more than this fraction... kept per step
    double rebuild_area_ratio = 0.3; // rebuild when a triangle falls below this fraction of its area at the last rebuild
    int max_steps = 100000;
    bool check_every_step = true;    // invariants + B-rep after every step (test mode)
};

struct EvolveReport {
    int steps = 0, rebuilds = 0, halvings = 0;
    std::vector<RepairEvent> events;
    std::vector<std::string> notes;
    std::vector<std::string> failures;   // invariant / B-rep violations (empty = ok)
};

// Move the mesh with the velocity field for `time`, detecting and repairing
// events and regenerating the mesh (with field transfer) as needed.
void evolve(Mesh &m, const Velocity &v, double time, const EvolveOptions &o, EvolveReport &rep);

}  // namespace suprem_core
