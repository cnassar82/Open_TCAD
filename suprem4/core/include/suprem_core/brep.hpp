// Boundary representation (B-rep) extracted from a Mesh.
//
// Every region's boundary is a set of closed loops (region on the left:
// outer loops counter-clockwise, holes clockwise).  Each loop is split into
// chains at junction points; a chain is a maximal run of boundary edges with
// the same left region and the same right side (another region, or one of the
// boundary codes EXPOSED / BACKSIDE / BC_OFFSET).  Interface chains therefore
// appear twice, once from each side, with reversed vertex order: that shared
// representation is what guarantees conforming meshes after regeneration.
#pragma once
#include <string>
#include <vector>

#include "suprem_core/mesh.hpp"

namespace suprem_core {

struct Loop {
    int region = 0;
    std::vector<int> point;   // closed: last connects to first
    double area = 0;          // signed (um^2), positive for outer loops
};

struct Chain {
    int left_region = 0;
    int right = 0;            // >= 0: region index; < 0: boundary code of the faces
    std::vector<int> point;   // open polyline, endpoints are junctions
};

struct Brep {
    std::vector<Loop> loop;
    std::vector<Chain> chain;
};

// Extract the B-rep; problems found during extraction are appended.
Brep extract_brep(const Mesh &m, std::vector<std::string> &problems);

// Validate: loops closed and simple (no two non-adjacent segments of any
// boundary intersect), per-region loop area equals the triangle area, every
// interface chain has an exact reversed twin.
std::vector<std::string> check_brep(const Mesh &m, const Brep &b, size_t max_reports = 20);

}  // namespace suprem_core
