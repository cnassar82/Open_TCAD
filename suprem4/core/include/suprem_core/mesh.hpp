// Semantic 2D mesh used by the geometry core.
//
// Mirrors SUPREM-IV's model: points carry one node per material that touches
// them (multi-valued fields at interfaces), triangles belong to regions, and
// each region has one material.  Coordinates are micrometres.
#pragma once
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "suprem_core/str_file.hpp"

namespace suprem_core {

// SUPREM-IV material numbers (include/material.h).
enum Material { GAS = 0, SIO2 = 1, SINI = 2, SI = 3, POLY = 4, OXNI = 5, AL = 6, PHRS = 7, GAAS = 8 };

// Neighbour codes for triangle faces without a neighbouring triangle.
constexpr int BC_OFFSET = -1024;           // reflecting / internal boundary
constexpr int EXPOSED = BC_OFFSET + 2;     // faces the ambient
constexpr int BACKSIDE = BC_OFFSET + 1;    // wafer backside

struct Point { double x = 0, y = 0; };

struct Triangle {
    int region = 0;              // 0-based
    int v[3] = {0, 0, 0};        // 0-based point indices (counter-clockwise)
    int nb[3] = {BC_OFFSET, BC_OFFSET, BC_OFFSET};  // face j is opposite v[j]: neighbour (0-based) or code
    int father = -1, offspring = -1;
};

struct Node {
    int point = 0;
    int material = 0;
    std::vector<double> value;   // per species slot
};

struct Mesh {
    std::vector<Point> point;
    std::vector<int> region_material;
    std::vector<Triangle> triangle;
    std::vector<int> species;          // impurity number per solution slot
    std::vector<Node> node;
    // Round-trip metadata
    std::string version;
    int dim = 2, nvrt = 3, nedg = 3;
    bool has_model = false;
    int orientation = 0;
    double last_temperature = 0;
    std::vector<int> implanted;
    std::vector<std::string> trailer;

    int material_of(int triangle_index) const { return region_material[triangle[triangle_index].region]; }
    // Node of a material at a point, or -1.
    int node_at(int point, int material) const;
    void rebuild_index();

  private:
    std::map<std::pair<int, int>, int> node_index_;
};

Mesh from_str(const StrFile &f);
StrFile to_str(const Mesh &m);

// Signed doubled area of triangle t (positive = counter-clockwise).
double orient2(const Mesh &m, int t);

// Structural invariants; returns human-readable violations (empty = valid).
// Checks: vertex indices, positive orientation, neighbour symmetry and shared
// faces, a node of the region's material at every vertex, node values sized
// to the species list, and no duplicate triangles.
std::vector<std::string> check_invariants(const Mesh &m, size_t max_reports = 20);

}  // namespace suprem_core
