// SUPREM-IV validity contract (M0 WP6a, doc/plan_m0_geometry_core.md 4b).
//
// Exact ports of the checks SUPREM applies to a structure it reloads with
// "structure infile=" (ig2_read + create_db) and to every trial motion:
//   orientation  no triangle flat in SUPREM's stored float cm coordinates
//                (orient2d_float_sign, src/geom/misc.c; "triangle is flat"
//                panic in geometry.c); negative ones are swapped on load by
//                sclk_tri (clktri.c), so the remaining rules use the swapped
//                order;
//   topology     each point-edge shared by at most two triangles (nxtel);
//   embedding    no two triangles overlap and no outer-boundary edge touches
//                another edge or vertex (mesh_embedding_check,
//                src/include/mesh_embedding.h; boundary = no neighbour);
//   loops        each region's border is one simple closed loop whose turn
//                sum is +-2 pi (skel_reg + ck_clock, src/refine/skel.c and
//                triang.c; a self-touching vertex or a crossing is a panic);
// Stricter than the reload (the core must still satisfy them, SUPREM loads
// files that violate them):
//   ordering     every triangle counter-clockwise as written;
//   components   each region one component without holes (skel_reg keeps only
//                the first loop, so later grid operations see a wrong region);
//   gas_added    a gas node at every exposed point (written, not added on load);
//   nodes        mtest1 (src/dbase/check.c): a gas node at every backside
//                point and at no unexposed point (flags from the face codes,
//                mk_ele; create_db adds missing ones at exposed points), one
//                node per material and point, every material node in a
//                triangle of its material, valid material codes.
// Coordinates are converted exactly as ig2_read does: (float)(x_um * 1e-4).
#pragma once
#include <string>
#include <vector>

#include "suprem_core/mesh.hpp"

namespace suprem_core {

struct SupremCheck {
    // counts per rule (0 = rule holds)
    int orientation = 0, topology = 0, embedding = 0, loops = 0, nodes = 0;   // reload fails
    int ordering = 0, components = 0, gas_added = 0;                         // strict only
    // the rules a reload into SUPREM enforces (panic or crash otherwise)
    bool reloadable() const { return !orientation && !topology && !embedding && !loops && !nodes; }
    bool ok() const { return reloadable() && !ordering && !components && !gas_added; }
    std::vector<std::string> messages;   // first reports, human readable
};

SupremCheck check_suprem(const Mesh &m, size_t max_reports = 20);

// SUPREM's exact orientation sign for float coordinates (cm).
int orient2d_float_sign(const float *a, const float *b, const float *c);
// The float cm coordinate SUPREM holds for a coordinate written in um.
inline float suprem_float(double um) { return (float)(um * 1e-4); }

}  // namespace suprem_core
