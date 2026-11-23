// WP4: field transfer from an old mesh to a new mesh of the same geometry.
//
//  - Per material: a new node (point, material) takes the P1 interpolant of
//    the old nodes of that material, evaluated in the old triangles of that
//    material (interface values stay per material: segregation jumps are
//    kept).  Points that exist in the old mesh with identical coordinates
//    copy the old node exactly.
//  - Gas nodes (ambient values) are copied from the nearest old gas node.
//  - Conserved species (total dopant concentrations): the inventory per
//    material, measured as the exact integral of the P1 field (lumped
//    area/3 volumes), is restored to the old one by the minimum change
//    weighted by each node's interpolation uncertainty (the spread of the
//    old values it was interpolated from; exact copies do not change), so the
//    correction lands where the old profile was under-resolved instead of
//    rescaling the whole material; values stay non-negative.  The active
//    concentration of a dopant changes in proportion to its total.
//    Why P1 and not SUPREM's mixed-Voronoi box volumes: on the sliver /
//    obtuse triangles of SUPREM meshes the box sum is a poor quadrature
//    (e.g. 25 % below the field integral in a pad oxide), so restoring it
//    would remove real dopant.  On the new Delaunay meshes both agree.  Everything else
//    (point defects, oxidant, velocities, stresses, potential) is
//    interpolated only.
#pragma once
#include <string>
#include <vector>

#include "suprem_core/mesh.hpp"

namespace suprem_core {

// SUPREM control volume of every node (index = node), mixed Voronoi areas
// as in dbase/geometry.c (or P1 lumping, area/3).
std::vector<double> node_volumes(const Mesh &m, bool p1_lumping = false);

// Inventory of species slot s per material: sum of value x volume.
std::vector<double> inventory(const Mesh &m, const std::vector<double> &vol, int slot, int nmat = 9);

bool conserved_species(int impurity);
// Total species whose factor an active species takes (-1 if not active).
int total_of_active(int impurity);

struct TransferReport {
    struct Row {
        int material, impurity;
        double old_inv, interpolated_inv, new_inv;
    };
    std::vector<Row> rows;            // conserved species and their actives, per material
    int exact = 0, interpolated = 0, gas = 0, outside = 0;
    double max_rel_interpolated = 0;  // inventory error before restoration
    double max_rel_final = 0;         // after restoration
    double max_rel_box = 0;           // change of the box-volume (SUPREM) inventory, after restoration
    std::vector<std::string> notes;
};

enum class Restore {
    WEIGHTED,   // minimum change weighted by the local interpolation spread (default)
    SCALE       // uniform factor per material
};
void transfer_fields(const Mesh &old, Mesh &neu, TransferReport &rep, bool restore = true, bool p1_lumping = true,
                     Restore mode = Restore::WEIGHTED);

// WP6b: fields after a Lagrangian update (update_mesh).  Nodes at kept points
// copy their old values exactly; inserted nodes are interpolated as in
// transfer_fields.  The conserved inventories (P1, per material) are then
// restored by one factor per material and species applied only to the nodes
// at changed points (inserted points and vertices of new triangles), so all
// other nodes keep their old values bit for bit.
// box: restore the box-volume (mixed Voronoi, SUPREM's own) inventory
// instead of the P1 lumped one; max_rel_final is then in that measure.
void update_fields(const Mesh &old, Mesh &neu, const std::vector<char> &changed_point, TransferReport &rep,
                   bool box = false);

// Round-trip profile error: interpolate the new mesh back onto the old
// nodes; per conserved species, the volume-weighted L1 error relative to
// the old inventory, max over materials (an interpolation-error measure).
double roundtrip_error(const Mesh &old, const Mesh &neu, int slot);

}  // namespace suprem_core
