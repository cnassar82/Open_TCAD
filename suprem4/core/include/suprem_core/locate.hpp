// Point location in a triangle mesh and the size field derived from it.
#pragma once
#include <array>
#include <vector>

#include "suprem_core/mesh.hpp"

namespace suprem_core {

// Bucket-grid point locator over the triangles of a mesh.
class TriangleLocator {
  public:
    // material >= 0: only triangles of that material are searched
    explicit TriangleLocator(const Mesh &m, int material = -1);
    // Triangle containing p and its barycentric coordinates.  When p lies
    // outside the mesh the nearest triangle is returned with the coordinates
    // clipped to it.  Returns -1 only for an empty mesh.
    int locate(const Point &p, double bary[3]) const;

  private:
    const Mesh &m_;
    double x0_ = 0, y0_ = 0, cell_ = 1;
    int nx_ = 1, ny_ = 1;
    std::vector<std::vector<int>> bucket_;
    double score(int t, const Point &p, double bary[3]) const;
};

// Isotropic mesh-size field of an existing mesh.  Per triangle the size is
// the area-equivalent edge length ("area"), the smallest altitude ("min")
// or the longest edge (LONGEST);
// vertex sizes are area-weighted means (degenerate triangles carry no
// weight) and are interpolated linearly, then scaled and clamped.
// Optionally tightened where the dopant profiles are steep.
class SizeField {
  public:
    enum Mode { AREA, FINEST, LONGEST };   // LONGEST: area-weighted longest edge
    // material >= 0: built from the triangles of that material only, and
    // evaluated in the nearest of them (the size a material had, also where
    // it has since moved)
    SizeField(const Mesh &m, Mode mode, double scale, double hmin, double hmax, int material = -1);
    double operator()(const Point &p) const;
    // Solution adaptation: where a conserved dopant varies by more than
    // `decades` (log10) across an old triangle, the size there shrinks in
    // proportion, but never below the old mesh's finest local size (smallest
    // altitude): the new mesh resolves at least what the old one did, and
    // is never refined beyond it.
    void adapt_to_solution(double decades);
    double hmin() const { return hmin_; }
    double hmax() const { return hmax_; }

  private:
    const Mesh &m_;
    TriangleLocator loc_;
    std::vector<double> size_;     // per point
    std::vector<double> finest_;   // per point: area-weighted smallest altitude
    Mode mode_;
    double scale_, hmin_, hmax_;
    int material_ = -1;
};

}  // namespace suprem_core
