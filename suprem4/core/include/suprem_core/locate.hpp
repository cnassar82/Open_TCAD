// Point location in a triangle mesh and the size field derived from it.
#pragma once
#include <array>
#include <vector>

#include "suprem_core/mesh.hpp"

namespace suprem_core {

// Bucket-grid point locator over the triangles of a mesh.
class TriangleLocator {
  public:
    explicit TriangleLocator(const Mesh &m);
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
// the area-equivalent edge length ("area") or the smallest altitude ("min");
// vertex sizes are area-weighted means (degenerate triangles carry no
// weight) and are interpolated linearly, then scaled and clamped.
class SizeField {
  public:
    enum Mode { AREA, FINEST };
    SizeField(const Mesh &m, Mode mode, double scale, double hmin, double hmax);
    double operator()(const Point &p) const;
    double hmin() const { return hmin_; }
    double hmax() const { return hmax_; }

  private:
    const Mesh &m_;
    TriangleLocator loc_;
    std::vector<double> size_;   // per point
    Mode mode_;
    double scale_, hmin_, hmax_;
};

}  // namespace suprem_core
