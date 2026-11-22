#include "suprem_core/locate.hpp"
#include "suprem_core/transfer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace suprem_core {

TriangleLocator::TriangleLocator(const Mesh &m, int material) : m_(m)
{
    if (m.triangle.empty()) return;
    double x1 = -1e300, y1 = -1e300;
    x0_ = y0_ = 1e300;
    for (const Point &p : m.point) {
        x0_ = std::min(x0_, p.x); y0_ = std::min(y0_, p.y);
        x1 = std::max(x1, p.x); y1 = std::max(y1, p.y);
    }
    double w = std::max(x1 - x0_, 1e-12), h = std::max(y1 - y0_, 1e-12);
    cell_ = std::sqrt(w * h / (double)m.triangle.size()) * 2;
    nx_ = std::max(1, std::min(4096, (int)std::ceil(w / cell_)));
    ny_ = std::max(1, std::min(4096, (int)std::ceil(h / cell_)));
    cell_ = std::max(w / nx_, h / ny_) * (1 + 1e-12);
    bucket_.assign((size_t)nx_ * ny_, {});
    for (int t = 0; t < (int)m.triangle.size(); ++t) {
        if (material >= 0 && m.material_of(t) != material) continue;
        double a = 1e300, b = 1e300, c = -1e300, d = -1e300;
        for (int v : m.triangle[t].v) {
            a = std::min(a, m.point[v].x); b = std::min(b, m.point[v].y);
            c = std::max(c, m.point[v].x); d = std::max(d, m.point[v].y);
        }
        int i0 = std::clamp((int)((a - x0_) / cell_), 0, nx_ - 1), i1 = std::clamp((int)((c - x0_) / cell_), 0, nx_ - 1);
        int j0 = std::clamp((int)((b - y0_) / cell_), 0, ny_ - 1), j1 = std::clamp((int)((d - y0_) / cell_), 0, ny_ - 1);
        for (int i = i0; i <= i1; ++i)
            for (int j = j0; j <= j1; ++j) bucket_[(size_t)j * nx_ + i].push_back(t);
    }
}

// Minimum barycentric coordinate (>= 0 inside); fills bary.
double TriangleLocator::score(int t, const Point &p, double bary[3]) const
{
    const Point &a = m_.point[m_.triangle[t].v[0]], &b = m_.point[m_.triangle[t].v[1]], &c = m_.point[m_.triangle[t].v[2]];
    double d = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    if (d == 0) return -1e300;
    bary[1] = ((p.x - a.x) * (c.y - a.y) - (p.y - a.y) * (c.x - a.x)) / d;
    bary[2] = ((b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x)) / d;
    bary[0] = 1 - bary[1] - bary[2];
    return std::min({bary[0], bary[1], bary[2]});
}

int TriangleLocator::locate(const Point &p, double bary[3]) const
{
    if (m_.triangle.empty()) return -1;
    int ci = std::clamp((int)std::floor((p.x - x0_) / cell_), 0, nx_ - 1);
    int cj = std::clamp((int)std::floor((p.y - y0_) / cell_), 0, ny_ - 1);
    int best = -1; double best_s = -1e300, b[3];
    for (int ring = 0; ring <= std::max(nx_, ny_); ++ring) {
        for (int i = ci - ring; i <= ci + ring; ++i)
            for (int j = cj - ring; j <= cj + ring; ++j) {
                if (i < 0 || j < 0 || i >= nx_ || j >= ny_) continue;
                if (std::max(std::abs(i - ci), std::abs(j - cj)) != ring) continue;
                for (int t : bucket_[(size_t)j * nx_ + i]) {
                    double s = score(t, p, b);
                    if (s > best_s) { best_s = s; best = t; std::copy(b, b + 3, bary); }
                }
            }
        if (best >= 0 && (best_s >= -1e-12 || ring >= 1)) break;
    }
    if (best < 0) return -1;
    if (best_s < 0) {   // outside: clip to the nearest triangle
        double s = 0;
        for (int k = 0; k < 3; ++k) { bary[k] = std::max(0.0, bary[k]); s += bary[k]; }
        for (int k = 0; k < 3; ++k) bary[k] = s > 0 ? bary[k] / s : 1.0 / 3;
    }
    return best;
}

SizeField::SizeField(const Mesh &m, Mode mode, double scale, double hmin, double hmax)
    : m_(m), loc_(m), size_(m.point.size(), 0), finest_(m.point.size(), 0), mode_(mode), scale_(scale), hmin_(hmin),
      hmax_(hmax)
{
    // Per triangle: area-equivalent size (side of the equilateral triangle of
    // equal area) and smallest altitude; per vertex: area-weighted means, so
    // degenerate (near-zero-area) triangles have no influence.
    std::vector<double> w(m.point.size(), 0);
    for (int t = 0; t < (int)m.triangle.size(); ++t) {
        const Triangle &tr = m.triangle[t];
        double area = 0.5 * std::fabs(orient2(m, t)), lmax = 0;
        for (int k = 0; k < 3; ++k) {
            const Point &p = m.point[tr.v[k]], &q = m.point[tr.v[(k + 1) % 3]];
            lmax = std::max(lmax, std::hypot(q.x - p.x, q.y - p.y));
        }
        double h_area = std::sqrt(4 * area / std::sqrt(3.0)), h_alt = lmax > 0 ? 2 * area / lmax : 0;
        for (int v : tr.v) {
            size_[v] += area * (mode == AREA ? h_area : h_alt);
            finest_[v] += area * h_alt;
            w[v] += area;
        }
    }
    for (size_t p = 0; p < size_.size(); ++p) {
        size_[p] = w[p] > 0 ? size_[p] / w[p] : hmax_;
        finest_[p] = w[p] > 0 ? finest_[p] / w[p] : hmax_;
    }
}

void SizeField::adapt_to_solution(double decades)
{
    if (decades <= 0) return;
    const Mesh &m = m_;
    std::vector<double> peak(m.species.size(), 0);
    for (const Node &n : m.node)
        for (size_t s = 0; s < m.species.size(); ++s) peak[s] = std::max(peak[s], n.value[s]);
    std::vector<double> sol(m.point.size(), 1e300);
    for (int t = 0; t < (int)m.triangle.size(); ++t) {
        const Triangle &tr = m.triangle[t];
        int mat = m.material_of(t);
        double var = 0;
        for (size_t s = 0; s < m.species.size(); ++s) {
            if (!conserved_species(m.species[s]) || peak[s] <= 0) continue;
            double floor = 1e-5 * peak[s], lo = 1e300, hi = -1e300;
            for (int v : tr.v) {
                int n = m.node_at(v, mat);
                if (n < 0) continue;
                double l = std::log10(std::max(m.node[n].value[s], floor));
                lo = std::min(lo, l); hi = std::max(hi, l);
            }
            var = std::max(var, hi - lo);
        }
        if (var <= decades) continue;
        double area = 0.5 * std::fabs(orient2(m, t));
        double h = std::sqrt(4 * area / std::sqrt(3.0)) * decades / var;
        for (int v : tr.v) sol[v] = std::min(sol[v], h);
    }
    for (size_t p = 0; p < size_.size(); ++p)
        if (sol[p] < size_[p]) size_[p] = std::max(finest_[p], sol[p]);
}

double SizeField::operator()(const Point &p) const
{
    double b[3];
    int t = loc_.locate(p, b);
    if (t < 0) return hmax_;
    double h = 0;
    for (int k = 0; k < 3; ++k) h += b[k] * size_[m_.triangle[t].v[k]];
    return std::clamp(scale_ * h, hmin_, hmax_);
}

}  // namespace suprem_core
