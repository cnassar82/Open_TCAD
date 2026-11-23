#include "suprem_core/evolve.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>

#include <CDT.h>
#include "clipper2/clipper.h"

#include "suprem_core/brep.hpp"
#include "suprem_core/locate.hpp"
#include "suprem_core/str_file.hpp"
#include "suprem_core/transfer.hpp"

namespace suprem_core {
namespace {

using namespace Clipper2Lib;
// Integer grid: 1e-6 um.  Coarser than SUPREM's float coordinates (~7e-8 um
// at device scale), so converting a result to SUPREM's floats can never move a
// vertex across an edge.
constexpr double S = 1e6;   // grid units per um

Point64 to64(const Point &p) { return Point64((int64_t)std::llround(p.x * S), (int64_t)std::llround(p.y * S)); }
Point from64(const Point64 &p) { return {p.x / S, p.y / S}; }

Paths64 shape_paths(const Shape &s)
{
    Paths64 out;
    for (const auto &l : s.loops) {
        Path64 p;
        for (const Point &q : l) p.push_back(to64(q));
        out.push_back(p);
    }
    return out;
}

Paths64 open_(const Paths64 &a, double r)
{
    if (a.empty()) return a;
    return InflatePaths(InflatePaths(a, -r, JoinType::Miter, EndType::Polygon, 10.0), r, JoinType::Miter,
                        EndType::Polygon, 10.0);
}

// Pieces that are real (wider than about d/8 somewhere), not offset-rounding
// slivers.  A real piece keeps all of its connected thin region within d/4,
// including the narrowest part (the tip of a closing wedge), which would
// otherwise stay behind and close by motion alone.
Paths64 significant(const Paths64 &a, double d)
{
    if (a.empty()) return a;
    Paths64 core = Intersect(open_(a, d / 16), a, FillRule::NonZero);
    if (core.empty()) return core;
    return Intersect(a, InflatePaths(core, d / 4, JoinType::Miter, EndType::Polygon, 10.0), FillRule::NonZero);
}

double area(const Paths64 &a) { return std::fabs(Area(a)) / (S * S); }

// Split a path set into connected pieces (outer + its holes each).
std::vector<Paths64> pieces(const Paths64 &a)
{
    std::vector<Paths64> out;
    if (a.empty()) return out;
    PolyTree64 tree;
    Clipper64 c;
    c.AddSubject(a);
    c.Execute(ClipType::Union, FillRule::NonZero, tree);
    std::vector<const PolyPath64 *> stack;
    for (size_t i = 0; i < tree.Count(); ++i) stack.push_back(tree.Child(i));
    while (!stack.empty()) {
        const PolyPath64 *outer = stack.back();
        stack.pop_back();
        Paths64 p{outer->Polygon()};
        for (size_t h = 0; h < outer->Count(); ++h) {
            p.push_back(outer->Child(h)->Polygon());
            for (size_t g = 0; g < outer->Child(h)->Count(); ++g) stack.push_back(outer->Child(h)->Child(g));
        }
        out.push_back(p);
    }
    return out;
}

Point any_point(const Paths64 &p)
{
    if (p.empty() || p[0].empty()) return {0, 0};
    return from64(p[0][0]);
}

// winding number of loop l around q (doubles)
int winding(const std::vector<Point> &l, const Point &q)
{
    int w = 0;
    for (size_t i = 0; i < l.size(); ++i) {
        const Point &a = l[i], &b = l[(i + 1) % l.size()];
        double c = (b.x - a.x) * (q.y - a.y) - (q.x - a.x) * (b.y - a.y);
        if (a.y <= q.y) { if (b.y > q.y && c > 0) ++w; }
        else if (b.y <= q.y && c < 0) --w;
    }
    return w;
}

uint64_t dkey(int a, int b) { return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b; }

double tri_area2(const Mesh &m, int t) { return orient2(m, t); }

}  // namespace

std::vector<Shape> mesh_shapes(const Mesh &m)
{
    std::vector<std::string> pr;
    Brep b = extract_brep(m, pr);
    std::vector<Shape> s(m.region_material.size());
    for (size_t r = 0; r < s.size(); ++r) s[r].material = m.region_material[r];
    for (const Loop &L : b.loop) {
        std::vector<Point> l;
        for (int p : L.point) l.push_back(m.point[p]);
        s[L.region].loops.push_back(l);
    }
    return s;
}

std::vector<Shape> mesh_shapes_union(const Mesh &m)
{
    std::map<int, Paths64> mat;
    for (size_t t = 0; t < m.triangle.size(); ++t) {
        const Triangle &tr = m.triangle[t];
        Path64 p{to64(m.point[tr.v[0]]), to64(m.point[tr.v[1]]), to64(m.point[tr.v[2]])};
        if (Area(p) < 0) std::reverse(p.begin(), p.end());
        if (Area(p) == 0) continue;
        mat[m.material_of((int)t)].push_back(p);
    }
    std::vector<Shape> out;
    for (auto &kv : mat) {
        Paths64 u = Union(kv.second, FillRule::NonZero);
        for (const Paths64 &piece : pieces(u)) {
            Shape s;
            s.material = kv.first;
            for (const Path64 &p : piece) {
                std::vector<Point> l;
                for (const Point64 &q : p) l.push_back(from64(q));
                if (l.size() >= 3) s.loops.push_back(l);
            }
            if (!s.loops.empty()) out.push_back(s);
        }
    }
    return out;
}

BoundaryCode boundary_codes_of(const Mesh &m)
{
    struct Seg { Point a, b; int code; };
    std::vector<Seg> segs;
    for (const Triangle &t : m.triangle)
        for (int j = 0; j < 3; ++j)
            if (t.nb[j] < 0 && t.nb[j] != EXPOSED)
                segs.push_back({m.point[t.v[(j + 1) % 3]], m.point[t.v[(j + 2) % 3]], t.nb[j]});
    return [segs](const Point &a, const Point &b) {
        const double tol = 1e-6;
        for (const Seg &s : segs) {
            double dx = s.b.x - s.a.x, dy = s.b.y - s.a.y, l = std::hypot(dx, dy);
            if (l == 0) continue;
            auto off = [&](const Point &p) {
                double t = ((p.x - s.a.x) * dx + (p.y - s.a.y) * dy) / (l * l);
                double d = std::fabs((p.x - s.a.x) * dy - (p.y - s.a.y) * dx) / l;
                return t > -tol / l && t < 1 + tol / l && d < tol;
            };
            if (off(a) && off(b)) return s.code;
        }
        return EXPOSED;
    };
}

std::vector<RepairEvent> repair_shapes(std::vector<Shape> &shapes, const RepairOptions &o)
{
    std::vector<RepairEvent> ev;
    std::map<int, Paths64> mat;
    for (const Shape &s : shapes) {
        Paths64 p = shape_paths(s);
        mat[s.material].insert(mat[s.material].end(), p.begin(), p.end());
    }
    for (auto &kv : mat) kv.second = Union(kv.second, FillRule::NonZero);
    const double d = o.delta * S;
    auto solid = [&]() {
        Paths64 all;
        for (auto &kv : mat) all.insert(all.end(), kv.second.begin(), kv.second.end());
        return Union(all, FillRule::NonZero);
    };
    Rect64 bb = GetBounds(solid());
    const int64_t margin = (int64_t)(20 * d) + 10;
    Path64 box{Point64(bb.left - margin, bb.top - margin), Point64(bb.right + margin, bb.top - margin),
               Point64(bb.right + margin, bb.bottom + margin), Point64(bb.left - margin, bb.bottom + margin)};
    auto shared = [&](const Paths64 &piece, const Paths64 &other) {
        if (other.empty()) return 0.0;
        return area(Intersect(InflatePaths(piece, d / 4, JoinType::Miter, EndType::Polygon, 10.0), other,
                              FillRule::NonZero));
    };
    auto rank = [&](int m) {
        auto it = std::find(o.fill_priority.begin(), o.fill_priority.end(), m);
        return it == o.fill_priority.end() ? (int)o.fill_priority.size() : (int)(it - o.fill_priority.begin());
    };
    // Each piece of space changes owner at most once per repair: a sliver
    // between two materials is thin for either owner, and would otherwise be
    // passed back and forth.
    Paths64 moved;
    for (int iter = 0; iter < 5; ++iter) {
        size_t before = ev.size();
        // thin parts of material regions
        auto protected_ = [&](int m) { return std::find(o.protect.begin(), o.protect.end(), m) != o.protect.end(); };
        for (auto &kv : mat) {
            if (protected_(kv.first)) continue;
            Paths64 thin = Difference(kv.second, open_(kv.second, d / 2), FillRule::NonZero);
            if (!moved.empty()) thin = Difference(thin, moved, FillRule::NonZero);
            thin = significant(thin, d);
            for (const Paths64 &piece : pieces(thin)) {
                Paths64 gas = Difference(Paths64{box}, solid(), FillRule::NonZero);
                int best = GAS;
                double bs = shared(piece, gas);
                for (auto &other : mat) {
                    if (other.first == kv.first) continue;
                    double sh = shared(piece, other.second);
                    if (sh > bs) { bs = sh; best = other.first; }
                }
                kv.second = Difference(kv.second, piece, FillRule::NonZero);
                if (best != GAS) mat[best] = Union(mat[best], piece, FillRule::NonZero);
                moved = Union(moved, piece, FillRule::NonZero);
                ev.push_back({"thin", kv.first, best, area(piece), any_point(piece)});
            }
        }
        // crumbs: whole connected components narrower than d/8 everywhere
        // (below the significance filter, e.g. remnants of a removed piece)
        for (auto &kv : mat) {
            for (const Paths64 &piece : pieces(kv.second)) {
                if (!open_(piece, d / 16).empty()) continue;
                if (protected_(kv.first) && area(piece) >= o.delta * o.delta) continue;
                Paths64 gas = Difference(Paths64{box}, solid(), FillRule::NonZero);
                int best = GAS;
                double bs = shared(piece, gas);
                for (auto &other : mat) {
                    if (other.first == kv.first) continue;
                    double sh = shared(piece, other.second);
                    if (sh > bs) { bs = sh; best = other.first; }
                }
                kv.second = Difference(kv.second, piece, FillRule::NonZero);
                if (best != GAS) mat[best] = Union(mat[best], piece, FillRule::NonZero);
                moved = Union(moved, piece, FillRule::NonZero);
                ev.push_back({"crumb", kv.first, best, area(piece), any_point(piece)});
            }
        }
        // ambient gaps and enclosed cavities
        Paths64 gas = Difference(Paths64{box}, solid(), FillRule::NonZero);
        std::vector<std::pair<Paths64, std::string>> fill;
        for (const Paths64 &piece : pieces(significant(Difference(gas, open_(gas, d / 2), FillRule::NonZero), d)))
            fill.push_back({piece, "gap"});
        for (const Paths64 &piece : pieces(gas)) {
            Rect64 pb = GetBounds(piece);
            bool outside = pb.left <= box[0].x || pb.top <= box[0].y || pb.right >= box[2].x || pb.bottom >= box[2].y;
            if (!outside) fill.push_back({piece, "cavity"});
        }
        for (auto &[piece, kind] : fill) {
            int best = GAS, br = 1 << 30;
            for (auto &other : mat)
                if (shared(piece, other.second) > 0 && rank(other.first) < br) { br = rank(other.first); best = other.first; }
            if (best == GAS) continue;
            mat[best] = Union(mat[best], piece, FillRule::NonZero);
            // pieces of the gap that were inside other materials' openings cannot overlap them: remove overlaps
            for (auto &other : mat)
                if (other.first != best) other.second = Difference(other.second, piece, FillRule::NonZero);
            moved = Union(moved, piece, FillRule::NonZero);
            ev.push_back({kind, GAS, best, area(piece), any_point(piece)});
        }
        if (ev.size() == before) break;
    }
    if (ev.empty()) return ev;   // shapes untouched (exact)
    shapes.clear();
    for (auto &kv : mat)
        for (const Paths64 &piece : pieces(kv.second)) {
            Shape s;
            s.material = kv.first;
            for (const Path64 &p : piece) {
                std::vector<Point> l;
                for (const Point64 &q : p) l.push_back(from64(q));
                if (l.size() >= 3) s.loops.push_back(l);
            }
            if (!s.loops.empty()) shapes.push_back(s);
        }
    return ev;
}

int clean_shapes(std::vector<Shape> &shapes, const RepairOptions &o)
{
    std::map<int, Paths64> mat;
    for (const Shape &sh : shapes) {
        Paths64 p = shape_paths(sh);
        mat[sh.material].insert(mat[sh.material].end(), p.begin(), p.end());
    }
    for (auto &kv : mat) kv.second = Union(kv.second, FillRule::NonZero);
    const double d = o.delta * S;
    auto solid = [&]() {
        Paths64 all;
        for (auto &kv : mat) all.insert(all.end(), kv.second.begin(), kv.second.end());
        return Union(all, FillRule::NonZero);
    };
    Rect64 bb = GetBounds(solid());
    const int64_t margin = (int64_t)(20 * d) + 10;
    Path64 box{Point64(bb.left - margin, bb.top - margin), Point64(bb.right + margin, bb.top - margin),
               Point64(bb.right + margin, bb.bottom + margin), Point64(bb.left - margin, bb.bottom + margin)};
    auto touches = [&](const Paths64 &piece, const Paths64 &other) {
        return !other.empty() && !Intersect(InflatePaths(piece, 2.0, JoinType::Miter, EndType::Polygon, 10.0), other,
                                            FillRule::NonZero).empty();
    };
    auto rank = [&](int m) {
        auto it = std::find(o.fill_priority.begin(), o.fill_priority.end(), m);
        return it == o.fill_priority.end() ? (int)o.fill_priority.size() : (int)(it - o.fill_priority.begin());
    };
    int n = 0;
    // ambient slits narrower than delta, including sub-grid ones that are no
    // event.  A slit is bounded by material along its length; the rounding
    // slivers of the opening along a slanted boundary have the open ambient
    // on one side and are left alone.
    Paths64 gas = Difference(Paths64{box}, solid(), FillRule::NonZero);
    Paths64 open_gas = open_(gas, d / 2), all = solid();
    auto contact = [&](const Paths64 &piece, const Paths64 &other) {
        return other.empty() ? 0.0 : area(Intersect(InflatePaths(piece, 2.0, JoinType::Miter, EndType::Polygon, 10.0), other,
                                                    FillRule::NonZero));
    };
    for (const Paths64 &piece : pieces(Difference(gas, open_gas, FillRule::NonZero))) {
        if (contact(piece, all) < 3 * contact(piece, open_gas)) continue;
        int best = GAS, br = 1 << 30;
        for (auto &other : mat)
            if (touches(piece, other.second) && rank(other.first) < br) { br = rank(other.first); best = other.first; }
        if (best == GAS) continue;
        mat[best] = Union(mat[best], piece, FillRule::NonZero);
        for (auto &other : mat)
            if (other.first != best) other.second = Difference(other.second, piece, FillRule::NonZero);
        ++n;
    }
    // needles: parts of any material (protected ones included) narrower than
    // delta/8; a film that thin is no film but geometric debris (e.g. a
    // folded-back spike), which SUPREM's float coordinates turn into a contact
    for (auto &kv : mat) {
        Paths64 raw = Difference(kv.second, open_(kv.second, d / 16), FillRule::NonZero);
        for (const Paths64 &piece : pieces(raw)) {
            // a needle sticks out (own material only at its base); a rounding
            // sliver lies along its own material and is left alone
            Paths64 rest = Difference(kv.second, piece, FillRule::NonZero);
            double own = contact(piece, rest), foreign = contact(piece, gas);
            for (auto &other : mat)
                if (other.first != kv.first) foreign += contact(piece, other.second);
            if (own >= 0.25 * foreign) continue;
            int best = GAS;
            double bs = contact(piece, gas);
            for (auto &other : mat) {
                if (other.first == kv.first) continue;
                double c = contact(piece, other.second);
                if (c > bs) { bs = c; best = other.first; }
            }
            kv.second = Difference(kv.second, piece, FillRule::NonZero);
            if (best != GAS) mat[best] = Union(mat[best], piece, FillRule::NonZero);
            ++n;
        }
    }
    if (!n) return 0;
    shapes.clear();
    for (auto &kv : mat)
        for (const Paths64 &piece : pieces(kv.second)) {
            Shape sh;
            sh.material = kv.first;
            for (const Path64 &p : piece) {
                std::vector<Point> l;
                for (const Point64 &q : p) l.push_back(from64(q));
                if (l.size() >= 3) sh.loops.push_back(l);
            }
            if (!sh.loops.empty()) shapes.push_back(sh);
        }
    return n;
}

// SUPREM regions are simply connected (one boundary loop): a shape with holes
// is cut horizontally through a hole, in the widest vertex-free band of the
// hole's height, until every piece is simple (both sides share the cut).
std::vector<Shape> simple_shapes(const std::vector<Shape> &in)
{
    std::vector<Shape> out;
    const int64_t big = (int64_t)4e15;
    for (const Shape &sh : in) {
        std::vector<Paths64> work(1, Union(shape_paths(sh), FillRule::NonZero));
        int guard = 0;
        while (!work.empty() && guard++ < 1000) {
            Paths64 cur = work.back();
            work.pop_back();
            for (const Paths64 &piece : pieces(cur)) {
                if (piece.size() == 1) {
                    Shape s;
                    s.material = sh.material;
                    std::vector<Point> l;
                    for (const Point64 &q : piece[0]) l.push_back(from64(q));
                    if (l.size() >= 3) { s.loops.push_back(l); out.push_back(s); }
                    continue;
                }
                Rect64 hb = GetBounds(piece[1]);
                std::vector<int64_t> ys{hb.top, hb.bottom};
                for (const Path64 &pp : piece)
                    for (const Point64 &q : pp)
                        if (q.y > hb.top && q.y < hb.bottom) ys.push_back(q.y);
                std::sort(ys.begin(), ys.end());
                int64_t gap = -1, yc = (hb.top + hb.bottom) / 2;
                for (size_t g = 1; g < ys.size(); ++g)
                    if (ys[g] - ys[g - 1] > gap) { gap = ys[g] - ys[g - 1]; yc = (ys[g] + ys[g - 1]) / 2; }
                Paths64 upper{Path64{Point64(-big, -big), Point64(big, -big), Point64(big, yc), Point64(-big, yc)}};
                Paths64 lower{Path64{Point64(-big, yc), Point64(big, yc), Point64(big, big), Point64(-big, big)}};
                work.push_back(Intersect(piece, upper, FillRule::NonZero));
                work.push_back(Intersect(piece, lower, FillRule::NonZero));
            }
        }
    }
    return out;
}

Mesh mesh_from_shapes(const std::vector<Shape> &shapes_in, const Mesh &reference, const BoundaryCode &code,
                      const SizeField *size, const RemeshOptions &o, RemeshReport &rep, std::vector<std::string> *notes)
{
    const std::vector<Shape> shapes = simple_shapes(shapes_in);
    // 1. planar straight-line graph on the integer grid
    std::map<std::pair<int64_t, int64_t>, int> vid;
    std::vector<Point64> V;
    auto vertex = [&](const Point64 &p) {
        auto it = vid.find({p.x, p.y});
        if (it != vid.end()) return it->second;
        vid[{p.x, p.y}] = (int)V.size();
        V.push_back(p);
        return (int)V.size() - 1;
    };
    std::vector<std::pair<int, int>> E;
    for (const Shape &s : shapes)
        for (const auto &l : s.loops)
            for (size_t i = 0; i < l.size(); ++i) {
                int a = vertex(to64(l[i])), b = vertex(to64(l[(i + 1) % l.size()]));
                if (a != b) E.push_back({std::min(a, b), std::max(a, b)});
            }
    // split edges at vertices lying exactly on them (T-junctions)
    std::vector<std::pair<int, int>> split;
    {
        std::vector<int> order(V.size());
        for (size_t i = 0; i < V.size(); ++i) order[i] = (int)i;
        std::sort(order.begin(), order.end(), [&](int a, int b) { return V[a].x < V[b].x; });
        std::vector<int64_t> xs;
        for (int i : order) xs.push_back(V[i].x);
        for (auto [a, b] : E) {
            int64_t x0 = std::min(V[a].x, V[b].x), x1 = std::max(V[a].x, V[b].x);
            int64_t y0 = std::min(V[a].y, V[b].y), y1 = std::max(V[a].y, V[b].y);
            std::vector<std::pair<__int128, int>> on;
            auto lo = std::lower_bound(xs.begin(), xs.end(), x0) - xs.begin();
            for (size_t k = lo; k < xs.size() && xs[k] <= x1; ++k) {
                int c = order[k];
                if (c == a || c == b || V[c].y < y0 || V[c].y > y1) continue;
                __int128 cr = (__int128)(V[b].x - V[a].x) * (V[c].y - V[a].y) - (__int128)(V[b].y - V[a].y) * (V[c].x - V[a].x);
                if (cr != 0) continue;
                __int128 t = (__int128)(V[c].x - V[a].x) * (V[b].x - V[a].x) + (__int128)(V[c].y - V[a].y) * (V[b].y - V[a].y);
                on.push_back({t, c});
            }
            std::sort(on.begin(), on.end());
            int prev = a;
            for (auto &q : on) { split.push_back({std::min(prev, q.second), std::max(prev, q.second)}); prev = q.second; }
            split.push_back({std::min(prev, b), std::max(prev, b)});
        }
    }
    std::sort(split.begin(), split.end());
    split.erase(std::unique(split.begin(), split.end()), split.end());
    // 2. CDT (crossings, which exact Booleans should not produce, are resolved)
    CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::TryResolve, 0.0);
    std::vector<CDT::V2d<double>> cv;
    for (const Point64 &p : V) cv.push_back(CDT::V2d<double>(p.x / S, p.y / S));
    cdt.insertVertices(cv);
    std::vector<CDT::Edge> ce;
    for (auto [a, b] : split) ce.push_back(CDT::Edge(a, b));
    cdt.insertEdges(ce);
    cdt.eraseSuperTriangle();
    if (notes && cdt.vertices.size() != V.size())
        notes->push_back(std::to_string(cdt.vertices.size() - V.size()) + " vertices added by the CDT at crossings");
    std::vector<Point> P;
    for (const auto &v : cdt.vertices) {
        Point q{v.x, v.y};
        if (o.float_cm_grid) { q.x = (double)(float)(1e-4 * q.x) * 1e4; q.y = (double)(float)(1e-4 * q.y) * 1e4; }
        P.push_back(q);
    }
    std::vector<std::array<int, 3>> T;
    for (const auto &t : cdt.triangles) {
        std::array<int, 3> v{(int)t.vertices[0], (int)t.vertices[1], (int)t.vertices[2]};
        const Point &a = P[v[0]], &b = P[v[1]], &c = P[v[2]];
        if ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) < 0) std::swap(v[1], v[2]);
        T.push_back(v);
    }
    // 3. labels: material of the shape containing the centroid (priority = shape order of materials)
    std::vector<int> label(T.size(), -1);
    struct Box { double x0, y0, x1, y1; };
    std::vector<Box> sb;
    for (const Shape &s : shapes) {
        Box b{1e300, 1e300, -1e300, -1e300};
        for (const auto &l : s.loops)
            for (const Point &q : l) { b.x0 = std::min(b.x0, q.x); b.y0 = std::min(b.y0, q.y); b.x1 = std::max(b.x1, q.x); b.y1 = std::max(b.y1, q.y); }
        sb.push_back(b);
    }
    for (size_t t = 0; t < T.size(); ++t) {
        Point c{(P[T[t][0]].x + P[T[t][1]].x + P[T[t][2]].x) / 3, (P[T[t][0]].y + P[T[t][1]].y + P[T[t][2]].y) / 3};
        for (size_t k = 0; k < shapes.size() && label[t] < 0; ++k) {
            if (c.x < sb[k].x0 || c.x > sb[k].x1 || c.y < sb[k].y0 || c.y > sb[k].y1) continue;
            int w = 0;
            for (const auto &l : shapes[k].loops) w += winding(l, c);
            if (w != 0) label[t] = (int)k;   // shape index (one region per shape)
        }
    }
    std::unordered_map<uint64_t, int> dir;
    for (size_t t = 0; t < T.size(); ++t)
        for (int k = 0; k < 3; ++k) dir[dkey(T[t][k], T[t][(k + 1) % 3])] = (int)t;
    // unlabelled components: ambient if they reach the hull, else filled (safety net; repair_shapes fills cavities)
    {
        std::vector<int> comp(T.size(), -1);
        for (size_t s0 = 0; s0 < T.size(); ++s0) {
            if (label[s0] >= 0 || comp[s0] >= 0) continue;
            std::vector<int> members{(int)s0}, stack{(int)s0};
            comp[s0] = (int)s0;
            bool hull = false;
            std::map<int, double> border;
            while (!stack.empty()) {
                int t = stack.back(); stack.pop_back();
                for (int k = 0; k < 3; ++k) {
                    int a = T[t][k], b = T[t][(k + 1) % 3];
                    auto it = dir.find(dkey(b, a));
                    if (it == dir.end()) { hull = true; continue; }
                    int n = it->second;
                    if (label[n] >= 0) { border[label[n]] += std::hypot(P[b].x - P[a].x, P[b].y - P[a].y); continue; }
                    if (comp[n] < 0) { comp[n] = (int)s0; members.push_back(n); stack.push_back(n); }
                }
            }
            if (hull || border.empty()) continue;
            int best = std::max_element(border.begin(), border.end(), [](auto &x, auto &y) { return x.second < y.second; })->first;
            for (int t : members) label[t] = best;
            if (notes) notes->push_back("mesh_from_shapes filled an uncovered interior pocket with material " + std::to_string(shapes[best].material));
        }
    }
    // 4. coarse mesh: one region per (simply connected) shape
    Mesh c;
    c.version = reference.version; c.dim = reference.dim; c.nvrt = reference.nvrt; c.nedg = reference.nedg;
    c.has_model = reference.has_model; c.orientation = reference.orientation;
    c.last_temperature = reference.last_temperature; c.implanted = reference.implanted;
    c.trailer = reference.trailer; c.species = reference.species;
    std::vector<int> local(P.size(), -1), keep;
    for (size_t t = 0; t < T.size(); ++t) if (label[t] >= 0) keep.push_back((int)t);
    std::vector<int> region(T.size(), -1);
    for (int s0 : keep) {
        if (region[s0] >= 0) continue;
        int r = (int)c.region_material.size();
        c.region_material.push_back(shapes[label[s0]].material);
        std::vector<int> stack{s0};
        region[s0] = r;
        while (!stack.empty()) {
            int t = stack.back(); stack.pop_back();
            for (int k = 0; k < 3; ++k) {
                auto it = dir.find(dkey(T[t][(k + 1) % 3], T[t][k]));
                if (it == dir.end()) continue;
                int n = it->second;
                if (label[n] == label[t] && region[n] < 0) { region[n] = r; stack.push_back(n); }
            }
        }
    }
    std::vector<int> tindex(T.size(), -1);
    for (int t : keep) {
        Triangle tr;
        tr.region = region[t];
        for (int k = 0; k < 3; ++k) {
            int p = T[t][k];
            if (local[p] < 0) { local[p] = (int)c.point.size(); c.point.push_back(P[p]); }
            tr.v[k] = local[p];
        }
        tindex[t] = (int)c.triangle.size();
        c.triangle.push_back(tr);
    }
    for (int t : keep) {
        Triangle &tr = c.triangle[tindex[t]];
        for (int j = 0; j < 3; ++j) {
            int a = T[t][(j + 1) % 3], b = T[t][(j + 2) % 3];
            auto it = dir.find(dkey(b, a));
            if (it != dir.end() && label[it->second] >= 0) tr.nb[j] = tindex[it->second];
            else tr.nb[j] = code(P[a], P[b]);
        }
    }
    std::set<std::pair<int, int>> nodes;
    for (const Triangle &t : c.triangle)
        for (int k = 0; k < 3; ++k) {
            nodes.insert({t.v[k], c.region_material[t.region]});
            if (t.nb[k] == EXPOSED || t.nb[k] == BACKSIDE) {
                nodes.insert({t.v[(k + 1) % 3], GAS});
                nodes.insert({t.v[(k + 2) % 3], GAS});
            }
        }
    for (auto &[p, m] : nodes) c.node.push_back({p, m, std::vector<double>(c.species.size(), 0.0)});
    c.rebuild_index();
    if (const char *dir = std::getenv("CORE_DEBUG_DIR")) write_str(to_str(c), std::string(dir) + "/coarse_last.str");
    try {
        return remesh(c, o, rep, size);
    } catch (...) {
        if (const char *dir = std::getenv("CORE_DEBUG_DIR")) write_str(to_str(c), std::string(dir) + "/coarse_fail.str");
        throw;
    }
}

void evolve(Mesh &m, const Velocity &v, double time, const EvolveOptions &o, EvolveReport &rep)
{
    const double step = std::min(o.max_step, o.repair.delta / 4);
    auto areas = [](const Mesh &x) {
        std::vector<double> a(x.triangle.size());
        for (size_t t = 0; t < a.size(); ++t) a[t] = tri_area2(x, (int)t);
        return a;
    };
    std::vector<double> ref = areas(m);
    // Target sizes come from the mesh as built at the last rebuild, not from
    // the deformed one: otherwise material compressed by the motion (e.g. a
    // consumed layer) would get finer at every rebuild without bound.
    Mesh built = m;
    auto rebuild = [&](std::vector<Shape> shapes) {
        BoundaryCode code = boundary_codes_of(m);
        SizeField h(built, o.remesh.size_mode == "min" ? SizeField::FINEST : SizeField::AREA, o.remesh.scale,
                    o.remesh.hmin, o.remesh.hmax);
        h.adapt_to_solution(o.remesh.solution_decades);
        RemeshReport rr;
        clean_shapes(shapes, o.repair);
        Mesh fresh = mesh_from_shapes(shapes, m, code, &h, o.remesh, rr, &rep.notes);
        TransferReport tr;
        transfer_fields(m, fresh, tr);
        m = std::move(fresh);
        built = m;
        ref = areas(m);
        rep.rebuilds++;
    };
    double t = 0;
    while (t < time * (1 - 1e-12)) {
        if (rep.steps >= o.max_steps) { rep.failures.push_back("step limit reached"); return; }
        std::vector<Point> V(m.point.size());
        double vmax = 0;
        for (size_t i = 0; i < V.size(); ++i) {
            V[i] = v(m.point[i], t);
            vmax = std::max(vmax, std::hypot(V[i].x, V[i].y));
        }
        double dt = time - t;
        if (vmax > 0) dt = std::min(dt, step / vmax);
        Mesh moved;
        bool ok = false;
        for (int h = 0; h < 40; ++h) {
            moved = m;
            for (size_t i = 0; i < V.size(); ++i) { moved.point[i].x += dt * V[i].x; moved.point[i].y += dt * V[i].y; }
            ok = true;
            for (size_t k = 0; k < moved.triangle.size() && ok; ++k)
                ok = tri_area2(moved, (int)k) > o.min_area_ratio * tri_area2(m, (int)k);
            if (ok) {   // boundaries must stay simple (gaps between fronts have no triangles)
                std::vector<std::string> pr;
                Brep b = extract_brep(moved, pr);
                ok = pr.empty() && check_brep(moved, b, 1).empty();
            }
            if (ok) break;
            dt *= 0.5;
            rep.halvings++;
        }
        if (!ok) {   // cannot move without degenerating: rebuild first
            std::vector<Shape> sh = mesh_shapes(m);
            auto ev = repair_shapes(sh, o.repair);
            rep.events.insert(rep.events.end(), ev.begin(), ev.end());
            rebuild(ev.empty() ? mesh_shapes(m) : sh);
            continue;
        }
        m = std::move(moved);
        t += dt;
        rep.steps++;
        if (std::getenv("CORE_EVOLVE_TRACE"))
            std::fprintf(stderr, "step %d t %.6g dt %.3g points %zu rebuilds %d halvings %d\n", rep.steps, t, dt,
                         m.point.size(), rep.rebuilds, rep.halvings);
        std::vector<Shape> sh = mesh_shapes(m);
        auto ev = repair_shapes(sh, o.repair);
        bool degraded = false;
        for (size_t k = 0; k < m.triangle.size() && !degraded; ++k)
            degraded = tri_area2(m, (int)k) < o.rebuild_area_ratio * ref[k];
        if (!ev.empty() || degraded) {
            rep.events.insert(rep.events.end(), ev.begin(), ev.end());
            rebuild(ev.empty() ? mesh_shapes(m) : sh);
        }
        if (o.check_every_step) {
            auto bad = check_invariants(m, 3);
            std::vector<std::string> pr;
            Brep b = extract_brep(m, pr);
            auto bb = check_brep(m, b, 3);
            bad.insert(bad.end(), pr.begin(), pr.end());
            bad.insert(bad.end(), bb.begin(), bb.end());
            for (auto &x : bad) rep.failures.push_back("step " + std::to_string(rep.steps) + ": " + x);
            if (!bad.empty()) return;
        }
    }
}

}  // namespace suprem_core
