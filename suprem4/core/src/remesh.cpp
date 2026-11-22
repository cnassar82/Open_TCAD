#include "suprem_core/remesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include <CDT.h>
#include "mmg/mmg2d/libmmg2d.h"

#include "suprem_core/brep.hpp"
#include "suprem_core/locate.hpp"

namespace suprem_core {
namespace {

constexpr double PI = 3.14159265358979323846;

uint64_t dkey(int a, int b) { return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b; }
uint64_t ukey(int a, int b) { return a < b ? dkey(a, b) : dkey(b, a); }

double cross(const Point &o, const Point &a, const Point &b)
{
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}
double dist(const Point &a, const Point &b) { return std::hypot(b.x - a.x, b.y - a.y); }
Point unit(double x, double y)
{
    double l = std::hypot(x, y);
    return l > 0 ? Point{x / l, y / l} : Point{0, 0};
}
// angle (radians) at a in triangle a b c
double angle_at(const Point &a, const Point &b, const Point &c)
{
    double ux = b.x - a.x, uy = b.y - a.y, vx = c.x - a.x, vy = c.y - a.y;
    return std::atan2(std::fabs(ux * vy - uy * vx), ux * vx + uy * vy);
}

// ---------------------------------------------------------------- segment grid
struct SegGrid {
    double x0 = 0, y0 = 0, cell = 1;
    int nx = 1, ny = 1;
    std::vector<std::vector<int>> cells;
    mutable std::vector<int> stamp;
    mutable int cur = 0;

    void build(const std::vector<std::pair<Point, Point>> &s)
    {
        double x1 = -1e300, y1 = -1e300, total = 0;
        x0 = y0 = 1e300;
        for (auto &e : s) {
            for (const Point *p : {&e.first, &e.second}) {
                x0 = std::min(x0, p->x); y0 = std::min(y0, p->y);
                x1 = std::max(x1, p->x); y1 = std::max(y1, p->y);
            }
            total += dist(e.first, e.second);
        }
        if (s.empty()) { x0 = y0 = 0; x1 = y1 = 1; }
        cell = std::max(2 * total / std::max<size_t>(1, s.size()), 1e-9);
        nx = std::max(1, std::min(2048, (int)std::ceil((x1 - x0) / cell)));
        ny = std::max(1, std::min(2048, (int)std::ceil((y1 - y0) / cell)));
        cell = std::max({(x1 - x0) / nx, (y1 - y0) / ny, 1e-9}) * (1 + 1e-12);
        cells.assign((size_t)nx * ny, {});
        stamp.assign(s.size(), -1);
        for (int i = 0; i < (int)s.size(); ++i) {
            const Point &a = s[i].first, &b = s[i].second;
            int i0, i1, j0, j1;
            range(std::min(a.x, b.x), std::min(a.y, b.y), std::max(a.x, b.x), std::max(a.y, b.y), i0, j0, i1, j1);
            for (int u = i0; u <= i1; ++u)
                for (int v = j0; v <= j1; ++v) cells[(size_t)v * nx + u].push_back(i);
        }
    }
    void range(double xa, double ya, double xb, double yb, int &i0, int &j0, int &i1, int &j1) const
    {
        i0 = std::clamp((int)std::floor((xa - x0) / cell), 0, nx - 1);
        i1 = std::clamp((int)std::floor((xb - x0) / cell), 0, nx - 1);
        j0 = std::clamp((int)std::floor((ya - y0) / cell), 0, ny - 1);
        j1 = std::clamp((int)std::floor((yb - y0) / cell), 0, ny - 1);
    }
    template <class F> void query(double xa, double ya, double xb, double yb, F f) const
    {
        ++cur;
        int i0, i1, j0, j1;
        range(std::min(xa, xb), std::min(ya, yb), std::max(xa, xb), std::max(ya, yb), i0, j0, i1, j1);
        for (int u = i0; u <= i1; ++u)
            for (int v = j0; v <= j1; ++v)
                for (int k : cells[(size_t)v * nx + u])
                    if (stamp[k] != cur) { stamp[k] = cur; f(k); }
    }
};

// proper or touching intersection of pq and rs
bool segments_touch(const Point &p, const Point &q, const Point &r, const Point &s)
{
    double d1 = cross(r, s, p), d2 = cross(r, s, q), d3 = cross(p, q, r), d4 = cross(p, q, s);
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) return true;
    auto on = [](const Point &a, const Point &b, const Point &c, double d) {
        return d == 0 && std::min(a.x, b.x) <= c.x && c.x <= std::max(a.x, b.x) && std::min(a.y, b.y) <= c.y &&
               c.y <= std::max(a.y, b.y);
    };
    return on(r, s, p, d1) || on(r, s, q, d2) || on(p, q, r, d3) || on(p, q, s, d4);
}

// ---------------------------------------------------------------- geometry
struct Canon {
    std::vector<int> pt;
    int left = 0, right = 0;   // right: region (> left) or boundary code
};
struct LoopRef {
    int region = 0;
    std::vector<std::pair<int, bool>> part;   // (canon, reversed)
};
struct Seg {
    int a, b, region, canon, idx;   // directed a->b with region on the left; canon segment idx
    bool rev;
};

struct Geo {
    std::vector<Point> P;
    std::vector<Canon> C;
    std::vector<LoopRef> L;

    // directed segments of one loop, in loop order (segment k starts at vertex k)
    std::vector<Seg> loop_segments(const LoopRef &l) const
    {
        std::vector<Seg> s;
        for (auto [c, rev] : l.part) {
            const auto &v = C[c].pt;
            if (!rev)
                for (int i = 0; i + 1 < (int)v.size(); ++i) s.push_back({v[i], v[i + 1], l.region, c, i, false});
            else
                for (int i = (int)v.size() - 1; i > 0; --i) s.push_back({v[i], v[i - 1], l.region, c, i - 1, true});
        }
        return s;
    }
};

Geo build_geo(const Mesh &m, const Brep &b)
{
    Geo g;
    g.P = m.point;
    g.L.resize(b.loop.size());
    for (size_t i = 0; i < b.loop.size(); ++i) g.L[i].region = b.loop[i].region;
    std::map<std::vector<int>, int> canon_of;
    for (const Chain &c : b.chain) {
        std::vector<int> pts = c.point;
        if (pts.size() > 2 && pts.front() == pts.back()) {   // closed chain: start at the smallest index
            pts.pop_back();
            std::rotate(pts.begin(), std::min_element(pts.begin(), pts.end()), pts.end());
            pts.push_back(pts.front());
        }
        bool rev = c.right >= 0 && c.left_region > c.right;
        if (rev) std::reverse(pts.begin(), pts.end());
        auto it = canon_of.find(pts);
        int id;
        if (it != canon_of.end()) {
            id = it->second;
        } else {
            id = (int)g.C.size();
            Canon k;
            k.pt = pts;
            k.left = rev ? c.right : c.left_region;
            k.right = rev ? c.left_region : c.right;
            g.C.push_back(k);
            canon_of[pts] = id;
        }
        g.L[c.loop].part.push_back({id, rev});
    }
    return g;
}

// Distance from x to segment ab.
double seg_dist(const Point &x, const Point &a, const Point &b)
{
    double dx = b.x - a.x, dy = b.y - a.y, l2 = dx * dx + dy * dy;
    double t = l2 > 0 ? std::clamp(((x.x - a.x) * dx + (x.y - a.y) * dy) / l2, 0.0, 1.0) : 0;
    return std::hypot(a.x + t * dx - x.x, a.y + t * dy - x.y);
}

// Indices (into g.C) of chains with a segment that crosses or touches a
// segment of another chain (or of itself) other than at shared endpoints.
std::set<int> crossing_chains(const Geo &g)
{
    std::vector<std::pair<Point, Point>> geom;
    std::vector<std::array<int, 3>> id;   // a, b, canon
    for (int c = 0; c < (int)g.C.size(); ++c)
        for (size_t i = 0; i + 1 < g.C[c].pt.size(); ++i) {
            geom.push_back({g.P[g.C[c].pt[i]], g.P[g.C[c].pt[i + 1]]});
            id.push_back({g.C[c].pt[i], g.C[c].pt[i + 1], c});
        }
    SegGrid grid;
    grid.build(geom);
    std::set<int> bad;
    for (int i = 0; i < (int)geom.size(); ++i) {
        const Point &p = geom[i].first, &q = geom[i].second;
        grid.query(p.x, p.y, q.x, q.y, [&](int j) {
            if (j <= i) return;
            if (id[i][0] == id[j][0] || id[i][0] == id[j][1] || id[i][1] == id[j][0] || id[i][1] == id[j][1]) return;
            if (segments_touch(p, q, geom[j].first, geom[j].second)) { bad.insert(id[i][2]); bad.insert(id[j][2]); }
        });
    }
    return bad;
}

// Remove chain vertices that are not needed: a vertex goes when the chord
// that replaces it stays within tol of every removed vertex and is not
// longer than the target size.  Chain end points and vertices shared with
// other chains stay.  Chains whose simplification creates an intersection
// are restored.  Returns the number of removed vertices.
int coarsen_chains(Geo &g, const SizeField &h, double tol, RemeshReport &rep)
{
    if (tol <= 0) return 0;
    std::unordered_map<int, int> uses;
    for (const Canon &c : g.C)
        for (size_t i = 0; i + 1 < c.pt.size() || i < c.pt.size(); ++i) uses[c.pt[i]]++;
    for (const Canon &c : g.C) { uses[c.pt.front()] += 2; uses[c.pt.back()] += 2; }
    std::vector<std::vector<int>> original;
    for (const Canon &c : g.C) original.push_back(c.pt);
    int removed = 0;
    for (Canon &c : g.C) {
        const auto &v = c.pt;
        std::vector<int> out{v[0]};
        size_t i = 0;
        while (i + 1 < v.size()) {
            size_t j = i + 1;
            while (j + 1 < v.size() && uses[v[j]] == 1) {
                const Point &a = g.P[v[i]], &b = g.P[v[j + 1]];
                Point mid{0.5 * (a.x + b.x), 0.5 * (a.y + b.y)};
                if (dist(a, b) > std::min({h(a), h(mid), h(b)})) break;
                bool ok = true;
                for (size_t k = i + 1; k <= j && ok; ++k) ok = seg_dist(g.P[v[k]], a, b) <= tol;
                if (!ok) break;
                ++j;
            }
            removed += (int)(j - i - 1);
            out.push_back(v[j]);
            i = j;
        }
        c.pt = out;
    }
    for (int round = 0; round < 10; ++round) {
        std::set<int> bad = crossing_chains(g);
        int restored = 0;
        for (int c : bad)
            if (g.C[c].pt != original[c]) {
                removed -= (int)(original[c].size() - g.C[c].pt.size());
                g.C[c].pt = original[c];
                ++restored;
            }
        if (!restored) {
            if (!bad.empty()) throw std::runtime_error("input chains intersect");
            break;
        }
        rep.notes.push_back("coarsening of " + std::to_string(restored) + " chains undone (would intersect)");
    }
    return removed;
}

void refine_chains(Geo &g, const SizeField &h)
{
    for (Canon &c : g.C) {
        std::vector<int> out{c.pt[0]};
        for (size_t i = 0; i + 1 < c.pt.size(); ++i) {
            Point a = g.P[c.pt[i]], b = g.P[c.pt[i + 1]];
            Point mid{0.5 * (a.x + b.x), 0.5 * (a.y + b.y)};
            double hh = std::min({h(a), h(mid), h(b)});
            int n = std::max(1, (int)std::lround(dist(a, b) / hh));
            for (int k = 1; k < n; ++k) {
                double t = (double)k / n;
                g.P.push_back({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)});
                out.push_back((int)g.P.size() - 1);
            }
            out.push_back(c.pt[i + 1]);
        }
        c.pt = out;
    }
}

// ---------------------------------------------------------------- thin films
struct Pair { int v, w, region; };

// One scan over the current region-loop vertices (a snapshot: points
// inserted by this scan are never projected again, so the construction
// terminates and cannot cascade in wedges).  From each vertex a ray is cast
// along the inward normal; a facing boundary of the same region closer than
// the film threshold makes the vertex a film vertex.  The hit is matched to
// an end vertex of the hit segment if that is within the matching tolerance,
// otherwise a vertex is inserted on the hit segment.  Every match is
// recorded as a pair (source vertex, target vertex).
void film_pass(Geo &g, const RemeshOptions &o, std::vector<Pair> &pairs)
{
    std::vector<std::vector<Seg>> LS;
    std::vector<Seg> all;
    for (const LoopRef &l : g.L) {
        LS.push_back(g.loop_segments(l));
        all.insert(all.end(), LS.back().begin(), LS.back().end());
    }
    std::vector<std::pair<Point, Point>> geom;
    for (const Seg &s : all) geom.push_back({g.P[s.a], g.P[s.b]});
    SegGrid grid;
    grid.build(geom);
    const double thr = o.film_threshold, cos_corner = std::cos(20 * PI / 180), cos_face = std::cos(30 * PI / 180);
    struct Ins { double p, sep; int source, region; };
    std::map<int, std::map<int, std::vector<Ins>>> ins;   // canon -> segment idx -> insertions
    for (const auto &segs : LS) {
        const int n = (int)segs.size();
        for (int k = 0; k < n; ++k) {
            const int v = segs[k].a, u = segs[(k + n - 1) % n].a, w = segs[k].b;
            const Point &pv = g.P[v];
            Point d1 = unit(pv.x - g.P[u].x, pv.y - g.P[u].y), d2 = unit(g.P[w].x - pv.x, g.P[w].y - pv.y);
            if (d1.x * d2.x + d1.y * d2.y < cos_corner) continue;
            Point nr = unit(-d1.y - d2.y, d1.x + d2.x);
            Point end{pv.x + thr * nr.x, pv.y + thr * nr.y};
            int best = -1; double bs = 2, bt = 0;
            grid.query(pv.x, pv.y, end.x, end.y, [&](int i) {
                const Seg &s = all[i];
                if (s.region != segs[k].region || s.a == v || s.b == v) return;
                const Point &a = g.P[s.a], &b = g.P[s.b];
                double rx = end.x - pv.x, ry = end.y - pv.y, dx = b.x - a.x, dy = b.y - a.y;
                double den = rx * dy - ry * dx;
                if (den == 0) return;
                double ax = a.x - pv.x, ay = a.y - pv.y;
                double sp = (ax * dy - ay * dx) / den, tp = (ax * ry - ay * rx) / den;
                if (sp <= 1e-9 || sp > 1 || tp < 0 || tp > 1) return;
                if (sp < bs) { bs = sp; best = i; bt = tp; }
            });
            if (best < 0) continue;
            const Seg &s = all[best];
            const Point &a = g.P[s.a], &b = g.P[s.b];
            Point ns = unit(-(b.y - a.y), b.x - a.x);
            if (ns.x * nr.x + ns.y * nr.y > -cos_face) continue;
            // Matching tolerance: small against the film thickness (bounded
            // skew) and against the vertex spacing on both sides (one-to-one).
            double thick = bs * thr, len = dist(a, b), off_a = bt * len, off_b = (1 - bt) * len;
            double src = std::min(dist(g.P[u], pv), dist(pv, g.P[w]));
            double tol = std::min({0.2 * thick, 0.3 * src, 0.3 * len});
            if (off_a <= tol && off_a <= off_b) pairs.push_back({v, s.a, segs[k].region});
            else if (off_b <= tol) pairs.push_back({v, s.b, segs[k].region});
            else ins[s.canon][s.idx].push_back({s.rev ? 1 - bt : bt, tol / len, v, segs[k].region});
        }
    }
    for (auto &[c, idxs] : ins) {
        Canon &K = g.C[c];
        std::vector<int> out{K.pt[0]};
        for (int i = 0; i + 1 < (int)K.pt.size(); ++i) {
            auto it = idxs.find(i);
            if (it != idxs.end()) {
                auto ps = it->second;
                std::sort(ps.begin(), ps.end(), [](const Ins &x, const Ins &y) { return x.p < y.p; });
                double last = 0;
                Point a = g.P[K.pt[i]], b = g.P[K.pt[i + 1]];
                for (const Ins &q : ps) {
                    if (q.p - last < q.sep || 1 - q.p < q.sep) continue;   // would make a sliver segment: unpaired
                    g.P.push_back({a.x + q.p * (b.x - a.x), a.y + q.p * (b.y - a.y)});
                    out.push_back((int)g.P.size() - 1);
                    pairs.push_back({q.source, (int)g.P.size() - 1, q.region});
                    last = q.p;
                }
            }
            out.push_back(K.pt[i + 1]);
        }
        K.pt = out;
    }
}

// Film pairs (deduplicated, mutually non-crossing, not crossing the boundary)
// and their interior layer points.
std::vector<Pair> film_pairs(Geo &g, const RemeshOptions &o, std::vector<int> &layer_points)
{
    std::vector<Pair> raw, pairs;
    film_pass(g, o, raw);
    std::unordered_set<uint64_t> seen;
    std::vector<std::pair<Point, Point>> bnd;
    for (const LoopRef &l : g.L)
        for (const Seg &s : g.loop_segments(l)) bnd.push_back({g.P[s.a], g.P[s.b]});
    std::vector<std::array<int, 2>> bnd_ids;
    for (const LoopRef &l : g.L)
        for (const Seg &s : g.loop_segments(l)) bnd_ids.push_back({s.a, s.b});
    SegGrid grid;
    grid.build(bnd);
    std::vector<std::pair<Point, Point>> kept_geom;
    for (const Pair &p : raw) {
        if (p.v == p.w || !seen.insert(ukey(p.v, p.w)).second) continue;
        const Point &a = g.P[p.v], &b = g.P[p.w];
        bool bad = false;
        grid.query(a.x, a.y, b.x, b.y, [&](int i) {
            if (bad) return;
            int s0 = bnd_ids[i][0], s1 = bnd_ids[i][1];
            if (s0 == p.v || s0 == p.w || s1 == p.v || s1 == p.w) return;
            if (segments_touch(a, b, bnd[i].first, bnd[i].second)) bad = true;
        });
        for (size_t j = 0; j < pairs.size() && !bad; ++j) {
            const Pair &q = pairs[j];
            if (q.v == p.v || q.v == p.w || q.w == p.v || q.w == p.w) continue;
            if (segments_touch(a, b, g.P[q.v], g.P[q.w])) bad = true;
        }
        if (!bad) pairs.push_back(p);
    }
    for (const Pair &p : pairs) {
        const Point &a0 = g.P[p.v], &b0 = g.P[p.w];
        int n = std::max(1, std::min(o.film_layers, (int)std::floor(dist(a0, b0) / o.film_min_layer + 1e-9)));
        for (int k = 1; k < n; ++k) {
            double t = (double)k / n;
            const Point &a = g.P[p.v], &b = g.P[p.w];
            g.P.push_back({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)});
            layer_points.push_back((int)g.P.size() - 1);
        }
    }
    return pairs;
}

// ---------------------------------------------------------------- working mesh
struct Work {
    std::vector<Point> P;
    std::vector<int> pid;                 // Geo point index, or -1 for points created by Mmg
    std::vector<std::array<int, 3>> T;    // counter-clockwise
    std::vector<int> R;                   // region
    std::vector<char> req;                // fixed film triangle
};

void make_ccw(Work &w, int t)
{
    auto &v = w.T[t];
    if (cross(w.P[v[0]], w.P[v[1]], w.P[v[2]]) < 0) std::swap(v[1], v[2]);
}

// edges on region boundaries (different regions on both sides, or one side)
std::vector<std::array<int, 2>> constrained_edges(const Work &w)
{
    std::unordered_map<uint64_t, int> first;
    std::vector<std::array<int, 2>> out;
    std::unordered_map<uint64_t, int> count;
    for (int t = 0; t < (int)w.T.size(); ++t)
        for (int k = 0; k < 3; ++k) {
            uint64_t key = ukey(w.T[t][k], w.T[t][(k + 1) % 3]);
            auto it = first.find(key);
            if (it == first.end()) first[key] = w.R[t];
            else if (it->second != w.R[t]) it->second = -1;
            count[key]++;
        }
    for (auto &kv : first)
        if (kv.second == -1 || count[kv.first] == 1) out.push_back({(int)(kv.first >> 32), (int)(kv.first & 0xffffffffu)});
    return out;
}

Work triangulate(const Geo &g, const std::vector<int> &layer, const std::vector<char> &film, const RemeshOptions &o,
                 RemeshReport &rep)
{
    std::vector<int> local(g.P.size(), -1);
    Work w;
    auto use = [&](int p) {
        if (local[p] < 0) { local[p] = (int)w.P.size(); w.P.push_back(g.P[p]); w.pid.push_back(p); }
        return local[p];
    };
    std::vector<CDT::Edge> edges;
    std::unordered_set<uint64_t> fixed;
    for (const Canon &c : g.C)
        for (size_t i = 0; i + 1 < c.pt.size(); ++i) {
            int a = use(c.pt[i]), b = use(c.pt[i + 1]);
            if (fixed.insert(ukey(a, b)).second) edges.push_back(CDT::Edge(a, b));
        }
    for (int p : layer) use(p);
    CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::NotAllowed, 0.0);
    std::vector<CDT::V2d<double>> V;
    for (const Point &p : w.P) V.push_back(CDT::V2d<double>(p.x, p.y));
    cdt.insertVertices(V);
    cdt.insertEdges(edges);
    cdt.eraseSuperTriangle();
    if (cdt.vertices.size() != w.P.size()) throw std::runtime_error("CDT changed the vertex set (duplicate points?)");
    std::vector<std::array<int, 3>> T;
    for (const auto &t : cdt.triangles) {
        std::array<int, 3> v{(int)t.vertices[0], (int)t.vertices[1], (int)t.vertices[2]};
        if (cross(w.P[v[0]], w.P[v[1]], w.P[v[2]]) < 0) std::swap(v[1], v[2]);
        T.push_back(v);
    }
    std::unordered_map<uint64_t, int> dir;   // directed edge -> triangle on its left
    for (int t = 0; t < (int)T.size(); ++t)
        for (int k = 0; k < 3; ++k) dir[dkey(T[t][k], T[t][(k + 1) % 3])] = t;
    std::vector<int> label(T.size(), -1), stack;
    for (const LoopRef &l : g.L)
        for (const Seg &s : g.loop_segments(l)) {
            auto it = dir.find(dkey(local[s.a], local[s.b]));
            if (it == dir.end()) throw std::runtime_error("boundary segment missing from CDT");
            int t = it->second;
            if (label[t] >= 0 && label[t] != l.region) throw std::runtime_error("region labelling conflict");
            if (label[t] < 0) { label[t] = l.region; stack.push_back(t); }
        }
    while (!stack.empty()) {
        int t = stack.back(); stack.pop_back();
        for (int k = 0; k < 3; ++k) {
            int a = T[t][k], b = T[t][(k + 1) % 3];
            if (fixed.count(ukey(a, b))) continue;
            auto it = dir.find(dkey(b, a));
            if (it == dir.end()) continue;
            int n = it->second;
            if (label[n] == label[t]) continue;
            if (label[n] >= 0) throw std::runtime_error("region labelling conflict (unclosed boundary)");
            label[n] = label[t];
            stack.push_back(n);
        }
    }
    const double maxang = o.film_max_angle * PI / 180;
    for (int t = 0; t < (int)T.size(); ++t) {
        if (label[t] < 0) continue;
        w.T.push_back(T[t]);
        w.R.push_back(label[t]);
        bool f = true;
        for (int v : T[t]) f = f && film[w.pid[v]];
        if (f) {
            const Point &a = w.P[T[t][0]], &b = w.P[T[t][1]], &c = w.P[T[t][2]];
            f = std::max({angle_at(a, b, c), angle_at(b, c, a), angle_at(c, a, b)}) <= maxang;
        }
        w.req.push_back(f);
        rep.film_triangles += f;
    }
    // drop points only used by outside triangles (cannot happen for valid input, kept for safety)
    return w;
}

// Exact-coordinate restoration: Mmg rescales coordinates internally, so
// output vertices are matched to input vertices within a small tolerance.
struct PointHash {
    double cell;
    std::unordered_map<uint64_t, std::vector<int>> h;
    const std::vector<Point> *P;
    uint64_t key(long i, long j) const { return ((uint64_t)(uint32_t)i << 32) | (uint32_t)j; }
    PointHash(const std::vector<Point> &pts, double c) : cell(c), P(&pts)
    {
        for (int i = 0; i < (int)pts.size(); ++i) h[key(std::lround(pts[i].x / cell), std::lround(pts[i].y / cell))].push_back(i);
    }
    int find(const Point &p, double tol) const
    {
        long i0 = std::lround(p.x / cell), j0 = std::lround(p.y / cell);
        int best = -1; double bd = tol;
        for (long i = i0 - 1; i <= i0 + 1; ++i)
            for (long j = j0 - 1; j <= j0 + 1; ++j) {
                auto it = h.find(key(i, j));
                if (it == h.end()) continue;
                for (int k : it->second) {
                    double d = dist(p, (*P)[k]);
                    if (d <= bd) { bd = d; best = k; }
                }
            }
        return best;
    }
};

bool run_mmg(Work &w, const SizeField &h, const RemeshOptions &o, RemeshReport &rep)
{
    MMG5_pMesh mesh = nullptr;
    MMG5_pSol met = nullptr;
    MMG2D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh, MMG5_ARG_ppMet, &met, MMG5_ARG_end);
    auto E = constrained_edges(w);
    const int np = (int)w.P.size(), nt = (int)w.T.size(), na = (int)E.size();
    bool ok = MMG2D_Set_meshSize(mesh, np, nt, 0, na) == 1;
    for (int i = 0; ok && i < np; ++i) ok = MMG2D_Set_vertex(mesh, w.P[i].x, w.P[i].y, 0, i + 1) == 1;
    for (int t = 0; ok && t < nt; ++t) {
        ok = MMG2D_Set_triangle(mesh, w.T[t][0] + 1, w.T[t][1] + 1, w.T[t][2] + 1, w.R[t] + 1, t + 1) == 1;
        if (ok && w.req[t]) ok = MMG2D_Set_requiredTriangle(mesh, t + 1) == 1;
    }
    std::vector<char> on_edge(np, 0);
    for (int e = 0; ok && e < na; ++e) {
        ok = MMG2D_Set_edge(mesh, E[e][0] + 1, E[e][1] + 1, 1, e + 1) == 1 && MMG2D_Set_requiredEdge(mesh, e + 1) == 1;
        on_edge[E[e][0]] = on_edge[E[e][1]] = 1;
    }
    for (int i = 0; ok && i < np; ++i)
        if (on_edge[i]) ok = MMG2D_Set_requiredVertex(mesh, i + 1) == 1;
    ok = ok && MMG2D_Set_solSize(mesh, met, MMG5_Vertex, np, MMG5_Scalar) == 1;
    for (int i = 0; ok && i < np; ++i) ok = MMG2D_Set_scalarSol(met, h(w.P[i]), i + 1) == 1;
    ok = ok && MMG2D_Set_iparameter(mesh, met, MMG2D_IPARAM_verbose, o.verbose > 1 ? 1 : -1) == 1;
    ok = ok && MMG2D_Set_iparameter(mesh, met, MMG2D_IPARAM_angle, 0) == 1;
    ok = ok && MMG2D_Set_iparameter(mesh, met, MMG2D_IPARAM_nosurf, 1) == 1;
    ok = ok && MMG2D_Set_dparameter(mesh, met, MMG2D_DPARAM_hmin, h.hmin()) == 1;
    ok = ok && MMG2D_Set_dparameter(mesh, met, MMG2D_DPARAM_hmax, h.hmax()) == 1;
    ok = ok && MMG2D_Set_dparameter(mesh, met, MMG2D_DPARAM_hgrad, o.hgrad) == 1;
    int status = ok ? MMG2D_mmg2dlib(mesh, met) : MMG5_STRONGFAILURE;
    if (status != MMG5_SUCCESS) {
        rep.mmg_failures++;
        rep.notes.push_back(std::string("Mmg2d ") + (status == MMG5_LOWFAILURE ? "low" : "strong") +
                            " failure; mesh kept from the previous stage");
        MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh, MMG5_ARG_ppMet, &met, MMG5_ARG_end);
        return false;
    }
    MMG5_int onp, ont, onq, ona;
    MMG2D_Get_meshSize(mesh, &onp, &ont, &onq, &ona);
    double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
    for (const Point &p : w.P) { x0 = std::min(x0, p.x); x1 = std::max(x1, p.x); y0 = std::min(y0, p.y); y1 = std::max(y1, p.y); }
    const double tol = 1e-9 * std::max(x1 - x0, y1 - y0);
    PointHash ph(w.P, std::max(tol * 4, 1e-12));
    Work out;
    for (MMG5_int i = 0; i < onp; ++i) {
        double x, y; MMG5_int ref; int corner, required;
        MMG2D_Get_vertex(mesh, &x, &y, &ref, &corner, &required);
        int k = ph.find({x, y}, tol);
        if (k >= 0) { out.P.push_back(w.P[k]); out.pid.push_back(w.pid[k]); }
        else { out.P.push_back({x, y}); out.pid.push_back(-1); }
    }
    for (MMG5_int t = 0; t < ont; ++t) {
        MMG5_int a, b, c, ref; int required;
        MMG2D_Get_triangle(mesh, &a, &b, &c, &ref, &required);
        out.T.push_back({(int)a - 1, (int)b - 1, (int)c - 1});
        out.R.push_back((int)ref - 1);
        out.req.push_back((char)required);
        make_ccw(out, (int)out.T.size() - 1);
    }
    MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh, MMG5_ARG_ppMet, &met, MMG5_ARG_end);
    // every required edge must have survived exactly
    std::unordered_set<uint64_t> have;
    for (auto &t : out.T)
        for (int k = 0; k < 3; ++k) {
            int a = out.pid[t[k]], b = out.pid[t[(k + 1) % 3]];
            if (a >= 0 && b >= 0) have.insert(ukey(a, b));
        }
    for (auto &e : E)
        if (!have.count(ukey(w.pid[e[0]], w.pid[e[1]]))) {
            rep.notes.push_back("Mmg2d lost a required edge; mesh kept from the previous stage");
            rep.mmg_failures++;
            return false;
        }
    w = std::move(out);
    return true;
}

// Lawson flips of non-Delaunay edges inside regions.
int delaunay_flips(Work &w)
{
    int total = 0;
    for (int pass = 0; pass < 200; ++pass) {
        std::unordered_map<uint64_t, std::pair<int, int>> dir;   // directed edge -> (triangle, local index of edge start)
        for (int t = 0; t < (int)w.T.size(); ++t)
            for (int k = 0; k < 3; ++k) dir[dkey(w.T[t][k], w.T[t][(k + 1) % 3])] = {t, k};
        std::vector<char> touched(w.T.size(), 0);
        int flips = 0;
        for (int t = 0; t < (int)w.T.size(); ++t) {
            if (touched[t]) continue;
            for (int k = 0; k < 3 && !touched[t]; ++k) {
                int a = w.T[t][k], b = w.T[t][(k + 1) % 3], c = w.T[t][(k + 2) % 3];
                auto it = dir.find(dkey(b, a));
                if (it == dir.end()) continue;
                int u = it->second.first;
                if (touched[u] || w.R[u] != w.R[t]) continue;
                int d = w.T[u][(it->second.second + 2) % 3];
                double ac = angle_at(w.P[c], w.P[a], w.P[b]), ad = angle_at(w.P[d], w.P[b], w.P[a]);
                if (ac + ad <= PI + 1e-9) continue;
                if (cross(w.P[a], w.P[d], w.P[c]) <= 0 || cross(w.P[d], w.P[b], w.P[c]) <= 0) continue;
                w.T[t] = {a, d, c};
                w.T[u] = {d, b, c};
                char rq = (char)(w.req[t] && w.req[u]);
                w.req[t] = w.req[u] = rq;
                touched[t] = touched[u] = 1;
                ++flips;
            }
        }
        total += flips;
        if (!flips) break;
    }
    return total;
}

// Region-boundary edges (as Geo point pairs) with an opposite angle above the limit.
std::vector<std::pair<int, int>> encroached(const Work &w, double limit_deg)
{
    std::unordered_set<uint64_t> con;
    for (auto &e : constrained_edges(w)) con.insert(ukey(e[0], e[1]));
    const double lim = limit_deg * PI / 180;
    std::set<std::pair<int, int>> out;
    for (int t = 0; t < (int)w.T.size(); ++t)
        for (int k = 0; k < 3; ++k) {
            int a = w.T[t][k], b = w.T[t][(k + 1) % 3], c = w.T[t][(k + 2) % 3];
            if (!con.count(ukey(a, b))) continue;
            if (angle_at(w.P[c], w.P[a], w.P[b]) <= lim) continue;
            int pa = w.pid[a], pb = w.pid[b];
            if (pa >= 0 && pb >= 0) out.insert({std::min(pa, pb), std::max(pa, pb)});
        }
    return {out.begin(), out.end()};
}

// One round of local encroachment repair (Ruppert): every region-boundary
// edge with an opposite angle above the limit is split at its midpoint in the
// working mesh (both adjacent triangles split) and in its chain.  Returns the
// number of splits.
int split_encroached(Work &w, Geo &g, double limit_deg)
{
    std::unordered_set<uint64_t> con;
    for (auto &e : constrained_edges(w)) con.insert(ukey(e[0], e[1]));
    std::unordered_map<uint64_t, std::pair<int, int>> dir;
    for (int t = 0; t < (int)w.T.size(); ++t)
        for (int k = 0; k < 3; ++k) dir[dkey(w.T[t][k], w.T[t][(k + 1) % 3])] = {t, k};
    const double lim = limit_deg * PI / 180;
    std::vector<char> touched(w.T.size(), 0);
    std::map<std::pair<int, int>, int> chain_split;   // (pa, pb) sorted -> new Geo point
    int n = 0;
    const int nt0 = (int)w.T.size();
    for (int t = 0; t < nt0; ++t)
        for (int k = 0; k < 3; ++k) {
            if (touched[t]) break;
            int a = w.T[t][k], b = w.T[t][(k + 1) % 3], c = w.T[t][(k + 2) % 3];
            if (!con.count(ukey(a, b)) || angle_at(w.P[c], w.P[a], w.P[b]) <= lim) continue;
            int pa = w.pid[a], pb = w.pid[b];
            if (pa < 0 || pb < 0) continue;
            auto other = dir.find(dkey(b, a));
            if (other != dir.end() && touched[other->second.first]) continue;
            Point mid{0.5 * (w.P[a].x + w.P[b].x), 0.5 * (w.P[a].y + w.P[b].y)};
            g.P.push_back(mid);
            int m = (int)w.P.size();
            w.P.push_back(mid);
            w.pid.push_back((int)g.P.size() - 1);
            chain_split[{std::min(pa, pb), std::max(pa, pb)}] = (int)g.P.size() - 1;
            auto split = [&](int tri, int x, int y, int z) {   // triangle (x, y, z) with edge x->y
                w.T[tri] = {x, m, z};
                w.T.push_back({m, y, z});
                w.R.push_back(w.R[tri]);
                w.req.push_back(w.req[tri]);
                touched[tri] = 1;
                touched.push_back(1);
            };
            split(t, a, b, c);
            if (other != dir.end()) {
                int u = other->second.first, j = other->second.second;
                split(u, w.T[u][j], w.T[u][(j + 1) % 3], w.T[u][(j + 2) % 3]);
            }
            ++n;
        }
    for (Canon &c : g.C) {
        std::vector<int> out{c.pt[0]};
        for (size_t i = 0; i + 1 < c.pt.size(); ++i) {
            auto it = chain_split.find({std::min(c.pt[i], c.pt[i + 1]), std::max(c.pt[i], c.pt[i + 1])});
            if (it != chain_split.end()) out.push_back(it->second);
            out.push_back(c.pt[i + 1]);
        }
        c.pt = out;
    }
    return n;
}

Mesh assemble(const Mesh &old, const Geo &g, const Work &w, RemeshReport &rep)
{
    Mesh m;
    m.version = old.version; m.dim = old.dim; m.nvrt = old.nvrt; m.nedg = old.nedg;
    m.has_model = old.has_model; m.orientation = old.orientation; m.last_temperature = old.last_temperature;
    m.implanted = old.implanted; m.trailer = old.trailer; m.species = old.species;
    m.region_material = old.region_material;
    m.point = w.P;
    std::unordered_map<uint64_t, int> code;   // Geo segment -> boundary code
    for (const Canon &c : g.C)
        if (c.right < 0)
            for (size_t i = 0; i + 1 < c.pt.size(); ++i) code[ukey(c.pt[i], c.pt[i + 1])] = c.right;
    std::unordered_map<uint64_t, int> dir;
    for (int t = 0; t < (int)w.T.size(); ++t)
        for (int k = 0; k < 3; ++k) dir[dkey(w.T[t][k], w.T[t][(k + 1) % 3])] = t;
    int missing = 0;
    for (int t = 0; t < (int)w.T.size(); ++t) {
        Triangle tr;
        tr.region = w.R[t];
        for (int k = 0; k < 3; ++k) tr.v[k] = w.T[t][k];
        for (int j = 0; j < 3; ++j) {   // face j is opposite v[j]: edge v[j+1] -> v[j+2]
            int a = tr.v[(j + 1) % 3], b = tr.v[(j + 2) % 3];
            auto it = dir.find(dkey(b, a));
            if (it != dir.end()) { tr.nb[j] = it->second; continue; }
            auto c = (w.pid[a] >= 0 && w.pid[b] >= 0) ? code.find(ukey(w.pid[a], w.pid[b])) : code.end();
            if (c == code.end()) { tr.nb[j] = BC_OFFSET; ++missing; }
            else tr.nb[j] = c->second;
        }
        m.triangle.push_back(tr);
    }
    if (missing) rep.notes.push_back(std::to_string(missing) + " boundary faces without a boundary code");
    std::set<std::pair<int, int>> nodes;
    for (const Triangle &t : m.triangle)
        for (int k = 0; k < 3; ++k) {
            nodes.insert({t.v[k], m.region_material[t.region]});
            if (t.nb[k] == EXPOSED || t.nb[k] == BACKSIDE) {   // SUPREM: gas node on surface and backside points
                nodes.insert({t.v[(k + 1) % 3], GAS});
                nodes.insert({t.v[(k + 2) % 3], GAS});
            }
        }
    for (auto &[p, mat] : nodes) m.node.push_back({p, mat, std::vector<double>(m.species.size(), 0.0)});
    m.rebuild_index();
    return m;
}

}  // namespace

// ---------------------------------------------------------------- public API

Quality measure_quality(const Mesh &m)
{
    Quality q;
    q.points = (int)m.point.size();
    q.triangles = (int)m.triangle.size();
    std::map<std::pair<uint64_t, int>, double> coupling;   // (edge, region) -> sum of cot / 2
    for (int t = 0; t < (int)m.triangle.size(); ++t) {
        const Triangle &tr = m.triangle[t];
        double mn = PI, mx = 0;
        for (int k = 0; k < 3; ++k) {
            const Point &a = m.point[tr.v[k]], &b = m.point[tr.v[(k + 1) % 3]], &c = m.point[tr.v[(k + 2) % 3]];
            double ang = angle_at(a, b, c);
            mn = std::min(mn, ang); mx = std::max(mx, ang);
            coupling[{ukey(tr.v[(k + 1) % 3], tr.v[(k + 2) % 3]), tr.region}] += 0.5 / std::tan(ang);
            int nb = tr.nb[k];
            bool region_boundary = nb < 0 || m.triangle[nb].region != tr.region;
            if (region_boundary && ang > 90.5 * PI / 180) q.obtuse_boundary++;
        }
        q.min_angle = std::min(q.min_angle, mn * 180 / PI);
        q.max_angle = std::max(q.max_angle, mx * 180 / PI);
        q.below_20 += mn < 20 * PI / 180;
        q.below_10 += mn < 10 * PI / 180;
        q.above_120 += mx > 120 * PI / 180;
    }
    // couplings are normalised per edge length: cot sum / 2 is dimensionless
    for (auto &kv : coupling) {
        q.negative_couplings += kv.second < -1e-9;
        q.significant_negative += kv.second < -0.01;
        q.worst_coupling = std::min(q.worst_coupling, kv.second);
    }
    return q;
}

std::string format_quality(const Quality &q)
{
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "points %d triangles %d | min angle %.2f max angle %.2f | <20deg %d <10deg %d >120deg %d | "
                  "negative couplings %d | obtuse boundary edges %d | couplings < -0.01: %d worst %.3g",
                  q.points, q.triangles, q.min_angle, q.max_angle, q.below_20, q.below_10, q.above_120,
                  q.negative_couplings, q.obtuse_boundary, q.significant_negative, q.worst_coupling);
    return buf;
}

Mesh remesh(const Mesh &m, const RemeshOptions &o, RemeshReport &rep)
{
    std::vector<std::string> problems;
    Brep b = extract_brep(m, problems);
    auto bad = check_brep(m, b);
    bad.insert(bad.begin(), problems.begin(), problems.end());
    if (!bad.empty()) throw std::runtime_error("input B-rep invalid: " + bad.front());
    SizeField h(m, o.size_mode == "min" ? SizeField::FINEST : SizeField::AREA, o.scale, o.hmin, o.hmax);
    h.adapt_to_solution(o.solution_decades);
    Geo g = build_geo(m, b);
    rep.coarsened = coarsen_chains(g, h, o.coarsen_tol, rep);
    refine_chains(g, h);
    std::vector<int> layer;
    std::vector<Pair> pairs = film_pairs(g, o, layer);
    rep.film_pairs = (int)pairs.size();
    rep.film_points = (int)layer.size();
    std::vector<char> film(g.P.size(), 0);
    for (const Pair &p : pairs) film[p.v] = film[p.w] = 1;
    for (int p : layer) film[p] = 1;
    Work w = triangulate(g, layer, film, o, rep);
    if (o.mmg)
        for (int pass = 0; pass < o.mmg_passes; ++pass)
            if (!run_mmg(w, h, o, rep)) break;
    if (o.flips) rep.flips += delaunay_flips(w);
    for (int round = 0; round < o.encroach_rounds; ++round) {
        int n = split_encroached(w, g, o.encroach_angle);
        if (!n) break;
        rep.encroach_splits += n;
        rep.flips += delaunay_flips(w);
    }
    if (auto enc = encroached(w, o.encroach_angle); !enc.empty())
        rep.notes.push_back(std::to_string(enc.size()) + " boundary edges still above the encroachment angle");
    rep.chain_points = 0;
    for (const Canon &c : g.C) rep.chain_points += (int)c.pt.size() - 1;
    return assemble(m, g, w, rep);
}

std::vector<std::string> compare_geometry(const Mesh &before, const Mesh &after, double tol, size_t max_reports)
{
    std::vector<std::string> bad;
    auto report = [&](const std::string &s) { if (bad.size() < max_reports) bad.push_back(s); };
    if (before.region_material != after.region_material) report("region list differs");
    std::vector<std::string> pr;
    Brep b0 = extract_brep(before, pr), b1 = extract_brep(after, pr);
    std::vector<double> a0(before.region_material.size(), 0), a1(after.region_material.size(), 0),
        perim(before.region_material.size(), 0);
    for (int t = 0; t < (int)before.triangle.size(); ++t) a0[before.triangle[t].region] += 0.5 * orient2(before, t);
    for (int t = 0; t < (int)after.triangle.size(); ++t)
        if (after.triangle[t].region < (int)a1.size()) a1[after.triangle[t].region] += 0.5 * orient2(after, t);
    auto segs = [](const Mesh &m, const Brep &b, std::vector<std::pair<Point, Point>> &seg,
                   std::set<std::pair<double, double>> &pts, std::vector<double> *per) {
        for (const Loop &L : b.loop)
            for (size_t i = 0; i < L.point.size(); ++i) {
                const Point &p = m.point[L.point[i]], &q = m.point[L.point[(i + 1) % L.point.size()]];
                pts.insert({p.x, p.y});
                seg.push_back({p, q});
                if (per) (*per)[L.region] += dist(p, q);
            }
    };
    std::vector<std::pair<Point, Point>> s0, s1;
    std::set<std::pair<double, double>> p0, p1;
    segs(before, b0, s0, p0, &perim);
    segs(after, b1, s1, p1, nullptr);
    for (size_t r = 0; r < std::min(a0.size(), a1.size()); ++r)
        if (std::fabs(a0[r] - a1[r]) > 1e-12 * std::max(1.0, std::fabs(a0[r])) + 1e-15 + tol * perim[r])
            report("region " + std::to_string(r + 1) + " area " + std::to_string(a0[r]) + " -> " + std::to_string(a1[r]));
    // two-sided Hausdorff distance between the boundaries (vertex to segment)
    auto check = [&](const std::set<std::pair<double, double>> &pts, const std::set<std::pair<double, double>> &exact,
                     const std::vector<std::pair<Point, Point>> &seg, const char *what) {
        SegGrid grid;
        grid.build(seg);
        const double lim = tol + 1e-10, box = lim + 1e-9;
        for (auto &p : pts) {
            if (exact.count(p)) continue;
            Point x{p.first, p.second};
            double best = 1e300;
            grid.query(x.x - box, x.y - box, x.x + box, x.y + box,
                       [&](int i) { best = std::min(best, seg_dist(x, seg[i].first, seg[i].second)); });
            if (best > lim) report(std::string(what) + " boundary vertex (" + std::to_string(x.x) + ", " +
                                   std::to_string(x.y) + ") is " + std::to_string(best) + " um off the other boundary");
        }
    };
    check(p1, p0, s0, "new");
    check(p0, p1, s1, "old");
    return bad;
}

std::vector<std::string> check_gas_nodes(const Mesh &m, size_t max_reports)
{
    std::vector<std::string> bad;
    std::set<int> exposed, gas;
    for (const Triangle &t : m.triangle)
        for (int j = 0; j < 3; ++j)
            if (t.nb[j] == EXPOSED || t.nb[j] == BACKSIDE) { exposed.insert(t.v[(j + 1) % 3]); exposed.insert(t.v[(j + 2) % 3]); }
    for (const Node &n : m.node)
        if (n.material == GAS) gas.insert(n.point);
    for (int p : exposed)
        if (!gas.count(p) && bad.size() < max_reports) bad.push_back("exposed/backside point " + std::to_string(p + 1) + " has no gas node");
    for (int p : gas)
        if (!exposed.count(p) && bad.size() < max_reports) bad.push_back("gas node at a point neither exposed nor backside " + std::to_string(p + 1));
    return bad;
}

}  // namespace suprem_core
