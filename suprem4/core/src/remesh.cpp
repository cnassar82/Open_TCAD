#include "suprem_core/remesh.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <functional>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
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

// Coordinate grid of the consumer.  SUPREM stores float cm: every point the
// remesher creates is put on that grid before it is triangulated, so the
// exact predicates see the coordinates SUPREM will hold.
thread_local bool g_float_cm = false;
Point snap(Point p)
{
    if (g_float_cm) {
        p.x = (double)(float)(1e-4 * p.x) * 1e4;
        p.y = (double)(float)(1e-4 * p.y) * 1e4;
    }
    return p;
}

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

double min_angle(const Point &a, const Point &b, const Point &c)
{
    return std::min({angle_at(a, b, c), angle_at(b, c, a), angle_at(c, a, b)});
}

// Largest minimum angle (radians) over all triangulations of the simple
// polygon P (counter-clockwise), by dynamic programming over its diagonals;
// -1 if P has more than max_n vertices or no valid triangulation.  Used for
// the cavity a removed point leaves (star polygon, typically 4-12 vertices).
double best_cavity_angle(const std::vector<Point> &P, size_t max_n = 24)
{
    const int n = (int)P.size();
    if (n < 3 || (size_t)n > max_n) return -1;
    auto proper = [](const Point &a, const Point &b, const Point &c, const Point &d) {
        double d1 = cross(a, b, c), d2 = cross(a, b, d), d3 = cross(c, d, a), d4 = cross(c, d, b);
        return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0)) && d1 != 0 && d2 != 0 && d3 != 0 && d4 != 0;
    };
    auto inside = [&](const Point &q) {   // even-odd rule
        bool in = false;
        for (int i = 0, j = n - 1; i < n; j = i++)
            if ((P[i].y > q.y) != (P[j].y > q.y) && q.x < (P[j].x - P[i].x) * (q.y - P[i].y) / (P[j].y - P[i].y) + P[i].x)
                in = !in;
        return in;
    };
    std::vector<char> diag(n * n, 0);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            if (j == i + 1 || (i == 0 && j == n - 1)) { diag[i * n + j] = 1; continue; }
            bool ok = inside({0.5 * (P[i].x + P[j].x), 0.5 * (P[i].y + P[j].y)});
            for (int k = 0; ok && k < n; ++k) {
                int l = (k + 1) % n;
                if (k == i || k == j || l == i || l == j) continue;
                if (proper(P[i], P[j], P[k], P[l])) ok = false;
            }
            diag[i * n + j] = ok;
        }
    std::vector<double> best(n * n, -1);
    for (int i = 0; i + 1 < n; ++i) best[i * n + i + 1] = 1e9;
    for (int len = 2; len < n; ++len)
        for (int i = 0; i + len < n; ++i) {
            const int j = i + len;
            if (!diag[i * n + j]) continue;
            double b = -1;
            for (int k = i + 1; k < j; ++k) {
                if (!diag[i * n + k] || !diag[k * n + j] || best[i * n + k] < 0 || best[k * n + j] < 0) continue;
                if (cross(P[i], P[k], P[j]) <= 0) continue;
                b = std::max(b, std::min({best[i * n + k], best[k * n + j], min_angle(P[i], P[k], P[j])}));
            }
            best[i * n + j] = b;
        }
    return best[n - 1];
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
        // clamp in double: a tiny cell (zero-length segments) would overflow the int cast
        nx = (int)std::clamp(std::ceil((x1 - x0) / cell), 1.0, 2048.0);
        ny = (int)std::clamp(std::ceil((y1 - y0) / cell), 1.0, 2048.0);
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

// After vertices were removed from chains: chains that now intersect, duplicate
// another chain's segment or collapse/flip their loop get their original vertex
// lists back (repeated until none is left).  Returns the number of vertices
// restored.
int restore_invalid_chains(Geo &g, std::vector<std::vector<int>> original, RemeshReport &rep, const char *what)
{
    int restored_vertices = 0;
    auto loop_area = [&](const LoopRef &l) {
        double a = 0;
        for (const Seg &sg : g.loop_segments(l)) a += 0.5 * (g.P[sg.a].x * g.P[sg.b].y - g.P[sg.b].x * g.P[sg.a].y);
        return a;
    };
    std::vector<double> area0;
    for (int c = 0; c < (int)g.C.size(); ++c) std::swap(g.C[c].pt, original[c]);
    for (const LoopRef &l : g.L) area0.push_back(loop_area(l));
    for (int c = 0; c < (int)g.C.size(); ++c) std::swap(g.C[c].pt, original[c]);
    for (int round = 0; round < 10; ++round) {
        std::set<int> bad = crossing_chains(g);
        // duplicate segments (two chains collapsed onto each other) ...
        std::map<std::pair<int, int>, int> seen;
        for (int c = 0; c < (int)g.C.size(); ++c)
            for (size_t i = 0; i + 1 < g.C[c].pt.size(); ++i) {
                std::pair<int, int> k{std::min(g.C[c].pt[i], g.C[c].pt[i + 1]), std::max(g.C[c].pt[i], g.C[c].pt[i + 1])};
                auto it = seen.find(k);
                if (it != seen.end()) { bad.insert(c); bad.insert(it->second); }
                else seen[k] = c;
            }
        // ... and loops whose orientation flipped or area collapsed (a region
        // narrower than the tolerance)
        for (size_t li = 0; li < g.L.size(); ++li) {
            double a = loop_area(g.L[li]);
            if (a * area0[li] <= 0 || std::fabs(a) < 0.5 * std::fabs(area0[li]))
                for (auto &[c, rev] : g.L[li].part) bad.insert(c);
        }
        int restored = 0;
        for (int c : bad)
            if (g.C[c].pt != original[c]) {
                restored_vertices += (int)(original[c].size() - g.C[c].pt.size());
                g.C[c].pt = original[c];
                ++restored;
            }
        if (!restored) {
            if (!bad.empty()) throw std::runtime_error("input chains intersect or loops are degenerate");
            break;
        }
        rep.notes.push_back(std::string(what) + " of " + std::to_string(restored) +
                            " chains undone (intersection, duplicate or collapsed loop)");
    }
    return restored_vertices;
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
    removed -= restore_invalid_chains(g, original, rep, "coarsening");
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
                g.P.push_back(snap({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)}));
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
                    g.P.push_back(snap({a.x + q.p * (b.x - a.x), a.y + q.p * (b.y - a.y)}));
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
            g.P.push_back(snap({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)}));
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

// extra: additional constraint edges (Geo point pairs) that do not bound
// regions, e.g. the edges of old triangles an update keeps as they are
Work triangulate(const Geo &g, const std::vector<int> &layer, const std::vector<char> &film, const RemeshOptions &o,
                 RemeshReport &rep, const std::vector<std::array<int, 2>> *extra = nullptr)
{
    std::vector<int> local(g.P.size(), -1);
    Work w;
    // points with identical coordinates (possible after snapping to the
    // consumer's coordinate grid) are one vertex
    std::map<std::pair<double, double>, int> at;
    auto use = [&](int p) {
        if (local[p] < 0) {
            auto it = at.find({g.P[p].x, g.P[p].y});
            if (it != at.end()) local[p] = it->second;
            else {
                local[p] = (int)w.P.size();
                at[{g.P[p].x, g.P[p].y}] = local[p];
                w.P.push_back(g.P[p]);
                w.pid.push_back(p);
            }
        }
        return local[p];
    };
    std::vector<CDT::Edge> edges;
    std::unordered_set<uint64_t> fixed;
    for (const Canon &c : g.C)
        for (size_t i = 0; i + 1 < c.pt.size(); ++i) {
            int a = use(c.pt[i]), b = use(c.pt[i + 1]);
            if (a != b && fixed.insert(ukey(a, b)).second) edges.push_back(CDT::Edge(a, b));
        }
    for (int p : layer) use(p);
    if (extra) {
        std::unordered_set<uint64_t> seen;
        for (auto &e : *extra) {
            int a = use(e[0]), b = use(e[1]);
            if (a != b && !fixed.count(ukey(a, b)) && seen.insert(ukey(a, b)).second) edges.push_back(CDT::Edge(a, b));
        }
    }
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
        else { out.P.push_back(snap({x, y})); out.pid.push_back(-1); }
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
using Triple = std::array<int, 3>;   // sorted Geo point indices of a triangle
Triple triple_of(const Work &w, int t)
{
    Triple k{w.pid[w.T[t][0]], w.pid[w.T[t][1]], w.pid[w.T[t][2]]};
    std::sort(k.begin(), k.end());
    return k;
}

// keep: triangles (Geo triples) whose shared edges are never flipped
int delaunay_flips(Work &w, const std::set<Triple> *keep = nullptr)
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
                if (keep && keep->count(triple_of(w, t)) && keep->count(triple_of(w, u))) continue;
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
// edge with an opposite angle above the limit is split at the foot of the
// perpendicular from that apex (inner 60 % of the edge, else midpoint) in the
// working mesh (both adjacent triangles split) and in its chain.  Returns the
// number of splits.
// Ruppert does not terminate near small input angles (the same edges are
// halved down to rounding); no split creates edges shorter than min_len.
// allow(pa, pb, apex, region): optional filter on the Geo points of the edge, its apex and the region of the encroached triangle
using SplitFilter = std::function<bool(int, int, int, int)>;
int split_encroached(Work &w, Geo &g, double limit_deg, double min_len, const SplitFilter &allow = nullptr)
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
            if (dist(w.P[a], w.P[b]) < 2 * min_len) continue;
            if (allow && !allow(pa, pb, w.pid[c], w.R[t])) continue;
            auto other = dir.find(dkey(b, a));
            if (other != dir.end() && touched[other->second.first]) continue;
            // split at the foot of the perpendicular from the encroaching
            // apex (clamped to the inner 60 %), not at the midpoint: the
            // apex then sees two near-right angles and is resolved by one
            // split.  A midpoint split leaves an apex close below a long
            // edge encroaching one half again, and the edge is halved down
            // to the apex's distance (a fan of needles on the other side:
            // LOCOS, silicon grid columns under the moving interface).
            const double ex = w.P[b].x - w.P[a].x, ey = w.P[b].y - w.P[a].y, L2 = ex * ex + ey * ey;
            double s = std::clamp(((w.P[c].x - w.P[a].x) * ex + (w.P[c].y - w.P[a].y) * ey) / L2, 0.2, 0.8);
            if (std::min(s, 1 - s) * std::sqrt(L2) < min_len) s = 0.5;
            Point mid = snap({w.P[a].x + s * ex, w.P[a].y + s * ey});
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

Mesh remesh(const Mesh &m, const RemeshOptions &o, RemeshReport &rep, const SizeField *size)
{
    struct Grid {
        bool old;
        explicit Grid(bool f) : old(g_float_cm) { g_float_cm = f; }
        ~Grid() { g_float_cm = old; }
    } grid_guard(o.float_cm_grid);
    std::vector<std::string> problems;
    Brep b = extract_brep(m, problems);
    auto bad = check_brep(m, b);
    bad.insert(bad.begin(), problems.begin(), problems.end());
    if (!bad.empty()) throw std::runtime_error("input B-rep invalid: " + bad.front());
    std::unique_ptr<SizeField> own;
    if (!size) {
        own = std::make_unique<SizeField>(m, o.size_mode == "min" ? SizeField::FINEST : SizeField::AREA, o.scale,
                                          o.hmin, o.hmax);
        own->adapt_to_solution(o.solution_decades);
        size = own.get();
    }
    const SizeField &h = *size;
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
        int n = split_encroached(w, g, o.encroach_angle, o.min_split_length);
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

namespace suprem_core {
namespace {

// Bucket hash of points that grows while points are added.
struct RowHash {
    double cell;
    std::unordered_map<int64_t, std::vector<int>> b;
    std::vector<Point> at;
    explicit RowHash(double c) : cell(c) {}
    int64_t key(long i, long j) const { return (int64_t)i * 4000037LL + j; }
    void add(int id, const Point &p)
    {
        if ((int)at.size() <= id) at.resize(id + 1);
        at[id] = p;
        b[key(std::lround(std::floor(p.x / cell)), std::lround(std::floor(p.y / cell)))].push_back(id);
    }
    template <class F> void query(const Point &x, double r, F f) const
    {
        long i0 = std::lround(std::floor((x.x - r) / cell)), i1 = std::lround(std::floor((x.x + r) / cell));
        long j0 = std::lround(std::floor((x.y - r) / cell)), j1 = std::lround(std::floor((x.y + r) / cell));
        for (long i = i0; i <= i1; ++i)
            for (long j = j0; j <= j1; ++j)
                if (auto it = b.find(key(i, j)); it != b.end())
                    for (int id : it->second) f(id);
    }
};

// M1 interface boundary layer (UpdateOptions::interface_h): the interfaces
// between a consumed material (silicon, poly) and oxide, with the normal into
// the consumed material, and the row structure below them.
class Interface {
  public:
    struct Near {
        double d = 0, hn = 0, ht = 0;
        Point n{0, 0}, t{0, 0};
        double metric(const Point &a, const Point &b) const
        {
            const double dx = b.x - a.x, dy = b.y - a.y;
            return std::hypot((dx * n.x + dy * n.y) / hn, (dx * t.x + dy * t.y) / ht);
        }
    };
    struct Vertex { int p; Point n; int region; double ht; double adv = 0; };   // adv: look-ahead clearance (update_mesh)

    Interface(const Geo &g, const Mesh &m, const UpdateOptions &o, std::function<double(const Point &, int)> fine)
        : m_(m), o_(o), fine_(std::move(fine)), loc_(m)
    {
        if (o.interface_h <= 0) return;
        cap_ = 5 * std::max(o.interface_depth, o.interface_h);
        auto consumed = [&](int mat) { return mat == SI || mat == POLY; };
        for (const Canon &c : g.C) {
            const int ml = m.region_material[c.left], mr = c.right >= 0 ? m.region_material[c.right] : -1;
            int region = -1;
            if (consumed(ml) && mr == SIO2) region = c.left;
            else if (consumed(mr) && ml == SIO2) region = c.right;
            if (region < 0 || c.pt.size() < 2) continue;
            const size_t first = seg_.size();
            size_t longest = first;
            for (size_t i = 0; i + 1 < c.pt.size(); ++i) {
                const Point a = g.P[c.pt[i]], b = g.P[c.pt[i + 1]];
                const double len = dist(a, b);
                if (len <= 0) continue;
                const Point t{(b.x - a.x) / len, (b.y - a.y) / len};
                const Point nl{-t.y, t.x};
                const double s = region == c.left ? 1 : -1;
                seg_.push_back({a, b, {s * nl.x, s * nl.y}, t, len, region, c.pt[i], c.pt[i + 1]});
                if (len > seg_[longest].len) longest = seg_.size() - 1;
            }
            if (seg_.size() == first) continue;
            // orientation check on the longest segment (the left/right convention is not trusted blindly)
            const S &L = seg_[longest];
            const double e = std::min(0.25 * L.len, 1e-4);
            const Point mid{0.5 * (L.a.x + L.b.x), 0.5 * (L.a.y + L.b.y)};
            if (region_at({mid.x + e * L.n.x, mid.y + e * L.n.y}) != region &&
                region_at({mid.x - e * L.n.x, mid.y - e * L.n.y}) == region)
                for (size_t k = first; k < seg_.size(); ++k) seg_[k].n = {-seg_[k].n.x, -seg_[k].n.y};
            // vertices: all but the chain ends (junctions), unless the chain is closed
            const bool closed = c.pt.front() == c.pt.back();
            const size_t ns = seg_.size() - first;
            if (!closed) {   // chain ends: junctions with the other boundaries of the region
                ends_.push_back({seg_[first].p, seg_[first].n, region, seg_[first].len});
                ends_.push_back({seg_.back().q, seg_.back().n, region, seg_.back().len});
            }
            for (size_t k = 0; k < ns; ++k) {
                const size_t prev = k > 0 ? first + k - 1 : (closed ? seg_.size() - 1 : SIZE_MAX);
                if (prev == SIZE_MAX) continue;
                const S &s0 = seg_[prev], &s1 = seg_[first + k];
                if (s0.q != s1.p) continue;   // a zero-length segment was skipped
                Point n{s0.n.x + s1.n.x, s0.n.y + s1.n.y};
                const double nn = std::hypot(n.x, n.y);
                if (nn < 1e-6) continue;
                vert_.push_back({s1.p, {n.x / nn, n.y / nn}, region, 0.5 * (s0.len + s1.len)});
            }
        }
        std::vector<std::pair<Point, Point>> s;
        for (const S &x : seg_) s.push_back({x.a, x.b});
        grid_.build(s);
    }
    bool active() const { return !seg_.empty(); }
    const std::vector<Vertex> &vertices() const { return vert_; }
    const std::vector<Vertex> &ends() const { return ends_; }
    // the old triangle at x (-1: outside)
    int old_triangle(const Point &x) const
    {
        double bary[3];
        return loc_.locate(x, bary);
    }
    // x inside one of the interface segments (a chain vertex created after
    // construction, e.g. by an encroachment split): its column direction
    bool vertex_at(const Point &x, Vertex &v) const
    {
        const double e = 1e-9;
        bool found = false;
        grid_.query(x.x - e, x.y - e, x.x + e, x.y + e, [&](int k) {
            const S &q = seg_[k];
            if (found || seg_dist(x, q.a, q.b) > e || dist(x, q.a) <= e || dist(x, q.b) <= e) return;
            v = {-1, q.n, q.region, std::max(dist(x, q.a), dist(x, q.b))};
            found = true;
        });
        return found;
    }

    int region_at(const Point &x) const
    {
        double bary[3];
        const int t = loc_.locate(x, bary);
        if (t < 0 || std::min({bary[0], bary[1], bary[2]}) < -1e-9) return -1;
        return m_.triangle[t].region;
    }
    // normal spacing at depth d below the interface (hf: the material's own size there)
    double hn(double d, double hf) const
    {
        const double h0 = o_.interface_h;
        if (d <= o_.interface_depth) return h0;
        return std::min(std::max(hf, h0), h0 + (o_.interface_growth - 1) * (d - o_.interface_depth));
    }
    // x in the layer zone of an interface of its own region
    bool zone(const Point &x, int region, Near &z) const
    {
        if (seg_.empty() || region < 0) return false;
        double best = 1e300;
        const S *bs = nullptr;
        grid_.query(x.x - cap_, x.y - cap_, x.x + cap_, x.y + cap_, [&](int k) {
            if (seg_[k].region != region) return;
            const double d = seg_dist(x, seg_[k].a, seg_[k].b);
            if (d < best) { best = d; bs = &seg_[k]; }
        });
        if (!bs || best >= cap_) return false;
        // the same depth as the columns reach (rows()), also where the
        // material is already as fine as the layer: every point of the zone
        // is a target or is removed, never left to the isotropic rule
        const double hf = fine_(x, region);
        if (best > end(hf)) return false;
        z.d = best;
        z.n = bs->n;
        z.t = bs->t;
        z.hn = hn(best, hf);
        z.ht = std::max(bs->len, z.hn);
        return true;
    }
    // end of the zone where the material's own size is hf: half a step
    // below the last row rows() places (same grading and rounding)
    double end(double hf) const
    {
        const double h0 = o_.interface_h;
        double d = o_.interface_depth, s = h0 * o_.interface_growth;
        while (true) {
            const double step = std::max(h0, std::round(s / h0) * h0);
            if (step >= hf || d + step >= cap_) return d + 0.5 * std::min(step, hf);
            d += step;
            s *= o_.interface_growth;
        }
    }
    // Lattice coordinate along the normal: every row lies on a level that is a
    // multiple of interface_h in space, so neighbouring columns share their
    // levels (rectangular cells) whenever and wherever they were created.
    double lattice_offset(const Point &x, const Point &n) const
    {
        const double h0 = o_.interface_h, s = (x.x * n.x + x.y * n.y) / h0;
        return std::fabs(s - std::round(s)) * h0;
    }
    // depths of the rows below vertex v (at: its position), all on the lattice
    // phase: depth (mod interface_h) of the column's existing levels; NaN:
    // the absolute lattice (a new column)
    std::vector<double> rows(const Vertex &v, const Point &at, double phase = std::nan("")) const
    {
        std::vector<double> r;
        const double h0 = o_.interface_h, sv = at.x * v.n.x + at.y * v.n.y;
        const double clear = o_.near * h0 + v.adv;   // first level clear of the interface (and of its advance within the step)
        double d = std::isnan(phase) ? std::ceil((sv + clear) / h0) * h0 - sv
                                     : phase + std::ceil((clear - phase) / h0) * h0;
        for (; d <= o_.interface_depth * (1 + 1e-9); d += h0) r.push_back(d);
        double s = h0 * o_.interface_growth;
        d = (r.empty() ? d - h0 : r.back());
        while (true) {
            const double step = std::max(h0, std::round(s / h0) * h0);
            d += step;
            if (d >= cap_ || step >= fine_({at.x + d * v.n.x, at.y + d * v.n.y}, v.region)) break;
            r.push_back(d);
            s *= o_.interface_growth;
        }
        return r;
    }

  private:
    struct S { Point a, b, n, t; double len; int region, p, q; };
    const Mesh &m_;
    const UpdateOptions &o_;
    std::function<double(const Point &, int)> fine_;
    TriangleLocator loc_;
    std::vector<S> seg_;
    std::vector<Vertex> vert_, ends_;
    SegGrid grid_;
    double cap_ = 0;
};

}  // namespace

Mesh update_mesh(const Mesh &m, const Mesh &reference, const UpdateOptions &o, UpdateReport &rep)
{
    // Size fields per material of the reference (each material keeps the
    // spacing it had, also where it has moved: growing oxide is measured
    // against the reference oxide, not the silicon it replaced); fallback:
    // all materials.
    struct Fields { std::unique_ptr<SizeField> fine, area, longest; };
    std::map<int, Fields> fields;
    {
        std::set<int> mats{-1};
        for (int t = 0; t < (int)reference.triangle.size(); ++t) mats.insert(reference.material_of(t));
        for (int mat : mats) {
            Fields &f = fields[mat];
            f.fine = std::make_unique<SizeField>(reference, SizeField::FINEST, o.scale, o.hmin, o.hmax, mat);
            f.area = std::make_unique<SizeField>(reference, SizeField::AREA, o.scale, o.hmin, o.hmax, mat);
            f.longest = std::make_unique<SizeField>(reference, SizeField::LONGEST, o.scale, o.hmin, o.hmax, mat);
            f.area->adapt_to_solution(o.solution_decades);
        }
    }
    auto F = [&](int region) -> const Fields & {
        auto it = region >= 0 ? fields.find(m.region_material[region]) : fields.end();
        return it != fields.end() ? it->second : fields.at(-1);
    };
    auto floor_of = [&](int region) {   // per-material lower bound (growing oxide: SUPREM's grid.oxide)
        if (region < 0 || o.min_size.empty()) return 0.0;
        auto it = o.min_size.find(m.region_material[region]);
        return it == o.min_size.end() ? 0.0 : it->second;
    };
    auto fine = [&](const Point &x, int region) { return std::max((*F(region).fine)(x), floor_of(region)); };
    auto area_h = [&](const Point &x, int region) { return std::max((*F(region).area)(x), floor_of(region)); };
    auto longest = [&](const Point &x, int region) { return std::max((*F(region).longest)(x), floor_of(region)); };
    // on a chain: the finer of its two sides
    auto cfine = [&](const Point &x, const Canon &c) { return std::min(fine(x, c.left), c.right >= 0 ? fine(x, c.right) : 1e300); };
    auto clongest = [&](const Point &x, const Canon &c) {
        return std::min(longest(x, c.left), c.right >= 0 ? longest(x, c.right) : 1e300);
    };
    auto clock0 = std::chrono::steady_clock::now();
    auto lap = [&](const char *what) {
        if (!o.timing) return;
        auto now = std::chrono::steady_clock::now();
        rep.notes.push_back(std::string("time ") + what + " " + std::to_string(std::chrono::duration<double>(now - clock0).count()) + " s");
        clock0 = now;
    };
    struct Grid {
        bool old;
        explicit Grid(bool f) : old(g_float_cm) { g_float_cm = f; }
        ~Grid() { g_float_cm = old; }
    } grid_guard(o.float_cm_grid);
    std::vector<std::string> problems;
    Brep b = extract_brep(m, problems);
    auto bad = check_brep(m, b);
    bad.insert(bad.begin(), problems.begin(), problems.end());
    if (!bad.empty()) throw std::runtime_error("input B-rep invalid: " + bad.front());
    const int n0 = (int)m.point.size();
    RemeshReport rr;
    Geo g = build_geo(m, b);
    std::vector<char> on_chain(n0, 0);
    for (const Canon &c : g.C)
        for (int p : c.pt) on_chain[p] = 1;

    // M1.2 look-ahead (UpdateOptions::step_dt): node velocities from SUPREM's
    // XVEL/YVEL solutions (cm/s, the last flow solve), per point and material.
    int sx = -1, sy = -1;
    for (size_t k = 0; k < m.species.size(); ++k) {
        if (m.species[k] == 8) sx = (int)k;
        if (m.species[k] == 9) sy = (int)k;
    }
    const bool ahead = o.step_dt > 0 && sx >= 0 && sy >= 0;
    std::unordered_map<uint64_t, Point> node_vel;   // (point, material) -> um/s
    double vmax = 0;
    if (ahead)
        for (const Node &nd : m.node) {
            const Point v{1e4 * nd.value[sx], 1e4 * nd.value[sy]};
            node_vel[(uint64_t)nd.point * 64 + (uint64_t)nd.material] = v;
            vmax = std::max(vmax, std::hypot(v.x, v.y));
        }
    auto vel = [&](int p, int mat) -> Point {
        if (p < 0 || p >= n0) return {0, 0};
        auto it = node_vel.find((uint64_t)p * 64 + (uint64_t)mat);
        return it == node_vel.end() ? Point{0, 0} : it->second;
    };
    // a point's own velocity (SUPREM moves an interface point with the
    // consumed material's node: the interface advances into it)
    auto pvel = [&](int p) -> Point {
        if (!ahead || p < 0 || p >= n0) return {0, 0};
        for (int mat : {SI, POLY, SIO2}) {
            auto it = node_vel.find((uint64_t)p * 64 + (uint64_t)mat);
            if (it != node_vel.end()) return it->second;
        }
        for (const auto &[k, v] : node_vel) if (k / 64 == (uint64_t)p) return v;   // (rare: other materials)
        return {0, 0};
    };
    // smallest distance of p from the chord (a, b) within the coming window
    // [0, 2 step_dt] under the point velocities (0 if p crosses the chord)
    auto chord_ahead = [&](int p, int a, int b) {
        const Point P = g.P[p], A = g.P[a], B = g.P[b];
        if (!ahead) return seg_dist(P, A, B);
        const Point vp = pvel(p), va = pvel(a), vb = pvel(b);
        double best = 1e300, s0 = 0;
        for (int k = 0; k <= 8; ++k) {
            const double t = 2 * o.step_dt * k / 8;
            const Point p1{P.x + t * vp.x, P.y + t * vp.y}, a1{A.x + t * va.x, A.y + t * va.y}, b1{B.x + t * vb.x, B.y + t * vb.y};
            const double ex = b1.x - a1.x, ey = b1.y - a1.y, l = std::hypot(ex, ey);
            const double h = l > 0 ? ((p1.x - a1.x) * ey - (p1.y - a1.y) * ex) / l : dist(p1, a1);
            if (k > 0 && h * s0 < 0) return 0.0;
            s0 = h;
            best = std::min(best, std::min(std::fabs(h), seg_dist(p1, a1, b1)));
        }
        return best;
    };
    // shortest length of segment (a, b) within the coming window
    auto len_ahead = [&](int a, int b) {
        const Point A = g.P[a], B = g.P[b];
        double best = dist(A, B);
        if (!ahead) return best;
        const Point va = pvel(a), vb = pvel(b);
        for (int k = 1; k <= 8; ++k) {
            const double t = 2 * o.step_dt * k / 8;
            best = std::min(best, dist({A.x + t * va.x, A.y + t * va.y}, {B.x + t * vb.x, B.y + t * vb.y}));
        }
        return best;
    };
    // clearance a column needs below an interface vertex so that SUPREM's
    // step limit (no element loses half its area within the step) does not
    // cut the step: twice the interface's advance into the material
    auto advance = [&](const Interface::Vertex &v) {
        if (!ahead || v.p < 0) return 0.0;
        const Point w = vel(v.p, m.region_material[v.region]);
        return 2 * std::max(0.0, w.x * v.n.x + w.y * v.n.y) * o.step_dt;
    };

    lap("setup");
    // 1. chain vertices on short segments (geometry kept within tol); not
    // on the side chains of the interface boundary layer (3-), whose
    // vertices inside the zone are its lattice levels
    std::unique_ptr<Interface> zone0;
    if (o.interface_h > 0) zone0 = std::make_unique<Interface>(g, m, o, [&](const Point &x, int region) { return fine(x, region); });
    auto in_zone0 = [&](const Point &x, const Canon &c) {
        if (!zone0 || !zone0->active()) return false;
        const int ml = m.region_material[c.left], mr = c.right >= 0 ? m.region_material[c.right] : -1;
        if (((ml == SI || ml == POLY) && mr == SIO2) || ((mr == SI || mr == POLY) && ml == SIO2)) return false;   // the interface itself
        Interface::Near z;
        return zone0->zone(x, c.left, z) || (c.right >= 0 && zone0->zone(x, c.right, z));
    };
    if (o.tol > 0 && o.short_edge > 0) {
        std::unordered_map<int, int> uses;
        for (const Canon &c : g.C)
            for (size_t i = 0; i + 1 < c.pt.size(); ++i) uses[c.pt[i]]++;
        for (const Canon &c : g.C) { uses[c.pt.front()] += 2; uses[c.pt.back()] += 2; }
        std::vector<std::vector<int>> original;
        for (const Canon &c : g.C) original.push_back(c.pt);
        int removed = 0;
        for (Canon &c : g.C) {
            const std::vector<int> v = c.pt;
            std::vector<int> out{v[0]}, skipped;
            for (size_t i = 1; i + 1 < v.size(); ++i) {
                const int p = v[i], prev = out.back(), next = v[i + 1];
                bool drop = uses[p] == 1 && !in_zone0(g.P[p], c) && dist(g.P[prev], g.P[next]) < clongest(g.P[p], c);
                // short now, or (look-ahead) within the coming step: a
                // junction sliding along a straight side towards p
                const bool short_seg =
                    std::min(len_ahead(prev, p), len_ahead(p, next)) < o.short_edge * cfine(g.P[p], c);
                // on a short segment (now or within the coming step), and
                // within tol of the chord now or (look-ahead) within the step:
                // a consumed corner flattening, whose cap triangle SUPREM
                // would otherwise halve step after step; the vertex is then
                // not more than a quarter of the size off the chord now.
                // Long segments are never merged by the look-ahead (a flat
                // chain is collinear at every vertex: that is resolution,
                // not redundancy).
                drop = drop && short_seg;
                for (size_t k = 0; drop && k <= skipped.size(); ++k) {
                    const int q = k < skipped.size() ? skipped[k] : p;
                    const double now = seg_dist(g.P[q], g.P[prev], g.P[next]);
                    drop = now <= o.tol || (ahead && now <= 0.25 * cfine(g.P[q], c) && chord_ahead(q, prev, next) <= o.tol);
                }
                if (drop) { skipped.push_back(p); ++removed; }
                else { out.push_back(p); skipped.clear(); }
            }
            out.push_back(v.back());
            c.pt = out;
        }
        removed -= restore_invalid_chains(g, original, rr, "short-edge removal");
        rep.removed_boundary = removed;
    }
    lap("short-edge removal");
    // 2. long chain segments
    for (Canon &c : g.C) {
        std::vector<int> out{c.pt[0]};
        for (size_t i = 0; i + 1 < c.pt.size(); ++i) {
            Point a = g.P[c.pt[i]], bb = g.P[c.pt[i + 1]];
            Point mid{0.5 * (a.x + bb.x), 0.5 * (a.y + bb.y)};
            double L = std::min({clongest(a, c), clongest(mid, c), clongest(bb, c)});
            int n = dist(a, bb) > o.long_edge * L ? (int)std::ceil(dist(a, bb) / L) : 1;
            for (int k = 1; k < n; ++k) {
                double t = (double)k / n;
                g.P.push_back(snap({a.x + t * (bb.x - a.x), a.y + t * (bb.y - a.y)}));
                out.push_back((int)g.P.size() - 1);
                rep.inserted_boundary++;
            }
            out.push_back(c.pt[i + 1]);
        }
        c.pt = out;
    }
    lap("chain refinement");
    // 3. interior points: removed near a chain or near an already kept interior point
    std::vector<int> point_region(n0, -1);
    for (const Triangle &t : m.triangle)
        for (int v : t.v) point_region[v] = t.region;
    // Flow quality on removal (3b): an interior point of a flow material that
    // the proximity rule would remove is kept when the cavity it leaves has
    // no triangulation with all angles >= flow_min_angle and none as good as
    // its star now.  Step 4 inserts by the same angle criterion, so the two
    // cannot undo each other (with the proximity rule alone ~250 points per
    // call were removed and re-inserted at the same place).  A sliver star
    // (point against a boundary) still goes: its cavity is better than it.
    std::vector<std::vector<int>> star_of;
    if (o.flow_min_angle > 0) {
        star_of.resize(n0);
        for (int t = 0; t < (int)m.triangle.size(); ++t) {
            const int mt = m.region_material[m.triangle[t].region];
            if (mt == SIO2 || mt == SINI || mt == OXNI)
                for (int v : m.triangle[t].v) star_of[v].push_back(t);
        }
    }
    auto needed_for_flow = [&](int p) {
        if (star_of.empty() || p < 0 || p >= n0 || star_of[p].size() < 3) return false;
        std::unordered_map<int, int> next;
        double star = 1e9;
        const int region = m.triangle[star_of[p][0]].region;
        for (int t : star_of[p]) {
            const Triangle &tr = m.triangle[t];
            if (tr.region != region) return false;
            int j = tr.v[0] == p ? 0 : tr.v[1] == p ? 1 : 2;
            int a = tr.v[(j + 1) % 3], b = tr.v[(j + 2) % 3];
            if (cross(m.point[p], m.point[a], m.point[b]) < 0) std::swap(a, b);
            if (!next.emplace(a, b).second) return false;
            star = std::min(star, min_angle(m.point[p], m.point[a], m.point[b]));
        }
        std::vector<Point> poly;
        int v = next.begin()->first;
        for (size_t k = 0; k < next.size(); ++k) {
            poly.push_back(m.point[v]);
            auto it = next.find(v);
            if (it == next.end()) return false;
            v = it->second;
        }
        if (v != next.begin()->first) return false;   // not a closed link (boundary point)
        const double cav = best_cavity_angle(poly);
        if (cav < 0) return false;
        return cav < o.flow_min_angle * PI / 180 && cav < star;
    };
    std::vector<int> layer;
    std::vector<char> il_flag;   // per g.P index: interface-layer point (kept or inserted)
    Interface iface(g, m, o, [&](const Point &x, int region) { return fine(x, region); });
    std::vector<int> side_points;   // chain vertices on lattice levels (side chains of the layer)
    // 3-. side chains of the boundary layer.  Where an interface ends on
    // another boundary of the consumed region (a reflecting side, a
    // silicon/nitride interface), that chain is the outermost column: within
    // the zone depth its vertices are put on the lattice levels, so the
    // cells next to it are as structured as the rest of the zone (else one
    // coarse boundary vertex fans out to many rows).  Only a straight part
    // (within tol) running along the column direction is changed (exact
    // geometry); vertices already on a level are kept.
    if (iface.active()) {
        std::unordered_map<int, int> uses;
        for (const Canon &c : g.C)
            for (size_t i = 0; i < c.pt.size(); ++i) uses[c.pt[i]] += (i == 0 || i + 1 == c.pt.size()) ? 2 : 1;
        auto is_interface = [&](const Canon &c) {
            const int ml = m.region_material[c.left], mr = c.right >= 0 ? m.region_material[c.right] : -1;
            return ((ml == SI || ml == POLY) && mr == SIO2) || ((mr == SI || mr == POLY) && ml == SIO2);
        };
        for (const Interface::Vertex &e : iface.ends()) {
            for (Canon &c : g.C) {
                if (is_interface(c) || (c.left != e.region && c.right != e.region) || c.pt.size() < 2) continue;
                const bool fwd = c.pt.front() == e.p;
                if (!fwd && c.pt.back() != e.p) continue;
                std::vector<int> v(c.pt);
                if (!fwd) std::reverse(v.begin(), v.end());
                const Point E = g.P[e.p];
                const double L1 = dist(E, g.P[v[1]]);
                if (L1 <= 0) continue;
                const Point u{(g.P[v[1]].x - E.x) / L1, (g.P[v[1]].y - E.y) / L1};
                if (u.x * e.n.x + u.y * e.n.y < std::cos(M_PI / 6)) continue;   // not along the column direction
                Interface::Vertex col{e.p, u, e.region, e.ht};
                col.adv = advance(e);
                const std::vector<double> lv = iface.rows(col, E);
                if (lv.empty()) continue;
                const double reach = lv.back() + 0.5 * iface.hn(lv.back(), fine(E, e.region));
                // straight part within reach: interior vertices on the line, used by this chain only
                size_t k = 1;   // v[1..k-1]: the part to rebuild; v[k]: first vertex kept beyond it
                while (k + 1 < v.size()) {
                    const Point &q = g.P[v[k]];
                    const double along = (q.x - E.x) * u.x + (q.y - E.y) * u.y;
                    const double off = std::fabs((q.x - E.x) * u.y - (q.y - E.y) * u.x);
                    if (along > reach || off > o.tol || uses[v[k]] != 1) break;
                    ++k;
                }
                const Point &Q = g.P[v[k]];
                const double stop = (Q.x - E.x) * u.x + (Q.y - E.y) * u.y;
                if ((Q.x - E.x) * u.y - (Q.y - E.y) * u.x > o.tol || -((Q.x - E.x) * u.y - (Q.y - E.y) * u.x) > o.tol) continue;
                std::vector<int> part{e.p};
                size_t j = 1;
                for (double d : lv) {
                    const double h = iface.hn(d, fine(E, e.region));
                    if (d > stop - o.near * h) break;
                    int id = -1;
                    for (; j < k; ++j) {   // an old vertex on this level stands for it
                        const Point &q = g.P[v[j]];
                        const double along = (q.x - E.x) * u.x + (q.y - E.y) * u.y;
                        if (along < d - 0.25 * h) continue;
                        if (along <= d + 0.25 * h) id = v[j++];
                        break;
                    }
                    if (id < 0) {
                        g.P.push_back(snap({E.x + d * u.x, E.y + d * u.y}));
                        id = (int)g.P.size() - 1;
                        rep.interface_inserted++;
                    }
                    part.push_back(id);
                    side_points.push_back(id);
                }
                part.insert(part.end(), v.begin() + k, v.end());
                if (!fwd) std::reverse(part.begin(), part.end());
                c.pt = part;
            }
        }
    }
    std::vector<std::pair<Point, Point>> seg;
    std::vector<std::pair<int, int>> seg_id;
    for (const Canon &c : g.C)
        for (size_t i = 0; i + 1 < c.pt.size(); ++i) {
            seg.push_back({g.P[c.pt[i]], g.P[c.pt[i + 1]]});
            seg_id.push_back({c.pt[i], c.pt[i + 1]});
        }
    SegGrid sg;
    sg.build(seg);
    auto near_chain = [&](const Point &x, double r) {
        bool hit = false;
        sg.query(x.x - r, x.y - r, x.x + r, x.y + r, [&](int k) {
            if (!hit && seg_dist(x, seg[k].first, seg[k].second) < r) hit = true;
        });
        return hit;
    };
    // An interior point p that a boundary of its material approaches: closer
    // than r plus twice the approach within the step (the same rule as the
    // columns: SUPREM's half-area step limit would otherwise cut the step to
    // a fraction, again at every step while the gap halves, SUPREM's own
    // ChooseKillNodes in the core's terms).
    auto overrun = [&](int p, const Point &x, int region, double r) {
        if (!ahead || vmax <= 0) return false;
        const int mat = m.region_material[region];
        const Point vp = vel(p, mat);
        const double R = r + 2 * (vmax + std::hypot(vp.x, vp.y)) * o.step_dt;
        bool hit = false;
        sg.query(x.x - R, x.y - R, x.x + R, x.y + R, [&](int k) {
            if (hit) return;
            const Point &a = seg[k].first, &b = seg[k].second;
            const double d = seg_dist(x, a, b);
            if (d >= R || d <= 0) return;
            const double ex = b.x - a.x, ey = b.y - a.y, l2 = ex * ex + ey * ey;
            const double u = l2 > 0 ? std::clamp(((x.x - a.x) * ex + (x.y - a.y) * ey) / l2, 0.0, 1.0) : 0;
            const Point c{a.x + u * ex, a.y + u * ey}, nrm{(x.x - c.x) / d, (x.y - c.y) / d};
            const Point va = vel(seg_id[k].first, mat), vb = vel(seg_id[k].second, mat);
            const double wx = (1 - u) * va.x + u * vb.x - vp.x, wy = (1 - u) * va.y + u * vb.y - vp.y;
            const double approach = wx * nrm.x + wy * nrm.y;   // > 0: the boundary closes in on p
            if (d < r + 2 * std::max(0.0, approach) * o.step_dt) hit = true;
        });
        return hit;
    };
    RowHash zh(0.01);   // accepted zone points (g.P ids)
    // The column below interface vertex v (at: its position): lattice levels
    // of rows(), stopping where the column leaves its material or meets
    // another boundary; levels closer than 0.6 spacings to accepted zone
    // points are skipped (columns converging on a concave interface).
    // match: an old point standing for the target (g.P id), or -1.
    auto column = [&](const Interface::Vertex &v, const Point &at,
                      const std::function<int(const Point &, const Interface::Near &, int)> &match, double phase) {
        std::vector<int> added;
        for (double d : iface.rows(v, at, phase)) {
            const Point x = snap({at.x + d * v.n.x, at.y + d * v.n.y});
            if (iface.region_at(x) != v.region) break;   // left the material
            Interface::Near z;
            z.d = d;
            z.n = v.n;
            z.t = {-v.n.y, v.n.x};
            z.hn = iface.hn(d, fine(x, v.region));
            z.ht = std::max(v.ht, z.hn);
            if (near_chain(x, o.near * z.hn)) break;   // another boundary
            const double R = std::max(z.ht, z.hn);
            bool clash = false;
            zh.query(x, 0.6 * R, [&](int q) {
                if (!clash && z.metric(x, zh.at[q]) < 0.6) clash = true;
            });
            if (clash) continue;
            int id = match ? match(x, z, v.region) : -1;
            if (id < 0) {
                g.P.push_back(x);
                id = (int)g.P.size() - 1;
                rep.interface_inserted++;
                added.push_back(id);
            }
            il_flag.resize(g.P.size(), 0);
            layer.push_back(id);
            il_flag[id] = 1;
            zh.add(id, g.P[id]);
            rep.interface_points++;
        }
        return added;
    };
    {
        std::vector<int> cand;
        std::vector<std::pair<Point, Point>> dots;
        for (int p = 0; p < n0; ++p)
            if (!on_chain[p]) { cand.push_back(p); dots.push_back({g.P[p], g.P[p]}); }
        SegGrid pg;
        pg.build(dots);
        std::vector<char> kept(cand.size(), 0), in_zone(cand.size(), 0);
        il_flag.assign(g.P.size(), 0);
        for (int p : side_points) { il_flag[p] = 1; zh.add(p, g.P[p]); }
        // 3a. interface boundary layer (structured zone).  Its point set is
        // exactly the columns below the interface vertices, on lattice levels
        // fixed in space (neighbouring columns share them and a level stays
        // where it is while the interface moves).  An old point within a
        // quarter spacing of a target stands for it (its values are kept);
        // other old points in the zone are removed.  So between calls only
        // the rows the interface eats and the rows added at the deep end
        // change.
        if (iface.active()) {
            RowHash oh(0.01);   // old zone points (cand indices)
            for (size_t i = 0; i < cand.size(); ++i) {
                Interface::Near z;
                if (iface.zone(g.P[cand[i]], point_region[cand[i]], z)) { in_zone[i] = 1; oh.add((int)i, g.P[cand[i]]); }
            }
            std::vector<char> matched(cand.size(), 0);
            // A kept old point may sit off its column (the interface moves
            // and turns under it).  Its offset along the interface is capped
            // by its depth, not by the column spacing: a point at depth d
            // leaning t towards a neighbouring interface segment makes the
            // angle opposite that segment 90 deg + atan(t / d), i.e. a
            // coupling of -t / 2d on the interface edge (LOCOS: 3.7 nm lean
            // at 2 nm depth gave -0.79).  t <= 0.1 d keeps it above -0.05.
            auto match = [&](const Point &x, const Interface::Near &z, int region) {
                int best = -1;
                double bm = 0.25;
                const double lean = 0.1 * z.d;
                oh.query(x, 0.25 * std::max(z.ht, z.hn), [&](int i) {
                    if (matched[i] || point_region[cand[i]] != region) return;
                    const Point &q = g.P[cand[i]];
                    if (std::fabs((q.x - x.x) * z.t.x + (q.y - x.y) * z.t.y) > lean) return;
                    const double mm = z.metric(x, q);
                    if (mm < bm) { bm = mm; best = i; }
                });
                if (best < 0) return -1;
                matched[best] = kept[best] = 1;
                return cand[best];
            };
            // The levels of a column stay where its old points are (the
            // absolute lattice x.n moves by |x| dn when the interface
            // rotates, e.g. in a bird's beak): phase from the old zone point
            // under the vertex nearest to the interface.
            auto phase_of = [&](const Interface::Vertex &v, const Point &at) {
                const double h0 = o.interface_h, reach = o.interface_depth;
                const Point t{-v.n.y, v.n.x};
                double best = std::nan(""), bd = 1e300;
                oh.query({at.x + 0.5 * reach * v.n.x, at.y + 0.5 * reach * v.n.y}, 0.5 * reach + 0.25 * v.ht, [&](int i) {
                    if (point_region[cand[i]] != v.region) return;
                    const Point &q = g.P[cand[i]];
                    const double dn = (q.x - at.x) * v.n.x + (q.y - at.y) * v.n.y, dt = (q.x - at.x) * t.x + (q.y - at.y) * t.y;
                    if (std::fabs(dt) > 0.25 * v.ht || dn < 0.5 * o.near * h0 || dn > reach || dn >= bd) return;
                    bd = dn;
                    best = dn - std::floor(dn / h0) * h0;
                });
                return best;
            };
            for (Interface::Vertex v : iface.vertices()) {
                v.adv = advance(v);
                column(v, g.P[v.p], match, phase_of(v, g.P[v.p]));
            }
            for (size_t i = 0; i < cand.size(); ++i)
                if (in_zone[i] && !matched[i]) rep.removed_interior++;
        }
        // 3b. outside the zone: interior points removed near a chain, a zone
        // point or an already kept interior point
        for (size_t i = 0; i < cand.size(); ++i) {
            if (in_zone[i]) continue;
            const Point &x = g.P[cand[i]];
            const int region = point_region[cand[i]];
            const double r = o.near * fine(x, region);
            const bool over = overrun(cand[i], x, region, r);
            bool close = !over && near_chain(x, r);
            if (!over && !close)
                zh.query(x, r, [&](int q) {
                    if (!close && dist(x, zh.at[q]) < r) close = true;
                });
            if (!over && !close)
                pg.query(x.x - r, x.y - r, x.x + r, x.y + r, [&](int k) {
                    if (!close && (size_t)k < i && kept[k] && !in_zone[k] && dist(x, dots[k].first) < r) close = true;
                });
            if (close && needed_for_flow(cand[i])) close = false;
            const bool ok = !over && !close;
            kept[i] = ok;
            if (ok) layer.push_back(cand[i]);
            else rep.removed_interior++;
        }
    }
    il_flag.resize(g.P.size(), 0);
    lap("interior removal");
    // 4a. old triangles kept as they are: all three points kept, no boundary
    // edge changed by steps 1-2, and no edge with a coupling below
    // keep_coupling (those are left to the Delaunay triangulation).  Their
    // edges become CDT constraints, so only the rest is re-triangulated.
    std::set<Triple> keep;
    if (o.keep_triangles) {
        std::vector<char> present(g.P.size(), 0);
        std::unordered_set<uint64_t> chain;
        for (const Canon &c : g.C)
            for (size_t i = 0; i < c.pt.size(); ++i) {
                present[c.pt[i]] = 1;
                if (i + 1 < c.pt.size()) chain.insert(ukey(c.pt[i], c.pt[i + 1]));
            }
        for (int p : layer) present[p] = 1;
        std::vector<char> ok(m.triangle.size(), 1);
        auto cot = [&](int c, int a, int b) {   // cotangent of the angle at c in (c, a, b)
            const Point &C = m.point[c], &A = m.point[a], &B = m.point[b];
            double ux = A.x - C.x, uy = A.y - C.y, vx = B.x - C.x, vy = B.y - C.y;
            double cr = std::fabs(ux * vy - uy * vx);
            return cr > 0 ? (ux * vx + uy * vy) / cr : -1e300;
        };
        for (size_t t = 0; t < m.triangle.size(); ++t) {
            const Triangle &tr = m.triangle[t];
            for (int j = 0; j < 3; ++j) {
                if (!present[tr.v[j]]) ok[t] = 0;
                int a = tr.v[(j + 1) % 3], b = tr.v[(j + 2) % 3], nb = tr.nb[j];
                bool boundary = nb < 0 || m.triangle[nb].region != tr.region;
                if (boundary) {
                    if (!chain.count(ukey(a, b))) ok[t] = 0;
                } else if ((int)t < nb) {
                    int d = -1;
                    for (int v : m.triangle[nb].v)
                        if (v != a && v != b) d = v;
                    if (d >= 0 && 0.5 * (cot(tr.v[j], a, b) + cot(d, a, b)) < o.keep_coupling) ok[t] = ok[nb] = 0;
                }
            }
        }
        for (size_t t = 0; t < m.triangle.size(); ++t)
            if (ok[t]) {
                Triple k{m.triangle[t].v[0], m.triangle[t].v[1], m.triangle[t].v[2]};
                std::sort(k.begin(), k.end());
                keep.insert(k);
            }
    }
    auto keep_edges = [&]() {
        std::vector<std::array<int, 2>> e;
        for (const Triple &k : keep) {
            e.push_back({k[0], k[1]});
            e.push_back({k[1], k[2]});
            e.push_back({k[0], k[2]});
        }
        return e;
    };
    // 4. triangulation; interior edges longer than long_edge x h get their midpoint
    RemeshOptions ro;
    Work w;
    // Steps 4-5 run again when the encroachment repair split an interface
    // of the boundary layer: the new vertex gets its column in this call
    // (else its triangles would fan out to the deep rows for one SUPREM
    // step, and the next call would interpolate the new column from them).
    const size_t chain_points0 = g.P.size();
    std::set<int> columned;   // chain vertices from encroachment splits that have their column
    for (int zone_round = 0;; ++zone_round) {
        for (int round = 0;; ++round) {
            std::vector<char> film(g.P.size(), 0);
            if (keep.empty()) {
                w = triangulate(g, layer, film, ro, rr);
            } else {
                auto ke = keep_edges();
                try {
                    w = triangulate(g, layer, film, ro, rr, &ke);
                } catch (const std::exception &e) {   // e.g. a new chain chord crosses a kept edge
                    rep.notes.push_back(std::string("kept triangles released (") + e.what() + ")");
                    keep.clear();
                    w = triangulate(g, layer, film, ro, rr);
                }
            }
            if (round >= o.rounds) break;
            std::unordered_set<uint64_t> con;
            for (auto &e : constrained_edges(w)) con.insert(ukey(e[0], e[1]));
            std::unordered_map<uint64_t, int> dirw;   // directed edge -> triangle on its left
            if (!keep.empty())
                for (int t = 0; t < (int)w.T.size(); ++t)
                    for (int k = 0; k < 3; ++k) dirw[dkey(w.T[t][k], w.T[t][(k + 1) % 3])] = t;
            int added = 0;
            std::unordered_set<uint64_t> split_edge;
            for (int t = 0; t < (int)w.T.size(); ++t) {
                const Point &A = w.P[w.T[t][0]], &B = w.P[w.T[t][1]], &C = w.P[w.T[t][2]];
                Point cen{(A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3};
                double area = 0.5 * std::fabs(cross(A, B, C)), hc = area_h(cen, w.R[t]);
                int k = 0;   // longest edge v[k] -> v[k+1]
                for (int j = 1; j < 3; ++j)
                    if (dist(w.P[w.T[t][j]], w.P[w.T[t][(j + 1) % 3]]) > dist(w.P[w.T[t][k]], w.P[w.T[t][(k + 1) % 3]])) k = j;
                int a = w.T[t][k], bb = w.T[t][(k + 1) % 3];
                // Flow quality (M1, LOCOS): the viscous/compress flow is a
                // finite-element problem on the oxide and nitride, and needle
                // elements there stall its Newton solve (stress-dependent
                // viscosity).  A coarse oxide on a finely divided boundary
                // (the silicon side sets the interface spacing) can only
                // connect by fans of needles, so a flow-material triangle
                // below flow_min_angle gets the midpoint of its longest
                // interior edge, with a clearance relative to its own size
                // instead of the material's floor: longest-edge bisection
                // with re-triangulation, i.e. Delaunay refinement graded from
                // the boundary to the material's own size.
                const int mat_t = m.region_material[w.R[t]];
                const bool flow_mat = mat_t == SIO2 || mat_t == SINI || mat_t == OXNI;
                bool quality = false;
                double qclear = 0;
                if (flow_mat && o.flow_min_angle > 0 && !con.count(ukey(a, bb))) {
                    const double lim = o.flow_min_angle * PI / 180;
                    double emin = 1e300;
                    for (int j = 0; j < 3; ++j) {
                        if (angle_at(w.P[w.T[t][j]], w.P[w.T[t][(j + 1) % 3]], w.P[w.T[t][(j + 2) % 3]]) < lim) quality = true;
                        emin = std::min(emin, dist(w.P[w.T[t][(j + 1) % 3]], w.P[w.T[t][(j + 2) % 3]]));
                    }
                    qclear = 0.5 * emin;
                    if (quality && 0.5 * dist(w.P[a], w.P[bb]) < std::max(o.min_split_length, o.hmin)) quality = false;
                }
                if (!quality && area <= o.grow * std::sqrt(3.0) / 4 * hc * hc) continue;
                Point x;
                if (!con.count(ukey(a, bb))) {
                    if (!split_edge.insert(ukey(a, bb)).second) continue;   // neighbour already splits it
                    x = {0.5 * (w.P[a].x + w.P[bb].x), 0.5 * (w.P[a].y + w.P[bb].y)};
                } else {
                    x = cen;
                }
                x = snap(x);
                if (Interface::Near z; iface.active() && iface.zone(x, w.R[t], z)) continue;   // the zone is structured
                const double rclear = quality ? std::min(qclear, o.near * fine(x, w.R[t])) : o.near * fine(x, w.R[t]);
                bool clear = true;
                for (int j = 0; j < 3; ++j) clear = clear && dist(x, w.P[w.T[t][j]]) >= rclear;
                clear = clear && !overrun(-1, x, w.R[t], rclear);   // not where a boundary arrives within the step
                if (clear && iface.active()) {   // as the removal rule (3b): not closer than near x fine to a zone point
                    const double r = rclear;
                    zh.query(x, r, [&](int q) {
                        if (clear && dist(x, zh.at[q]) < r) clear = false;
                    });
                }
                if (!clear) continue;
                if (!keep.empty()) {   // the triangles the point goes into are re-triangulated
                    keep.erase(triple_of(w, t));
                    if (!con.count(ukey(a, bb)))
                        if (auto it = dirw.find(dkey(bb, a)); it != dirw.end()) keep.erase(triple_of(w, it->second));
                }
                g.P.push_back(x);
                layer.push_back((int)g.P.size() - 1);
                ++added;
            }
            rep.inserted_interior += added;
            if (!added) break;
        }
        // The kept edges are CDT constraints, so an edge between a kept and a new
        // triangle can be non-Delaunay: a negative box coupling, which costs SUPREM
        // dose (WP7: flat oxidation lost 5e-4 of its boron with them, 5e-5 without).
        // Flip those edges; the kept triangle on such an edge is re-triangulated.
        // Edges between two kept triangles keep their (>= keep_coupling) quality.
        if (!keep.empty()) rep.flips += delaunay_flips(w, &keep);
        lap("triangulation+insertion");
        // 5. encroachment repair (as WP3)
        // Only encroachment the update introduced is repaired: an edge whose
        // apex already faced it in the old mesh keeps the old quality (thin
        // films of a SUPREM grid would otherwise be refined down to their
        // thickness), no split creates a segment that step 1 would remove again,
        // and the number of splits is bounded.
        {
            std::set<std::array<int, 3>> old_faces;   // (min, max, apex) of old triangles
            for (const Triangle &t : m.triangle)
                for (int k = 0; k < 3; ++k) {
                    int a = t.v[(k + 1) % 3], bb = t.v[(k + 2) % 3];
                    old_faces.insert({std::min(a, bb), std::max(a, bb), t.v[k]});
                }
            const int budget = std::max(100, (int)(o.encroach_budget * n0));
            SplitFilter allow = [&](int pa, int pb, int apex, int region) {
                if (pa < n0 && pb < n0 && apex >= 0 && apex < n0 && old_faces.count({std::min(pa, pb), std::max(pa, pb), apex}))
                    return false;
                // interface-layer rows are anisotropic on purpose: their flat
                // triangles must not refine the interface down to the row spacing
                if (apex >= 0 && apex < (int)il_flag.size() && il_flag[apex]) return false;
                Point mid{0.5 * (g.P[pa].x + g.P[pb].x), 0.5 * (g.P[pa].y + g.P[pb].y)};
                return 0.5 * dist(g.P[pa], g.P[pb]) >= o.short_edge * fine(mid, region);   // never what step 1 would remove
            };
            for (int round = 0; round < o.encroach_rounds && rep.encroach_splits < budget; ++round) {
                int n = split_encroached(w, g, o.encroach_angle, std::max(o.min_split_length, o.hmin), allow);
                if (!n) break;
                rep.encroach_splits += n;
                rep.flips += delaunay_flips(w, keep.empty() ? nullptr : &keep);
            }
            if (rep.encroach_splits >= budget) rep.notes.push_back("encroachment split budget reached");
        }
        if (!iface.active() || zone_round >= 3) break;
        std::vector<int> fresh_pts;
        for (const Canon &c : g.C)
            for (int p : c.pt) {
                if ((size_t)p < chain_points0 || (p < (int)il_flag.size() && il_flag[p]) || !columned.insert(p).second) continue;
                Interface::Vertex v;
                if (!iface.vertex_at(g.P[p], v)) continue;
                v.p = p;
                for (int id : column(v, g.P[p], nullptr, std::nan(""))) fresh_pts.push_back(id);
                for (auto it = keep.begin(); it != keep.end();) {   // kept triangles on the split edge
                    const Triple &k = *it;
                    bool on = false;
                    for (int e = 0; e < 3 && !on; ++e) on = seg_dist(g.P[p], m.point[k[e]], m.point[k[(e + 1) % 3]]) < 1e-9;
                    it = on ? keep.erase(it) : std::next(it);
                }
            }
        if (fresh_pts.empty()) break;
        rep.notes.push_back("boundary layer: " + std::to_string(fresh_pts.size()) + " points for interface vertices from encroachment splits");
        for (int id : fresh_pts) {   // old triangles that now hold a point are re-triangulated
            const int t = iface.old_triangle(g.P[id]);
            if (t < 0) continue;
            Triple k{m.triangle[t].v[0], m.triangle[t].v[1], m.triangle[t].v[2]};
            std::sort(k.begin(), k.end());
            keep.erase(k);
        }
    }
    if (auto enc = encroached(w, o.encroach_angle); !enc.empty())
        rep.notes.push_back(std::to_string(enc.size()) + " boundary edges still above the encroachment angle");
    lap("encroachment");
    Mesh out = assemble(m, g, w, rr);
    rep.notes.insert(rep.notes.end(), rr.notes.begin(), rr.notes.end());

    lap("assemble");
    // 6. bookkeeping: origin of every point, points whose neighbourhood changed
    rep.origin.assign(out.point.size(), -1);
    rep.changed.assign(out.point.size(), 0);
    for (size_t i = 0; i < out.point.size(); ++i) {
        int p = w.pid[i];
        if (p >= 0 && p < n0) rep.origin[i] = p; else rep.changed[i] = 1;
    }
    std::set<std::array<int, 4>> old_tri;
    for (const Triangle &t : m.triangle) {
        std::array<int, 4> k{t.v[0], t.v[1], t.v[2], 0};
        std::sort(k.begin(), k.begin() + 3);
        k[3] = t.region;
        old_tri.insert(k);
    }
    for (const Triangle &t : out.triangle) {
        std::array<int, 4> k{rep.origin[t.v[0]], rep.origin[t.v[1]], rep.origin[t.v[2]], 0};
        std::sort(k.begin(), k.begin() + 3);
        k[3] = t.region;
        if (k[0] >= 0 && old_tri.count(k)) continue;
        rep.new_triangles++;
        for (int v : t.v) rep.changed[v] = 1;
    }
    for (const Node &n : out.node)   // a kept point that gained a material
        if (rep.origin[n.point] >= 0 && m.node_at(rep.origin[n.point], n.material) < 0) rep.changed[n.point] = 1;
    rep.kept = 0;
    for (int p : rep.origin) rep.kept += p >= 0;
    lap("bookkeeping");
    return out;
}

}  // namespace suprem_core
