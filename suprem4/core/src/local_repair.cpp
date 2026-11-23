// WP6c: geometric events applied to the mesh locally (see evolve.hpp).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>

#include <CDT.h>

#include "suprem_core/evolve.hpp"
#include "suprem_core/mesh.hpp"
#include "suprem_core/suprem_check.hpp"

namespace suprem_core {

namespace {

double crs(const Point &o, const Point &a, const Point &b) { return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x); }

uint64_t ekey(int a, int b)
{
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

struct Box {
    double x0 = 1e300, y0 = 1e300, x1 = -1e300, y1 = -1e300;
    void add(const Point &p) { x0 = std::min(x0, p.x); y0 = std::min(y0, p.y); x1 = std::max(x1, p.x); y1 = std::max(y1, p.y); }
    bool has(const Point &p) const { return p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1; }
    bool overlaps(const Box &b) const { return b.x0 <= x1 && b.x1 >= x0 && b.y0 <= y1 && b.y1 >= y0; }
    Box grown(double d) const { return {x0 - d, y0 - d, x1 + d, y1 + d}; }
};

int winding(const std::vector<std::vector<Point>> &loops, const Point &p)
{
    int w = 0;
    for (const auto &l : loops)
        for (size_t i = 0; i < l.size(); ++i) {
            const Point &a = l[i], &b = l[(i + 1) % l.size()];
            if (a.y <= p.y) { if (b.y > p.y && crs(a, b, p) > 0) ++w; }
            else if (b.y <= p.y && crs(a, b, p) < 0) --w;
        }
    return w;
}

bool segments_cross(const Point &a, const Point &b, const Point &c, const Point &d)   // proper crossing
{
    double d1 = crs(a, b, c), d2 = crs(a, b, d), d3 = crs(c, d, a), d4 = crs(c, d, b);
    return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}

// Components of a region other than its largest: those smaller than tiny go
// to the neighbour region with the longest shared border (or the ambient).
// Returns the number of triangles moved.
int clean_components(std::vector<Triangle> &tris, const Mesh &m, double tiny, std::vector<char> &touched)
{
    const auto &P = m.point;
    std::unordered_map<uint64_t, std::vector<int>> et;
    for (int t = 0; t < (int)tris.size(); ++t)
        for (int k = 0; k < 3; ++k) et[ekey(tris[t].v[k], tris[t].v[(k + 1) % 3])].push_back(t);
    std::vector<int> comp(tris.size(), -1);
    std::vector<std::vector<int>> comps;
    for (int t0 = 0; t0 < (int)tris.size(); ++t0) {
        if (comp[t0] >= 0) continue;
        comps.push_back({});
        std::vector<int> st{t0};
        comp[t0] = (int)comps.size() - 1;
        while (!st.empty()) {
            int t = st.back(); st.pop_back();
            comps.back().push_back(t);
            for (int k = 0; k < 3; ++k)
                for (int u : et[ekey(tris[t].v[k], tris[t].v[(k + 1) % 3])])
                    if (comp[u] < 0 && tris[u].region == tris[t].region) { comp[u] = comp[t]; st.push_back(u); }
        }
    }
    auto area = [&](int t) { return 0.5 * std::fabs(crs(P[tris[t].v[0]], P[tris[t].v[1]], P[tris[t].v[2]])); };
    std::map<int, int> largest;   // region -> component
    std::vector<double> ca(comps.size(), 0);
    for (size_t c = 0; c < comps.size(); ++c) {
        for (int t : comps[c]) ca[c] += area(t);
        int r = tris[comps[c][0]].region;
        if (!largest.count(r) || ca[c] > ca[largest[r]]) largest[r] = (int)c;
    }
    std::vector<char> erase(tris.size(), 0);
    int moved = 0;
    for (size_t c = 0; c < comps.size(); ++c) {
        int r = tris[comps[c][0]].region;
        if (largest[r] == (int)c || ca[c] >= tiny) continue;
        std::map<int, int> votes;
        for (int t : comps[c])
            for (int k = 0; k < 3; ++k)
                for (int u : et[ekey(tris[t].v[k], tris[t].v[(k + 1) % 3])])
                    if (tris[u].region != r) votes[tris[u].region]++;
        int to = -1, best = 0;
        for (auto &[q, n] : votes)
            if (n > best) { best = n; to = q; }
        for (int t : comps[c]) {
            for (int v : tris[t].v) touched[v] = 1;
            if (to >= 0) tris[t].region = to; else erase[t] = 1;
            ++moved;
        }
    }
    std::vector<Triangle> keep;
    for (size_t t = 0; t < tris.size(); ++t)
        if (!erase[t]) keep.push_back(tris[t]);
    tris.swap(keep);
    return moved;
}

// A region whose triangles around a vertex form more than one fan touches
// itself there, which SUPREM does not accept.  Resolved locally: the smaller
// fan goes to its neighbour material when it is tiny, otherwise the ambient
// wedge at the vertex is closed by one triangle.  Returns false when a pinch
// cannot be resolved.
bool resolve_pinches(std::vector<Triangle> &tris, const Mesh &m, double tiny, std::vector<char> &touched, LocalRepairReport &rep)
{
    const auto &P = m.point;
    for (int iter = 0; iter < 100; ++iter) {
        std::unordered_map<int, std::vector<int>> star;
        std::unordered_map<uint64_t, std::vector<int>> et;
        for (int t = 0; t < (int)tris.size(); ++t)
            for (int k = 0; k < 3; ++k) {
                star[tris[t].v[k]].push_back(t);
                et[ekey(tris[t].v[k], tris[t].v[(k + 1) % 3])].push_back(t);
            }
        auto area = [&](int t) { return 0.5 * std::fabs(crs(P[tris[t].v[0]], P[tris[t].v[1]], P[tris[t].v[2]])); };
        int pv = -1, pr = -1;
        std::vector<std::vector<int>> comps;
        for (auto &[v, st] : star) {
            std::map<int, std::vector<int>> byreg;
            for (int t : st) byreg[tris[t].region].push_back(t);
            for (auto &[r, ts] : byreg) {
                if (ts.size() < 2) continue;
                std::map<int, int> par;
                for (int t : ts) par[t] = t;
                std::function<int(int)> find = [&](int x) { return par[x] == x ? x : par[x] = find(par[x]); };
                for (int t : ts)
                    for (int k = 0; k < 3; ++k) {
                        int a = tris[t].v[k], b = tris[t].v[(k + 1) % 3];
                        if (a != v && b != v) continue;
                        for (int u : et[ekey(a, b)])
                            if (u != t && par.count(u)) par[find(u)] = find(t);
                    }
                std::map<int, std::vector<int>> cc;
                for (int t : ts) cc[find(t)].push_back(t);
                if (cc.size() > 1) {
                    pv = v; pr = r;
                    for (auto &kv : cc) comps.push_back(kv.second);
                    break;
                }
            }
            if (pv >= 0) break;
        }
        if (pv < 0) return true;
        std::sort(comps.begin(), comps.end(), [&](const std::vector<int> &a, const std::vector<int> &b) {
            double sa = 0, sb = 0;
            for (int t : a) sa += area(t);
            for (int t : b) sb += area(t);
            return sa < sb;
        });
        double small = 0;
        for (int t : comps[0]) small += area(t);
        if (small < tiny) {   // a tiny fan: to the neighbour material (or the ambient)
            std::map<int, int> votes;
            for (int t : comps[0])
                for (int k = 0; k < 3; ++k)
                    for (int u : et[ekey(tris[t].v[k], tris[t].v[(k + 1) % 3])])
                        if (u != t && tris[u].region != pr) votes[tris[u].region]++;
            int to = -1, best = 0;
            for (auto &[r, n] : votes)
                if (n > best) { best = n; to = r; }
            std::vector<int> sel = comps[0];
            std::sort(sel.rbegin(), sel.rend());
            for (int t : sel) {
                for (int v : tris[t].v) touched[v] = 1;
                if (to >= 0) tris[t].region = to;
                else tris.erase(tris.begin() + t);
            }
            rep.notes.push_back("pinch: tiny fan given away");
            continue;
        }
        // close the ambient wedge: a triangle (v, a, b) between two boundary edges at v
        std::vector<int> bn;
        for (int t : star[pv])
            for (int k = 0; k < 3; ++k) {
                int a = tris[t].v[k], b = tris[t].v[(k + 1) % 3];
                if ((a == pv || b == pv) && et[ekey(a, b)].size() == 1) bn.push_back(a == pv ? b : a);
            }
        std::set<int> near_t;
        for (int t : star[pv]) near_t.insert(t);
        for (int x : bn)
            for (int t : star[x]) near_t.insert(t);
        int ba = -1, bb = -1;
        double barea = 1e300;
        for (int a : bn)
            for (int b : bn) {
                if (a == b || crs(P[pv], P[a], P[b]) <= 0) continue;
                const Point &V = P[pv], &A = P[a], &B = P[b];
                double ar = 0.5 * crs(V, A, B);
                if (ar >= barea) continue;
                Point g{(V.x + A.x + B.x) / 3, (V.y + A.y + B.y) / 3};
                bool ok = true;
                for (int t : near_t) {
                    const Point &p0 = P[tris[t].v[0]], &p1 = P[tris[t].v[1]], &p2 = P[tris[t].v[2]];
                    if (crs(p0, p1, g) > 0 && crs(p1, p2, g) > 0 && crs(p2, p0, g) > 0) { ok = false; break; }
                    for (int k = 0; k < 3 && ok; ++k) {
                        int x = tris[t].v[k], y = tris[t].v[(k + 1) % 3];
                        if (segments_cross(A, B, P[x], P[y]) || segments_cross(V, A, P[x], P[y]) || segments_cross(V, B, P[x], P[y]))
                            ok = false;
                        if (x != pv && x != a && x != b && crs(V, A, P[x]) > 0 && crs(A, B, P[x]) > 0 && crs(B, V, P[x]) > 0) ok = false;
                    }
                    if (!ok) break;
                }
                if (ok) { barea = ar; ba = a; bb = b; }
            }
        if (ba < 0) {
            rep.notes.push_back("pinch at (" + std::to_string(P[pv].x) + ", " + std::to_string(P[pv].y) + ") cannot be closed");
            return false;
        }
        Triangle tr;
        tr.v[0] = pv; tr.v[1] = ba; tr.v[2] = bb;
        tr.region = pr;
        tris.push_back(tr);
        touched[pv] = touched[ba] = touched[bb] = 1;
        rep.filled++;
        rep.notes.push_back("pinch: ambient wedge closed");
    }
    rep.notes.push_back("pinches not resolved in 100 rounds");
    return false;
}

// Boundary edges are split at the given points that lie on them (within
// tol, more than end_tol from both ends), on SUPREM's float cm grid: the
// mouths of a gap piece, where the opening boundary meets the faces.
// Returns the number of points added (appended to m.point).
int split_boundary_at(Mesh &m, const std::vector<Point> &qs, double tol, double end_tol)
{
    int added = 0;
    for (Point q : qs) {
        q.x = (double)(float)(1e-4 * q.x) * 1e4;
        q.y = (double)(float)(1e-4 * q.y) * 1e4;
        int bt = -1, bj = -1;
        double bd = tol, bs = 0;
        for (int t = 0; t < (int)m.triangle.size(); ++t)
            for (int j = 0; j < 3; ++j) {
                if (m.triangle[t].nb[j] >= 0) continue;
                const Point &a = m.point[m.triangle[t].v[(j + 1) % 3]], &b = m.point[m.triangle[t].v[(j + 2) % 3]];
                double ex = b.x - a.x, ey = b.y - a.y, l2 = ex * ex + ey * ey;
                if (l2 <= 0) continue;
                double sp = ((q.x - a.x) * ex + (q.y - a.y) * ey) / l2;
                if (sp <= 0 || sp >= 1) continue;
                double d = std::hypot(a.x + sp * ex - q.x, a.y + sp * ey - q.y), l = std::sqrt(l2);
                if (d < bd && sp * l > end_tol && (1 - sp) * l > end_tol) { bd = d; bt = t; bj = j; bs = sp; }
            }
        if (bt < 0) continue;
        const Triangle t0 = m.triangle[bt];
        int a = t0.v[(bj + 1) % 3], b = t0.v[(bj + 2) % 3], code = t0.nb[bj];
        {   // on the float grid the point must not cross to the other side of anything (sub-angstrom gaps)
            int c = t0.v[bj];
            bool bad = crs(m.point[c], m.point[a], q) <= 0 || crs(m.point[c], q, m.point[b]) <= 0;
            for (int t = 0; t < (int)m.triangle.size() && !bad; ++t)
                for (int j = 0; j < 3 && !bad; ++j) {
                    int x = m.triangle[t].v[(j + 1) % 3], y = m.triangle[t].v[(j + 2) % 3];
                    bad = segments_cross(m.point[a], q, m.point[x], m.point[y]) || segments_cross(q, m.point[b], m.point[x], m.point[y]) ||
                          segments_cross(m.point[c], q, m.point[x], m.point[y]);
                }
            if (bad) continue;
        }
        int p = (int)m.point.size();
        m.point.push_back(q);
        Triangle t1 = t0, t2 = t0;
        t1.v[(bj + 2) % 3] = p;   // (c, a, p)
        t2.v[(bj + 1) % 3] = p;   // (c, p, b)
        m.triangle[bt] = t1;
        m.triangle.push_back(t2);
        // neighbours: rebuilt from shared faces, the two halves keep the parent's code
        std::map<uint64_t, int> codes;
        for (const Triangle &t : m.triangle)
            for (int j = 0; j < 3; ++j)
                if (t.nb[j] < 0) codes[ekey(t.v[(j + 1) % 3], t.v[(j + 2) % 3])] = t.nb[j];
        codes[ekey(a, p)] = code;
        codes[ekey(p, b)] = code;
        std::unordered_map<uint64_t, int> dir;
        for (int t = 0; t < (int)m.triangle.size(); ++t)
            for (int k = 0; k < 3; ++k) dir[((uint64_t)(uint32_t)m.triangle[t].v[k] << 32) | (uint32_t)m.triangle[t].v[(k + 1) % 3]] = t;
        for (Triangle &t : m.triangle)
            for (int j = 0; j < 3; ++j) {
                int x = t.v[(j + 1) % 3], y = t.v[(j + 2) % 3];
                auto it = dir.find(((uint64_t)(uint32_t)y << 32) | (uint32_t)x);
                if (it != dir.end()) { t.nb[j] = it->second; continue; }
                auto c = codes.find(ekey(x, y));
                t.nb[j] = c == codes.end() ? BC_OFFSET : c->second;
            }
        // nodes of the new point: interpolated along the edge
        std::vector<int> mats{m.region_material[t0.region]};
        if (code == EXPOSED || code == BACKSIDE) mats.push_back(GAS);
        for (int mat : mats) {
            int na = m.node_at(a, mat), nb = m.node_at(b, mat);
            std::vector<double> v(m.species.size(), 0.0);
            if (na >= 0 && nb >= 0)
                for (size_t k = 0; k < v.size(); ++k) v[k] = (1 - bs) * m.node[na].value[k] + bs * m.node[nb].value[k];
            m.node.push_back({p, mat, v});
        }
        m.rebuild_index();
        ++added;
    }
    return added;
}

int orient_f(const Mesh &m, int a, int b, int c)   // SUPREM's float predicate
{
    float A[2] = {suprem_float(m.point[a].x), suprem_float(m.point[a].y)};
    float B[2] = {suprem_float(m.point[b].x), suprem_float(m.point[b].y)};
    float C[2] = {suprem_float(m.point[c].x), suprem_float(m.point[c].y)};
    return orient2d_float_sign(A, B, C);
}

void relink(Mesh &m, const std::map<uint64_t, int> &codes)   // neighbours from shared faces; boundary codes by point pair
{
    std::unordered_map<uint64_t, int> dir;
    for (int t = 0; t < (int)m.triangle.size(); ++t)
        for (int k = 0; k < 3; ++k) dir[((uint64_t)(uint32_t)m.triangle[t].v[k] << 32) | (uint32_t)m.triangle[t].v[(k + 1) % 3]] = t;
    for (Triangle &t : m.triangle)
        for (int j = 0; j < 3; ++j) {
            int x = t.v[(j + 1) % 3], y = t.v[(j + 2) % 3];
            auto it = dir.find(((uint64_t)(uint32_t)y << 32) | (uint32_t)x);
            if (it != dir.end()) { t.nb[j] = it->second; continue; }
            auto c = codes.find(ekey(x, y));
            t.nb[j] = c == codes.end() ? EXPOSED : c->second;
        }
}

// Contacts: surfaces closer than eps (a gap the motion has closed below the
// shapes' resolution; SUPREM's collision limit leaves pairs closer than
// delta/2 to the core).  Zipper: an exposed vertex within eps of a facing
// exposed edge (of another part, or of its own surface folded back) becomes a vertex of that edge (its triangle
// is split there; no point moves), and an ambient wedge between two exposed
// edges at a vertex whose far ends are closer than eps is closed by a
// triangle.  Nodes of new (point, material) pairs are copies; touched marks
// the vertices involved.  Returns the number of operations.
int close_contacts(Mesh &m, double eps, std::vector<char> &touched)
{
    int ops = 0;
    for (int iter = 0; iter < 500; ++iter) {
        std::map<uint64_t, int> codes;
        std::unordered_map<uint64_t, std::vector<std::pair<int, int>>> et;   // edge -> (triangle, face)
        for (int t = 0; t < (int)m.triangle.size(); ++t)
            for (int j = 0; j < 3; ++j) {
                const Triangle &tr = m.triangle[t];
                uint64_t k = ekey(tr.v[(j + 1) % 3], tr.v[(j + 2) % 3]);
                et[k].push_back({t, j});
                if (tr.nb[j] < 0) codes[k] = tr.nb[j];
            }
        struct Ex { int a, b, t, c; };
        std::vector<Ex> ex;
        std::unordered_map<int, std::vector<int>> at;   // exposed vertex -> exposed edges
        for (auto &[k, v] : et) {
            if (v.size() != 1) continue;
            const Triangle &tr = m.triangle[v[0].first];
            int j = v[0].second;
            if (tr.nb[j] != EXPOSED) continue;
            ex.push_back({tr.v[(j + 1) % 3], tr.v[(j + 2) % 3], v[0].first, tr.v[j]});
            at[ex.back().a].push_back((int)ex.size() - 1);
            at[ex.back().b].push_back((int)ex.size() - 1);
        }
        // spatial hash of exposed edges
        const double cell = std::max(eps * 64, 1e-6);
        std::unordered_map<uint64_t, std::vector<int>> grid;
        auto ck = [&](double x, double y) {
            return ((uint64_t)(uint32_t)(int64_t)std::floor(x / cell) << 32) | (uint32_t)(int64_t)std::floor(y / cell);
        };
        for (int e = 0; e < (int)ex.size(); ++e) {
            const Point &a = m.point[ex[e].a], &b = m.point[ex[e].b];
            for (double x = std::floor((std::min(a.x, b.x) - eps) / cell); x <= std::floor((std::max(a.x, b.x) + eps) / cell); ++x)
                for (double y = std::floor((std::min(a.y, b.y) - eps) / cell); y <= std::floor((std::max(a.y, b.y) + eps) / cell); ++y)
                    grid[((uint64_t)(uint32_t)(int64_t)x << 32) | (uint32_t)(int64_t)y].push_back(e);
        }
        // best zip: vertex p, edge e (interior projection, p on the ambient side, not adjacent)
        double bd = eps;
        int bp = -1, be = -1;
        for (auto &[p, es] : at) {
            const Point &P = m.point[p];
            auto it = grid.find(ck(P.x, P.y));
            if (it == grid.end()) continue;
            for (int e : it->second) {
                const Ex &E = ex[e];
                if (E.a == p || E.b == p) continue;   // (p joined to an end of e: a fold, zipped as well)
                const Point &A = m.point[E.a], &B = m.point[E.b], &C = m.point[E.c];
                double vx = B.x - A.x, vy = B.y - A.y, l2 = vx * vx + vy * vy;
                double u = ((P.x - A.x) * vx + (P.y - A.y) * vy) / l2;
                if (u <= 0 || u >= 1) continue;
                double d = std::hypot(A.x + u * vx - P.x, A.y + u * vy - P.y);
                if (d >= bd || crs(A, B, P) * crs(A, B, C) > 0) continue;   // not facing: p on the material side
                if (orient_f(m, E.a, p, E.c) <= 0 || orient_f(m, p, E.b, E.c) <= 0) continue;
                bd = d; bp = p; be = e;
            }
        }
        if (bp >= 0) {
            const Ex E = ex[be];
            Triangle t0 = m.triangle[E.t], t1 = t0, t2 = t0;
            int j = -1;
            for (int k = 0; k < 3; ++k)
                if (t0.v[k] == E.c) j = k;
            t1.v[(j + 2) % 3] = bp;   // (c, a, p)
            t2.v[(j + 1) % 3] = bp;   // (c, p, b)
            m.triangle[E.t] = t1;
            m.triangle.push_back(t2);
            codes.erase(ekey(E.a, E.b));
            relink(m, codes);
            int mat = m.region_material[t0.region];
            if (m.node_at(bp, mat) < 0) {
                int on = -1;
                for (int mm = 1; mm < 16 && on < 0; ++mm) on = m.node_at(bp, mm);
                m.node.push_back({bp, mat, on >= 0 ? m.node[on].value : std::vector<double>(m.species.size(), 0.0)});
            }
            touched[bp] = touched[E.a] = touched[E.b] = touched[E.c] = 1;
            ++ops;
        } else {
            // wedge: exposed edges (v, p) and (v, q) of different triangles, |pq| < eps, (v, p, q) in the ambient
            int wv = -1, wp = -1, wq = -1, wt = -1;
            double wd = eps;
            for (auto &[v, es] : at)
                for (int e1 : es)
                    for (int e2 : es) {
                        if (e1 == e2) continue;
                        // e1 = (v, p) with the material on its left (a -> b with c on the left): p = e1.b when v = e1.a
                        const Ex &X = ex[e1], &Y = ex[e2];
                        if (X.a != v || Y.b != v) continue;   // ambient wedge between X (v -> p) and Y (q -> v)
                        int p = X.b, q = Y.a;
                        if (p == q) continue;
                        double d = std::hypot(m.point[p].x - m.point[q].x, m.point[p].y - m.point[q].y);
                        if (d >= wd) continue;
                        if (orient_f(m, v, q, p) <= 0) continue;   // (v, q, p) counter-clockwise, in the ambient wedge
                        wd = d; wv = v; wp = p; wq = q; wt = X.t;
                    }
            if (wv < 0) break;
            Triangle tr;
            tr.v[0] = wv; tr.v[1] = wq; tr.v[2] = wp;
            tr.region = m.triangle[wt].region;
            m.triangle.push_back(tr);
            relink(m, codes);
            int mat = m.region_material[tr.region];
            for (int v : tr.v)
                if (m.node_at(v, mat) < 0) {
                    int on = -1;
                    for (int mm = 1; mm < 16 && on < 0; ++mm) on = m.node_at(v, mm);
                    m.node.push_back({v, mat, on >= 0 ? m.node[on].value : std::vector<double>(m.species.size(), 0.0)});
                }
            touched[wv] = touched[wp] = touched[wq] = 1;
            ++ops;
        }
        m.rebuild_index();
    }
    return ops;
}

}  // namespace

bool local_event_repair(Mesh &m_in, const std::vector<RepairEvent> &ev, LocalRepairReport &rep, double delta)
{
    Mesh m = m_in;
    const int n_in = (int)m_in.point.size();
    {   // gap mouths: points where the pieces meet the faces become boundary vertices
        std::vector<Point> qs;
        for (const RepairEvent &e : ev)
            if (e.material == GAS && e.kind == "gap")
                for (auto &l : e.piece) qs.insert(qs.end(), l.begin(), l.end());
        rep.split = split_boundary_at(m, qs, 3e-6, 1e-5);
    }
    std::vector<char> pre_touched(m.point.size(), 0);
    rep.contacts = close_contacts(m, delta / 2, pre_touched);
    const int nt = (int)m.triangle.size();
    std::vector<Point> cen(nt);
    for (int t = 0; t < nt; ++t) {
        const Triangle &tr = m.triangle[t];
        cen[t] = {(m.point[tr.v[0]].x + m.point[tr.v[1]].x + m.point[tr.v[2]].x) / 3,
                  (m.point[tr.v[0]].y + m.point[tr.v[1]].y + m.point[tr.v[2]].y) / 3};
    }
    std::unordered_map<uint64_t, std::vector<int>> edge_tri;
    for (int t = 0; t < nt; ++t)
        for (int j = 0; j < 3; ++j) edge_tri[ekey(m.triangle[t].v[(j + 1) % 3], m.triangle[t].v[(j + 2) % 3])].push_back(t);
    auto region_of_material = [&](int mat, const std::vector<int> &near) {   // a region of mat next to the triangles `near`
        for (int t : near)
            for (int j = 0; j < 3; ++j) {
                const Triangle &tr = m.triangle[t];
                for (int u : edge_tri[ekey(tr.v[(j + 1) % 3], tr.v[(j + 2) % 3])])
                    if (u != t && m.material_of(u) == mat) return m.triangle[u].region;
            }
        return -1;
    };

    std::vector<int> new_region(nt, -2);   // -2 unchanged, -1 removed, else region
    {   // needles: triangles with two or more edges on their region's border and a smallest altitude below
        // delta/8 (cusps and folds left by the motion, below the shapes' resolution) go to the material
        // across those edges, thinnest first
        std::vector<int> reg(nt);
        for (int t = 0; t < nt; ++t) reg[t] = m.triangle[t].region;
        std::vector<double> alt(nt);
        for (int t = 0; t < nt; ++t) {
            const Triangle &tr = m.triangle[t];
            double lmax = 0;
            for (int k = 0; k < 3; ++k) {
                const Point &a = m.point[tr.v[k]], &b = m.point[tr.v[(k + 1) % 3]];
                lmax = std::max(lmax, std::hypot(b.x - a.x, b.y - a.y));
            }
            alt[t] = lmax > 0 ? std::fabs(crs(m.point[tr.v[0]], m.point[tr.v[1]], m.point[tr.v[2]])) / lmax : 0;
        }
        const double thin = delta / 8;
        for (int round = 0; round < 1000; ++round) {
            int best = -1, to = -2;
            for (int t = 0; t < nt; ++t) {
                if (reg[t] < 0 || alt[t] >= thin || (best >= 0 && alt[t] >= alt[best])) continue;
                // borders by material (two regions of one material count as one side)
                auto mat_of = [&](int region) { return region < 0 ? -1 : m.region_material[region]; };
                const int own = mat_of(reg[t]);
                std::map<int, int> across_n;                       // material (-1: ambient / domain boundary) -> edges
                std::map<int, std::pair<double, int>> longest;     // material -> (longest shared edge, its region)
                int border = 0;
                for (int j = 0; j < 3; ++j) {
                    const Triangle &tr = m.triangle[t];
                    int a = tr.v[(j + 1) % 3], b = tr.v[(j + 2) % 3];
                    int other = -1;
                    for (int u : edge_tri[ekey(a, b)])
                        if (u != t && reg[u] >= 0) other = u;
                    int r = other >= 0 ? reg[other] : -1, q = mat_of(r);
                    if (q == own) continue;
                    ++border;
                    across_n[q]++;
                    double l = std::hypot(m.point[b].x - m.point[a].x, m.point[b].y - m.point[a].y);
                    if (l > longest[q].first) longest[q] = {l, r};
                }
                if (border < 2) continue;
                // the move must reduce the triangle's border edges (more edges with the target material
                // than with its own): every move shortens the total border, so the loop terminates
                int r = -2;
                double len = -1;
                for (auto &[q, n] : across_n)
                    if (n > 3 - border && longest[q].first > len) { len = longest[q].first; r = longest[q].second; }
                if (r == -2) continue;
                if (r == -1 && m.triangle[t].nb[0] != EXPOSED && m.triangle[t].nb[1] != EXPOSED && m.triangle[t].nb[2] != EXPOSED)
                    continue;   // never remove material at a coded (non-ambient) boundary
                best = t;
                to = r;
            }
            if (best < 0) break;
            reg[best] = to;
            new_region[best] = to;
            ++rep.needles;
        }
    }
    std::vector<std::array<int, 3>> added;
    std::vector<int> added_region;
    for (const RepairEvent &e : ev) {
        Box pb;
        for (auto &l : e.piece)
            for (auto &p : l) pb.add(p);
        if (e.material != GAS) {   // thin / crumb: relabel or remove
            std::vector<int> sel;
            for (int t = 0; t < nt; ++t)
                if (m.material_of(t) == e.material && new_region[t] == -2 && pb.has(cen[t]) && winding(e.piece, cen[t]) != 0)
                    sel.push_back(t);
            if (sel.empty()) {
                rep.notes.push_back(e.kind + " piece contains no triangle centroid (below the mesh scale)");
                continue;
            }
            int r = -1;
            if (e.to != GAS) {
                r = region_of_material(e.to, sel);
                if (r < 0) {   // grow the selection's neighbourhood once: through shared vertices
                    std::set<int> vs;
                    for (int t : sel)
                        for (int v : m.triangle[t].v) vs.insert(v);
                    for (int t = 0; t < nt && r < 0; ++t)
                        if (m.material_of(t) == e.to)
                            for (int v : m.triangle[t].v)
                                if (vs.count(v)) { r = m.triangle[t].region; break; }
                }
                if (r < 0) {
                    rep.notes.push_back(e.kind + " piece: no region of material " + std::to_string(e.to) + " next to it");
                    return false;
                }
            }
            for (int t : sel) new_region[t] = r;
            (r < 0 ? rep.removed : rep.relabeled) += (int)sel.size();
        } else {   // gap / cavity: fill from a local constrained triangulation of the boundary points
            const double size = std::max(pb.x1 - pb.x0, pb.y1 - pb.y0);
            Box lb = pb.grown(0.5 * size + 1e-9);
            std::vector<std::array<int, 2>> bedges;   // boundary edges near the piece
            for (int t = 0; t < nt; ++t)
                for (int j = 0; j < 3; ++j) {
                    if (m.triangle[t].nb[j] >= 0) continue;
                    int a = m.triangle[t].v[(j + 1) % 3], b = m.triangle[t].v[(j + 2) % 3];
                    Box eb;
                    eb.add(m.point[a]);
                    eb.add(m.point[b]);
                    if (eb.overlaps(lb)) bedges.push_back({a, b});
                }
            std::map<int, int> local;
            std::vector<int> glob;
            std::vector<CDT::V2d<double>> V;
            auto use = [&](int p) {
                auto it = local.find(p);
                if (it != local.end()) return it->second;
                int k = (int)glob.size();
                local[p] = k;
                glob.push_back(p);
                V.push_back(CDT::V2d<double>(m.point[p].x, m.point[p].y));
                return k;
            };
            std::vector<CDT::Edge> E;
            for (auto &be : bedges) E.push_back(CDT::Edge(use(be[0]), use(be[1])));
            if (V.size() < 3) {   // e.g. a hole of the union of overlapping triangles: not one of the mesh
                rep.notes.push_back(e.kind + " piece has no mesh boundary around it");
                continue;
            }
            std::vector<std::array<int, 3>> T;
            try {
                CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::NotAllowed, 0.0);
                cdt.insertVertices(V);
                cdt.insertEdges(E);
                cdt.eraseSuperTriangle();
                if (cdt.vertices.size() != V.size()) throw std::runtime_error("duplicate boundary points");
                for (const auto &t : cdt.triangles) T.push_back({glob[t.vertices[0]], glob[t.vertices[1]], glob[t.vertices[2]]});
            } catch (const std::exception &x) {
                rep.notes.push_back(e.kind + " piece: local triangulation failed (" + x.what() + ")");
                return false;
            }
            std::vector<int> near;   // mesh triangles around the piece (overlap test)
            for (int t = 0; t < nt; ++t) {
                Box tb;
                for (int v : m.triangle[t].v) tb.add(m.point[v]);
                if (tb.overlaps(lb)) near.push_back(t);
            }
            std::vector<std::array<int, 3>> fill;
            int n_out = 0, n_in_mesh = 0, n_cross = 0;
            for (auto v : T) {
                const Point &A = m.point[v[0]], &B = m.point[v[1]], &C = m.point[v[2]];
                if (crs(A, B, C) < 0) std::swap(v[1], v[2]);
                const Point &a = m.point[v[0]], &b = m.point[v[1]], &c = m.point[v[2]];
                if (crs(a, b, c) <= 0) continue;
                Point g{(a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3};
                if (winding(e.piece, g) == 0) { ++n_out; continue; }
                bool ok = true;   // must lie in the ambient: no mesh triangle contains its centroid, no edge crosses
                for (int t : near) {
                    const Triangle &tr = m.triangle[t];
                    const Point &p0 = m.point[tr.v[0]], &p1 = m.point[tr.v[1]], &p2 = m.point[tr.v[2]];
                    double s = crs(p0, p1, p2) > 0 ? 1 : -1;
                    if (s * crs(p0, p1, g) > 0 && s * crs(p1, p2, g) > 0 && s * crs(p2, p0, g) > 0) { ok = false; ++n_in_mesh; break; }
                }
                for (int k = 0; k < 3 && ok; ++k)
                    for (auto &be : bedges)
                        if (segments_cross(m.point[v[k]], m.point[v[(k + 1) % 3]], m.point[be[0]], m.point[be[1]])) { ok = false; ++n_cross; break; }
                if (ok) fill.push_back(v);
            }
            if (fill.empty()) {
                rep.notes.push_back(e.kind + " piece: no fill triangle (" + std::to_string(T.size()) + " local, " + std::to_string(n_out) +
                                    " outside the piece, " + std::to_string(n_in_mesh) + " in the mesh, " + std::to_string(n_cross) +
                                    " crossing; " + std::to_string(V.size()) + " points, " + std::to_string(E.size()) + " edges)");
                continue;
            }
            int r = -1;
            for (auto &v : fill)
                for (int k = 0; k < 3 && r < 0; ++k)
                    for (int u : edge_tri[ekey(v[k], v[(k + 1) % 3])])
                        if (m.material_of(u) == e.to && new_region[u] == -2) { r = m.triangle[u].region; break; }
            if (r < 0) {
                rep.notes.push_back(e.kind + " piece: no region of material " + std::to_string(e.to) + " next to the fill");
                return false;
            }
            for (auto &v : fill) { added.push_back(v); added_region.push_back(r); }
            rep.filled += (int)fill.size();
        }
    }
    if (!rep.relabeled && !rep.removed && !rep.filled && !rep.needles && !rep.contacts) {   // nothing to do
        rep.split = 0;
        rep.origin.resize(n_in);
        for (int i = 0; i < n_in; ++i) rep.origin[i] = i;
        rep.touched.assign(n_in, 0);
        return true;
    }

    // assemble: kept / relabelled triangles and the fill, unused points dropped
    std::map<uint64_t, int> code;   // old boundary codes by point pair
    for (int t = 0; t < nt; ++t)
        for (int j = 0; j < 3; ++j)
            if (m.triangle[t].nb[j] < 0) code[ekey(m.triangle[t].v[(j + 1) % 3], m.triangle[t].v[(j + 2) % 3])] = m.triangle[t].nb[j];
    std::vector<Triangle> tris;
    std::vector<char> changed_tri;
    std::vector<char> touched_old(m.point.size(), 0);
    for (size_t p = n_in; p < m.point.size(); ++p) touched_old[p] = 1;   // split points
    for (size_t p = 0; p < pre_touched.size(); ++p) touched_old[p] |= pre_touched[p];
    for (int t = 0; t < nt; ++t) {
        if (new_region[t] == -1) {
            for (int v : m.triangle[t].v) touched_old[v] = 1;
            continue;
        }
        Triangle tr = m.triangle[t];
        if (new_region[t] >= 0) {
            tr.region = new_region[t];
            for (int v : tr.v) touched_old[v] = 1;
        }
        tris.push_back(tr);
    }
    for (size_t i = 0; i < added.size(); ++i) {
        Triangle tr;
        for (int k = 0; k < 3; ++k) { tr.v[k] = added[i][k]; touched_old[added[i][k]] = 1; }
        tr.region = added_region[i];
        tris.push_back(tr);
    }
    for (int round = 0; round < 3; ++round) {
        bool moved = clean_components(tris, m, delta * delta, touched_old) > 0;
        if (!resolve_pinches(tris, m, delta * delta, touched_old, rep)) return false;
        if (!moved) break;
    }
    std::vector<int> renum(m.point.size(), -1);
    Mesh out = m;
    out.point.clear();
    rep.origin.clear();
    for (Triangle &tr : tris)
        for (int &v : tr.v) {
            if (renum[v] < 0) { renum[v] = (int)out.point.size(); out.point.push_back(m.point[v]); rep.origin.push_back(v); }
            v = renum[v];
        }
    rep.dropped_points = (int)(m.point.size() - out.point.size());
    rep.touched.assign(out.point.size(), 0);
    for (size_t i = 0; i < out.point.size(); ++i) rep.touched[i] = touched_old[rep.origin[i]];
    // neighbours: shared faces; old boundary codes kept, new boundary faces the ambient
    std::unordered_map<uint64_t, int> dir;
    for (int t = 0; t < (int)tris.size(); ++t)
        for (int k = 0; k < 3; ++k) dir[((uint64_t)(uint32_t)tris[t].v[k] << 32) | (uint32_t)tris[t].v[(k + 1) % 3]] = t;
    for (Triangle &tr : tris)
        for (int j = 0; j < 3; ++j) {
            int a = tr.v[(j + 1) % 3], b = tr.v[(j + 2) % 3];
            auto it = dir.find(((uint64_t)(uint32_t)b << 32) | (uint32_t)a);
            if (it != dir.end()) { tr.nb[j] = it->second; continue; }
            auto c = code.find(ekey(rep.origin[a], rep.origin[b]));
            tr.nb[j] = c == code.end() ? EXPOSED : c->second;
        }
    out.triangle = tris;
    // nodes: one per (point, material), gas nodes on exposed / backside points
    std::set<std::pair<int, int>> want;
    for (const Triangle &t : out.triangle)
        for (int k = 0; k < 3; ++k) {
            want.insert({t.v[k], out.region_material[t.region]});
            if (t.nb[k] == EXPOSED || t.nb[k] == BACKSIDE) {
                want.insert({t.v[(k + 1) % 3], GAS});
                want.insert({t.v[(k + 2) % 3], GAS});
            }
        }
    out.node.clear();
    for (auto &[p, mat] : want) {
        int on = m.node_at(rep.origin[p], mat);
        if (on < 0)   // a material new at this point: copy any of its nodes
            for (int mm = 0; mm < 16 && on < 0; ++mm)
                if (mm != GAS) on = m.node_at(rep.origin[p], mm);
        out.node.push_back({p, mat, on >= 0 ? m.node[on].value : std::vector<double>(m.species.size(), 0.0)});
    }
    out.rebuild_index();
    for (int &o : rep.origin)
        if (o >= n_in) o = -1;
    m_in = std::move(out);
    return true;
}

bool apply_events_locally(Mesh &m, const RepairOptions &o, LocalRepairReport &rep, std::vector<RepairEvent> *events,
                          int rounds)
{
    rep = LocalRepairReport();
    rep.origin.resize(m.point.size());
    for (size_t i = 0; i < m.point.size(); ++i) rep.origin[i] = (int)i;
    rep.touched.assign(m.point.size(), 0);
    for (int r = 0; r < rounds; ++r) {
        std::vector<Shape> shapes = mesh_shapes_union(m);
        std::vector<RepairEvent> ev = repair_shapes(shapes, o);
        if (events) events->insert(events->end(), ev.begin(), ev.end());
        Mesh trial = m;
        LocalRepairReport one;
        if (!local_event_repair(trial, ev, one, o.delta)) {
            rep.notes.insert(rep.notes.end(), one.notes.begin(), one.notes.end());
            return false;
        }
        rep.relabeled += one.relabeled; rep.removed += one.removed; rep.filled += one.filled;
        rep.dropped_points += one.dropped_points;
        rep.split += one.split;
        rep.notes.insert(rep.notes.end(), one.notes.begin(), one.notes.end());
        rep.needles += one.needles;
        rep.contacts += one.contacts;
        if (!one.relabeled && !one.removed && !one.filled && !one.needles && !one.contacts) return true;   // the rest is below the mesh scale
        std::vector<int> origin(trial.point.size());
        std::vector<char> touched(trial.point.size());
        for (size_t i = 0; i < trial.point.size(); ++i) {
            int o = one.origin[i];
            origin[i] = o < 0 ? -1 : rep.origin[o];
            touched[i] = o < 0 || rep.touched[o] || one.touched[i];
        }
        rep.origin.swap(origin);
        rep.touched.swap(touched);
        m = std::move(trial);
    }
    rep.notes.push_back("events left after " + std::to_string(rounds) + " local rounds");
    return true;
}

}  // namespace suprem_core
