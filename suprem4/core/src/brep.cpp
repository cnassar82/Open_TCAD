#include "suprem_core/brep.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>

namespace suprem_core {
namespace {

struct Edge { int a, b, right; };   // directed, region on the left

double cross(const Point &o, const Point &a, const Point &b)
{
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// Proper or touching intersection of segments pq and rs (shared endpoints by
// index are excluded by the caller).
bool segments_touch(const Point &p, const Point &q, const Point &r, const Point &s)
{
    double d1 = cross(r, s, p), d2 = cross(r, s, q), d3 = cross(p, q, r), d4 = cross(p, q, s);
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) return true;
    auto on = [](const Point &a, const Point &b, const Point &c, double d) {
        return d == 0 && std::min(a.x, b.x) <= c.x && c.x <= std::max(a.x, b.x) &&
               std::min(a.y, b.y) <= c.y && c.y <= std::max(a.y, b.y);
    };
    return on(r, s, p, d1) || on(r, s, q, d2) || on(p, q, r, d3) || on(p, q, s, d4);
}

}  // namespace

Brep extract_brep(const Mesh &m, std::vector<std::string> &problems)
{
    Brep b;
    const int nr = (int)m.region_material.size();
    std::vector<std::vector<Edge>> edges(nr);
    for (int t = 0; t < (int)m.triangle.size(); ++t) {
        const Triangle &tr = m.triangle[t];
        for (int j = 0; j < 3; ++j) {
            int nb = tr.nb[j];
            if (nb >= 0 && m.triangle[nb].region == tr.region) continue;
            int a = tr.v[(j + 1) % 3], c = tr.v[(j + 2) % 3];
            // counter-clockwise triangle: face (v[j+1] -> v[j+2]) has the triangle on its left
            edges[tr.region].push_back({a, c, nb >= 0 ? m.triangle[nb].region : nb});
        }
    }
    for (int r = 0; r < nr; ++r) {
        auto &E = edges[r];
        std::multimap<int, int> from;
        for (int i = 0; i < (int)E.size(); ++i) from.insert({E[i].a, i});
        std::vector<char> used(E.size(), 0);
        for (int start = 0; start < (int)E.size(); ++start) {
            if (used[start]) continue;
            Loop L; L.region = r;
            std::vector<int> right;
            int cur = start;
            for (size_t guard = 0; guard <= E.size(); ++guard) {
                used[cur] = 1;
                L.point.push_back(E[cur].a);
                right.push_back(E[cur].right);
                if (E[cur].b == E[start].a) { cur = -1; break; }
                // next unused edge leaving E[cur].b; at a pinch take the one that turns most to the left
                int next = -1; double best = -10;
                auto range = from.equal_range(E[cur].b);
                for (auto it = range.first; it != range.second; ++it) {
                    int k = it->second;
                    if (used[k]) continue;
                    const Point &p = m.point[E[cur].a], &q = m.point[E[cur].b], &s = m.point[E[k].b];
                    double ang = std::atan2(cross(p, q, s), (q.x - p.x) * (s.x - q.x) + (q.y - p.y) * (s.y - q.y));
                    if (ang > best) { best = ang; next = k; }
                }
                if (next < 0) break;
                cur = next;
            }
            if (cur != -1) { problems.push_back("open boundary loop in region " + std::to_string(r + 1)); continue; }
            for (size_t i = 0; i < L.point.size(); ++i) {
                const Point &p = m.point[L.point[i]], &q = m.point[L.point[(i + 1) % L.point.size()]];
                L.area += 0.5 * (p.x * q.y - q.x * p.y);
            }
            // chains: rotate the loop to start where the right side changes,
            // then group consecutive edges with the same right side
            size_t n = L.point.size(), s0 = 0;
            while (s0 < n && right[s0] == right[(s0 + n - 1) % n]) ++s0;
            if (s0 == n) {      // one right side all round: closed chain, first point repeated
                Chain c{r, right[0], L.point};
                c.point.push_back(L.point[0]);
                b.chain.push_back(c);
            } else {
                std::vector<int> P(n), R(n);
                for (size_t k = 0; k < n; ++k) { P[k] = L.point[(s0 + k) % n]; R[k] = right[(s0 + k) % n]; }
                for (size_t i = 0; i < n;) {
                    Chain c{r, R[i], {P[i]}};
                    size_t j = i;
                    while (j < n && R[j] == R[i]) { c.point.push_back(P[(j + 1) % n]); ++j; }
                    b.chain.push_back(c);
                    i = j;
                }
            }
            b.loop.push_back(L);
        }
    }
    return b;
}

std::vector<std::string> check_brep(const Mesh &m, const Brep &b, size_t max_reports)
{
    std::vector<std::string> bad;
    auto report = [&](const std::string &s) { if (bad.size() < max_reports) bad.push_back(s); };
    // areas per region
    std::vector<double> tri_area(m.region_material.size(), 0), loop_area(m.region_material.size(), 0);
    for (int t = 0; t < (int)m.triangle.size(); ++t) tri_area[m.triangle[t].region] += 0.5 * orient2(m, t);
    for (const Loop &L : b.loop) loop_area[L.region] += L.area;
    for (size_t r = 0; r < tri_area.size(); ++r) {
        double scale = std::max(1.0, std::fabs(tri_area[r]));
        if (std::fabs(tri_area[r] - loop_area[r]) > 1e-9 * scale)
            report("region " + std::to_string(r + 1) + ": loop area " + std::to_string(loop_area[r]) +
                   " != triangle area " + std::to_string(tri_area[r]));
    }
    // interface chains have reversed twins
    std::map<std::pair<int, int>, std::vector<const Chain *>> by_pair;
    for (const Chain &c : b.chain) if (c.right >= 0) by_pair[{c.left_region, c.right}].push_back(&c);
    for (const Chain &c : b.chain) {
        if (c.right < 0) continue;
        bool found = false;
        for (const Chain *o : by_pair[{c.right, c.left_region}]) {
            std::vector<int> rev(o->point.rbegin(), o->point.rend());
            if (rev == c.point) { found = true; break; }
        }
        if (!found) report("interface chain " + std::to_string(c.left_region + 1) + "|" + std::to_string(c.right + 1) +
                           " (" + std::to_string(c.point.size()) + " points) has no exact twin");
    }
    // simplicity: no two boundary segments of the whole mesh cross or touch
    // except at shared endpoints (each undirected segment taken once)
    struct Seg { int a, b; double x0, x1, y0, y1; };
    std::vector<Seg> seg;
    std::set<std::pair<int, int>> taken;
    for (const Loop &L : b.loop)
        for (size_t i = 0; i < L.point.size(); ++i) {
            int a = L.point[i], c = L.point[(i + 1) % L.point.size()];
            if (!taken.insert({std::min(a, c), std::max(a, c)}).second) continue;
            const Point &p = m.point[a], &q = m.point[c];
            seg.push_back({a, c, std::min(p.x, q.x), std::max(p.x, q.x), std::min(p.y, q.y), std::max(p.y, q.y)});
        }
    std::sort(seg.begin(), seg.end(), [](const Seg &u, const Seg &v) { return u.x0 < v.x0; });
    for (size_t i = 0; i < seg.size(); ++i)
        for (size_t j = i + 1; j < seg.size() && seg[j].x0 <= seg[i].x1; ++j) {
            const Seg &u = seg[i], &v = seg[j];
            if (v.y0 > u.y1 || v.y1 < u.y0) continue;
            if (u.a == v.a || u.a == v.b || u.b == v.a || u.b == v.b) continue;
            if (segments_touch(m.point[u.a], m.point[u.b], m.point[v.a], m.point[v.b]))
                report("boundary segments " + std::to_string(u.a + 1) + "-" + std::to_string(u.b + 1) + " and " +
                       std::to_string(v.a + 1) + "-" + std::to_string(v.b + 1) + " intersect");
        }
    return bad;
}

}  // namespace suprem_core
