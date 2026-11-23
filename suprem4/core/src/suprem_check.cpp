// SUPREM-IV validity contract; see suprem_core/suprem_check.hpp.
#include "suprem_core/suprem_check.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

namespace suprem_core {

// Verbatim port of src/geom/misc.c.
int orient2d_float_sign(const float *a, const float *b, const float *c)
{
    double terms[6], expansion[6], next[6], q, sum, bv, error;
    int i, j, n = 0, m;
    terms[0] = (double)a[0] * b[1]; terms[1] = -(double)a[0] * c[1];
    terms[2] = (double)b[0] * c[1]; terms[3] = -(double)b[0] * a[1];
    terms[4] = (double)c[0] * a[1]; terms[5] = -(double)c[0] * b[1];
    for (i = 0; i < 6; i++) {
        q = terms[i]; m = 0;
        for (j = 0; j < n; j++) {
            sum = q + expansion[j]; bv = sum - q;
            error = (q - (sum - bv)) + (expansion[j] - bv);
            if (error != 0) next[m++] = error;
            q = sum;
        }
        if (q != 0 || m == 0) next[m++] = q;
        n = m;
        for (j = 0; j < n; j++) expansion[j] = next[j];
    }
    return expansion[n - 1] > 0 ? 1 : expansion[n - 1] < 0 ? -1 : 0;
}

namespace {

using XY = std::array<float, 2>;

struct Elem { int v[3]; int boundary; int id; };

// --- mesh_embedding.h pair predicates, verbatim ---------------------------
struct Emb {
    const std::vector<XY> &xy;
    int o(int a, int b, int c) const { return orient2d_float_sign(xy[a].data(), xy[b].data(), xy[c].data()); }
    bool separated(const int *a, const int *b) const
    {
        for (int i = 0; i < 3; i++) {
            int j;
            for (j = 0; j < 3; j++) if (o(a[i], a[(i + 1) % 3], b[j]) > 0) break;
            if (j == 3) return true;
        }
        return false;
    }
    bool on_segment(int p, int a, int b) const
    {
        if (o(a, b, p)) return false;
        for (int k = 0; k < 2; k++)
            if (xy[p][k] < std::fmin(xy[a][k], xy[b][k]) || xy[p][k] > std::fmax(xy[a][k], xy[b][k])) return false;
        return true;
    }
    bool boundary_contact(int a, int b, int c, int d) const
    {
        int ac = o(a, b, c), ad = o(a, b, d), ca = o(c, d, a), cb = o(c, d, b);
        if (ac * ad < 0 && ca * cb < 0) return true;
        if (!ac && !ad && !ca && !cb) {
            int axis = xy[a][0] != xy[b][0] ? 0 : 1;
            if (std::fmin(std::fmax(xy[a][axis], xy[b][axis]), std::fmax(xy[c][axis], xy[d][axis])) >
                std::fmax(std::fmin(xy[a][axis], xy[b][axis]), std::fmin(xy[c][axis], xy[d][axis]))) return true;
        }
        return (a != c && a != d && on_segment(a, c, d)) || (b != c && b != d && on_segment(b, c, d)) ||
               (c != a && c != b && on_segment(c, a, b)) || (d != a && d != b && on_segment(d, a, b));
    }
    // true = valid pair
    bool pair(const Elem &a, const Elem &b) const
    {
        if (!separated(a.v, b.v) && !separated(b.v, a.v)) return false;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                if ((a.boundary & (1 << i)) || (b.boundary & (1 << j)))
                    if (boundary_contact(a.v[i], a.v[(i + 1) % 3], b.v[j], b.v[(j + 1) % 3])) return false;
        return true;
    }
};

// SUPREM's intang (src/dbase/geom.c) on float coordinates.
double intang(const XY &p1, const XY &p2, const XY &p3)
{
    double dx1 = (double)p2[0] - p1[0], dy1 = (double)p2[1] - p1[1];
    double dx2 = (double)p3[0] - p2[0], dy2 = (double)p3[1] - p2[1];
    if ((dx1 == 0 && dy1 == 0) || (dx2 == 0 && dy2 == 0)) return 3.40282347e+38;   // MAXFLOAT
    return M_PI - std::atan2(dx1 * dy2 - dx2 * dy1, dx1 * dx2 + dy1 * dy2);
}

std::string um(const Mesh &m, int p)
{
    std::ostringstream s;
    s.precision(9);
    s << "(" << m.point[p].x << "," << m.point[p].y << ")";
    return s.str();
}

}  // namespace

SupremCheck check_suprem(const Mesh &m, size_t max_reports)
{
    SupremCheck r;
    auto report = [&](int &counter, const std::string &s) {
        ++counter;
        if (r.messages.size() < max_reports) r.messages.push_back(s);
    };
    const int np = (int)m.point.size(), nt = (int)m.triangle.size(), nreg = (int)m.region_material.size();
    std::vector<XY> xy(np);
    for (int p = 0; p < np; ++p) xy[p] = {suprem_float(m.point[p].x), suprem_float(m.point[p].y)};

    // --- nodes (mtest1) -----------------------------------------------------
    for (int i = 0; i < nreg; ++i) {
        int mat = m.region_material[i];
        if (mat != SI && mat != POLY && mat != SINI && mat != SIO2 && mat != OXNI && mat != PHRS && mat != AL && mat != GAAS)
            report(r.nodes, "region " + std::to_string(i + 1) + " has bad material code " + std::to_string(mat));
    }
    for (int t = 0; t < nt; ++t) {
        const Triangle &tr = m.triangle[t];
        if (tr.region < 0 || tr.region >= nreg) { report(r.nodes, "triangle " + std::to_string(t + 1) + " has a bad region"); return r; }
        for (int i = 0; i < 3; ++i)
            if (tr.v[i] < 0 || tr.v[i] >= np) { report(r.nodes, "triangle " + std::to_string(t + 1) + " has a bad vertex"); return r; }
        if (tr.v[0] == tr.v[1] || tr.v[1] == tr.v[2] || tr.v[0] == tr.v[2])
            report(r.nodes, "triangle " + std::to_string(t + 1) + " has duplicate nodes");
    }
    std::set<std::pair<int, int>> node_seen, material_at;   // (point, material)
    for (const Triangle &tr : m.triangle)
        for (int i = 0; i < 3; ++i) material_at.insert({tr.v[i], m.region_material[tr.region]});
    std::vector<char> surf(np, 0), back(np, 0), gas(np, 0);
    for (const Triangle &tr : m.triangle)
        for (int j = 0; j < 3; ++j) {
            if (tr.nb[j] == EXPOSED) surf[tr.v[(j + 1) % 3]] = surf[tr.v[(j + 2) % 3]] = 1;
            if (tr.nb[j] == BACKSIDE) back[tr.v[(j + 1) % 3]] = back[tr.v[(j + 2) % 3]] = 1;
        }
    for (size_t n = 0; n < m.node.size(); ++n) {
        const Node &nd = m.node[n];
        if (nd.point < 0 || nd.point >= np) { report(r.nodes, "node " + std::to_string(n) + " has a bad point"); continue; }
        if (!node_seen.insert({nd.point, nd.material}).second)
            report(r.nodes, "two nodes of material " + std::to_string(nd.material) + " at point " + std::to_string(nd.point + 1));
        if (nd.material == GAS) gas[nd.point] = 1;
        else if (!material_at.count({nd.point, nd.material}))
            report(r.nodes, "material node " + std::to_string(n) + " belongs to no triangles " + um(m, nd.point));
    }
    for (const auto &pm : material_at)
        if (!node_seen.count(pm))
            report(r.nodes, "no node of material " + std::to_string(pm.second) + " at point " + std::to_string(pm.first + 1));
    for (int p = 0; p < np; ++p) {
        // create_db adds a missing gas node at a surface point, not at a backside-only one
        if (gas[p]) {
            if (!surf[p] && !back[p]) report(r.nodes, "unexposed point " + std::to_string(p + 1) + " has gas node " + um(m, p));
        } else if (surf[p]) {
            report(r.gas_added, "exposed point " + std::to_string(p + 1) + " has no gas node (added on load) " + um(m, p));
        } else if (back[p]) {
            report(r.nodes, "backside point " + std::to_string(p + 1) + " has no gas node " + um(m, p));
        }
    }

    // --- orientation (sclk_tri swaps negative triangles; flat ones panic) -------
    std::vector<std::array<int, 3>> vert(nt);   // vertex order SUPREM holds after loading
    for (int t = 0; t < nt; ++t) {
        const Triangle &tr = m.triangle[t];
        int s = orient2d_float_sign(xy[tr.v[0]].data(), xy[tr.v[1]].data(), xy[tr.v[2]].data());
        vert[t] = {tr.v[0], tr.v[1], tr.v[2]};
        if (s < 0) {
            std::swap(vert[t][1], vert[t][2]);
            report(r.ordering, "triangle " + std::to_string(t + 1) + " clockwise as written (SUPREM swaps it) at " + um(m, tr.v[0]));
        }
        if (s == 0)
            report(r.orientation, "triangle " + std::to_string(t + 1) + " flat in float cm at " + um(m, tr.v[0]) +
                                      um(m, tr.v[1]) + um(m, tr.v[2]));
    }

    // --- neighbours from shared points (nxtel) ----------------------------------
    std::map<std::pair<int, int>, std::vector<int>> edge_tris;   // sorted point pair -> triangles
    for (int t = 0; t < nt; ++t)
        for (int j = 0; j < 3; ++j) {
            int a = m.triangle[t].v[j], b = m.triangle[t].v[(j + 1) % 3];
            edge_tris[{std::min(a, b), std::max(a, b)}].push_back(t);
        }
    for (const auto &e : edge_tris)
        if (e.second.size() > 2)
            report(r.topology, "edge " + um(m, e.first.first) + "-" + um(m, e.first.second) + " shared by " +
                                   std::to_string(e.second.size()) + " triangles");

    // --- embedding (only meaningful without flat triangles) ----------------------
    if (!r.orientation) {
        std::vector<Elem> el(nt);
        for (int t = 0; t < nt; ++t) {
            el[t].id = t; el[t].boundary = 0;
            for (int j = 0; j < 3; ++j) {
                el[t].v[j] = vert[t][j];
                int a = vert[t][j], b = vert[t][(j + 1) % 3];
                if (edge_tris[{std::min(a, b), std::max(a, b)}].size() < 2) el[t].boundary |= 1 << j;
            }
        }
        // broad phase: sweep in x over inclusive float boxes (mesh_embedding.h tests touching boxes)
        struct Box { float lo[2], hi[2]; int t; };
        std::vector<Box> box(nt);
        for (int t = 0; t < nt; ++t) {
            Box &b = box[t]; b.t = t;
            for (int k = 0; k < 2; ++k) {
                b.lo[k] = std::min({xy[el[t].v[0]][k], xy[el[t].v[1]][k], xy[el[t].v[2]][k]});
                b.hi[k] = std::max({xy[el[t].v[0]][k], xy[el[t].v[1]][k], xy[el[t].v[2]][k]});
            }
        }
        std::sort(box.begin(), box.end(), [](const Box &a, const Box &b) { return a.lo[0] < b.lo[0]; });
        Emb emb{xy};
        std::set<int> bad_tris;
        for (int i = 0; i < nt; ++i)
            for (int k = i + 1; k < nt && box[k].lo[0] <= box[i].hi[0]; ++k) {
                if (box[k].hi[1] < box[i].lo[1] || box[i].hi[1] < box[k].lo[1]) continue;
                if (!emb.pair(el[box[i].t], el[box[k].t])) {
                    int a = box[i].t, b = box[k].t;
                    if (bad_tris.insert(a).second | bad_tris.insert(b).second)
                        report(r.embedding, "triangles " + std::to_string(a + 1) + " and " + std::to_string(b + 1) +
                                                " overlap or touch at boundary near " + um(m, m.triangle[a].v[0]));
                }
            }
    }

    // --- region border loops (skel_reg + ck_clock) ---------------------------
    // Border edge of region r: a point-edge of a triangle in r whose other
    // triangle (if any) is in another region.  (Different materials have
    // different nodes in SUPREM, so their shared edge is a border of both.)
    std::vector<std::vector<std::pair<int, int>>> border(nreg);   // directed, CCW as in the triangle
    for (int t = 0; t < nt; ++t)
        for (int j = 0; j < 3; ++j) {
            int a = vert[t][j], b = vert[t][(j + 1) % 3], reg = m.triangle[t].region;
            const auto &ts = edge_tris[{std::min(a, b), std::max(a, b)}];
            bool inner = false;
            for (int o : ts) if (o != t && m.triangle[o].region == reg) inner = true;
            if (!inner) border[reg].push_back({a, b});
        }
    for (int reg = 0; reg < nreg; ++reg) {
        const auto &be = border[reg];
        std::string id = "region " + std::to_string(reg + 1) + " (material " + std::to_string(m.region_material[reg]) + ")";
        if (be.empty()) { report(r.loops, id + " has no triangles"); continue; }
        std::map<int, int> next, degree;
        bool bad = false;
        for (const auto &e : be) {
            degree[e.first]++; degree[e.second]++;
            if (!next.insert({e.first, e.second}).second) bad = true;
        }
        for (const auto &d : degree)
            if (d.second != 2) {
                report(r.loops, id + " border touches itself at point " + std::to_string(d.first + 1) + " " + um(m, d.first));
                bad = true;
                break;
            }
        if (bad) continue;
        // every border point has one incoming and one outgoing edge: the border is a set of
        // disjoint simple cycles.  skel_reg walks only the first; ck_clock checks its turn sum.
        std::set<int> unvisited;
        for (const auto &e : be) unvisited.insert(e.first);
        int outer = 0, holes = 0;
        while (!unvisited.empty()) {
            std::vector<int> loop;
            int start = *unvisited.begin(), p = start;
            do {
                loop.push_back(p);
                unvisited.erase(p);
                p = next.at(p);
            } while (p != start && loop.size() <= be.size());
            double tang = 0, area = 0;
            const size_t n = loop.size();
            for (size_t i = 0; i < n; ++i) {
                tang += M_PI - intang(xy[loop[(i + n - 1) % n]], xy[loop[i]], xy[loop[(i + 1) % n]]);
                const XY &a = xy[loop[i]], &b = xy[loop[(i + 1) % n]];
                area += (double)a[0] * b[1] - (double)b[0] * a[1];
            }
            if (!(std::fabs(tang - 2 * M_PI) < 1e-6 || std::fabs(tang + 2 * M_PI) < 1e-6))
                report(r.loops, id + " turn sum " + std::to_string(tang) + " (crosses itself or malformed) near " +
                                    um(m, loop[0]));
            (area > 0 ? outer : holes)++;
        }
        if (outer != 1 || holes)
            report(r.components, id + " has " + std::to_string(outer) + " components and " + std::to_string(holes) +
                                     " holes (SUPREM loads it but its skeleton holds only one loop)");
    }
    return r;
}

}  // namespace suprem_core
