#include "suprem_core/mesh.hpp"

#include <algorithm>
#include <array>
#include <set>
#include <sstream>
#include <stdexcept>

namespace suprem_core {

int Mesh::node_at(int p, int material) const
{
    auto it = node_index_.find({p, material});
    return it == node_index_.end() ? -1 : it->second;
}

void Mesh::rebuild_index()
{
    node_index_.clear();
    for (size_t n = 0; n < node.size(); ++n) node_index_[{node[n].point, node[n].material}] = (int)n;
}

Mesh from_str(const StrFile &f)
{
    if (f.dim != 2 || f.nvrt != 3) throw std::runtime_error("only 2D triangle meshes are supported");
    Mesh m;
    m.version = f.version; m.dim = f.dim; m.nvrt = f.nvrt; m.nedg = f.nedg;
    m.has_model = f.has_model; m.orientation = f.orientation; m.last_temperature = f.last_temperature;
    m.implanted = f.implanted; m.trailer = f.trailer; m.species = f.species;
    m.region_material = f.region_material;
    for (const auto &c : f.coord) m.point.push_back({c[0], c[1]});
    for (const StrTriangle &s : f.triangle) {
        Triangle t;
        t.region = s.region - 1;
        for (int i = 0; i < 3; ++i) {
            t.v[i] = s.vertex[i] - 1;
            t.nb[i] = s.code[i] >= 0 ? s.code[i] - 1 : s.code[i];   // tcode(): neighbours are 1-based
        }
        t.father = s.father; t.offspring = s.offspring;
        m.triangle.push_back(t);
    }
    for (const StrNode &s : f.node) m.node.push_back({s.point, s.material, s.value});
    m.rebuild_index();
    return m;
}

StrFile to_str(const Mesh &m)
{
    StrFile f;
    f.version = m.version; f.dim = m.dim; f.nvrt = m.nvrt; f.nedg = m.nedg;
    f.has_model = m.has_model; f.orientation = m.orientation; f.last_temperature = m.last_temperature;
    f.implanted = m.implanted; f.trailer = m.trailer; f.species = m.species;
    f.region_material = m.region_material;
    for (const Point &p : m.point) f.coord.push_back({p.x, p.y});
    for (const Triangle &t : m.triangle) {
        StrTriangle s;
        s.region = t.region + 1;
        for (int i = 0; i < 3; ++i) {
            s.vertex[i] = t.v[i] + 1;
            s.code[i] = t.nb[i] >= 0 ? t.nb[i] + 1 : t.nb[i];
        }
        s.father = t.father; s.offspring = t.offspring;
        f.triangle.push_back(s);
    }
    for (const Node &n : m.node) f.node.push_back({n.point, n.material, n.value});
    return f;
}

double orient2(const Mesh &m, int t)
{
    const Point &a = m.point[m.triangle[t].v[0]], &b = m.point[m.triangle[t].v[1]], &c = m.point[m.triangle[t].v[2]];
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

std::vector<std::string> check_invariants(const Mesh &m, size_t max_reports)
{
    std::vector<std::string> bad;
    auto report = [&](const std::string &s) { if (bad.size() < max_reports) bad.push_back(s); };
    const int np = (int)m.point.size(), nt = (int)m.triangle.size();
    std::set<std::array<int, 3>> seen;
    for (int t = 0; t < nt; ++t) {
        const Triangle &tr = m.triangle[t];
        std::ostringstream id; id << "triangle " << t + 1;
        if (tr.region < 0 || tr.region >= (int)m.region_material.size()) { report(id.str() + ": bad region"); continue; }
        bool ok = true;
        for (int i = 0; i < 3; ++i) if (tr.v[i] < 0 || tr.v[i] >= np) { report(id.str() + ": bad vertex"); ok = false; }
        if (!ok) continue;
        if (orient2(m, t) <= 0) report(id.str() + ": non-positive orientation");
        std::array<int, 3> key = {tr.v[0], tr.v[1], tr.v[2]};
        std::sort(key.begin(), key.end());
        if (!seen.insert(key).second) report(id.str() + ": duplicate of another triangle");
        int mat = m.region_material[tr.region];
        for (int i = 0; i < 3; ++i) {
            int n = m.node_at(tr.v[i], mat);
            if (n < 0) report(id.str() + ": no node of material " + std::to_string(mat) + " at point " + std::to_string(tr.v[i] + 1));
        }
        for (int j = 0; j < 3; ++j) {
            int nb = tr.nb[j];
            if (nb < 0) continue;
            if (nb >= nt) { report(id.str() + ": neighbour out of range"); continue; }
            // face j joins v[j+1], v[j+2]; the neighbour must point back across the same face
            int a = tr.v[(j + 1) % 3], b = tr.v[(j + 2) % 3], back = -1;
            const Triangle &o = m.triangle[nb];
            for (int k = 0; k < 3; ++k)
                if (o.nb[k] == t) {
                    int c = o.v[(k + 1) % 3], d = o.v[(k + 2) % 3];
                    if ((c == a && d == b) || (c == b && d == a)) back = k;
                }
            if (back < 0) report(id.str() + ": neighbour " + std::to_string(nb + 1) + " does not share face " + std::to_string(j));
        }
    }
    for (size_t n = 0; n < m.node.size(); ++n) {
        if (m.node[n].point < 0 || m.node[n].point >= np) report("node " + std::to_string(n) + ": bad point");
        if (m.node[n].value.size() != m.species.size()) report("node " + std::to_string(n) + ": value count");
    }
    return bad;
}

}  // namespace suprem_core
