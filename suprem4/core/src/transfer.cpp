#include "suprem_core/transfer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <unordered_map>

#include "suprem_core/locate.hpp"

namespace suprem_core {
namespace {

struct PtKey {
    uint64_t x, y;
    bool operator==(const PtKey &o) const { return x == o.x && y == o.y; }
};
struct PtHash {
    size_t operator()(const PtKey &k) const { return std::hash<uint64_t>()(k.x * 0x9e3779b97f4a7c15ULL ^ k.y); }
};
PtKey key(const Point &p)
{
    PtKey k;
    std::memcpy(&k.x, &p.x, 8);
    std::memcpy(&k.y, &p.y, 8);
    return k;
}

}  // namespace

bool conserved_species(int imp)
{
    return (imp >= 2 && imp <= 5) || imp == 13 || imp == 18 || imp == 25 || (imp >= 31 && imp <= 47 && imp % 2 == 1);
}

int total_of_active(int imp)
{
    if (imp >= 20 && imp <= 23) return imp - 18;
    if (imp >= 32 && imp <= 48 && imp % 2 == 0) return imp - 1;
    return -1;
}

std::vector<double> node_volumes(const Mesh &m, bool p1)
{
    std::vector<double> vol(m.node.size(), 0.0);
    for (int t = 0; t < (int)m.triangle.size(); ++t) {
        const Triangle &tr = m.triangle[t];
        const Point *P[3] = {&m.point[tr.v[0]], &m.point[tr.v[1]], &m.point[tr.v[2]]};
        double s = 0.5 * std::fabs(orient2(m, t)), area[3];
        if (p1 || s == 0) {
            area[0] = area[1] = area[2] = s / 3;
        } else {
            double c[3], d2[3];   // half cotangent of the angle at vertex k, squared length of the opposite edge
            for (int k = 0; k < 3; ++k) {
                const Point &a = *P[k], &b = *P[(k + 1) % 3], &e = *P[(k + 2) % 3];
                double ux = b.x - a.x, uy = b.y - a.y, vx = e.x - a.x, vy = e.y - a.y;
                c[k] = 0.5 * (ux * vx + uy * vy) / std::fabs(ux * vy - uy * vx);
                d2[k] = (e.x - b.x) * (e.x - b.x) + (e.y - b.y) * (e.y - b.y);
            }
            if (c[0] >= 0 && c[1] >= 0 && c[2] >= 0) {   // circumcentric Voronoi pieces
                double a0 = 0.25 * d2[0] * c[0], a1 = 0.25 * d2[1] * c[1], a2 = 0.25 * d2[2] * c[2];
                area[0] = a1 + a2; area[1] = a2 + a0; area[2] = a0 + a1;
            } else {                                     // mixed Voronoi (Meyer et al. 2003)
                for (int k = 0; k < 3; ++k) area[k] = c[k] < 0 ? 0.5 * s : 0.25 * s;
            }
        }
        int mat = m.material_of(t);
        for (int k = 0; k < 3; ++k) {
            int n = m.node_at(tr.v[k], mat);
            if (n >= 0) vol[n] += area[k];
        }
    }
    return vol;
}

std::vector<double> inventory(const Mesh &m, const std::vector<double> &vol, int slot, int nmat)
{
    std::vector<double> inv(nmat, 0.0);
    for (size_t n = 0; n < m.node.size(); ++n)
        if (m.node[n].material < nmat) inv[m.node[n].material] += m.node[n].value[slot] * vol[n];
    return inv;
}

void transfer_fields(const Mesh &old, Mesh &neu, TransferReport &rep, bool restore, bool p1, Restore mode)
{
    const size_t ns = old.species.size();
    neu.species = old.species;
    std::unordered_map<PtKey, int, PtHash> exact;
    for (int p = 0; p < (int)old.point.size(); ++p) exact[key(old.point[p])] = p;
    std::map<int, std::unique_ptr<TriangleLocator>> loc;
    std::vector<int> old_gas;
    for (int n = 0; n < (int)old.node.size(); ++n)
        if (old.node[n].material == GAS) old_gas.push_back(n);
    // per new node and species: spread of the old values it was interpolated
    // from (zero for exact copies) = weight of the inventory correction
    std::vector<std::vector<double>> spread(neu.node.size(), std::vector<double>(ns, 0.0));
    for (size_t ni = 0; ni < neu.node.size(); ++ni) {
        Node &nd = neu.node[ni];
        nd.value.assign(ns, 0.0);
        const Point &p = neu.point[nd.point];
        if (nd.material == GAS) {
            int best = -1; double bd = 1e300;
            for (int n : old_gas) {
                const Point &q = old.point[old.node[n].point];
                double d = (q.x - p.x) * (q.x - p.x) + (q.y - p.y) * (q.y - p.y);
                if (d < bd) { bd = d; best = n; }
            }
            if (best >= 0) nd.value = old.node[best].value;
            rep.gas++;
            continue;
        }
        auto it = exact.find(key(p));
        if (it != exact.end()) {
            int on = old.node_at(it->second, nd.material);
            if (on >= 0) { nd.value = old.node[on].value; rep.exact++; continue; }
        }
        auto &L = loc[nd.material];
        if (!L) L = std::make_unique<TriangleLocator>(old, nd.material);
        double b[3];
        int t = L->locate(p, b);
        if (t < 0) {
            rep.notes.push_back("material " + std::to_string(nd.material) + " absent from the old mesh");
            continue;
        }
        // clipped coordinates (point slightly outside the old material) are reported
        {
            const Triangle &tr = old.triangle[t];
            double raw[3];
            const Point &a = old.point[tr.v[0]], &bb = old.point[tr.v[1]], &c = old.point[tr.v[2]];
            double d = (bb.x - a.x) * (c.y - a.y) - (bb.y - a.y) * (c.x - a.x);
            raw[1] = ((p.x - a.x) * (c.y - a.y) - (p.y - a.y) * (c.x - a.x)) / d;
            raw[2] = ((bb.x - a.x) * (p.y - a.y) - (bb.y - a.y) * (p.x - a.x)) / d;
            raw[0] = 1 - raw[1] - raw[2];
            if (std::min({raw[0], raw[1], raw[2]}) < -1e-6) rep.outside++;
        }
        std::vector<double> lo(ns, 1e300), hi(ns, -1e300);
        for (int k = 0; k < 3; ++k) {
            int on = old.node_at(old.triangle[t].v[k], nd.material);
            if (on < 0) continue;   // cannot happen for a valid mesh
            for (size_t s = 0; s < ns; ++s) {
                double v = old.node[on].value[s];
                nd.value[s] += b[k] * v;
                lo[s] = std::min(lo[s], v); hi[s] = std::max(hi[s], v);
            }
        }
        for (size_t s = 0; s < ns; ++s) spread[ni][s] = hi[s] > lo[s] ? hi[s] - lo[s] : 0.0;
        rep.interpolated++;
    }
    std::vector<double> v0 = node_volumes(old, p1), v1 = node_volumes(neu, p1);
    std::vector<double> box0 = node_volumes(old, false), box1;
    for (size_t s = 0; s < ns; ++s) {
        if (!conserved_species(old.species[s])) continue;
        if (box1.empty()) box1 = node_volumes(neu, false);
        int active = -1;
        for (size_t a = 0; a < ns; ++a)
            if (total_of_active(old.species[a]) == old.species[s]) active = (int)a;
        auto i0 = inventory(old, v0, (int)s), i1 = inventory(neu, v1, (int)s);
        for (int mat = 0; mat < (int)i0.size(); ++mat) {
            if (i0[mat] == 0 && i1[mat] == 0) continue;
            std::vector<double> before;   // totals before correction (for the active species)
            std::vector<size_t> idx;
            for (size_t n = 0; n < neu.node.size(); ++n)
                if (neu.node[n].material == mat) { idx.push_back(n); before.push_back(neu.node[n].value[s]); }
            double cur = i1[mat];
            // nothing to restore when the interpolated inventory already agrees
            // to rounding (sums of ~1e5 terms: 1e-12 relative)
            if (restore && i0[mat] > 0 && cur > 0 && std::fabs(cur - i0[mat]) > 1e-12 * i0[mat]) {
                if (mode == Restore::SCALE) {
                    for (size_t n : idx) neu.node[n].value[s] *= i0[mat] / cur;
                } else {
                    // minimum weighted change: c_i += lambda * w_i, w_i = spread of
                    // the old values the node was interpolated from; nodes that
                    // would go negative are pinned at zero and the rest re-solved
                    std::vector<char> pinned(idx.size(), 0);
                    for (int it = 0; it < 50; ++it) {
                        double sum = 0, sw = 0;
                        for (size_t k = 0; k < idx.size(); ++k) {
                            sum += neu.node[idx[k]].value[s] * v1[idx[k]];
                            if (!pinned[k]) sw += spread[idx[k]][s] * v1[idx[k]];
                        }
                        double resid = i0[mat] - sum;
                        if (std::fabs(resid) <= 1e-12 * i0[mat]) break;
                        if (sw <= 0) {   // nothing uncertain left: uniform scale
                            double now = i0[mat] - resid;
                            if (now > 0)
                                for (size_t n : idx) neu.node[n].value[s] *= i0[mat] / now;
                            break;
                        }
                        double lam = resid / sw;
                        bool clipped = false;
                        for (size_t k = 0; k < idx.size(); ++k) {
                            if (pinned[k]) continue;
                            double &c = neu.node[idx[k]].value[s];
                            c += lam * spread[idx[k]][s];
                            if (c < 0) { c = 0; pinned[k] = 1; clipped = true; }
                        }
                        if (!clipped) break;
                    }
                }
                if (active >= 0)
                    for (size_t k = 0; k < idx.size(); ++k) {
                        Node &nd = neu.node[idx[k]];
                        if (before[k] > 0) nd.value[active] *= nd.value[s] / before[k];
                    }
            }
            double fin = inventory(neu, v1, (int)s)[mat];
            {
                double b0 = inventory(old, box0, (int)s)[mat], b1 = inventory(neu, box1, (int)s)[mat];
                if (b0 != 0) rep.max_rel_box = std::max(rep.max_rel_box, std::fabs(b1 - b0) / std::fabs(b0));
            }
            double rel = std::fabs(i1[mat] - i0[mat]) / std::max(std::fabs(i0[mat]), 1e-300);
            double frel = std::fabs(fin - i0[mat]) / std::max(std::fabs(i0[mat]), 1e-300);
            rep.max_rel_interpolated = std::max(rep.max_rel_interpolated, rel);
            rep.max_rel_final = std::max(rep.max_rel_final, frel);
            rep.rows.push_back({mat, old.species[s], i0[mat], i1[mat], fin});
        }
    }
}

double roundtrip_error(const Mesh &old, const Mesh &neu, int slot)
{
    Mesh back = old;
    TransferReport r;
    transfer_fields(neu, back, r, false);
    std::vector<double> vol = node_volumes(old);
    std::map<int, std::pair<double, double>> acc;   // material -> (error, total)
    for (size_t n = 0; n < old.node.size(); ++n) {
        if (old.node[n].material == GAS) continue;
        auto &a = acc[old.node[n].material];
        a.first += std::fabs(back.node[n].value[slot] - old.node[n].value[slot]) * vol[n];
        a.second += std::fabs(old.node[n].value[slot]) * vol[n];
    }
    if (std::getenv("CORE_TRANSFER_DEBUG")) {   // largest weighted round-trip errors
        std::vector<std::pair<double, size_t>> e;
        for (size_t n = 0; n < old.node.size(); ++n)
            if (old.node[n].material != GAS)
                e.push_back({-std::fabs(back.node[n].value[slot] - old.node[n].value[slot]) * vol[n], n});
        std::sort(e.begin(), e.end());
        for (size_t i = 0; i < 15 && i < e.size(); ++i) {
            const Node &nd = old.node[e[i].second];
            std::printf("  RT node %zu pt %d (%.5f, %.5f) mat %d old %.4e back %.4e vol %.3e\n", e[i].second,
                        nd.point + 1, old.point[nd.point].x, old.point[nd.point].y, nd.material, nd.value[slot],
                        back.node[e[i].second].value[slot], vol[e[i].second]);
        }
    }
    double worst = 0;
    for (auto &kv : acc)
        if (kv.second.second > 0) worst = std::max(worst, kv.second.first / kv.second.second);
    return worst;
}

}  // namespace suprem_core
