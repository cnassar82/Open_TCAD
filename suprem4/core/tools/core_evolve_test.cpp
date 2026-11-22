// core_evolve_test REFERENCE.str [OUTDIR]
// WP5 synthetic moving-boundary tests on a 1 x 1 um domain (y down, bottom =
// backside, sides reflecting, top exposed); REFERENCE.str only supplies the
// file metadata.  Every step is checked (invariants + B-rep).
//   shrink  poly square in oxide consumed from all sides (interface moving
//           inward; poly interior and far oxide still): collapses to a thin
//           remnant that is removed; one oxide region fills the domain
//   merge   silicon slot closing from both sides: top closes first (gap
//           event), the trapped pocket below is a cavity event; one silicon
//           region, no ambient below the top surface
//   slit    slit through an oxide layer on silicon closing by lateral oxide
//           motion: filled with oxide (priority over silicon); one oxide region
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <string>

#include "suprem_core/brep.hpp"
#include "suprem_core/evolve.hpp"
#include "suprem_core/remesh.hpp"
#include "suprem_core/str_file.hpp"
#include "suprem_core/transfer.hpp"

using namespace suprem_core;

namespace {

std::vector<Point> rect(double x0, double y0, double x1, double y1, bool reverse = false)
{
    std::vector<Point> r{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
    if (reverse) std::reverse(r.begin(), r.end());
    return r;
}

int domain_code(const Point &a, const Point &b)
{
    const double e = 1e-9;
    if (std::fabs(a.y - 1) < e && std::fabs(b.y - 1) < e) return BACKSIDE;
    if ((std::fabs(a.x) < e && std::fabs(b.x) < e) || (std::fabs(a.x - 1) < e && std::fabs(b.x - 1) < e)) return BC_OFFSET;
    return EXPOSED;
}

double ramp(double d, double full, double zero)   // 1 for d <= full, 0 for d >= zero
{
    if (d <= full) return 1;
    if (d >= zero) return 0;
    return (zero - d) / (zero - full);
}

struct Case {
    std::string name;
    std::vector<Shape> shapes;
    Velocity v;
    double time;
    std::function<std::vector<std::string>(const Mesh &)> verify;
};

std::map<int, int> regions_per_material(const Mesh &m)
{
    std::map<int, int> c;
    for (int mat : m.region_material) c[mat]++;
    return c;
}

double total_area(const Mesh &m)
{
    double a = 0;
    for (size_t t = 0; t < m.triangle.size(); ++t) a += 0.5 * orient2(m, (int)t);
    return a;
}

double max_exposed_y(const Mesh &m)
{
    double y = -1e300;
    for (const Triangle &t : m.triangle)
        for (int j = 0; j < 3; ++j)
            if (t.nb[j] == EXPOSED) y = std::max({y, m.point[t.v[(j + 1) % 3]].y, m.point[t.v[(j + 2) % 3]].y});
    return y;
}

}  // namespace

int main(int argc, char **argv)
{
    if (argc < 2) { std::fprintf(stderr, "usage: core_evolve_test REFERENCE.str [OUTDIR]\n"); return 2; }
    std::string outdir = argc > 2 ? argv[2] : "";
    Mesh ref = from_str(read_str(argv[1]));
    ref.species = {5};   // boron, to exercise the transfer
    std::vector<Case> cases;
    {   // shrinking square: poly square [0.45, 0.55]^2 in oxide, consumed from all
        // sides (interface moves inward at unit speed); only a band around
        // the moving interface moves, the poly interior and far oxide are still
        Shape ox{SIO2, {rect(0, 0, 1, 1), rect(0.45, 0.45, 0.55, 0.55, true)}};
        Shape po{POLY, {rect(0.45, 0.45, 0.55, 0.55)}};
        // homogeneous contraction of the poly (interface at unit speed),
        // continued smoothly into the oxide and decaying to zero 0.03 out
        Velocity v = [](const Point &p, double t) {
            double dx = p.x - 0.5, dy = p.y - 0.5, a = std::max(1e-6, 0.05 - t);
            double r = std::max(std::fabs(dx), std::fabs(dy));   // L-infinity distance
            if (r < 1e-12) return Point{0, 0};
            double w = r <= a ? r / a : ramp(r - a, 0, 0.03);      // speed along the ray, 1 at the interface
            return Point{-w * dx / r, -w * dy / r};
        };
        cases.push_back({"shrink", {ox, po}, v, 0.055, [](const Mesh &m) {
            std::vector<std::string> bad;
            auto c = regions_per_material(m);
            if (c.count(POLY)) bad.push_back("poly left: " + std::to_string(c[POLY]) + " regions");
            if (c[SIO2] != 1) bad.push_back("oxide regions " + std::to_string(c[SIO2]));
            if (std::fabs(total_area(m) - 1) > 1e-6) bad.push_back("area " + std::to_string(total_area(m)));
            return bad;
        }});
    }
    {   // merging fronts: silicon slot x in [0.45, 0.55], y in [0, 0.5]
        Shape si{SI, {{{0, 0}, {0.45, 0}, {0.45, 0.5}, {0.55, 0.5}, {0.55, 0}, {1, 0}, {1, 1}, {0, 1}}}};
        Velocity v = [](const Point &p, double) {
            double d = std::fabs(p.x - 0.5), s = ramp(d, 0.1, 0.3) * ramp(p.y, 0.4, 0.5);
            return Point{(p.x < 0.5 ? 1 : -1) * s, 0};
        };
        cases.push_back({"merge", {si}, v, 0.052, [](const Mesh &m) {
            std::vector<std::string> bad;
            auto c = regions_per_material(m);
            if (c[SI] != 1 || c.size() != 1) bad.push_back("expected one silicon region");
            if (max_exposed_y(m) > 1e-3) bad.push_back("ambient below the top surface, y = " + std::to_string(max_exposed_y(m)));
            return bad;
        }});
    }
    {   // slit in an oxide layer [0,1] x [0,0.1] over silicon
        Shape si{SI, {rect(0, 0.1, 1, 1)}};
        Shape oxl{SIO2, {rect(0, 0, 0.49, 0.1)}}, oxr{SIO2, {rect(0.51, 0, 1, 0.1)}};
        Velocity v = [](const Point &p, double) {
            double d = std::fabs(p.x - 0.5), s = ramp(d, 0.05, 0.2) * ramp(p.y, 0.08, 0.1);
            return Point{(p.x < 0.5 ? 1 : -1) * s, 0};
        };
        cases.push_back({"slit", {si, oxl, oxr}, v, 0.0105, [](const Mesh &m) {
            std::vector<std::string> bad;
            auto c = regions_per_material(m);
            if (c[SIO2] != 1) bad.push_back("oxide regions " + std::to_string(c[SIO2]));
            if (c[SI] != 1) bad.push_back("silicon regions " + std::to_string(c[SI]));
            if (max_exposed_y(m) > 1e-3) bad.push_back("ambient below the top surface, y = " + std::to_string(max_exposed_y(m)));
            if (std::fabs(total_area(m) - 1) > 1e-3) bad.push_back("area " + std::to_string(total_area(m)));
            return bad;
        }});
    }
    int failed = 0;
    const char *only = std::getenv("CORE_EVOLVE_ONLY");
    for (Case &c : cases) {
        if (only && c.name != only) continue;
        auto t0 = std::chrono::steady_clock::now();
        EvolveOptions o;
        o.remesh.hmax = 0.02;
        RemeshReport rr;
        std::vector<std::string> notes;
        Mesh m = mesh_from_shapes(c.shapes, ref, domain_code, nullptr, o.remesh, rr, &notes);
        for (Node &n : m.node) if (n.material != GAS) n.value[0] = 1e17;
        if (!outdir.empty()) write_str(to_str(m), outdir + "/" + c.name + "_start.str");
        EvolveReport rep;
        std::vector<std::string> bad;
        try {
            evolve(m, c.v, c.time, o, rep);
        } catch (const std::exception &e) {
            bad.push_back(std::string("exception: ") + e.what());
        }
        bad.insert(bad.end(), rep.failures.begin(), rep.failures.end());
        if (bad.empty()) {
            auto inv = check_invariants(m);
            bad.insert(bad.end(), inv.begin(), inv.end());
            auto gas = check_gas_nodes(m);
            bad.insert(bad.end(), gas.begin(), gas.end());
            auto v = c.verify(m);
            bad.insert(bad.end(), v.begin(), v.end());
        }
        double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::map<std::string, int> kinds;
        for (auto &e : rep.events) kinds[e.kind]++;
        std::printf("%-7s %s | %.1f s | steps %d rebuilds %d halvings %d | points %zu triangles %zu | events:",
                    c.name.c_str(), bad.empty() ? "PASS" : "FAIL", sec, rep.steps, rep.rebuilds, rep.halvings,
                    m.point.size(), m.triangle.size());
        for (auto &k : kinds) std::printf(" %s %d", k.first.c_str(), k.second);
        std::printf("\n");
        for (auto &e : rep.events)
            std::printf("    event %-6s mat %d -> %d area %.3e at (%.4f, %.4f)\n", e.kind.c_str(), e.material, e.to,
                        e.area, e.at.x, e.at.y);
        for (auto &n : rep.notes) std::printf("    NOTE %s\n", n.c_str());
        for (auto &b : bad) std::printf("    FAIL %s\n", b.c_str());
        if (!outdir.empty()) write_str(to_str(m), outdir + "/" + c.name + "_end.str");
        failed += !bad.empty();
    }
    std::printf("%zu/%zu passed\n", cases.size() - failed, cases.size());
    return failed ? 1 : 0;
}
