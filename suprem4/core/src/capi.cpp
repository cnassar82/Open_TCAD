// C entry point for SUPREM (WP6).  SUPREM writes its current structure,
// calls this, and reloads the result when a new mesh was written.
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>

#include "suprem_core/brep.hpp"
#include "suprem_core/evolve.hpp"
#include "suprem_core/locate.hpp"
#include "suprem_core/mesh.hpp"
#include "suprem_core/remesh.hpp"
#include "suprem_core/str_file.hpp"
#include "suprem_core/transfer.hpp"

using namespace suprem_core;

namespace {
double env(const char *name, double dflt)
{
    const char *v = std::getenv(name);
    return v && *v ? std::atof(v) : dflt;
}
}  // namespace

// Returns 1 when a repaired / regenerated mesh was written to out_path, 0
// when the structure needs nothing, < 0 on error (msg holds the reason).
// Rebuild when: a geometric event of size delta_um exists, or the worst
// box-method coupling is below min_coupling, or force != 0.
extern "C" int suprem_core_oxidation_check(const char *in_path, const char *out_path, double delta_um,
                                           double min_coupling, int force, char *msg, int msglen)
{
    auto say = [&](const std::string &s) {
        if (msg && msglen > 0) std::snprintf(msg, (size_t)msglen, "%s", s.c_str());
    };
    try {
        Mesh m = from_str(read_str(in_path));
        // SUPREM's own motion can leave crossing boundaries (a collision):
        // then the shapes come from the union of the triangles per material
        bool valid;
        {
            std::vector<std::string> pr;
            Brep b = extract_brep(m, pr);
            valid = pr.empty() && check_brep(m, b, 1).empty();
        }
        std::vector<Shape> shapes = valid ? mesh_shapes(m) : mesh_shapes_union(m);
        RepairOptions ro;
        ro.delta = delta_um;
        std::vector<RepairEvent> ev = repair_shapes(shapes, ro);
        Quality q = measure_quality(m);
        bool degraded = q.worst_coupling < min_coupling;
        if (valid && ev.empty() && !degraded && !force) {
            say("ok");
            return 0;
        }
        RemeshOptions o;
        o.hmax = env("SUPREM4_CORE_HMAX", o.hmax);
        o.scale = env("SUPREM4_CORE_SCALE", o.scale);
        o.float_cm_grid = true;
        // target sizes from the mesh the core built last (in this process), not
        // from SUPREM's deformed mesh: otherwise sizes ratchet down every call
        static Mesh built;
        static bool have_built = false;
        SizeField h(have_built ? built : m, SizeField::AREA, o.scale, o.hmin, o.hmax);
        h.adapt_to_solution(o.solution_decades);
        RemeshReport rr;
        std::vector<std::string> notes;
        // no event: remesh the mesh itself (exact geometry); events: rebuild
        // from the repaired polygons
        Mesh fresh = (valid && ev.empty()) ? remesh(m, o, rr, &h)
                                           : mesh_from_shapes(shapes, m, boundary_codes_of(m), &h, o, rr, &notes);
        TransferReport tr;
        transfer_fields(m, fresh, tr);
        // SUPREM stores coordinates as float (cm): use exactly the values it
        // will hold and require every triangle to stay positively oriented
        for (Point &p : fresh.point) {
            p.x = (double)(float)(1e-4 * p.x) * 1e4;
            p.y = (double)(float)(1e-4 * p.y) * 1e4;
        }
        for (size_t t = 0; t < fresh.triangle.size(); ++t) {
            const Triangle &tr = fresh.triangle[t];
            long double ax = (float)(1e-4 * fresh.point[tr.v[0]].x), ay = (float)(1e-4 * fresh.point[tr.v[0]].y);
            long double bx = (float)(1e-4 * fresh.point[tr.v[1]].x), by = (float)(1e-4 * fresh.point[tr.v[1]].y);
            long double cx = (float)(1e-4 * fresh.point[tr.v[2]].x), cy = (float)(1e-4 * fresh.point[tr.v[2]].y);
            if ((bx - ax) * (cy - ay) - (by - ay) * (cx - ax) <= 0) {
                char b[200];
                std::snprintf(b, sizeof b, "triangle %zu not positive in SUPREM float coordinates at (%.7f, %.7f)", t,
                              fresh.point[tr.v[0]].x, fresh.point[tr.v[0]].y);
                say(b);
                return -3;
            }
        }
        auto bad = check_invariants(fresh, 3);
        auto gas = check_gas_nodes(fresh, 3);
        bad.insert(bad.end(), gas.begin(), gas.end());
        if (!bad.empty()) {
            say("rebuilt mesh invalid: " + bad.front());
            return -2;
        }
        write_str(to_str(fresh), out_path);
        built = fresh;
        have_built = true;
        char buf[512];
        std::snprintf(buf, sizeof buf, "%s%zu events, worst coupling %.3g -> %.3g, points %zu -> %zu, inventory error %.2g",
                      valid ? "" : "crossing boundaries repaired, ", ev.size(), q.worst_coupling, measure_quality(fresh).worst_coupling, m.point.size(),
                      fresh.point.size(), tr.max_rel_final);
        std::string s = buf;
        for (auto &e : ev) {
            std::snprintf(buf, sizeof buf, "; %s %d->%d %.2e um2 at (%.4f, %.4f)", e.kind.c_str(), e.material, e.to,
                          e.area, e.at.x, e.at.y);
            s += buf;
        }
        say(s);
        return 1;
    } catch (const std::exception &e) {
        say(std::string("error: ") + e.what());
        return -1;
    }
}
