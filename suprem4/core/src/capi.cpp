// C entry point for SUPREM (WP6).  SUPREM writes its current structure at the
// start of every oxidation step, calls this, and reloads the result when a
// new mesh was written.
//
// Division of work: SUPREM computes velocities and moves points inside its
// time step (the motion is part of its moving-mesh diffusion solve); the core
// owns the topology.  SUPREM's own node insertion / removal / loop removal is
// switched off while the core is active, and the core regenerates the mesh
// when
//   - the geometry has an event (thin / crumb / gap / cavity of size delta),
//   - SUPREM's motion left crossing boundaries,
//   - a triangle was compressed below 30 % of the area it had when the core
//     built it (identified by its three points), or has grown so that its
//     longest edge exceeds twice the target size there,
//   - (first call) the worst box-method coupling is below min_coupling.
// The output is made valid in SUPREM's own representation: coordinates are
// rounded to the float cm values SUPREM will hold and any triangle that is
// flat or inverted there is repaired by edge flips evaluated on those floats,
// and the result must pass check_suprem (suprem_check.hpp) before it is written.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>

#include "suprem_core/brep.hpp"
#include "suprem_core/evolve.hpp"
#include "suprem_core/locate.hpp"
#include "suprem_core/mesh.hpp"
#include "suprem_core/remesh.hpp"
#include "suprem_core/str_file.hpp"
#include "suprem_core/suprem_check.hpp"
#include "suprem_core/transfer.hpp"

using namespace suprem_core;

namespace {

double env(const char *name, double dflt)
{
    const char *v = std::getenv(name);
    return v && *v ? std::atof(v) : dflt;
}

using Key = std::array<int, 3>;
Key key_of(const Triangle &t)
{
    Key k{t.v[0], t.v[1], t.v[2]};
    std::sort(k.begin(), k.end());
    return k;
}

// exact orientation sign in SUPREM's float cm coordinates (its own predicate)
int orient_float(const Mesh &m, int a, int b, int c)
{
    float fa[2] = {suprem_float(m.point[a].x), suprem_float(m.point[a].y)};
    float fb[2] = {suprem_float(m.point[b].x), suprem_float(m.point[b].y)};
    float fc[2] = {suprem_float(m.point[c].x), suprem_float(m.point[c].y)};
    return orient2d_float_sign(fa, fb, fc);
}

uint64_t dkey(int a, int b) { return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b; }

void rebuild_neighbours(Mesh &m, const std::map<std::pair<int, int>, int> &code)
{
    std::unordered_map<uint64_t, int> dir;
    for (int t = 0; t < (int)m.triangle.size(); ++t)
        for (int k = 0; k < 3; ++k) dir[dkey(m.triangle[t].v[k], m.triangle[t].v[(k + 1) % 3])] = t;
    for (Triangle &t : m.triangle)
        for (int j = 0; j < 3; ++j) {
            int a = t.v[(j + 1) % 3], b = t.v[(j + 2) % 3];
            auto it = dir.find(dkey(b, a));
            if (it != dir.end()) { t.nb[j] = it->second; continue; }
            auto c = code.find({std::min(a, b), std::max(a, b)});
            t.nb[j] = c == code.end() ? BC_OFFSET : c->second;
        }
}

// Round to SUPREM's floats, then flip interior edges of triangles that are
// not positive there.  Returns the number of triangles still not positive.
int float_repair(Mesh &m, int &flips)
{
    for (Point &p : m.point) {
        p.x = (double)(float)(1e-4 * p.x) * 1e4;
        p.y = (double)(float)(1e-4 * p.y) * 1e4;
    }
    std::map<std::pair<int, int>, int> code;
    for (const Triangle &t : m.triangle)
        for (int j = 0; j < 3; ++j)
            if (t.nb[j] < 0) {
                int a = t.v[(j + 1) % 3], b = t.v[(j + 2) % 3];
                code[{std::min(a, b), std::max(a, b)}] = t.nb[j];
            }
    auto bad_count = [&]() {
        int n = 0;
        for (const Triangle &t : m.triangle) n += orient_float(m, t.v[0], t.v[1], t.v[2]) <= 0;
        return n;
    };
    for (int pass = 0; pass < 200; ++pass) {
        bool changed = false;
        for (int t = 0; t < (int)m.triangle.size() && !changed; ++t) {
            const Triangle tr = m.triangle[t];
            if (orient_float(m, tr.v[0], tr.v[1], tr.v[2]) > 0) continue;
            for (int j = 0; j < 3 && !changed; ++j) {
                int u = tr.nb[j];
                if (u < 0 || m.triangle[u].region != tr.region) continue;
                int c = tr.v[j], a = tr.v[(j + 1) % 3], b = tr.v[(j + 2) % 3];
                const Triangle tu = m.triangle[u];
                int d = -1;
                for (int k = 0; k < 3; ++k)
                    if (tu.v[k] != a && tu.v[k] != b) d = tu.v[k];
                if (d < 0) continue;
                if (orient_float(m, c, a, d) <= 0 || orient_float(m, c, d, b) <= 0) continue;
                Triangle n1 = tr, n2 = tu;
                n1.v[0] = c; n1.v[1] = a; n1.v[2] = d;
                n2.v[0] = c; n2.v[1] = d; n2.v[2] = b;
                m.triangle[t] = n1;
                m.triangle[u] = n2;
                rebuild_neighbours(m, code);
                ++flips;
                changed = true;
            }
        }
        if (!changed) break;
    }
    return bad_count();
}

// state kept between calls (the core runs inside SUPREM's process)
struct State {
    bool have_reference = false;
    Mesh reference;                    // grid at the start of oxidation: target sizes
    std::unique_ptr<SizeField> target; // size field of the reference (no solution part)
    std::map<Key, double> built_area;  // triangles of the last mesh the core wrote
};
State &state()
{
    static State s;
    return s;
}

}  // namespace

// Returns 1 when a new mesh was written to out_path, 0 when nothing is
// needed, < 0 on error (msg holds the reason).
extern "C" int suprem_core_oxidation_check(const char *in_path, const char *out_path, double delta_um,
                                           double min_coupling, int force, char *msg, int msglen)
{
    auto say = [&](const std::string &s) {
        if (msg && msglen > 0) std::snprintf(msg, (size_t)msglen, "%s", s.c_str());
    };
    try {
        State &st = state();
        if (const char *dump = std::getenv("SUPREM4_CORE_DUMP"); dump && *dump) {   // test data: every input
            static int n = 0;
            char name[1024];
            std::snprintf(name, sizeof name, "%s/core_in_%04d.str", dump, n++);
            write_str(read_str(in_path), name);
        }
        Mesh m = from_str(read_str(in_path));
        if (!st.have_reference) {
            st.reference = m;
            st.have_reference = true;
            RemeshOptions d;
            st.target = std::make_unique<SizeField>(st.reference, SizeField::AREA, env("SUPREM4_CORE_SCALE", d.scale), d.hmin,
                                                    env("SUPREM4_CORE_HMAX", d.hmax));
        }
        bool valid;
        {
            std::vector<std::string> pr;
            Brep b = extract_brep(m, pr);
            valid = pr.empty() && check_brep(m, b, 1).empty();
        }
        std::vector<Shape> shapes = mesh_shapes_union(m);
        RepairOptions ro;
        ro.delta = delta_um;
        std::vector<RepairEvent> ev = repair_shapes(shapes, ro);
        std::string why;
        if (!valid) why += "crossing boundaries; ";
        if (!ev.empty()) why += std::to_string(ev.size()) + " events; ";
        if (!st.built_area.empty()) {
            int squeezed = 0, stretched = 0;
            for (size_t t = 0; t < m.triangle.size(); ++t) {
                auto it = st.built_area.find(key_of(m.triangle[t]));
                if (it == st.built_area.end()) continue;
                double a = std::fabs(orient2(m, (int)t)) * 0.5;
                squeezed += a < 0.3 * it->second;
                if (a > 2.0 * it->second) {   // grown since built and now too coarse for the target size?
                    const Triangle &tr = m.triangle[t];
                    double lmax = 0;
                    Point c{0, 0};
                    for (int k = 0; k < 3; ++k) {
                        const Point &p = m.point[tr.v[k]], &q = m.point[tr.v[(k + 1) % 3]];
                        lmax = std::max(lmax, std::hypot(q.x - p.x, q.y - p.y));
                        c.x += p.x / 3; c.y += p.y / 3;
                    }
                    stretched += lmax > 2.0 * (*st.target)(c);
                }
            }
            if (squeezed) why += std::to_string(squeezed) + " compressed triangles; ";
            if (stretched) why += std::to_string(stretched) + " stretched triangles; ";
        } else {
            // first call: the core owns the topology from the start of the
            // oxidation (SUPREM's own node insertion/removal is off)
            why += "first step; ";
            (void)min_coupling;
        }
        if (force) why += "forced; ";
        if (why.empty()) {
            say("ok");
            return 0;
        }
        RemeshOptions o;
        o.hmax = env("SUPREM4_CORE_HMAX", o.hmax);
        o.scale = env("SUPREM4_CORE_SCALE", o.scale);
        o.float_cm_grid = true;
        // features below the event size are unresolved by definition: boundary
        // simplification at delta/10 removes the sub-angstrom debris of the
        // morphology operations (geometry kept to 1 A for delta = 1 nm)
        o.coarsen_tol = env("SUPREM4_CORE_COARSEN_TOL", delta_um / 10);
        SizeField h(st.reference, SizeField::AREA, o.scale, o.hmin, o.hmax);
        h.adapt_to_solution(env("SUPREM4_CORE_SOLUTION_DECADES", o.solution_decades));
        RemeshReport rr;
        std::vector<std::string> notes;
        int cleaned = clean_shapes(shapes, ro);
        Mesh fresh = mesh_from_shapes(shapes, m, boundary_codes_of(m), &h, o, rr, &notes);
        TransferReport tr;
        transfer_fields(m, fresh, tr);
        int flips = 0, still = float_repair(fresh, flips);
        if (still) {
            say(why + "rebuilt mesh has " + std::to_string(still) + " triangles not positive in SUPREM float coordinates");
            return -3;
        }
        auto bad = check_invariants(fresh, 3);
        auto gas = check_gas_nodes(fresh, 3);
        bad.insert(bad.end(), gas.begin(), gas.end());
        if (!bad.empty()) {
            say(why + "rebuilt mesh invalid: " + bad.front());
            return -2;
        }
        // WP6a: never hand SUPREM a mesh it cannot reload or move
        SupremCheck sc = check_suprem(fresh, 3);
        if (!sc.ok()) {
            say(why + "rebuilt mesh violates the SUPREM contract: " + (sc.messages.empty() ? std::string("?") : sc.messages.front()));
            return -4;
        }
        write_str(to_str(fresh), out_path);
        st.built_area.clear();
        for (size_t t = 0; t < fresh.triangle.size(); ++t)
            st.built_area[key_of(fresh.triangle[t])] = std::fabs(orient2(fresh, (int)t)) * 0.5;
        char buf[512];
        std::snprintf(buf, sizeof buf, "%sslivers cleaned %d, worst coupling -> %.3g, points %zu -> %zu, float flips %d, inventory error %.2g",
                      why.c_str(), cleaned, measure_quality(fresh).worst_coupling, m.point.size(), fresh.point.size(), flips,
                      tr.max_rel_final);
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
