// C entry point for SUPREM (WP6).  SUPREM writes its current structure at the
// start of every oxidation step, calls this, and reloads the result when a
// new mesh was written.
//
// Division of work: SUPREM computes velocities and moves points inside its
// time step (the motion is part of its moving-mesh diffusion solve); the core
// owns the topology.  SUPREM's own node insertion / removal / loop removal is
// switched off while the core is active.  Policy (WP6d):
//   - every call: the Lagrangian update (update_mesh / update_fields, 6b):
//     points and their values are kept, points are removed / inserted only
//     where the mesh is too fine / too coarse for the size of the grid at the
//     start of the oxidation, and the inventories are restored on the changed
//     points only; nothing is written when nothing changed;
//   - full regeneration (WP3) from the repaired shapes when the geometry has
//     an event (thin / crumb / gap / cavity of size delta; 6c will make this
//     local), SUPREM's motion left crossing boundaries, SUPREM asked for it
//     (force: collision prediction), or the update fails;
//   - SUPREM4_CORE_REGENERATE=1: the old policy (regeneration on the first
//     call and when triangles were compressed below 30 %), for comparison.
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
    int calls = 0;
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
        }
        ++st.calls;
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
        (void)min_coupling;
        const bool regen_only = env("SUPREM4_CORE_REGENERATE", 0) != 0;   // old policy (comparison)
        std::string why;
        if (!valid) why += "crossing boundaries; ";
        if (regen_only && !ev.empty()) why += std::to_string(ev.size()) + " events; ";
        if (regen_only && force) why += "forced; ";
        if (regen_only && st.calls == 1) why += "first step; ";
        if (regen_only && st.calls > 1) {
            int squeezed = 0;
            for (size_t t = 0; t < m.triangle.size(); ++t) {
                auto it = st.built_area.find(key_of(m.triangle[t]));
                squeezed += it != st.built_area.end() && std::fabs(orient2(m, (int)t)) * 0.5 < 0.3 * it->second;
            }
            if (squeezed) why += std::to_string(squeezed) + " compressed triangles; ";
        }

        char buf[1024];
        // validity in SUPREM's representation, then write; returns the check failure or ""
        auto finish = [&](Mesh &fresh, int &flips) -> std::string {
            int still = float_repair(fresh, flips);
            if (still) return std::to_string(still) + " triangles not positive in SUPREM float coordinates";
            auto bad = check_invariants(fresh, 3);
            auto gas = check_gas_nodes(fresh, 3);
            bad.insert(bad.end(), gas.begin(), gas.end());
            if (!bad.empty()) return "invalid: " + bad.front();
            SupremCheck sc = check_suprem(fresh, 3);   // WP6a: never hand SUPREM a mesh it cannot reload or move
            if (!sc.ok()) return "violates the SUPREM contract: " + (sc.messages.empty() ? std::string("?") : sc.messages.front());
            write_str(to_str(fresh), out_path);
            return "";
        };
        // WP6c: events applied locally to the mesh (no point created or moved)
        Mesh geom;
        LocalRepairReport lr;
        std::string local;
        bool have_geom = false;
        if (why.empty() && !regen_only) {   // also without events: needles below the shapes' resolution
            geom = m;
            bool ok = apply_events_locally(geom, ro, lr);
            if (ok && !lr.relabeled && !lr.removed && !lr.filled && !lr.needles && !lr.contacts) {
                // nothing changed: the plain update
            } else if (ok) {
                std::vector<std::string> pr;
                Brep b = extract_brep(geom, pr);
                if (pr.empty() && check_brep(geom, b, 1).empty()) {
                    have_geom = true;
                    std::snprintf(buf, sizeof buf, "%zu events local (relabelled %d, removed %d, filled %d, needles %d, contacts %d, split %d); ",
                                  ev.size(), lr.relabeled, lr.removed, lr.filled, lr.needles, lr.contacts, lr.split);
                    local = buf;
                } else {
                    why += std::to_string(ev.size()) + " events, local repair left an invalid B-rep; ";
                }
            } else {
                why += std::to_string(ev.size()) + " events, local repair failed (" +
                       (lr.notes.empty() ? std::string("?") : lr.notes.back()) + "); ";
            }
        }
        // WP6d: the Lagrangian update (6b) on every regular step
        std::string update_failure;
        if (why.empty() && !regen_only) {
            const Mesh &g0 = have_geom ? geom : m;
            UpdateOptions uo;
            uo.hmax = env("SUPREM4_CORE_HMAX", uo.hmax);
            uo.scale = env("SUPREM4_CORE_SCALE", uo.scale);
            uo.solution_decades = env("SUPREM4_CORE_SOLUTION_DECADES", uo.solution_decades);
            uo.float_cm_grid = true;
            if (force) uo.keep_triangles = false;   // SUPREM predicts squeezed elements: full Delaunay triangulation
            try {
                UpdateReport ur;
                Mesh fresh = update_mesh(g0, st.reference, uo, ur);
                if (have_geom)   // vertices of the locally repaired triangles
                    for (size_t i = 0; i < fresh.point.size(); ++i)
                        if (ur.origin[i] >= 0 && lr.touched[ur.origin[i]]) ur.changed[i] = 1;
                int nchanged = 0;
                for (char c : ur.changed) nchanged += c;
                if (!have_geom && nchanged == 0 && fresh.point.size() == m.point.size()) {
                    say("ok");
                    return 0;
                }
                TransferReport tr;
                update_fields(m, fresh, ur.changed, tr, true);   // from the original mesh: conserve what SUPREM conserves
                int flips = 0;
                update_failure = finish(fresh, flips);
                if (update_failure.empty()) {
                    std::snprintf(buf, sizeof buf,
                                  "%s%supdate: points %zu -> %zu (removed %d+%d, inserted %d+%d, encroach %d), new triangles %d, "
                                  "worst coupling %.3g, float flips %d, inventory error %.2g, box change %.2g",
                                  local.c_str(), force ? "forced, full triangulation; " : "", m.point.size(), fresh.point.size(), ur.removed_interior, ur.removed_boundary,
                                  ur.inserted_interior, ur.inserted_boundary, ur.encroach_splits, ur.new_triangles,
                                  measure_quality(fresh).worst_coupling, flips, tr.max_rel_final, tr.max_rel_box);
                    say(buf);
                    return 1;
                }
            } catch (const std::exception &e) {
                update_failure = e.what();
            }
            why = "update failed (" + update_failure + "); ";
        }

        // fallback / events: full regeneration (WP3) from the repaired shapes
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
        int flips = 0;
        std::string f = finish(fresh, flips);
        if (!f.empty()) {
            say(why + "rebuilt mesh " + f);
            return -2;
        }
        st.built_area.clear();
        if (regen_only)
            for (size_t t = 0; t < fresh.triangle.size(); ++t)
                st.built_area[key_of(fresh.triangle[t])] = std::fabs(orient2(fresh, (int)t)) * 0.5;
        std::snprintf(buf, sizeof buf, "%sregenerated: slivers cleaned %d, worst coupling -> %.3g, points %zu -> %zu, float flips %d, inventory error %.2g",
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
