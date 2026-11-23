// core_update: WP6b Lagrangian update of a SUPREM .str (remesh mode "update").
//   core_update IN.str OUT.str [--reference REF.str] [--scale S] [--hmax H]
//               [--near F] [--short F] [--tol T] [--long F] [--grow F]
//               [--solution-decades D] [--float] [--box]
// --retriangulate: full CDT (no old triangles kept); --keep-coupling C.
// --events DELTA: geometric events of size DELTA (um) applied locally first (WP6c, as the hook).
// --box: restore the box-volume inventory (as the hook) instead of P1.
// The size fields come from REF.str (default: IN.str itself), as in
// the oxidation hook (grid at the start of the oxidation).  Prints the
// counts of kept / removed / inserted points, the quality of both meshes and
// the checks:
//   invariants, brep, geometry (B-rep within tol), gas-nodes, suprem (WP6a),
//   retained  every node at a point outside the changed set has exactly its
//             old values (bit for bit),
//   inventory conserved inventories (P1, per material) restored to 1e-10.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <memory>
#include <string>
#include <vector>

#include "suprem_core/brep.hpp"
#include "suprem_core/evolve.hpp"
#include "suprem_core/locate.hpp"
#include "suprem_core/mesh.hpp"
#include "suprem_core/remesh.hpp"
#include "suprem_core/str_file.hpp"
#include "suprem_core/suprem_check.hpp"
#include "suprem_core/transfer.hpp"

using namespace suprem_core;

int main(int argc, char **argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: core_update IN.str OUT.str [options]\n"); return 2; }
    UpdateOptions o;
    std::string ref_path;
    bool box = false;
    double events = 0;
    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> const char * {
            if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", a.c_str()); std::exit(2); }
            return argv[++i];
        };
        if (a == "--reference") ref_path = next();
        else if (a == "--scale") o.scale = std::atof(next());
        else if (a == "--hmax") o.hmax = std::atof(next());
        else if (a == "--grow") o.grow = std::atof(next());
        else if (a == "--solution-decades") o.solution_decades = std::atof(next());
        else if (a == "--near") o.near = std::atof(next());
        else if (a == "--short") o.short_edge = std::atof(next());
        else if (a == "--tol") o.tol = std::atof(next());
        else if (a == "--long") o.long_edge = std::atof(next());
        else if (a == "--float") o.float_cm_grid = true;
        else if (a == "--timing") o.timing = true;
        else if (a == "--box") box = true;
        else if (a == "--events") events = std::atof(next());
        else if (a == "--retriangulate") o.keep_triangles = false;
        else if (a == "--keep-coupling") o.keep_coupling = std::atof(next());
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    try {
        Mesh in = from_str(read_str(argv[1]));
        Mesh ref = ref_path.empty() ? in : from_str(read_str(ref_path.c_str()));
        std::printf("input  %s\n", format_quality(measure_quality(in)).c_str());
        auto t0 = std::chrono::steady_clock::now();
        Mesh geom = in;
        LocalRepairReport lr;
        if (events > 0) {
            RepairOptions ro;
            ro.delta = events;
            std::vector<RepairEvent> ev;
            bool ok = apply_events_locally(geom, ro, lr, &ev);
            for (auto &e : ev)
                std::printf("event %s %d->%d %.3g um2 at (%.5f, %.5f)\n", e.kind.c_str(), e.material, e.to, e.area, e.at.x, e.at.y);
            std::printf("local repair %s: relabelled %d removed %d filled %d, points dropped %d split %d needles %d contacts %d\n", ok ? "ok" : "FAILED",
                        lr.relabeled, lr.removed, lr.filled, lr.dropped_points, lr.split, lr.needles, lr.contacts);
            for (auto &n : lr.notes) std::printf("NOTE %s\n", n.c_str());
            if (!ok) return 1;
            for (auto &x : check_invariants(geom, 5)) std::printf("repaired invariants %s\n", x.c_str());
            std::vector<std::string> pr;
            Brep b = extract_brep(geom, pr);
            for (auto &x : pr) std::printf("repaired brep %s\n", x.c_str());
            for (auto &x : check_brep(geom, b, 5)) std::printf("repaired brep %s\n", x.c_str());
            std::vector<Shape> after = mesh_shapes_union(geom);
            std::printf("events after the local repair %zu\n", repair_shapes(after, ro).size());
        }
        UpdateReport rep;
        Mesh out = update_mesh(geom, ref, o, rep);
        if (events > 0 && !lr.touched.empty())
            for (size_t i = 0; i < out.point.size(); ++i)
                if (rep.origin[i] >= 0 && lr.touched[rep.origin[i]]) rep.changed[i] = 1;
        auto t1 = std::chrono::steady_clock::now();
        TransferReport tr;
        update_fields(in, out, rep.changed, tr, box);
        if (o.timing) std::printf("time fields %.2f s\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count());
        double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        int nchanged = 0;
        for (char c : rep.changed) nchanged += c;
        std::printf("output %s\n", format_quality(measure_quality(out)).c_str());
        std::printf("update %.2f s | points %zu -> %zu: kept %d, removed interior %d boundary %d, inserted interior %d "
                    "boundary %d encroach %d | new triangles %d of %zu | changed points %d\n",
                    sec, in.point.size(), out.point.size(), rep.kept, rep.removed_interior, rep.removed_boundary,
                    rep.inserted_interior, rep.inserted_boundary, rep.encroach_splits, rep.new_triangles,
                    out.triangle.size(), nchanged);
        std::printf("fields | nodes exact %d interpolated %d gas %d | inventory error before restoration %.3g, final %.3g "
                    "(box-volume inventory change %.3g)\n",
                    tr.exact, tr.interpolated, tr.gas, tr.max_rel_interpolated, tr.max_rel_final, tr.max_rel_box);
        for (auto &n : rep.notes) std::printf("NOTE %s\n", n.c_str());
        for (auto &n : tr.notes) std::printf("NOTE %s\n", n.c_str());
        int failed = 0;
        auto section = [&](const char *name, const std::vector<std::string> &v) {
            for (auto &x : v) std::printf("%s %s\n", name, x.c_str());
            std::printf("%s %s\n", name, v.empty() ? "ok" : "FAILED");
            failed += !v.empty();
        };
        section("invariants", check_invariants(out));
        std::vector<std::string> pr;
        Brep b = extract_brep(out, pr);
        auto bb = check_brep(out, b);
        bb.insert(bb.begin(), pr.begin(), pr.end());
        section("brep", bb);
        section("geometry", compare_geometry(geom, out, o.tol));   // the events changed the input on purpose
        section("gas-nodes", check_gas_nodes(out));
        {
            SupremCheck sc = check_suprem(out, 20);
            section("suprem", sc.ok() ? std::vector<std::string>{} : sc.messages);
        }
        {
            std::vector<std::string> r;
            int checked = 0;
            for (const Node &n : out.node) {
                if (rep.changed[n.point]) continue;
                int op = rep.origin[n.point];
                if (op >= 0 && !lr.origin.empty()) op = lr.origin[op];
                int on = op < 0 ? -1 : in.node_at(op, n.material);
                ++checked;
                if (on < 0 || in.node[on].value != n.value) {
                    if (r.size() < 10)
                        r.push_back("node at point " + std::to_string(n.point + 1) + " material " + std::to_string(n.material) +
                                    (on < 0 ? " has no old node" : " changed"));
                }
            }
            std::printf("retained %d nodes outside the changed set compared\n", checked);
            section("retained", r);
        }
        section("inventory", tr.max_rel_final <= 1e-10 ? std::vector<std::string>{}
                                                        : std::vector<std::string>{"final inventory error " +
                                                                                   std::to_string(tr.max_rel_final)});
        write_str(to_str(out), argv[2]);
        return failed ? 1 : 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 3;
    }
}
