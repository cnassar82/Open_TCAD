// core_remesh: regenerate the mesh of a SUPREM .str from its B-rep (WP3).
//   core_remesh IN.str OUT.str [--size area|min] [--scale S] [--hmin H] [--hmax H]
//               [--hgrad G] [--film-threshold T] [--film-layers N] [--no-mmg]
//               [--no-flips] [--passes N] [--encroach-rounds N] [-v]
// Prints quality of the input and output meshes, the WP3 checks
// (invariants, B-rep, geometry preservation, gas nodes) and the WP4 field
// transfer report (inventories per material and conserved species, round-
// trip interpolation error, identity transfer).  --no-transfer leaves the
// node values zero; --no-restore skips the inventory restoration.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "suprem_core/brep.hpp"
#include "suprem_core/mesh.hpp"
#include "suprem_core/remesh.hpp"
#include "suprem_core/str_file.hpp"
#include "suprem_core/transfer.hpp"

using namespace suprem_core;

int main(int argc, char **argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: core_remesh IN.str OUT.str [options]\n"); return 2; }
    RemeshOptions o;
    int worst = 0;
    bool worst_by_max = false, do_transfer = true, restore = true;
    Restore rmode = Restore::WEIGHTED;
    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> const char * {
            if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", a.c_str()); std::exit(2); }
            return argv[++i];
        };
        if (a == "--size") o.size_mode = next();
        else if (a == "--scale") o.scale = std::atof(next());
        else if (a == "--hmin") o.hmin = std::atof(next());
        else if (a == "--hmax") o.hmax = std::atof(next());
        else if (a == "--coarsen-tol") o.coarsen_tol = std::atof(next());
        else if (a == "--solution-decades") o.solution_decades = std::atof(next());
        else if (a == "--hgrad") o.hgrad = std::atof(next());
        else if (a == "--film-threshold") o.film_threshold = std::atof(next());
        else if (a == "--film-layers") o.film_layers = std::atoi(next());
        else if (a == "--passes") o.mmg_passes = std::atoi(next());
        else if (a == "--encroach-rounds") o.encroach_rounds = std::atoi(next());
        else if (a == "--no-mmg") o.mmg = false;
        else if (a == "--no-flips") o.flips = false;
        else if (a == "-v") o.verbose++;
        else if (a == "--no-transfer") do_transfer = false;
        else if (a == "--no-restore") restore = false;
        else if (a == "--restore-scale") rmode = Restore::SCALE;
        else if (a == "--worst") worst = std::atoi(next());
        else if (a == "--worst-max") { worst = std::atoi(next()); worst_by_max = true; }
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    try {
        Mesh in = from_str(read_str(argv[1]));
        std::printf("input  %s\n", format_quality(measure_quality(in)).c_str());
        auto t0 = std::chrono::steady_clock::now();
        RemeshReport rep;
        Mesh out = remesh(in, o, rep);
        double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("output %s\n", format_quality(measure_quality(out)).c_str());
        std::printf("remesh %.2f s | coarsened %d | chain segments %d | film pairs %d layer points %d fixed film triangles %d | "
                    "encroach splits %d | flips %d | mmg failures %d\n",
                    sec, rep.coarsened, rep.chain_points, rep.film_pairs, rep.film_points, rep.film_triangles, rep.encroach_splits,
                    rep.flips, rep.mmg_failures);
        for (auto &n : rep.notes) std::printf("NOTE %s\n", n.c_str());
        if (worst > 0) {   // worst triangles by min angle: location, material, edge lengths, boundary edges
            std::vector<std::pair<double, int>> q;
            for (int t = 0; t < (int)out.triangle.size(); ++t) {
                double mn = 180, mx = 0;
                for (int k = 0; k < 3; ++k) {
                    const Point &a = out.point[out.triangle[t].v[k]], &b = out.point[out.triangle[t].v[(k + 1) % 3]],
                                &c = out.point[out.triangle[t].v[(k + 2) % 3]];
                    double ux = b.x - a.x, uy = b.y - a.y, vx = c.x - a.x, vy = c.y - a.y;
                    double ang = std::atan2(std::fabs(ux * vy - uy * vx), ux * vx + uy * vy) * 180 / M_PI;
                    mn = std::min(mn, ang); mx = std::max(mx, ang);
                }
                q.push_back({worst_by_max ? -mx : mn, t});
            }
            std::sort(q.begin(), q.end());
            for (int i = 0; i < worst && i < (int)q.size(); ++i) {
                const Triangle &t = out.triangle[q[i].second];
                std::printf("WORST %.3f deg mat %d at (%.5f, %.5f) edges", std::fabs(q[i].first), out.material_of(q[i].second),
                            out.point[t.v[0]].x, out.point[t.v[0]].y);
                for (int j = 0; j < 3; ++j) {
                    const Point &a = out.point[t.v[(j + 1) % 3]], &b = out.point[t.v[(j + 2) % 3]];
                    bool bnd = t.nb[j] < 0 || out.triangle[t.nb[j]].region != t.region;
                    std::printf(" %.2e%s", std::hypot(b.x - a.x, b.y - a.y), bnd ? "B" : "");
                }
                std::printf("\n");
            }
        }
        int failed = 0;
        std::vector<std::string> transfer_bad;
        if (do_transfer) {
            // identity: transferring onto the unchanged mesh must reproduce it exactly
            Mesh same = in;
            for (Node &n : same.node) std::fill(n.value.begin(), n.value.end(), -1.0);
            TransferReport ri;
            transfer_fields(in, same, ri, restore, true, rmode);
            for (size_t n = 0; n < in.node.size(); ++n)
                if (same.node[n].value != in.node[n].value) {
                    std::string d;
                    for (size_t k = 0; k < in.species.size(); ++k)
                        if (same.node[n].value[k] != in.node[n].value[k])
                            d += " imp " + std::to_string(in.species[k]) + ": " + std::to_string(in.node[n].value[k]) +
                                 " -> " + std::to_string(same.node[n].value[k]);
                    transfer_bad.push_back("identity transfer changed node " + std::to_string(n) + " (mat " +
                                           std::to_string(in.node[n].material) + ")" + d);
                    break;
                }
            auto t1 = std::chrono::steady_clock::now();
            TransferReport tr;
            transfer_fields(in, out, tr, restore, true, rmode);
            double tsec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
            std::printf("transfer %.2f s | nodes exact %d interpolated %d gas %d clipped %d | inventory error "
                        "interpolated %.3g final %.3g (box-volume inventory change %.3g)\n", tsec, tr.exact, tr.interpolated,
                        tr.gas, tr.outside, tr.max_rel_interpolated, tr.max_rel_final, tr.max_rel_box);
            for (auto &r : tr.rows)
                if (o.verbose || std::fabs(r.interpolated_inv - r.old_inv) > 1e-3 * std::fabs(r.old_inv))
                    std::printf("  inventory mat %d imp %d: %.6e -> %.6e (interp %+.2e) -> %.6e\n", r.material,
                                r.impurity, r.old_inv, r.interpolated_inv, (r.interpolated_inv - r.old_inv) / r.old_inv,
                                r.new_inv);
            double rt = 0;
            for (size_t s = 0; s < in.species.size(); ++s)
                if (conserved_species(in.species[s])) rt = std::max(rt, roundtrip_error(in, out, (int)s));
            std::printf("transfer roundtrip L1 error %.3g\n", rt);
            for (auto &n : tr.notes) std::printf("NOTE %s\n", n.c_str());
            if (restore && tr.max_rel_final > 1e-6)
                transfer_bad.push_back("inventory error after restoration " + std::to_string(tr.max_rel_final));
        }
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
        section("geometry", compare_geometry(in, out, o.coarsen_tol));
        section("gas-nodes", check_gas_nodes(out));
        if (do_transfer) section("transfer", transfer_bad);
        write_str(to_str(out), argv[2]);
        return failed ? 1 : 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 3;
    }
}
