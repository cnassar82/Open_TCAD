// core_str: read a SUPREM .str file into the geometry core and report.
//   core_str IN.str                  invariants + summary
//   core_str IN.str --write OUT.str  round trip through the semantic Mesh
#include <cstdio>
#include <cstring>
#include <exception>
#include <map>

#include "suprem_core/brep.hpp"
#include "suprem_core/mesh.hpp"
#include "suprem_core/remesh.hpp"
#include "suprem_core/str_file.hpp"
#include "suprem_core/transfer.hpp"

int main(int argc, char **argv)
{
    if (argc < 2) { std::fprintf(stderr, "usage: core_str IN.str [--write OUT.str]\n"); return 2; }
    const char *out = nullptr;
    for (int i = 2; i < argc; ++i)
        if (!std::strcmp(argv[i], "--write") && i + 1 < argc) out = argv[++i];
    try {
        suprem_core::Mesh m = suprem_core::from_str(suprem_core::read_str(argv[1]));
        std::map<int, int> tri_per_mat;
        for (size_t t = 0; t < m.triangle.size(); ++t) tri_per_mat[m.material_of((int)t)]++;
        std::printf("points %zu triangles %zu regions %zu nodes %zu species %zu\n", m.point.size(),
                    m.triangle.size(), m.region_material.size(), m.node.size(), m.species.size());
        std::printf("triangles by material:");
        for (auto &kv : tri_per_mat) std::printf(" %d:%d", kv.first, kv.second);
        std::printf("\n");
        auto bad = suprem_core::check_invariants(m);
        for (auto &b : bad) std::printf("INVARIANT %s\n", b.c_str());
        for (auto &g : suprem_core::check_gas_nodes(m)) bad.push_back("gas: " + g);
        for (auto &b2 : bad) if (b2.rfind("gas: ", 0) == 0) std::printf("INVARIANT %s\n", b2.c_str());
        std::printf("invariants %s\n", bad.empty() ? "ok" : "FAILED");
        std::vector<std::string> problems;
        suprem_core::Brep b = suprem_core::extract_brep(m, problems);
        auto bb = suprem_core::check_brep(m, b);
        bb.insert(bb.begin(), problems.begin(), problems.end());
        int interface = 0, exposed = 0, back = 0, side = 0;
        for (auto &c : b.chain) {
            if (c.right >= 0) interface++;
            else if (c.right == suprem_core::EXPOSED) exposed++;
            else if (c.right == suprem_core::BACKSIDE) back++;
            else side++;
        }
        std::printf("brep loops %zu chains %zu (interface %d exposed %d backside %d other %d)\n",
                    b.loop.size(), b.chain.size(), interface, exposed, back, side);
        for (auto &x : bb) std::printf("BREP %s\n", x.c_str());
        std::printf("brep %s\n", bb.empty() ? "ok" : "FAILED");
        if (!bb.empty()) bad.push_back("brep");
        auto vol = suprem_core::node_volumes(m);
        for (size_t s = 0; s < m.species.size(); ++s) {
            if (!suprem_core::conserved_species(m.species[s])) continue;
            auto inv = suprem_core::inventory(m, vol, (int)s);
            std::printf("inventory imp %d:", m.species[s]);
            for (size_t mat = 0; mat < inv.size(); ++mat)
                if (inv[mat] != 0) std::printf(" mat%zu=%.9e", mat, inv[mat]);
            std::printf("\n");
        }
        if (out) suprem_core::write_str(suprem_core::to_str(m), out);
        return bad.empty() ? 0 : 1;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 3;
    }
}
