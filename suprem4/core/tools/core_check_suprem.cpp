// core_check_suprem: apply SUPREM-IV's reload/motion validity contract
// (suprem_core/suprem_check.hpp) to .str files.
//   core_check_suprem A.str [B.str ...]
// One line per file: "OK", "WARN" (SUPREM reloads it, but a strict rule
// fails) or "FAIL" (SUPREM would panic or crash) with rule counts, then the
// first reports.
// Exit status: 0 all valid, 1 some invalid, 2 read error.
#include <cstdio>
#include <exception>

#include "suprem_core/str_file.hpp"
#include "suprem_core/suprem_check.hpp"

int main(int argc, char **argv)
{
    if (argc < 2) { std::fprintf(stderr, "usage: core_check_suprem A.str [B.str ...]\n"); return 2; }
    int status = 0;
    for (int i = 1; i < argc; ++i) {
        try {
            suprem_core::Mesh m = suprem_core::from_str(suprem_core::read_str(argv[i]));
            m.rebuild_index();
            auto c = suprem_core::check_suprem(m, 5);
            std::printf("%s %s orientation %d topology %d embedding %d loops %d nodes %d | ordering %d components %d gas_added %d\n",
                        c.ok() ? "OK  " : c.reloadable() ? "WARN" : "FAIL",
                        argv[i], c.orientation, c.topology, c.embedding, c.loops, c.nodes, c.ordering, c.components, c.gas_added);
            for (auto &s : c.messages) std::printf("    %s\n", s.c_str());
            if (!c.ok()) status = 1;
        } catch (const std::exception &e) {
            std::printf("ERR  %s %s\n", argv[i], e.what());
            return 2;
        }
    }
    return status;
}
