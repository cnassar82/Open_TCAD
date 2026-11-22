#include "suprem_core/str_file.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace suprem_core {
namespace {

[[noreturn]] void fail(const std::string &path, int line, const std::string &what)
{
    throw std::runtime_error(path + ":" + std::to_string(line) + ": " + what);
}

double to_double(const std::string &s, const std::string &path, int line)
{
    char *end = nullptr;
    double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || *end) fail(path, line, "bad number '" + s + "'");
    return v;
}

int to_int(const std::string &s, const std::string &path, int line)
{
    char *end = nullptr;
    long v = std::strtol(s.c_str(), &end, 10);
    if (end == s.c_str() || *end) fail(path, line, "bad integer '" + s + "'");
    return (int)v;
}

std::string fmt_g17(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

}  // namespace

StrFile read_str(const std::string &path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    StrFile f;
    std::string text;
    int line = 0;
    bool seen_I = false;
    while (std::getline(in, text)) {
        ++line;
        if (seen_I) { f.trailer.push_back(text); continue; }
        std::istringstream ss(text);
        std::string tag;
        if (!(ss >> tag)) continue;
        std::vector<std::string> tok;
        for (std::string t; ss >> t;) tok.push_back(t);
        if (tag == "v") {
            f.version = text.size() > 2 ? text.substr(2) : "";
        } else if (tag == "D") {
            if (tok.size() != 3) fail(path, line, "D record");
            f.dim = to_int(tok[0], path, line);
            f.nvrt = to_int(tok[1], path, line);
            f.nedg = to_int(tok[2], path, line);
        } else if (tag == "c") {
            if ((int)tok.size() != f.dim + 2) fail(path, line, "c record");
            int id = to_int(tok[0], path, line);
            if (id != (int)f.coord.size() + 1) fail(path, line, "points out of order");
            std::vector<double> c;
            for (int i = 0; i < f.dim; ++i) c.push_back(to_double(tok[1 + i], path, line));
            f.coord.push_back(c);
        } else if (tag == "r") {
            if (tok.size() != 2) fail(path, line, "r record");
            int id = to_int(tok[0], path, line);
            if (id != (int)f.region_material.size() + 1) fail(path, line, "regions out of order");
            f.region_material.push_back(to_int(tok[1], path, line));
        } else if (tag == "t") {
            if ((int)tok.size() != 2 + f.nvrt + f.nedg + 2) fail(path, line, "t record");
            StrTriangle t;
            int id = to_int(tok[0], path, line);
            if (id != (int)f.triangle.size() + 1) fail(path, line, "triangles out of order");
            t.region = to_int(tok[1], path, line);
            for (int i = 0; i < 3; ++i) t.vertex[i] = to_int(tok[2 + i], path, line);
            for (int i = 0; i < 3; ++i) t.code[i] = to_int(tok[5 + i], path, line);
            t.father = to_int(tok[8], path, line);
            t.offspring = to_int(tok[9], path, line);
            f.triangle.push_back(t);
        } else if (tag == "M") {
            if (tok.size() != 2) fail(path, line, "M record");
            f.has_model = true;
            f.orientation = to_int(tok[0], path, line);
            f.last_temperature = to_double(tok[1], path, line);
        } else if (tag == "s") {
            int n = to_int(tok.at(0), path, line);
            if ((int)tok.size() != n + 1) fail(path, line, "s record");
            for (int i = 0; i < n; ++i) f.species.push_back(to_int(tok[1 + i], path, line));
        } else if (tag == "n") {
            if (tok.size() != 2 + f.species.size()) fail(path, line, "n record");
            StrNode n;
            n.point = to_int(tok[0], path, line);
            n.material = to_int(tok[1], path, line);
            for (size_t i = 0; i < f.species.size(); ++i) n.value.push_back(to_double(tok[2 + i], path, line));
            f.node.push_back(n);
        } else if (tag == "I") {
            int n = to_int(tok.at(0), path, line);
            if ((int)tok.size() != n + 1) fail(path, line, "I record");
            for (int i = 0; i < n; ++i) f.implanted.push_back(to_int(tok[1 + i], path, line));
            seen_I = true;
        } else {
            fail(path, line, "unknown record '" + tag + "'");
        }
    }
    return f;
}

void write_str(const StrFile &f, const std::string &path)
{
    std::FILE *out = std::fopen(path.c_str(), "w");
    if (!out) throw std::runtime_error("cannot write " + path);
    std::fprintf(out, "v %s\n", f.version.c_str());
    std::fprintf(out, "D %d %d %d\n", f.dim, f.nvrt, f.nedg);
    for (size_t p = 0; p < f.coord.size(); ++p) {
        std::fprintf(out, "c %zu", p + 1);
        for (double c : f.coord[p]) std::fprintf(out, " %s", fmt_g17(c).c_str());
        std::fprintf(out, "  0\n");
    }
    for (size_t r = 0; r < f.region_material.size(); ++r)
        std::fprintf(out, "r %zu   %d\n", r + 1, f.region_material[r]);
    for (size_t t = 0; t < f.triangle.size(); ++t) {
        const StrTriangle &tr = f.triangle[t];
        std::fprintf(out, "t %zu %d ", t + 1, tr.region);
        for (int v : tr.vertex) std::fprintf(out, "%d ", v);
        for (int c : tr.code) std::fprintf(out, "%d ", c);
        std::fprintf(out, "%d %d ", tr.father, tr.offspring);
        std::fprintf(out, "\n");
    }
    if (f.has_model) std::fprintf(out, "M %d %e\n", f.orientation, f.last_temperature);
    std::fprintf(out, "s %zu  ", f.species.size());
    for (int s : f.species) std::fprintf(out, " %d", s);
    std::fprintf(out, "\n");
    for (const StrNode &n : f.node) {
        std::fprintf(out, "n %d   %d  ", n.point, n.material);
        for (double v : n.value) std::fprintf(out, " %s", fmt_g17(v).c_str());
        std::fprintf(out, "\n");
    }
    std::fprintf(out, "I %zu  ", f.implanted.size());
    for (int i : f.implanted) std::fprintf(out, " %d", i);
    std::fprintf(out, "\n");
    for (const std::string &l : f.trailer) std::fprintf(out, "%s\n", l.c_str());
    std::fclose(out);
}

}  // namespace suprem_core
