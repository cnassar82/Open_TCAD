// Exact reader/writer for SUPREM-IV.GS structure files (.str, "ig2" format).
//
// The record layout follows mesh/ig2_meshio.c (ig2_write).  Writing a file
// that was read reproduces it byte for byte, which is the WP1 round-trip test.
#pragma once
#include <string>
#include <vector>

namespace suprem_core {

struct StrTriangle {
    int region = 0;          // 1-based region index as written
    int vertex[3] = {0, 0, 0};   // 1-based point indices
    int code[3] = {0, 0, 0};     // neighbour triangle (1-based) or boundary code (<0)
    int father = -1, offspring = -1;
};

struct StrNode {
    int point = 0;           // 0-based point index (as written by ig2_write)
    int material = 0;
    std::vector<double> value;   // one per species in StrFile::species
};

struct StrFile {
    std::string version;                 // text after "v "
    int dim = 2, nvrt = 3, nedg = 3;     // "D" record
    std::vector<std::vector<double>> coord;  // per point, dim values (um), as written
    std::vector<int> region_material;    // per region (0-based index = id-1)
    std::vector<StrTriangle> triangle;
    bool has_model = false;
    int orientation = 0;                 // "M" record
    double last_temperature = 0;
    std::vector<int> species;            // "s" record: impurity number per solution slot
    std::vector<StrNode> node;
    std::vector<int> implanted;          // "I" record
    std::vector<std::string> trailer;    // later records (process surface), kept verbatim
};

// Throws std::runtime_error with the line number on malformed input.
StrFile read_str(const std::string &path);
void write_str(const StrFile &f, const std::string &path);

}  // namespace suprem_core
