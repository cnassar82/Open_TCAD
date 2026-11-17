#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/etch_corner.c" \
    "$root/src/refine/sp_edge.c" "$root/src/refine/dbase.c" \
    "${SUPREM_SKEL_SOURCE:-$root/src/refine/skel.c}" \
    "$root/src/refine/bfill.c" \
    "$root/src/refine/triang.c" "$root/src/dbase/point.c" \
    "$root/src/dbase/node.c" "$root/src/dbase/edge.c" \
    "$root/src/dbase/element.c" "$root/src/dbase/list.c" \
    "$root/src/dbase/remove.c" \
    "$root/src/dbase/mesh_patch.c" "$root/src/dbase/mesh_patch_remove.c" \
    "$root/src/dbase/locate.c" "$root/src/dbase/clktri.c" \
    "$root/src/dbase/geometry.c" "$root/src/dbase/geom.c" \
    "$root/src/geom/misc.c" \
    -Wl,--gc-sections -lm -o "$test_dir/etch_corner"
"$test_dir/etch_corner"
