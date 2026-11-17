#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/rectangular_ears.c" \
    "$root/src/refine/rect_tri.c" "$root/src/refine/triheur.c" \
    "$root/src/refine/tridata.c" "$root/src/refine/dbase.c" \
    "$root/src/refine/triang.c" "$root/src/dbase/point.c" \
    "$root/src/dbase/node.c" "$root/src/dbase/edge.c" \
    "$root/src/dbase/element.c" "$root/src/dbase/list.c" \
    "$root/src/dbase/geom.c" "$root/src/geom/misc.c" \
    -Wl,--gc-sections -lm -o "$test_dir/rectangular_ears"
"$test_dir/rectangular_ears"
