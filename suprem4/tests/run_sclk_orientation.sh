#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections -fdata-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/sclk_orientation.c" \
    "${SUPREM_CLKTRI_SOURCE:-$root/src/dbase/clktri.c}" \
    "${SUPREM_GEOMETRY_SOURCE:-$root/src/dbase/geometry.c}" "$root/src/geom/misc.c" \
    -Wl,--gc-sections -lm -o "$test_dir/sclk_orientation"
"$test_dir/sclk_orientation"
