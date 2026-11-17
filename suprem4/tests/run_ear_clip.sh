#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/ear_clip.c" \
    "$root/src/refine/triheur.c" "$root/src/dbase/geom.c" \
    "$root/src/refine/triang.c" \
    "$root/src/geom/misc.c" \
    -Wl,--gc-sections -lm -o "$test_dir/ear_clip"
"$test_dir/ear_clip"
