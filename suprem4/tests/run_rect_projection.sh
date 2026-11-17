#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections -fdata-sections \
 -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
 -I "$root/src/include" "$root/tests/rect_projection.c" \
 "${SUPREM_RECT_SOURCE:-$root/src/refine/rect_tri.c}" "$root/src/dbase/geom.c" "$root/src/geom/misc.c" \
 -Wl,--gc-sections -lm -o "$dir/test"
"$dir/test"
