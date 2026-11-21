#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
gcc -std=gnu89 -O1 -g -fcommon -ffunction-sections -fdata-sections \
 -fsanitize=address,undefined -fno-pie -no-pie -Wl,--gc-sections \
 -I "$root/src/include" "$root/tests/mesh_embedding.c" "$root/src/geom/misc.c" \
 -lm -o "$work/test"
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$work/test"
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$work/test" 424
