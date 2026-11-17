#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections -fdata-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -DSUPREM4_TEST_INTERFACE -I "$root/src/include" \
    "$root/tests/etch_interface.c" "$root/src/refine/etch.c" "$root/src/geom/misc.c" \
    -Wl,--gc-sections -lm -o "$test_dir/etch_interface"
UBSAN_OPTIONS=halt_on_error=1 "$test_dir/etch_interface"
