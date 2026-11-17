#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections -fdata-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -DSUPREM_MESH_PATCH_TEST -Dmalloc=mesh_patch_test_malloc -Dcalloc=mesh_patch_test_calloc \
    -I "$root/src/include" "$root/tests/mesh_patch_commit.c" \
    "$root/src/dbase/mesh_patch.c" "$root/src/dbase/mesh_patch_remove.c" \
    "$root/src/geom/misc.c" -Wl,--gc-sections -lm -o "$dir/test"
UBSAN_OPTIONS=halt_on_error=1 "$dir/test"
