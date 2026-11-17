#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/assembly_boundary.c" \
    "$root/src/diffuse/prepare.c" "$root/src/refine/dbase.c" \
    -Wl,--gc-sections -lm -o "$test_dir/assembly_boundary"
"$test_dir/assembly_boundary"
