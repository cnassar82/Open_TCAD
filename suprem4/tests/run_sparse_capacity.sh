#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections -DDEVICE \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/sparse_capacity.c" \
    "$root/src/math/generate.c" "$root/src/math/symb.c" \
    -Wl,--gc-sections -lm -o "$test_dir/sparse_capacity"
"$test_dir/sparse_capacity"
