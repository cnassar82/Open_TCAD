#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/new_layer_corner.c" \
    "$root/src/dbase/new_layer.c" "$root/src/dbase/locate.c" \
    -Wl,--gc-sections -lm -o "$test_dir/new_layer_corner"
"$test_dir/new_layer_corner"
