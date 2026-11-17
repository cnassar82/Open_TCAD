#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections -fdata-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/reflecting_normal.c" \
    "${SUPREM_NEW_LAYER_SOURCE:-$root/src/dbase/new_layer.c}" \
    -Wl,--gc-sections -lm -o "$test_dir/reflecting_normal"
"$test_dir/reflecting_normal"
