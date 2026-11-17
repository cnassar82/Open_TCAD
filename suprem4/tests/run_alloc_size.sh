#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/alloc_size.c" -o "$test_dir/alloc_size"
"$test_dir/alloc_size"
