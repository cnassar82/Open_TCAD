#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
ordering_flags=
ordering_libs=
if [ "${SUPREM_TEST_AMD:-0}" = 1 ]; then
    ordering_flags=-DSUPREM_USE_AMD
    ordering_libs=-lamd
fi
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon $ordering_flags \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/min_fill_pool.c" \
    $ordering_libs -lm -o "$test_dir/min_fill_pool"
"$test_dir/min_fill_pool"
