#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/key_reader.c" \
    "${KEY_READ_SOURCE:-$root/src/misc/read.c}" "$root/src/misc/reader.c" \
    -o "$test_dir/key_reader"
UBSAN_OPTIONS=halt_on_error=1 "$test_dir/key_reader" "$test_dir/key.bin"
