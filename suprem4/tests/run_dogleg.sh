#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
gcc -std=gnu89 -O1 -g -fsanitize=address,undefined -fno-pie -no-pie \
 -I "$root/src/include" "$root/tests/dogleg.c" -lm -o "$work/test"
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$work/test"
