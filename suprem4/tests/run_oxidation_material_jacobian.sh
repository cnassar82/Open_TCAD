#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
src=${SUPREM_TEST_SOURCE:-"$root/src"}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
gcc -std=gnu89 -O1 -g -fcommon -ffunction-sections -fdata-sections -Wl,--gc-sections \
 -fsanitize=address,undefined -fno-pie -no-pie -I "$src/include" \
 "$root/tests/oxidation_material_jacobian.c" "$src/oxide/Oxidant.c" -lm -o "$work/test"
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$work/test"
