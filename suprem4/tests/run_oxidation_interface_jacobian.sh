#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
src=${SUPREM_TEST_SOURCE:-$root/src}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
gcc -std=gnu89 -O1 -g -fcommon -ffunction-sections -fdata-sections \
 -fsanitize=address,undefined -fno-pie -no-pie -Wl,--gc-sections \
 -I "$src/include" "$root/tests/oxidation_interface_jacobian.c" \
 "$src/finel/FEassmb.c" "$src/oxide/triox.c" -lm -o "$work/test"
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$work/test"
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$work/test" two-oxide-sides
