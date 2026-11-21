#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
src=${SUPREM_TEST_SOURCE:-$root/src}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
gcc -std=gnu89 -O1 -g -fcommon -ffunction-sections -fdata-sections \
 -fsanitize=address,undefined -fno-pie -no-pie -Wl,--gc-sections \
 -I "$src/include" -I "$root/src/include" "$root/tests/charged_flux.c" "$src/diffuse/block_set.c" \
 -lm -o "$work/test"
SUPREM4_CHARGED_FLUX=sg ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$work/test"
