#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
gcc -std=gnu89 -O1 -g -fsanitize=address,undefined -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/directional_etch.c" -lm -o "$work/test"
UBSAN_OPTIONS=halt_on_error=1 "$work/test"
