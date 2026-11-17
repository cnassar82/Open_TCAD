#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -O1 -g -fsanitize=address,undefined -fno-pie -no-pie \
 -I "$root/src/include" "$root/tests/pivoted_solver.c" -lumfpack -lm -o "$work/test"
OPENBLAS_NUM_THREADS=1 "$work/test"
