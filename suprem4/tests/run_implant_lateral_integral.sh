#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -O1 -g -fcommon -fsanitize=address,undefined -fno-pie -no-pie \
 -I "$root/src/include" "$root/tests/implant_lateral_integral.c" \
 "$root/src/implant/imp_vert.c" "$root/src/implant/imp_lat.c" -lm -o "$work/test"
"$work/test"
