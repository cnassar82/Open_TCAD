#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/mesh_trial.c" \
    "$root/src/dbase/mesh_trial.c" "${MESH_TRIAL_FLIP_SOURCE:-$root/src/dbase/flip.c}" \
    "$root/src/dbase/list.c" "$root/src/dbase/locate.c" \
    "$root/src/dbase/clktri.c" "$root/src/dbase/geometry.c" \
    "$root/src/geom/misc.c" \
    -Wl,--gc-sections -lm -o "$test_dir/mesh_trial"
"$test_dir/mesh_trial"
