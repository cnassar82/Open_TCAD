#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -I "$root/src/include" "$root/tests/implant_tables.c" \
    "$root/src/implant/fgtmom.c" "$root/src/implant/imp_setprs.c" \
    -lm -o "$test_dir/implant_tables"
unset IMPFILE
"$test_dir/implant_tables" "$root/data/sup4gs.imp"
