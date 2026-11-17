#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -DSUPREM_DOMAIN_TEST -Dmalloc=domain_test_malloc -Dcalloc=domain_test_calloc \
    -I "$root/src/include" "$root/tests/material_domain.c" \
    "$root/src/dbase/material_domain.c" -Wl,--gc-sections -lm -o "$test_dir/domain"
"$test_dir/domain"
"${CC:-gcc}" -std=gnu89 -g -O1 -fcommon -ffunction-sections \
    -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
    -Dmalloc=domain_test_malloc -Dcalloc=domain_test_calloc \
    -I "$root/src/include" "$root/tests/material_domain_db.c" \
    "$root/src/dbase/material_domain.c" "$root/src/dbase/material_domain_db.c" \
    -Wl,--gc-sections -lm -o "$test_dir/domain_db"
"$test_dir/domain_db"
