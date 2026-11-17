#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT HUP INT TERM
c++ -std=c++17 -g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie \
 -I "$root/src/third_party/clipper2/include" "$root/tests/domain_offset.cpp" \
 "${SUPREM_DOMAIN_OFFSET_SOURCE:-$root/src/refine/domain_offset.cpp}" "$root/src/third_party/clipper2/src/clipper.engine.cpp" \
 "$root/src/third_party/clipper2/src/clipper.offset.cpp" "$root/src/third_party/clipper2/src/clipper.rectclip.cpp" \
 -o "$dir/test"
"$dir/test" "$root/tests/fixtures/rit_pmos_offset_surface.dat"
