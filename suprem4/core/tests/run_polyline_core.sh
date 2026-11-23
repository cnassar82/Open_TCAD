#!/bin/bash
# Run the six poly-line re-oxidation decks (doc/poly_oxidation/test_decks) with
# the geometry core active, in parallel, and summarise how far each got.
#
# usage: run_polyline_core.sh SUPREM_BUILD_DIR OUT_DIR [timeout_s]
#   SUPREM_BUILD_DIR  a native64 build linked with --geometry-core
#                     (contains src/suprem4 and data/)
#   OUT_DIR           one sub-directory per deck (log, before.str, after.str or step_failure.str)
# Extra SUPREM4_CORE_* variables in the environment are passed through.
set -u
BD=$(cd "$1" && pwd); OUT=$2; TMO=${3:-3600}
B=$BD/src/suprem4; D=$BD/data
T=$(cd "$(dirname "$0")/../../doc/poly_oxidation/test_decks" && pwd)
export SUP4KEYFILE=$D/suprem.uk KEYFILE=$D/suprem.uk SUP4MODELRC=$D/modelrc SUP4IMPDATA=$D/sup4gs.imp IMPFILE=$D/sup4gs.imp
export SUPREM4_GEOMETRY_CORE=1
mkdir -p "$OUT"
run() {
  local d=$OUT/$1; mkdir -p "$d"; cp "$T/poly_line_$1.in" "$d/"
  ( cd "$d" && ulimit -v 4000000 && timeout "$TMO" "$B" < "poly_line_$1.in" > log 2>&1 )
}
for k in base coarse dry1000_120 narrow025 wet900_20 wide080; do run $k & done
wait
printf '%-12s %10s %7s %6s %6s  %s\n' deck t_reached_s points regen reject result
for k in base coarse dry1000_120 narrow025 wet900_20 wide080; do
  L=$OUT/$k/log
  t=$(grep Solving "$L" | tail -n 1 | awk '{print $6}'); np=$(grep Solving "$L" | tail -n 1 | awk '{print $NF}')
  r=$( [ -f "$OUT/$k/after.str" ] && echo COMPLETED || (grep -h "panic\|ERROR" "$L" | head -n 1 | cut -c1-60) )
  printf '%-12s %10s %7s %6s %6s  %s\n' $k "${t%,}" "$np" "$(grep -c regenerated "$L")" "$(grep -c continuing "$L")" "$r"
done
