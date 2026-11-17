#!/bin/sh
# Launch the accepted build with the model/key/implant files it was tested with.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
candidate="$root/build/validated-native64"

if [ ! -x "$candidate/src/suprem4" ]; then
    echo "SUPREM: validated 64-bit executable is missing: $candidate/src/suprem4" >&2
    exit 127
fi
for name in suprem.uk modelrc sup4gs.imp; do
    if [ ! -r "$candidate/data/$name" ]; then
        echo "SUPREM: validated runtime data is missing: $candidate/data/$name" >&2
        exit 127
    fi
done

SUPREM4_HOME=$root
SUP4KEYFILE="$candidate/data/suprem.uk"
# The legacy key reader checks KEYFILE before the path supplied by main().
KEYFILE=$SUP4KEYFILE
SUP4MODELRC="$candidate/data/modelrc"
SUP4IMPDATA="$candidate/data/sup4gs.imp"
SUP4MANDIR="$root/help"
export SUPREM4_HOME SUP4KEYFILE KEYFILE SUP4MODELRC SUP4IMPDATA SUP4MANDIR

exec "$candidate/src/suprem4" "$@"
