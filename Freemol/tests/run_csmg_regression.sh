#!/bin/sh
# CSMG regression test: run every test input that has a stored reference
# output and compare numerically (see numdiff.py).
# Run from anywhere after building: Freemol/tests/run_csmg_regression.sh
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(dirname "$HERE")
EXE="$ROOT/bin/CSMG.exe"
REF="$ROOT/data/CSMG/outs/osx"
IN="$ROOT/data/CSMG/tests"
TMP=${TMPDIR:-/tmp}/csmg_regression.$$
mkdir -p "$TMP"

[ -x "$EXE" ] || { echo "CSMG.exe not found, build first"; exit 2; }

pass=0; fail=0
for ref in "$REF"/*.out; do
    name=$(basename "$ref" .out)
    # references were produced reading stdin; stderr (e.g. gfortran STOP
    # and IEEE notes) is kept apart and not compared
    # run inside $TMP: CSMG writes a CSMD.bufferInput scratch file in the cwd
    (cd "$TMP" && "$EXE" < "$IN/$name.mld" > "$name.out" 2> "$name.err")
    if python3 "$HERE/numdiff.py" "$ref" "$TMP/$name.out" > "$TMP/$name.diff"; then
        echo "PASS  $name"; pass=$((pass + 1))
    else
        echo "FAIL  $name"; sed 's/^/      /' "$TMP/$name.diff"; fail=$((fail + 1))
    fi
done
echo "CSMG regression: $pass passed, $fail failed"
rm -rf "$TMP"
[ "$fail" -eq 0 ]
