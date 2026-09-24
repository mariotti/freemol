#!/bin/bash
# CSMG regression tests: replay every reference case under data/CSMG/outs/
# and compare (via numdiff.py) against the input under data/CSMG/tests/
# with the same base name.
#
# Two reference sets, both under data/CSMG/outs/:
#   osx/       8 cases, generated on macOS in 2016 -- independent,
#              cross-decade, cross-platform validation.
#   generated/ the other 51 CSMG test inputs (out of 59 total), which had
#              no reference output at all until this directory was added.
#              These were captured from this codebase's own (correct,
#              CI-verified) output, so they only catch *future*
#              regressions from that point on -- they don't carry the
#              same independent cross-platform/cross-decade weight as
#              osx/. Still valuable: without them, 51 of 59 CSMG inputs
#              had zero regression coverage.
#
# CSMG writes a CSMD.bufferInput scratch file into the cwd, so each case runs
# in its own temp dir. stderr is captured separately and NOT compared: for
# example gfortran prints "STOP modcsm_nmopt exceeding maximum iterations"
# for G_C2, which is expected behavior, not a failure.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FREEMOL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$FREEMOL_DIR/bin/CSMG.exe"
OUTS_DIRS=("$FREEMOL_DIR/data/CSMG/outs/osx" "$FREEMOL_DIR/data/CSMG/outs/generated")
TESTS_DIR="$FREEMOL_DIR/data/CSMG/tests"

if [[ ! -x "$BIN" ]]; then
    echo "error: $BIN not found (run build_all.sh first)" >&2
    exit 1
fi

total=0
failed=0

for OUTS_DIR in "${OUTS_DIRS[@]}"; do
for ref in "$OUTS_DIR"/*.out; do
    name="$(basename "$ref" .out)"
    mld="$TESTS_DIR/$name.mld"
    total=$((total + 1))

    if [[ ! -f "$mld" ]]; then
        echo "FAIL $name (no matching $mld)"
        failed=$((failed + 1))
        continue
    fi

    workdir="$(mktemp -d)"
    (cd "$workdir" && "$BIN" < "$mld" > stdout.txt 2> stderr.txt)

    if python3 "$SCRIPT_DIR/numdiff.py" "$ref" "$workdir/stdout.txt" > "$workdir/numdiff.log" 2>&1; then
        echo "PASS $name"
    else
        echo "FAIL $name"
        sed 's/^/    /' "$workdir/numdiff.log"
        failed=$((failed + 1))
    fi
    rm -rf "$workdir"
done
done

echo "----"
echo "CSMG regression: $((total - failed))/$total passed"
[[ $failed -eq 0 ]]
