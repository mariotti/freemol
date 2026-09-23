#!/bin/bash
# Smoke tests for tools that have no reference outputs to regress against.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FREEMOL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$FREEMOL_DIR/bin"
DATA="$FREEMOL_DIR/data"

fail=0

check() {
    if "$@"; then
        echo "PASS: $*"
    else
        echo "FAIL: $*"
        fail=1
    fi
}

# --- fit1Dpol -----------------------------------------------------------
# fit.inp is an exact y=x^2 fit: it must converge to f2=1.0, f0~=0.
fit1dpol_converges() {
    local workdir out f0 f2
    workdir="$(mktemp -d)"
    "$BIN/fit1Dpol.exe" -i "$DATA/fit1Dpol/examples/fit.inp" -o "$workdir/out" || true
    out="$workdir/out"
    if [[ ! -f "$out" ]]; then
        echo "  no output file produced"
        rm -rf "$workdir"
        return 1
    fi
    f0="$(awk '/^ Par: f0 /{v=$3} END{print v}' "$out")"
    f2="$(awk '/^ Par: f2 /{v=$3} END{print v}' "$out")"
    rm -rf "$workdir"
    if [[ -z "$f0" || -z "$f2" ]]; then
        echo "  could not find final f0/f2 values"
        return 1
    fi
    python3 - "$f0" "$f2" <<'EOF'
import sys
f0, f2 = float(sys.argv[1]), float(sys.argv[2])
ok = abs(f2 - 1.0) < 1e-4 and abs(f0) < 1e-4
if not ok:
    print(f"  f0={f0} f2={f2} out of tolerance")
sys.exit(0 if ok else 1)
EOF
}
check fit1dpol_converges

fit1dpol_runs() {
    local workdir="$1"
    "$BIN/fit1Dpol.exe" -i "$DATA/fit1Dpol/examples/$2" -o "$workdir/out"
}
for f in fit_debug.inp fit_debug_02.inp; do
    wd="$(mktemp -d)"
    check fit1dpol_runs "$wd" "$f"
    rm -rf "$wd"
done

# --- ch4sym2cart ----------------------------------------------------------
# CH4, r=1.09, symmetric stretch S1=0.08 -> every bond length +0.01 == 1.1000,
# and the tetrahedral H-C-H angle (109.47 deg) is preserved.
ch4sym2cart_geometry() {
    local workdir stdout
    workdir="$(mktemp -d)"
    stdout="$workdir/stdout"
    "$BIN/ch4sym2cart.exe" -i "$DATA/ch4sym2cart/tests/ch4_s1.inp" -o "$workdir/echo" > "$stdout"
    if ! grep -q "NEW Coordinates H1" "$stdout"; then
        echo "  stdout has no 'NEW Coordinates H1' line"
        rm -rf "$workdir"
        return 1
    fi
    python3 - "$stdout" <<'EOF'
import re
import sys
import math

text = open(sys.argv[1]).read()
coords = {}
for m in re.finditer(r"NEW Coordinates H(\d)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", text):
    coords[int(m.group(1))] = tuple(float(m.group(i)) for i in (2, 3, 4))

ok = True
if 1 not in coords or coords[1] != (1.1, 0.0, 0.0):
    print(f"  H1 = {coords.get(1)}, expected (1.1, 0.0, 0.0)")
    ok = False

for i in (2, 3, 4):
    if i not in coords:
        print(f"  H{i} coordinates not found")
        ok = False
        continue
    dist = math.sqrt(sum(c * c for c in coords[i]))
    if abs(dist - 1.1) > 1e-3:
        print(f"  H{i} distance from origin = {dist}, expected 1.1")
        ok = False

if 1 in coords and 2 in coords:
    h1, h2 = coords[1], coords[2]
    d1 = math.sqrt(sum(c * c for c in h1))
    d2 = math.sqrt(sum(c * c for c in h2))
    cos_angle = sum(a * b for a, b in zip(h1, h2)) / (d1 * d2)
    angle = math.degrees(math.acos(cos_angle))
    if abs(angle - 109.47) > 0.1:
        print(f"  H1-C-H2 angle = {angle}, expected ~109.47")
        ok = False

sys.exit(0 if ok else 1)
EOF
    local rc=$?
    rm -rf "$workdir"
    return $rc
}
check ch4sym2cart_geometry

# --- Frimol -----------------------------------------------------------
check "$BIN/Frimol.exe"

echo "----"
if [[ $fail -eq 0 ]]; then
    echo "smoke tests: all passed"
else
    echo "smoke tests: FAILURES"
fi
exit $fail
