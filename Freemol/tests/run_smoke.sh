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

# Individual symmetry-displacement dimensions beyond S1. ch4sym2cart's
# vzmatcrd formulas (ch4sym2cart.F90:200-236) show r1-r4 (bond lengths)
# depend only on S1/S2z/S2a/S2b, and a12-a34 (angles) only on
# S2x/S2y/S4x/S4y/S4z/Sr -- a DIFFERENT grouping than XY4Coord's below,
# despite the identical field names/order (confirmed empirically before
# writing this: the two programs do not share a coordinate convention).
# Verified: each "radial" dimension changes bond lengths unevenly while
# leaving all angles at the tetrahedral 109.4712; each "angular" dimension
# leaves all bonds at the reference 1.0900 while perturbing angles.
ch4sym2cart_displacement() {
    local name="$1" group="$2"  # group: radial | angular
    local workdir stdout
    workdir="$(mktemp -d)"
    stdout="$workdir/stdout"
    "$BIN/ch4sym2cart.exe" -i "$DATA/ch4sym2cart/tests/${name}.inp" -o "$workdir/echo" > "$stdout" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "  exit $rc"
        cat "$stdout"
        rm -rf "$workdir"
        return 1
    fi
    python3 - "$stdout" "$group" <<'EOF'
import re
import sys
import math
import itertools

text = open(sys.argv[1]).read()
group = sys.argv[2]
coords = {}
for m in re.finditer(r"NEW Coordinates H(\d)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", text):
    coords[int(m.group(1))] = tuple(float(m.group(i)) for i in (2, 3, 4))
if len(coords) != 4:
    print(f"  expected 4 NEW Coordinates lines, found {len(coords)}")
    sys.exit(1)

dists = {i: math.sqrt(sum(c * c for c in coords[i])) for i in coords}
angles = {}
for i, j in itertools.combinations(range(1, 5), 2):
    a, b = coords[i], coords[j]
    cosv = sum(x * y for x, y in zip(a, b)) / (dists[i] * dists[j])
    angles[(i, j)] = math.degrees(math.acos(max(-1.0, min(1.0, cosv))))

bonds_changed = any(abs(d - 1.09) > 1e-3 for d in dists.values())
angles_changed = any(abs(a - 109.4712) > 1e-2 for a in angles.values())

ok = True
if group == "radial":
    if not bonds_changed:
        print("  expected bond lengths to change, none did")
        ok = False
    if angles_changed:
        print(f"  expected all angles ~109.4712, got {angles}")
        ok = False
else:
    if bonds_changed:
        print(f"  expected all bonds ~1.0900, got {dists}")
        ok = False
    if not angles_changed:
        print("  expected some angle to change, none did")
        ok = False

sys.exit(0 if ok else 1)
EOF
    local rc2=$?
    rm -rf "$workdir"
    return $rc2
}
check ch4sym2cart_displacement s2x_only angular
check ch4sym2cart_displacement s2y_only angular
check ch4sym2cart_displacement s2z_only radial
check ch4sym2cart_displacement s2a_only radial
check ch4sym2cart_displacement s2b_only radial
check ch4sym2cart_displacement s4x_only angular
check ch4sym2cart_displacement s4y_only angular
check ch4sym2cart_displacement s4z_only angular
check ch4sym2cart_displacement sr_only angular

# --- XY4PolySphere ----------------------------------------------------
# ch4_poly.inp: tetrahedral CH4 at r=1.10 (bonds all +0.01 from the r=1.09
# reference geometry) -- must produce 4x1.1000 bonds and 6x109.4712 angles,
# and the Metpot4 value at r=1.10 must be larger than at the r=1.09
# reference (both are printed by a single run).
xy4polysphere_geometry() {
    local workdir out
    workdir="$(mktemp -d)"
    out="$workdir/out"
    "$BIN/XY4PolySphere.exe" -i "$DATA/XY4PolySphere/tests/ch4_poly.inp" -o "$out" > "$workdir/stdout" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "  exit $rc"
        cat "$workdir/stdout"
        rm -rf "$workdir"
        return 1
    fi
    python3 - "$out" <<'EOF'
import re
import sys

text = open(sys.argv[1]).read()

def row_after(marker):
    idx = text.index(marker)
    rest = text[idx + len(marker):]
    m = re.search(r"((?:\s+-?\d+\.\d+){11})", rest)
    return [float(x) for x in m.group(1).split()]

ref = row_after("[x-out-result] Bonds, Angles DEGREE and Metpot4 data:")
new = row_after("[x-out-result] Bonds and Angles DEGREE")

bonds, angles, metpot_new, metpot_ref = new[0:4], new[4:10], new[10], ref[10]

ok = True
for b in bonds:
    if abs(b - 1.1000) > 1e-3:
        print(f"  bond {b} != 1.1000")
        ok = False
for a in angles:
    if abs(a - 109.4712) > 1e-3:
        print(f"  angle {a} != 109.4712")
        ok = False
if not (metpot_new > metpot_ref):
    print(f"  Metpot4 at r=1.10 ({metpot_new}) not > r=1.09 ({metpot_ref})")
    ok = False

sys.exit(0 if ok else 1)
EOF
    local rc2=$?
    rm -rf "$workdir"
    return $rc2
}
check xy4polysphere_geometry

# -R <seed> must be reproducible: same seed, same -r results file, twice.
xy4polysphere_random_reproducible() {
    local workdir
    workdir="$(mktemp -d)"
    "$BIN/XY4PolySphere.exe" -i "$DATA/XY4PolySphere/tests/ch4_poly_random.inp" \
        -o "$workdir/out1" -R 12345 -r "$workdir/res1" > /dev/null 2>&1
    "$BIN/XY4PolySphere.exe" -i "$DATA/XY4PolySphere/tests/ch4_poly_random.inp" \
        -o "$workdir/out2" -R 12345 -r "$workdir/res2" > /dev/null 2>&1
    if diff -q "$workdir/res1" "$workdir/res2" > /dev/null 2>&1; then
        rm -rf "$workdir"
        return 0
    else
        diff "$workdir/res1" "$workdir/res2"
        rm -rf "$workdir"
        return 1
    fi
}
check xy4polysphere_random_reproducible

# --- XY4Coord -----------------------------------------------------------
# Self-check: output up to "Generate Redundancies" must show "Check OK at
# Sr input value" and must NOT show "Error in Cartesian routine" (errors
# *after* that point are expected -- unphysical Sr branches get rejected).
xy4coord_selfcheck() {
    local fixture="$1"
    local workdir prefix
    workdir="$(mktemp -d)"
    "$BIN/XY4coord.exe" -i "$fixture" -o "$workdir/out" > "$workdir/stdout" 2>&1
    prefix="$(sed '/Generate Redundancies/q' "$workdir/stdout")"
    rm -rf "$workdir"
    if ! grep -q "Check OK at Sr input value" <<< "$prefix"; then
        echo "  no 'Check OK at Sr input value' before Generate Redundancies"
        return 1
    fi
    if grep -q "Error in Cartesian routine" <<< "$prefix"; then
        echo "  unexpected 'Error in Cartesian routine' before Generate Redundancies"
        return 1
    fi
    return 0
}
check xy4coord_selfcheck "$DATA/XY4Coord/tests/ch4_s1.inp"
check xy4coord_selfcheck "$DATA/XY4Coord/tests/equilibrium.inp"

# Individual symmetry-displacement dimensions beyond S1. Confirmed
# empirically that XY4Coord's Sdr group (S1, S2x, S2y, S2z) all pass the
# self-check (bonds change unevenly, angles stay at 109.4712), while every
# dimension in the Sda group (S2a, S2b, S4x, S4y, S4z, Sr) currently fails
# -- via one of two distinct paths. NOT the same grouping as ch4sym2cart's
# above, despite the identical field names/order.
check xy4coord_selfcheck "$DATA/XY4Coord/tests/s2x_only.inp"
check xy4coord_selfcheck "$DATA/XY4Coord/tests/s2y_only.inp"
check xy4coord_selfcheck "$DATA/XY4Coord/tests/s2z_only.inp"

# Sda-group (angular) displacements: S2a, S2b, S4x, S4y, S4z all leave a
# second-order residual the direct Sr=0 input can't satisfy on its own
# (expected -- see PRECISION_NOTES.md section 3), so unlike the
# self-checks above, these are *not* required to pass before "Generate
# Redundancies"; they're required to pass via it, at the Sr eval_sr()
# actually finds. (Sr alone, pushing all six angles the same direction,
# stays a separate, expected-by-design failure -- no fixture for it here.)
xy4coord_redundancy_selfcheck() {
    local fixture="$1"
    local workdir
    workdir="$(mktemp -d)"
    "$BIN/XY4coord.exe" -i "$fixture" -o "$workdir/out" > "$workdir/stdout" 2>&1
    local stdout
    stdout="$(cat "$workdir/stdout")"
    rm -rf "$workdir"
    if ! grep -q "Check OK at Sr solution" <<< "$stdout"; then
        echo "  no 'Check OK at Sr solution' found (redundancy solve didn't find a valid Sr)"
        return 1
    fi
    return 0
}
check xy4coord_redundancy_selfcheck "$DATA/XY4Coord/tests/s2a_only.inp"
check xy4coord_redundancy_selfcheck "$DATA/XY4Coord/tests/s2b_only.inp"
check xy4coord_redundancy_selfcheck "$DATA/XY4Coord/tests/s4x_only.inp"
check xy4coord_redundancy_selfcheck "$DATA/XY4Coord/tests/s4y_only.inp"
check xy4coord_redundancy_selfcheck "$DATA/XY4Coord/tests/s4z_only.inp"

# --- Frimol -----------------------------------------------------------
check "$BIN/Frimol.exe"

echo "----"
if [[ $fail -eq 0 ]]; then
    echo "smoke tests: all passed"
else
    echo "smoke tests: FAILURES"
fi
exit $fail
