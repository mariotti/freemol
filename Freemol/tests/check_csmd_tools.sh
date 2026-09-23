#!/bin/bash
# Check tools/CSMD/genGTOcoeff.m against the values used in the README CSMG
# example: O (nc=8, ac=0, radius ard(8)) and H (nc=1, ac=0, radius ard(1)).
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FREEMOL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CSMD_DIR="$FREEMOL_DIR/tools/CSMD"

OCTAVE=octave-cli
command -v "$OCTAVE" >/dev/null 2>&1 || OCTAVE=octave

result="$(cd "$CSMD_DIR" && "$OCTAVE" --no-gui --eval "
source('loadRadii.m');
[g,z] = genGTOcoeff(8,0,ard(8));
printf('%.6f %.6f\n', g, z);
[g,z] = genGTOcoeff(1,0,ard(1));
printf('%.6f %.6f\n', g, z);
")"

echo "$result"

python3 - "$result" <<'EOF'
import sys

expected = {
    "O": (2.883843, 1.591247),
    "H": (3.838395, 7.701635),
}
tol = 1e-6

lines = [l for l in sys.argv[1].strip().splitlines() if l.strip()]
if len(lines) != 2:
    print(f"expected 2 result lines, got {len(lines)}")
    sys.exit(1)

ok = True
for label, line in zip(("O", "H"), lines):
    gnorm, zeta = (float(x) for x in line.split())
    exp_g, exp_z = expected[label]
    if abs(gnorm - exp_g) > tol or abs(zeta - exp_z) > tol:
        print(f"{label}: got gnorm={gnorm} zeta={zeta}, expected gnorm={exp_g} zeta={exp_z}")
        ok = False

sys.exit(0 if ok else 1)
EOF
rc=$?

if [[ $rc -eq 0 ]]; then
    echo "check_csmd_tools: PASS"
else
    echo "check_csmd_tools: FAIL"
fi
exit $rc
