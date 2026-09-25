#!/usr/bin/env python3
"""Validate a fix for XY4Coord's broken Sr redundancy solve, against the
real built binary -- not a reimplementation trusted on its own.

Background (see PRECISION_NOTES.md section 3): XY4Coord's angular inputs
(S2a, S2b, S4x, S4y, S4z) are linear symmetry-coordinate displacements
that leave a second-order residual do_checks() rejects, unless the
redundant coordinate Sr is set correctly to compensate. XY4Coord's own
solver for Sr (eval_sr/qtcrt) is broken: it solves for a *cosine-based*
Sr from the (unread -- see below) Wang-Carrington paper, while the rest
of the program (get_ra/get_symc) uses a *radian-displacement-based* Sr,
so the physical root is never found.

This script does NOT use that paper (wang-carrington-2003.pdf is not in
this repo). Instead it derives Sr from first principles: 4 unit
bond-direction vectors are realizable in 3D only if their geometry
closes -- the same closure condition do_checks() already checks via its
Gamma-Sum formula (XY4coord.F90:577-589,610-639). Given Sda(1:5) fixed,
get_ra's own formula (XY4coord.F90:390-395) makes every angle an affine
function of Sr, so the Gamma-Sum residual is a well-defined scalar
function of one variable -- root-findable near Sr=0 by bracketing
outward and bisecting, exactly what a correct eval_sr() should do.

For each fixture, the found Sr is written directly into a scratch input
file (bypassing eval_sr entirely) and run through the real
bin/XY4coord.exe, asserting do_checks() itself ("Check OK at Sr input
value") accepts it -- this is what actually matters, not agreement with
a second Python implementation.

NOTE, found empirically while using this: for S4x/S4y/S4z-type inputs,
the true closure residual sits right at the ~1e-15 scale (float64's own
precision floor for this formula -- confirmed by tightening the
bisection well past any plausible benefit and watching the residual stop
improving). At that scale, THIS script's independent Python
recomputation of Sr and XY4Coord's own internal Fortran bisection
(implemented in eval_sr() using the identical formula, see
XY4coord.F90) can each round to a value a few ULPs apart -- enough to
land on opposite sides of do_checks()'s zero-tolerance boundary. So this
script can report a "FAIL" here even when the real fix (eval_sr's own
internal solve, exercised end-to-end via
Freemol/tests/run_smoke.sh's xy4coord_selfcheck on the s4x_only/
s4y_only/s4z_only fixtures) passes cleanly -- that end-to-end run is the
actual regression test; this script is a derivation/validation aid, not
a CI gate, and isn't wired into any CI job.

Usage: validate_xy4coord_redundancy.py [freemol_dir]
"""
from __future__ import annotations

import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

FR2 = math.sqrt(2.0)
FR3 = math.sqrt(3.0)
FR6 = math.sqrt(6.0)
YXY = math.acos(-1.0 / 3.0)  # reference (equilibrium) angle, all six pairs
TWO_PI = 2.0 * math.pi

MOLECULE = """[molecule] nrec=5 format=nscxyz
c  1 6 0.0 0.0 0.0
h  2 1  0.62931179  0.62931179  0.62931179
h  3 1  0.62931179 -0.62931179 -0.62931179
h  4 1 -0.62931179  0.62931179 -0.62931179
h  5 1 -0.62931179 -0.62931179  0.62931179
"""

# angular-only fixtures under test: name -> (S2a, S2b, S4x, S4y, S4z)
FIXTURES = {
    "equilibrium": (0.0, 0.0, 0.0, 0.0, 0.0),
    "s2a_only": (0.03, 0.0, 0.0, 0.0, 0.0),
    "s2b_only": (0.0, 0.03, 0.0, 0.0, 0.0),
    "s4x_only": (0.0, 0.0, 0.03, 0.0, 0.0),
    "s4y_only": (0.0, 0.0, 0.0, 0.03, 0.0),
    "s4z_only": (0.0, 0.0, 0.0, 0.0, 0.03),
}


def da_of(sda1_5, sr):
    s1, s2, s3, s4, s5 = sda1_5
    return [
        sr / FR6 + s1 / FR3 - s5 / FR2,
        sr / FR6 - s1 / FR3 + s2 / FR2 - s3 / FR2,
        sr / FR6 - s1 / FR3 - s2 / FR2 - s4 / FR2,
        sr / FR6 - s1 / FR3 - s2 / FR2 + s4 / FR2,
        sr / FR6 - s1 / FR3 + s2 / FR2 + s3 / FR2,
        sr / FR6 + s1 / FR3 + s5 / FR2,
    ]


def gamma_sums(sda1_5, sr):
    """The four do_checks() Gamma-Sum values (XY4coord.F90:577-589,
    610-639), for the given trial Sr. Returns [Y1, Y2, Y3, Y4] sums,
    each ideally == 2*pi at a valid (closing) geometry."""
    da = da_of(sda1_5, sr)
    va = [YXY + d for d in da]
    c = [math.cos(v) for v in va]
    s = [math.sin(v) for v in va]

    def agam(i, j, k):
        # cos(va[i]) - cos(va[j])*cos(va[k]) over sin(va[j])*sin(va[k])
        return (c[i] - c[j] * c[k]) / (s[j] * s[k])

    # index mapping (0-based): 0=a12 1=a13 2=a14 3=a23 4=a24 5=a34
    agam4 = agam(3, 0, 1)
    agam5 = agam(4, 0, 2)
    agam6 = agam(5, 1, 2)
    agam2 = agam(1, 0, 3)
    agam3 = agam(2, 0, 4)
    agam12 = agam(5, 3, 4)
    agam1 = agam(0, 3, 4)
    agam9 = agam(2, 1, 5)
    agam11 = agam(4, 3, 5)
    agam7 = agam(0, 2, 4)
    agam8 = agam(1, 2, 5)
    agam10 = agam(3, 4, 5)

    def sacos(x):
        return math.acos(max(-1.0, min(1.0, x)))

    y1 = sacos(agam4) + sacos(agam5) + sacos(agam6)
    y2 = sacos(agam2) + sacos(agam3) + sacos(agam12)
    y3 = sacos(agam1) + sacos(agam9) + sacos(agam11)
    y4 = sacos(agam7) + sacos(agam8) + sacos(agam10)
    return [y1, y2, y3, y4]


def max_excess(sda1_5, sr):
    """do_checks() fails iff this is > 0 for any atom; root-find target."""
    return max(g - TWO_PI for g in gamma_sums(sda1_5, sr))


def find_sr(sda1_5, window=0.5, steps=2000):
    """Bracket the root of max_excess nearest Sr=0, then bisect."""
    f = lambda sr: max_excess(sda1_5, sr)
    if abs(f(0.0)) < 1e-13:
        return 0.0
    xs = [i * (window / steps) for i in range(-steps, steps + 1)]
    best = None
    for i in range(len(xs) - 1):
        a, b = xs[i], xs[i + 1]
        fa, fb = f(a), f(b)
        if fa == 0.0:
            return a
        if fa * fb < 0.0:
            if best is None or abs(a) < abs(best[0]):
                best = (a, b)
    if best is None:
        raise RuntimeError(f"no sign change found for max_excess within +/-{window}")
    a, b = best
    fa, fb = f(a), f(b)
    for _ in range(200):
        m = 0.5 * (a + b)
        fm = f(m)
        if abs(fm) < 1e-15 or (b - a) < 1e-16:
            return m
        if fa * fm < 0.0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def run_xy4coord(exe, sda1_5, sr):
    symmline = "0.0 0.0 0.0 0.0 {} {} {} {} {} {!r}".format(*sda1_5, sr)
    text = MOLECULE + "[x-xy4-symmcoord]\n" + symmline + "\n"
    with tempfile.TemporaryDirectory(prefix="xy4coord-redundancy-") as d:
        infile = Path(d) / "input.txt"
        outfile = Path(d) / "output.txt"
        infile.write_text(text)
        proc = subprocess.run(
            [str(exe), "-i", str(infile), "-o", str(outfile)],
            cwd=d,
            capture_output=True,
            text=True,
            timeout=30,
        )
        return proc.stdout


def main(argv):
    freemol_dir = (
        Path(argv[1]).resolve() if len(argv) > 1 else Path(__file__).resolve().parents[1]
    )
    exe = freemol_dir / "bin" / "XY4coord.exe"
    if not exe.is_file():
        print(f"error: {exe} not built. Build freemol first.", file=sys.stderr)
        return 2

    failures = []
    for name, sda1_5 in FIXTURES.items():
        sr = find_sr(sda1_5)
        residual = max_excess(sda1_5, sr)
        stdout = run_xy4coord(exe, sda1_5, sr)
        ok = "Check OK at Sr input value" in stdout
        status = "PASS" if ok else "FAIL"
        print(
            f"{status}: {name:12s} Sr={sr: .10f}  max Gamma-Sum excess={residual: .3e}"
        )
        if not ok:
            failures.append(name)
            print(f"      full-precision Sr = {sr!r}")
            for line in stdout.splitlines():
                if (
                    "Check OK" in line
                    or "Check not Passed" in line
                    or "Error in Cartesian" in line
                    or "Gamma Sum for atom" in line
                    or "Evaluation of Angles" in line
                    or line.startswith("#MESERRO:        ")
                ):
                    print(f"      {line}")

    if failures:
        print(f"\n{len(failures)} fixture(s) did not reach 'Check OK': {failures}")
        return 1
    print("\nall fixtures: XY4coord.exe accepts the derived Sr directly")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
