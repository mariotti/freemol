#!/usr/bin/env python3
"""Compare two Freemol outputs token by token.

Text tokens must match exactly; numeric tokens (incl. Fortran D exponents)
must agree within an absolute + relative tolerance, so harmless round-off
between compilers/platforms (e.g. -0.275D-15 vs -0.278D-15) does not fail.
Only lines matching SELECT are compared (default: the CSM value trace and
results). Optimizer parameters and final coordinates are skipped on purpose:
when the CSM is invariant along a parameter (e.g. a rotation about a
symmetry axis) the simplex drifts along that flat direction and round-off
decides where it ends up, without changing the CSM itself.
Usage: numdiff.py reference.out new.out [abs_tol] [rel_tol]
"""
import re
import sys

SELECT = re.compile(r'CSM|Minimization at step')
NUM = re.compile(r'^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$')


def to_float(tok):
    return float(tok.replace('D', 'E').replace('d', 'e'))


def main():
    ref, new = sys.argv[1], sys.argv[2]
    atol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10
    rtol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-8
    a = [l for l in open(ref).read().splitlines() if SELECT.search(l)]
    b = [l for l in open(new).read().splitlines() if SELECT.search(l)]
    if len(a) != len(b):
        print(f"selected line count differs: {len(a)} vs {len(b)}")
        return 1
    bad = 0
    for n, (la, lb) in enumerate(zip(a, b), 1):
        ta, tb = la.split(), lb.split()
        ok = len(ta) == len(tb)
        if ok:
            for x, y in zip(ta, tb):
                if NUM.match(x) and NUM.match(y):
                    fx, fy = to_float(x), to_float(y)
                    if abs(fx - fy) > atol + rtol * abs(fx):
                        ok = False
                        break
                elif x != y:
                    ok = False
                    break
        if not ok:
            bad += 1
            if bad <= 5:
                print(f"line {n}:\n  ref: {la}\n  new: {lb}")
    if bad:
        print(f"{bad} differing line(s)")
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
