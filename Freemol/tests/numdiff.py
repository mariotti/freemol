#!/usr/bin/env python3
"""Compare two CSMG outputs, tolerant of compiler-dependent optimizer noise.

Usage: numdiff.py ref new [atol] [rtol]

Only lines matching the regex 'CSM|Minimization at step' are compared: the
Minuit simplex optimizer's intermediate parameters and final coordinates
legitimately differ across compilers/platforms (e.g. G_C3h has a flat
rotation direction the simplex drifts along with no effect on the fit), so a
full-text diff is not usable. The CSM value itself, however, must match
regardless of compiler: it is what regression testing CSMG actually cares
about.

Each matched line is compared token by token: numeric tokens (including
Fortran D-exponent notation like 0.123D+01) are parsed as floats and compared
within atol/rtol; everything else is compared as exact text.
"""
import re
import sys

LINE_RE = re.compile(r"CSM|Minimization at step")
NUMBER_RE = re.compile(
    r"^[+-]?(\d+\.\d*|\.\d+|\d+)([DdEe][+-]?\d+)?$"
)


def parse_float(token):
    if NUMBER_RE.match(token) is None:
        return None
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def matched_lines(path):
    with open(path) as f:
        return [line.rstrip("\n") for line in f if LINE_RE.search(line)]


def close_enough(a, b, atol, rtol):
    return abs(a - b) <= atol + rtol * abs(b)


def compare(ref_path, new_path, atol, rtol):
    ref_lines = matched_lines(ref_path)
    new_lines = matched_lines(new_path)

    diffs = []
    if len(ref_lines) != len(new_lines):
        diffs.append(
            f"line count mismatch: ref has {len(ref_lines)} matching lines, "
            f"new has {len(new_lines)}"
        )

    for i, (rline, nline) in enumerate(zip(ref_lines, new_lines), start=1):
        rtoks = rline.split()
        ntoks = nline.split()
        if len(rtoks) != len(ntoks):
            diffs.append(f"line {i}: token count differs\n  ref: {rline}\n  new: {nline}")
            continue
        for rtok, ntok in zip(rtoks, ntoks):
            rnum = parse_float(rtok)
            nnum = parse_float(ntok)
            if rnum is not None and nnum is not None:
                if not close_enough(rnum, nnum, atol, rtol):
                    diffs.append(
                        f"line {i}: {rtok!r} != {ntok!r} "
                        f"(diff {abs(rnum - nnum):.3g} > atol+rtol*|ref|)\n"
                        f"  ref: {rline}\n  new: {nline}"
                    )
            elif rtok != ntok:
                diffs.append(f"line {i}: {rtok!r} != {ntok!r}\n  ref: {rline}\n  new: {nline}")

    return diffs


def main(argv):
    if len(argv) < 3 or len(argv) > 5:
        print(__doc__, file=sys.stderr)
        return 2
    ref_path, new_path = argv[1], argv[2]
    atol = float(argv[3]) if len(argv) > 3 else 1e-10
    rtol = float(argv[4]) if len(argv) > 4 else 1e-8

    diffs = compare(ref_path, new_path, atol, rtol)
    if diffs:
        for d in diffs:
            print(d, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
