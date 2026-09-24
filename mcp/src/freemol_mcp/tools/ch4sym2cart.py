"""Wraps ch4sym2cart: a CH4 symmetric-coordinate displacement -> Cartesian.

Only the program's "NEW Coordinates H1-H4" output is exposed. Its other
diagnostic block ("Unchenged Coordina...") is documented in README.md's
Known Issues as geometrically wrong (a missing degrees->radians
conversion) and is deliberately left out of this tool's response.
"""

from __future__ import annotations

import re

from .. import binary
from ..citation import Citation
from ..molecule import CH4_REFERENCE, Atom, format_molecule_section

CITATION = Citation(
    routine="ch4sym2cart (program body)",
    file="Freemol/programs/ch4sym2cart/ch4sym2cart.F90",
    line_start=241,
    line_end=333,
    note=(
        "Only the 'NEW Coordinates H1-H4' block is exposed here. The "
        "program's other diagnostic output ('Unchenged Coordina...', "
        "lines ~337-344) is known-incorrect -- see README.md's Known "
        "Issues -- and is intentionally not surfaced by this tool."
    ),
)

_COORD_RE = re.compile(
    r"NEW Coordinates H(\d)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)"
)


def ch4sym2cart_apply_displacement(
    s1: float,
    s2x: float,
    s2y: float,
    s2z: float,
    s2a: float,
    s2b: float,
    s4x: float,
    s4y: float,
    s4z: float,
    sr: float = 0.0,
    reference_molecule: tuple[Atom, ...] | None = None,
) -> dict:
    """Apply a symmetric-coordinate displacement to a CH4(-like) reference
    geometry and return the displaced Cartesian coordinates of H1-H4, via
    ch4sym2cart.

    s1: symmetric stretch (all 4 bonds change equally: new bond = old
    bond + s1/4... in practice s1=0.08 -> each bond +0.01, per the
    program's own worked example). s2x/s2y/s2z, s2a/s2b, s4x/s4y/s4z, sr:
    the remaining symmetry-adapted displacement coordinates, passed
    through to [x-ch4-symmcoord] in that order (same convention as
    xy4coord_apply_displacement's Sdr/Sda).

    The carbon atom is not itself computed or printed by ch4sym2cart --
    by construction it stays at the origin (H1's bond is placed along the
    +x axis), so it's reported here as a fixed (0, 0, 0) for convenience,
    not because the Fortran code confirms it.

    reference_molecule (optional): the CH4(-like) reference geometry (5
    atoms: X then 4 Y). Defaults to the CH4 reference used throughout
    freemol's own test fixtures (r=1.09).
    """
    atoms = reference_molecule or CH4_REFERENCE
    molecule_section = format_molecule_section(atoms)
    symm_line = f"{s1} {s2x} {s2y} {s2z} {s2a} {s2b} {s4x} {s4y} {s4z} {sr}"
    input_text = (
        f"{molecule_section}\n[x-ch4-symmcoord]\n{symm_line}\n"
    )

    result = binary.run("ch4sym2cart", input_text)
    if result.returncode != 0:
        return {
            "error": "ch4sym2cart.exe exited non-zero",
            "returncode": result.returncode,
            "stderr": result.stderr,
            "citation": CITATION.as_dict(),
        }

    coords = {}
    for m in _COORD_RE.finditer(result.stdout):
        idx = int(m.group(1))
        coords[f"H{idx}"] = [float(m.group(i)) for i in (2, 3, 4)]

    if len(coords) != 4:
        return {
            "error": "could not find all 4 NEW Coordinates lines in stdout",
            "stdout": result.stdout,
            "citation": CITATION.as_dict(),
        }

    return {
        "cartesian": {
            "C": [0.0, 0.0, 0.0],
            **coords,
        },
        "citation": CITATION.as_dict(),
    }
