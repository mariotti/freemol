"""Wraps XY4PolySphere's poly2cart: polyspherical -> Cartesian."""

from __future__ import annotations

import re

from .. import binary
from ..citation import Citation
from ..molecule import CH4_REFERENCE, Atom, format_molecule_section

CITATION = Citation(
    routine="poly2cart",
    file="Freemol/programs/XY4PolySphere/XY4PolySphere.F90",
    line_start=1580,
    line_end=1696,
    note=(
        "Called from the [x-xy4-polyspherical] input section handler "
        "(XY4PolySphere.F90:403-460)."
    ),
)

_ROW_RE = re.compile(r"((?:\s+-?\d+\.\d+){11})")


def _row_after(marker: str, text: str) -> list[float]:
    idx = text.index(marker)
    m = _ROW_RE.search(text[idx + len(marker) :])
    if not m:
        raise ValueError(f"could not find a data row after {marker!r}")
    return [float(x) for x in m.group(1).split()]


def _parse_xyz(text: str) -> list[dict]:
    idx = text.index("[x-out-result] XYZ format")
    lines = text[idx:].splitlines()[1:]  # drop the marker line
    n = int(lines[0].strip())
    atoms = []
    for line in lines[2 : 2 + n]:  # lines[1] is blank
        parts = line.split()
        atoms.append(
            {
                "element": parts[0],
                "x": float(parts[1]),
                "y": float(parts[2]),
                "z": float(parts[3]),
            }
        )
    return atoms


def xy4polysphere_to_cartesian(
    r1: float,
    r2: float,
    r3: float,
    r4: float,
    th3: float,
    th2: float,
    th1: float,
    phi2: float,
    phi1: float,
    reference_molecule: tuple[Atom, ...] | None = None,
) -> dict:
    """Convert polyspherical coordinates to Cartesian for an XY4-type
    (e.g. methane-like) molecule, via poly2cart.

    Parameter names and order mirror the Fortran routine's own read order
    exactly (plr(1:4), pla(1:5) = r1,r2,r3,r4,th3,th2,th1,phi2,phi1) --
    see the citation. r1-r4 are the four X-Y bond lengths; th1-th3 are
    polar angles and phi1-phi2 azimuthal angles, all in degrees, with phi
    in (0, 180]. Convention, from the code's own comment: Y4 sits on the
    z axis, Y3 at phi=180, Y2 at +phi2, Y1 at -phi1.

    reference_molecule (optional): 5 atoms (X then 4 Y) used only to seed
    the orientation of the *output* Cartesian frame -- the bond lengths
    and angles in the result do not depend on it (verified empirically:
    two very different reference molecules give the same
    bonds/angles/Metpot4 for the same polyspherical input, only a
    different rotation of the Cartesian result). Defaults to the CH4
    reference geometry used throughout freemol's own test fixtures.
    """
    atoms = reference_molecule or CH4_REFERENCE
    molecule_section = format_molecule_section(atoms)
    poly_line = f"{r1} {r2} {r3} {r4} {th3} {th2} {th1} {phi2} {phi1}"
    input_text = f"{molecule_section}\n[x-xy4-polyspherical]\n{poly_line}\n"

    result = binary.run("XY4PolySphere", input_text)
    if result.returncode != 0:
        return {
            "error": "XY4PolySphere.exe exited non-zero",
            "returncode": result.returncode,
            "stderr": result.stderr,
            "citation": CITATION.as_dict(),
        }

    output = result.output_file
    reference_row = _row_after(
        "[x-out-result] Bonds, Angles DEGREE and Metpot4 data:", output
    )
    result_row = _row_after("[x-out-result] Bonds and Angles DEGREE", output)
    cartesian = _parse_xyz(output)

    return {
        "bonds": result_row[0:4],
        "angles_degrees": {
            "a12": result_row[4],
            "a13": result_row[5],
            "a14": result_row[6],
            "a23": result_row[7],
            "a24": result_row[8],
            "a34": result_row[9],
        },
        "metpot4": result_row[10],
        "reference_geometry_metpot4": reference_row[10],
        "cartesian": cartesian,
        "citation": CITATION.as_dict(),
    }
