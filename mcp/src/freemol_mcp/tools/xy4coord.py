"""Wraps XY4Coord: a symmetric-coordinate displacement -> Cartesian, for
general XY4-type molecules (see ch4sym2cart for the CH4-specific tool).
"""

from __future__ import annotations

import re

from .. import binary
from ..citation import Citation
from ..molecule import CH4_REFERENCE, Atom, format_molecule_section

CITATION = Citation(
    routine="get_ra + get_cart",
    file="Freemol/programs/XY4Coord/XY4coord.F90",
    line_start=244,
    line_end=273,
    note=(
        "Self-check at XY4coord.F90:756-763. See README.md's Known "
        "Issues: a pure angle displacement (s2a only, no bond-length "
        "change) can still fail this self-check -- root cause traced "
        "but not fixed. This tool surfaces that as an error rather than "
        "returning a result the program itself flagged as wrong."
    ),
)

_ROW_RE = re.compile(r"((?:\s+-?\d+\.\d+){10})")


def _row_after(marker: str, text: str) -> list[float]:
    idx = text.index(marker)
    m = _ROW_RE.search(text[idx + len(marker) :])
    if not m:
        raise ValueError(f"could not find a data row after {marker!r}")
    return [float(x) for x in m.group(1).split()]


def _parse_xyz_block(text: str) -> list[dict]:
    idx = text.index("[x-xyz]")
    lines = text[idx:].splitlines()
    n = int(lines[0].split()[1])
    atoms = []
    for line in lines[1 : 1 + n]:
        parts = line.split()
        atoms.append(
            {
                "label": parts[0],
                "x": float(parts[1]),
                "y": float(parts[2]),
                "z": float(parts[3]),
            }
        )
    return atoms


def xy4coord_apply_displacement(
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
    """Apply a symmetric-coordinate displacement to an XY4-type reference
    geometry and return the displaced Cartesian coordinates, via
    XY4Coord's get_ra/get_cart path with the requested Sr taken directly
    (not the "Generate Redundancies" solution search that follows it in
    the program -- that's not exposed by this tool).

    s1, s2x/s2y/s2z, s2a/s2b, s4x/s4y/s4z, sr: the 10
    [x-xy4-symmcoord] values in the program's own order (S1, S2x, S2y,
    S2z, S2a, S2b, S4x, S4y, S4z, Sr). s1=0.08 with everything else 0, on
    the default CH4 reference, stretches every bond by +0.01 (same
    example used throughout freemol's own README/tests).

    reference_molecule (optional): 5 atoms (X then 4 Y). Defaults to the
    CH4 reference geometry used throughout freemol's own test fixtures.
    """
    atoms = reference_molecule or CH4_REFERENCE
    molecule_section = format_molecule_section(atoms)
    symm_line = f"{s1} {s2x} {s2y} {s2z} {s2a} {s2b} {s4x} {s4y} {s4z} {sr}"
    input_text = f"{molecule_section}\n[x-xy4-symmcoord]\n{symm_line}\n"

    result = binary.run("XY4coord", input_text)  # note lowercase "c": the
    # program directory is XY4Coord but the built binary is XY4coord.exe.
    if result.returncode != 0:
        return {
            "error": "XY4Coord.exe exited non-zero",
            "returncode": result.returncode,
            "stderr": result.stderr,
            "citation": CITATION.as_dict(),
        }

    # Only look at the part of stdout before "Generate Redundancies":
    # errors after that point are about unphysical Sr branches, unrelated
    # to whether *this* requested displacement is valid.
    stdout = result.stdout
    prefix = stdout.split("Generate Redundancies", 1)[0]

    if "Error in Cartesian routine" in prefix:
        return {
            "error": (
                "XY4Coord's self-check failed for this displacement "
                "(it printed 'Error in Cartesian routine' -- the "
                "requested internal coordinates and the coordinates it "
                "actually built from them didn't match within "
                "tolerance). This can happen for pure-angle "
                "displacements; see the citation note."
            ),
            "stdout": prefix,
            "citation": CITATION.as_dict(),
        }
    if "Check OK at Sr input value" not in prefix:
        return {
            "error": "did not find 'Check OK at Sr input value' in output",
            "stdout": prefix,
            "citation": CITATION.as_dict(),
        }

    bonds_angles = _row_after(
        "Displaced Coordinates from SymCoord: Bonds and Angles (GRAD).", prefix
    )
    cartesian = _parse_xyz_block(prefix)

    return {
        "bonds": bonds_angles[0:4],
        "angles_degrees": {
            "a12": bonds_angles[4],
            "a13": bonds_angles[5],
            "a14": bonds_angles[6],
            "a23": bonds_angles[7],
            "a24": bonds_angles[8],
            "a34": bonds_angles[9],
        },
        "cartesian": cartesian,
        "citation": CITATION.as_dict(),
    }
