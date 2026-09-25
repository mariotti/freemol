"""Wraps XY4Coord: a symmetric-coordinate displacement -> Cartesian, for
general XY4-type molecules (see ch4sym2cart for the CH4-specific tool).
"""

from __future__ import annotations

import re

from .. import binary
from ..citation import Citation
from ..molecule import CH4_REFERENCE, Atom, format_molecule_section

CITATION = Citation(
    routine="get_ra + get_cart, plus eval_sr's redundancy solve when the "
    "direct displacement needs it",
    file="Freemol/programs/XY4Coord/XY4coord.F90",
    line_start=244,
    line_end=273,
    note=(
        "Pure-angle displacements (s2a/s2b/s4x/s4y/s4z) leave a "
        "second-order residual the direct Sr=0 input can't satisfy on "
        "its own -- expected, see PRECISION_NOTES.md section 3. "
        "XY4Coord's own redundancy solve (eval_sr(), "
        "XY4coord.F90:1361-1417) finds the Sr that resolves it; this "
        "tool falls back to that solution when the direct check fails, "
        "and reports the Sr value it found. Raising Sr alone (with no "
        "other displacement) is geometrically impossible and still "
        "returns an error -- see README.md's Known issues."
    ),
)

_ROW_RE = re.compile(r"((?:\s+-?\d+\.\d+){10})")
_SR_VALUE_RE = re.compile(r"Sr Value:\s*(-?\d+\.\d+)")


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
    when that's already self-consistent (the common case for radial
    displacements), or via XY4Coord's own "Generate Redundancies" ->
    eval_sr() solve when it isn't (pure-angle displacements leave a
    second-order residual the direct Sr can't satisfy on its own -- see
    the citation note). Either way, the response's "sr_used" field says
    which Sr the returned geometry actually corresponds to.

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

    stdout = result.stdout

    # Fast path: the direct displacement (requested Sr taken as-is)
    # already self-consistent -- the common case for radial (S1/S2x/S2y/
    # S2z) displacements, and for Sr=0 with nothing else displaced.
    prefix = stdout.split("Generate Redundancies", 1)[0]
    if "Check OK at Sr input value" in prefix:
        return _parse_solution_block(prefix, requested_sr=sr, sr_found=sr)

    # Fallback: the direct displacement leaves a second-order residual
    # (expected for pure-angle displacements -- see the citation note).
    # XY4Coord's own redundancy solve tries corrected Sr candidates
    # after "Generate Redundancies", computed from s2a/s2b/s4x/s4y/s4z
    # alone -- it does NOT know or care what Sr was originally
    # requested, so if the caller asked for a specific nonzero Sr and
    # this path is reached, the returned geometry may correspond to a
    # different Sr than requested. _parse_solution_block flags that
    # explicitly rather than silently substituting one displacement for
    # another.
    for block in stdout.split("# [XY4C] Checking Solution: ")[1:]:
        if "Check OK at Sr solution" not in block:
            continue
        m = _SR_VALUE_RE.search(block)
        sr_found = float(m.group(1)) if m else None
        return _parse_solution_block(block, requested_sr=sr, sr_found=sr_found)

    return {
        "error": (
            "XY4Coord could not find a self-consistent geometry for this "
            "displacement, either directly or via its Sr redundancy "
            "solve. If this is a pure Sr displacement (raising all six "
            "angles together, nothing else set), that's expected -- see "
            "README.md's Known issues."
        ),
        "stdout": stdout,
        "citation": CITATION.as_dict(),
    }


def _parse_solution_block(
    text: str, requested_sr: float, sr_found: float | None
) -> dict:
    bonds_angles = _row_after(
        "Displaced Coordinates from SymCoord: Bonds and Angles (GRAD).", text
    )
    cartesian = _parse_xyz_block(text)

    result = {
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
        "requested_sr": requested_sr,
        "sr_used": sr_found,
        "citation": CITATION.as_dict(),
    }
    # Only worth flagging if the caller explicitly asked for a nonzero Sr
    # (the default, 0.0, just means "let the redundancy solve pick one" --
    # that's the normal case for pure-angle displacements, not a mismatch
    # worth a note).
    if (
        requested_sr != 0.0
        and sr_found is not None
        and abs(sr_found - requested_sr) > 1.0e-6
    ):
        result["note"] = (
            f"The requested sr ({requested_sr!r}) did not give a "
            "self-consistent geometry on its own; this result uses the "
            f"Sr XY4Coord's own redundancy solve found instead "
            f"({sr_found!r}), computed from s2a/s2b/s4x/s4y/s4z alone "
            "-- it may not be what you expect if you specifically "
            "wanted that Sr value applied."
        )
    return result
