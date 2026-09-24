"""Shared [molecule] section formatting -- all three tools need this."""

from __future__ import annotations

from dataclasses import dataclass

# Atomic number lookup for the handful of elements freemol's test fixtures
# actually use. Extend if a tool needs more.
_ATOMIC_NUMBER = {"H": 1, "C": 6, "N": 7, "O": 8}


@dataclass(frozen=True)
class Atom:
    element: str
    x: float
    y: float
    z: float


# The CH4 reference geometry used throughout Freemol/data/*/tests/*.inp
# (r = 1.09 Angstrom, tetrahedral).
CH4_REFERENCE = (
    Atom("C", 0.0, 0.0, 0.0),
    Atom("H", 0.62931179, 0.62931179, 0.62931179),
    Atom("H", 0.62931179, -0.62931179, -0.62931179),
    Atom("H", -0.62931179, 0.62931179, -0.62931179),
    Atom("H", -0.62931179, -0.62931179, 0.62931179),
)


def format_molecule_section(atoms: tuple[Atom, ...]) -> str:
    lines = [f"[molecule] nrec={len(atoms)} format=nscxyz"]
    for i, atom in enumerate(atoms, start=1):
        z = _ATOMIC_NUMBER.get(atom.element.upper())
        if z is None:
            raise ValueError(
                f"unsupported element {atom.element!r}; known: "
                f"{sorted(_ATOMIC_NUMBER)}"
            )
        lines.append(
            f"{atom.element.lower():<2s} {i} {z} "
            f"{atom.x:.8f} {atom.y:.8f} {atom.z:.8f}"
        )
    return "\n".join(lines)
