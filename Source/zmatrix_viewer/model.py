"""Data models for interactive Z-matrix visualization."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ZMatrixAtom:
    """One row from a numeric CSPToolbox Z-matrix."""

    row_index: int
    line_number: int
    label: str
    element: str
    bond_to: int | None
    bond_length: float | None
    angle_to: int | None
    angle_degrees: float | None
    dihedral_to: int | None
    dihedral_degrees: float | None

    @property
    def has_dihedral(self) -> bool:
        return (
            self.bond_to is not None
            and self.angle_to is not None
            and self.dihedral_to is not None
            and self.dihedral_degrees is not None
        )


@dataclass(frozen=True)
class ZMatrixDocument:
    """Parsed Z-matrix file before Cartesian reconstruction."""

    title: str
    atoms: tuple[ZMatrixAtom, ...]
    explicit_bonds: frozenset[tuple[int, int]]
    source_name: str | None = None


@dataclass(frozen=True)
class ViewerAtom:
    """An atom ready to be serialized to the browser viewer."""

    index: int
    label: str
    element: str
    coordinates: tuple[float, float, float]
    color: str
    display_radius: float


@dataclass(frozen=True)
class ViewerBond:
    """A bond or construction link ready for display."""

    left: int
    right: int
    kind: str


@dataclass(frozen=True)
class ViewerDihedral:
    """A selectable dihedral row from the source Z-matrix."""

    id: str
    row_index: int
    atom_indices: tuple[int, int, int, int]
    atom_labels: tuple[str, str, str, str]
    value_degrees: float
    kind: str
    links: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class ViewerMolecule:
    """Complete molecule payload used by the static HTML viewer."""

    title: str
    atoms: tuple[ViewerAtom, ...]
    bonds: tuple[ViewerBond, ...]
    dihedrals: tuple[ViewerDihedral, ...]
    warnings: tuple[str, ...]
    source_name: str | None = None

