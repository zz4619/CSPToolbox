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


@dataclass(frozen=True)
class ViewerCoordinate:
    """An internal coordinate with a stable bnd/ang/dih row identity."""

    name: str
    kind: str
    atom_indices: tuple[int, ...]
    atom_labels: tuple[str, ...]
    value: float
    unit: str
    independent: bool = False
    reference_value: float | None = None
    lower: float | None = None
    upper: float | None = None


@dataclass(frozen=True)
class ViewerScene:
    """One selectable workflow; coordinates are prepared in Python, never refitted in JS."""

    title: str
    molecule: ViewerMolecule
    reference: ViewerMolecule | None = None
    full_cell: ViewerMolecule | None = None
    cell: tuple[tuple[float, float, float], ...] = ()
    coordinates: tuple[ViewerCoordinate, ...] = ()
    mapping: tuple[tuple[str, str], ...] = ()
    rmsd: float | None = None
    mapping_method: str = ""
    notices: tuple[str, ...] = ()
    source_paths: tuple[str, ...] = ()
    source_hashes: tuple[tuple[str, str], ...] = ()
    mapped_zmatrix: str | None = None


@dataclass(frozen=True)
class ViewerDocument:
    """A portable collection of scenes, using one renderer and payload schema."""

    scenes: tuple[ViewerScene, ...]
    title: str = "CSPToolbox molecular viewer"
    schema_version: int = 1
