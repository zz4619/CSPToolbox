"""Plain data records, defaults and type aliases shared by the crystal package."""

from __future__ import annotations

from dataclasses import dataclass


# (a, b, c) in Å and (alpha, beta, gamma) in degrees.
CellParameters = tuple[float, float, float, float, float, float]
Vector3 = tuple[float, float, float]

# A pair is bonded when its distance <= scale * (sum of ASE covalent radii);
# pairs involving hydrogen use at least HYDROGEN_COVALENT_SCALE.
DEFAULT_COVALENT_SCALE = 1.20
HYDROGEN_COVALENT_SCALE = 1.30

# Edge of the cubic box (Å) for gas-phase VASP cells.
DEFAULT_GAS_PHASE_BOX_LENGTH = 20.0

# Z-matrix angle references within this many degrees of 0/180 are flagged.
DEFAULT_LINEAR_THRESHOLD = 15.0

# Symmetry copies of a CIF site closer than this (Å) are merged.
DEFAULT_CIF_SITE_MERGE_TOLERANCE = 1e-2

# Preferred order of elements when choosing Z-matrix roots and branches.
ELEMENT_PRIORITY = {
    "C": 0,
    "O": 1,
    "N": 2,
    "S": 3,
    "F": 4,
    "CL": 5,
    "BR": 6,
    "H": 99,
}


@dataclass(frozen=True)
class AtomRecord:
    label: str
    element: str
    coordinates: Vector3


@dataclass(frozen=True)
class BondEdge:
    neighbor: int
    shift: tuple[int, int, int]
    distance: float
    ratio: float


@dataclass(frozen=True)
class MoleculeGroup:
    signature: str
    representative_molecule: list[AtomRecord]
    duplicate_molecules: list[list[AtomRecord]]


@dataclass(frozen=True)
class CifExpansionReport:
    """Metadata from expanding CIF atom sites to an explicit unit cell."""

    name: str
    space_group: str
    expanded_atom_count: int
    partial_occupancies: list[tuple[str, float]]
    ase_expands_to_unit_cell: bool
    ase_atom_count: int
    raw_atom_row_count: int
    ase_matches_manual: bool
    ase_comparison_message: str | None
    site_merge_tolerance: float = DEFAULT_CIF_SITE_MERGE_TOLERANCE
    duplicate_sites_merged: int = 0


@dataclass(frozen=True)
class ZMatrixEntry:
    label: str
    element: str
    bond_to: int | None
    bond_length: float | None
    angle_to: int | None
    angle_degrees: float | None
    dihedral_to: int | None
    dihedral_degrees: float | None


@dataclass(frozen=True)
class ZMatrixRepresentation:
    molecule_index: int
    entries: list[ZMatrixEntry]
    ordered_atom_labels: list[str]
    warnings: list[str]


def zmatrix_topology_signature(
    zmatrix: ZMatrixRepresentation,
) -> tuple[tuple[str, str, int | None, int | None, int | None], ...]:
    """Return the geometry-independent identity of a Z-matrix definition."""

    if len(zmatrix.entries) != len(zmatrix.ordered_atom_labels):
        raise ValueError(
            "Z-matrix entry count does not match its ordered atom-label count."
        )
    return tuple(
        (
            label,
            entry.element,
            entry.bond_to,
            entry.angle_to,
            entry.dihedral_to,
        )
        for label, entry in zip(
            zmatrix.ordered_atom_labels,
            zmatrix.entries,
        )
    )


@dataclass(frozen=True)
class SpaceGroupDetection:
    symbol: str
    number: int
    hall_symbol: str
    hall_number: int
    equivalent_atoms: tuple[int, ...]
    shelx_latt_value: int
    symmetry_operations: tuple[str, ...]


@dataclass(frozen=True)
class CSORMSymmetrySanityCheck:
    supported: bool
    latt_value: int
    original_operations: tuple[str, ...]
    parsed_operations: tuple[str, ...]
    unsupported_operations: tuple[str, ...]
    unsupported_reasons: tuple[str, ...]
