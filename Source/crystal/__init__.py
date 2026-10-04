"""Crystal-structure core: the CrystalStructure class and the algorithms behind it.

Modules, in dependency order:

- ``records``         data records (AtomRecord, ZMatrixEntry, ...), defaults, type aliases
- ``elements``        element-symbol normalisation and ordering
- ``lattice``         cell parameters, fractional and Cartesian coordinates
- ``symmetry``        spglib detection, symmetry operations, SHELX LATT/SYMM records
- ``fileio``          helpers shared by the file formats
- ``cif``, ``pdb``, ``res``   file formats
- ``connectivity``    bond graphs, molecules, periodic unwrapping
- ``zmatrix_builder`` Z-matrix construction and shared-topology evaluation
- ``plotting``        unit-cell PNG images
- ``structure``       the CrystalStructure class (public API; delegates to the above)

``Source.crystal_structure`` re-exports this API for existing imports.
"""

from .records import (
    DEFAULT_CIF_SITE_MERGE_TOLERANCE,
    DEFAULT_COVALENT_SCALE,
    DEFAULT_GAS_PHASE_BOX_LENGTH,
    DEFAULT_LINEAR_THRESHOLD,
    ELEMENT_PRIORITY,
    HYDROGEN_COVALENT_SCALE,
    AtomRecord,
    BondEdge,
    CellParameters,
    CifExpansionReport,
    CSORMSymmetrySanityCheck,
    MoleculeGroup,
    SpaceGroupDetection,
    Vector3,
    ZMatrixEntry,
    ZMatrixRepresentation,
    zmatrix_topology_signature,
)
from .structure import CrystalStructure

__all__ = [
    "AtomRecord",
    "BondEdge",
    "CSORMSymmetrySanityCheck",
    "CellParameters",
    "CifExpansionReport",
    "CrystalStructure",
    "DEFAULT_CIF_SITE_MERGE_TOLERANCE",
    "DEFAULT_COVALENT_SCALE",
    "DEFAULT_GAS_PHASE_BOX_LENGTH",
    "DEFAULT_LINEAR_THRESHOLD",
    "ELEMENT_PRIORITY",
    "HYDROGEN_COVALENT_SCALE",
    "MoleculeGroup",
    "SpaceGroupDetection",
    "Vector3",
    "ZMatrixEntry",
    "ZMatrixRepresentation",
    "zmatrix_topology_signature",
]
