"""The CrystalStructure container.

One in-memory representation for the periodic structures used across the CSP
workflow: atom labels, element types, Cartesian coordinates in Å, cell
parameters (a, b, c, alpha, beta, gamma), a Hermann-Mauguin space-group label and,
for asymmetric units, the SHELX LATT/SYMM records that regenerate the cell.

The class is the public API. The algorithms live in the topic modules of this
package (file formats, symmetry, connectivity, Z-matrices, plotting); methods
here validate input and delegate.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.cell import Cell

from .cif import expand_cif_unit_cell, read_cif, write_cif
from .connectivity import (
    build_connectivity,
    build_connectivity_from_template,
    chemical_signature,
    component_records,
    connected_components,
    local_graph,
    molecule_aware_representative_indices,
)
from .elements import normalize_element_symbol
from .fileio import normalize_format
from .pdb import read_pdb, write_pdb
from .plotting import write_unit_cell_molecule_image
from .records import (
    DEFAULT_CIF_SITE_MERGE_TOLERANCE,
    DEFAULT_COVALENT_SCALE,
    DEFAULT_GAS_PHASE_BOX_LENGTH,
    DEFAULT_LINEAR_THRESHOLD,
    AtomRecord,
    BondEdge,
    CellParameters,
    CifExpansionReport,
    CSORMSymmetrySanityCheck,
    MoleculeGroup,
    SpaceGroupDetection,
    Vector3,
    ZMatrixRepresentation,
)
from .res import read_res, write_res
from .symmetry import (
    expand_shelx_atoms,
    normalize_csorm_symmetry_operation,
    normalize_space_group_symbol,
    shelx_symmetry_records,
    shelx_symmetry_records_from_symmetry,
    symmetry_dataset,
)
from .zmatrix_builder import apply_zmatrix_topology, generate_zmatrices


@dataclass
class CrystalStructure:
    """A periodic structure: either an explicit unit cell or an asymmetric unit.

    ``explict_unit_cell`` (historical spelling, kept because it is written into
    CIF/PDB/RES files) is True when ``atoms`` already fill the cell. Otherwise the
    atoms are an asymmetric unit, and ``shelx_latt_value`` plus
    ``symmetry_operations`` (or, failing those, ``space_group``/``hall_number``)
    regenerate the cell.
    """

    atoms: list[AtomRecord]
    cell_parameters: CellParameters
    lattice_matrix: tuple[Vector3, Vector3, Vector3] | None = None
    space_group: str = "P 1"
    hall_number: int | None = None
    name: str = "CrystalStructure"
    explict_unit_cell: bool = False
    shelx_latt_value: int | None = None
    symmetry_operations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.atoms:
            raise ValueError("CrystalStructure requires at least one atom.")
        if len(self.cell_parameters) != 6:
            raise ValueError("cell_parameters must be (a, b, c, alpha, beta, gamma).")
        if self.lattice_matrix is not None:
            if len(self.lattice_matrix) != 3 or any(len(vector) != 3 for vector in self.lattice_matrix):
                raise ValueError("lattice_matrix must be a 3x3 matrix when provided.")
            self.lattice_matrix = tuple(
                tuple(float(value) for value in vector) for vector in self.lattice_matrix
            )
        self.atoms = [
            AtomRecord(
                label=atom.label,
                element=normalize_element_symbol(atom.element or atom.label),
                coordinates=atom.coordinates,
            )
            for atom in self.atoms
        ]
        self.explict_unit_cell = bool(self.explict_unit_cell)
        self.symmetry_operations = tuple(str(operation).strip() for operation in self.symmetry_operations)

    # ------------------------------------------------------------------ cell
    @property
    def explicit_unit_cell(self) -> bool:
        """Correctly spelled alias of ``explict_unit_cell``."""

        return self.explict_unit_cell

    @property
    def cell(self) -> Cell:
        if self.lattice_matrix is not None:
            return Cell(np.array(self.lattice_matrix, dtype=float))
        return Cell.fromcellpar(self.cell_parameters)

    def fractional_coordinates(self, cartesian: Vector3) -> Vector3:
        scaled = self.cell.scaled_positions(np.array([cartesian], dtype=float))[0]
        return tuple(float(value) for value in scaled)

    def cartesian_coordinates(self, fractional: Vector3) -> Vector3:
        cart = np.dot(np.array(fractional, dtype=float), self.cell.array)
        return tuple(float(value) for value in cart)

    def to_ase_atoms(self) -> Atoms:
        """Return the atoms as a periodic ASE ``Atoms`` object (no symmetry expansion)."""

        return Atoms(
            symbols=[atom.element for atom in self.atoms],
            positions=np.array([atom.coordinates for atom in self.atoms], dtype=float),
            cell=self.cell.array,
            pbc=True,
        )

    # ------------------------------------------------------------- file I/O
    @classmethod
    def from_file(cls, path: str | Path, fmt: str | None = None) -> "CrystalStructure":
        """Read a CIF, PDB or RES file as written (no symmetry expansion)."""

        file_path = Path(path)
        file_format = normalize_format(file_path, fmt)
        if file_format == "cif":
            return read_cif(cls, file_path)
        if file_format == "pdb":
            return read_pdb(cls, file_path)
        if file_format == "res":
            return read_res(cls, file_path)
        raise ValueError(f"Unsupported format: {file_format}")

    @classmethod
    def expand_cif_to_unit_cell(
        cls,
        path: str | Path,
        *,
        site_merge_tolerance: float = DEFAULT_CIF_SITE_MERGE_TOLERANCE,
    ) -> "CrystalStructure":
        """Read a CIF and apply its symmetry operations to fill the unit cell."""

        structure, _ = expand_cif_unit_cell(cls, path, site_merge_tolerance=site_merge_tolerance)
        return structure

    @classmethod
    def from_cif_unit_cell(
        cls,
        path: str | Path,
        *,
        site_merge_tolerance: float = DEFAULT_CIF_SITE_MERGE_TOLERANCE,
    ) -> "CrystalStructure":
        """Alias of :meth:`expand_cif_to_unit_cell`."""

        return cls.expand_cif_to_unit_cell(path, site_merge_tolerance=site_merge_tolerance)

    @classmethod
    def inspect_cif_unit_cell_expansion(
        cls,
        path: str | Path,
        *,
        site_merge_tolerance: float = DEFAULT_CIF_SITE_MERGE_TOLERANCE,
    ) -> CifExpansionReport:
        """Expand a CIF and report merges, partial occupancies and the ASE cross-check."""

        _, report = expand_cif_unit_cell(cls, path, site_merge_tolerance=site_merge_tolerance)
        return report

    def to_file(
        self,
        path: str | Path,
        fmt: str | None = None,
        *,
        rounding: bool = False,
    ) -> None:
        """Write CIF, PDB or RES.

        With ``rounding``, RES SYMM translations within 0.01 of 0, 1/4, 1/3, 1/2,
        2/3 or 3/4 are snapped to that value.
        """

        file_path = Path(path)
        file_format = normalize_format(file_path, fmt)
        if file_format == "cif":
            write_cif(self, file_path)
            return
        if file_format == "pdb":
            write_pdb(self, file_path)
            return
        if file_format == "res":
            write_res(self, file_path, rounding=rounding)
            return
        raise ValueError(f"Unsupported format: {file_format}")

    # ------------------------------------------------------------- symmetry
    def detect_space_group_symmetry(
        self,
        *,
        symprec: float = 0.05,
        angle_tolerance: float = -1.0,
    ) -> SpaceGroupDetection:
        dataset = symmetry_dataset(self, symprec=symprec, angle_tolerance=angle_tolerance)
        symbol = normalize_space_group_symbol(str(dataset.international))
        shelx_latt_value, symmetry_operations = shelx_symmetry_records_from_symmetry(
            space_group_symbol=symbol,
            rotations=np.asarray(dataset.rotations, dtype=int),
            translations=np.asarray(dataset.translations, dtype=float),
        )
        return SpaceGroupDetection(
            symbol=symbol,
            number=int(dataset.number),
            hall_symbol=str(dataset.hall),
            hall_number=int(dataset.hall_number),
            equivalent_atoms=tuple(int(value) for value in dataset.equivalent_atoms),
            shelx_latt_value=shelx_latt_value,
            symmetry_operations=tuple(symmetry_operations),
        )

    def reduce_to_asymmetric_unit(
        self,
        *,
        symprec: float = 0.05,
        angle_tolerance: float = -1.0,
    ) -> "CrystalStructure":
        """Detect the space group and keep one atom per orbit, molecules kept whole."""

        detection = self.detect_space_group_symmetry(
            symprec=symprec,
            angle_tolerance=angle_tolerance,
        )
        representative_indices = molecule_aware_representative_indices(
            self, detection.equivalent_atoms
        )

        reduced_atoms = [self.atoms[index] for index in representative_indices]
        return CrystalStructure(
            atoms=reduced_atoms,
            cell_parameters=self.cell_parameters,
            lattice_matrix=self.lattice_matrix,
            space_group=detection.symbol,
            hall_number=detection.hall_number,
            name=self.name,
            explict_unit_cell=False,
            shelx_latt_value=detection.shelx_latt_value,
            symmetry_operations=detection.symmetry_operations,
        )

    def expand_to_explicit_unit_cell(self) -> "CrystalStructure":
        """Apply LATT/SYMM (or the space group's operations) to fill the unit cell."""

        if self.explict_unit_cell:
            return CrystalStructure(
                atoms=self.atoms,
                cell_parameters=self.cell_parameters,
                lattice_matrix=self.lattice_matrix,
                space_group=self.space_group,
                hall_number=self.hall_number,
                name=self.name,
                explict_unit_cell=True,
                shelx_latt_value=self.shelx_latt_value,
                symmetry_operations=self.symmetry_operations,
            )

        if self.shelx_latt_value is not None:
            latt_value = self.shelx_latt_value
            symmetry_operations = list(self.symmetry_operations)
        else:
            latt_value, symmetry_operations = shelx_symmetry_records(
                self.space_group,
                hall_number=self.hall_number,
            )

        expanded_atoms = expand_shelx_atoms(
            atoms=self.atoms,
            cell=self.cell.array,
            latt_value=latt_value,
            symmetry_operations=symmetry_operations,
        )
        return CrystalStructure(
            atoms=expanded_atoms,
            cell_parameters=self.cell_parameters,
            lattice_matrix=self.lattice_matrix,
            space_group=self.space_group,
            hall_number=self.hall_number,
            name=self.name,
            explict_unit_cell=True,
        )

    def csorm_symmetry_sanity_check(self) -> CSORMSymmetrySanityCheck:
        """Check that CSO-RM's symmetry parser can read every SYMM operation."""

        if self.shelx_latt_value is not None:
            latt_value = self.shelx_latt_value
            symmetry_operations = tuple(self.symmetry_operations)
        else:
            latt_value, derived_operations = shelx_symmetry_records(
                self.space_group,
                hall_number=self.hall_number,
            )
            symmetry_operations = tuple(derived_operations)

        parsed_operations: list[str] = []
        unsupported_operations: list[str] = []
        unsupported_reasons: list[str] = []

        for operation in symmetry_operations:
            try:
                parsed_operations.append(normalize_csorm_symmetry_operation(operation))
            except ValueError as exc:
                unsupported_operations.append(operation)
                unsupported_reasons.append(str(exc))

        return CSORMSymmetrySanityCheck(
            supported=not unsupported_operations,
            latt_value=latt_value,
            original_operations=symmetry_operations,
            parsed_operations=tuple(parsed_operations),
            unsupported_operations=tuple(unsupported_operations),
            unsupported_reasons=tuple(unsupported_reasons),
        )

    def expand_to_csorm_explicit_unit_cell(self) -> "CrystalStructure":
        """Expand the cell exactly as CSO-RM would, failing on unsupported operations."""

        if self.explict_unit_cell:
            return CrystalStructure(
                atoms=self.atoms,
                cell_parameters=self.cell_parameters,
                lattice_matrix=self.lattice_matrix,
                space_group=self.space_group,
                hall_number=self.hall_number,
                name=self.name,
                explict_unit_cell=True,
                shelx_latt_value=self.shelx_latt_value,
                symmetry_operations=self.symmetry_operations,
            )

        sanity = self.csorm_symmetry_sanity_check()
        if not sanity.supported:
            detail = "; ".join(
                f"{operation}: {reason}"
                for operation, reason in zip(sanity.unsupported_operations, sanity.unsupported_reasons)
            )
            raise ValueError(
                "CSORM symmetry parser cannot reproduce one or more SYMM operations"
                + (f" ({detail})" if detail else "")
            )

        expanded_atoms = expand_shelx_atoms(
            atoms=self.atoms,
            cell=self.cell.array,
            latt_value=sanity.latt_value,
            symmetry_operations=sanity.parsed_operations,
        )
        return CrystalStructure(
            atoms=expanded_atoms,
            cell_parameters=self.cell_parameters,
            lattice_matrix=self.lattice_matrix,
            space_group=self.space_group,
            hall_number=self.hall_number,
            name=self.name,
            explict_unit_cell=True,
            shelx_latt_value=sanity.latt_value,
            symmetry_operations=sanity.parsed_operations,
        )

    # ------------------------------------------------------------ molecules
    def detect_molecules(
        self,
        covalent_scale: float = DEFAULT_COVALENT_SCALE,
    ) -> list[list[AtomRecord]]:
        """Detect molecular components in the unit cell.

        The bond criterion matches the workflow in `CSP-personal/6_Zmatgen`:
        covalent radii with default scale factor 1.20, but any pair involving
        hydrogen is relaxed to 1.30. Periodic connectivity is evaluated through
        ASE's neighbor list, and hydrogen is restricted to its shortest single
        covalent contact.

        This method does not deduplicate chemically identical molecules. If the
        unit cell contains multiple copies of the same molecule, each connected
        component is returned separately.
        """

        adjacency = build_connectivity(self, covalent_scale)
        components = connected_components(adjacency)
        return [
            component_records(self, component, adjacency, unwrap=True)
            for component in components
        ]

    def detect_molecules_with_template(
        self,
        template_structure: "CrystalStructure",
        covalent_scale: float = DEFAULT_COVALENT_SCALE,
    ) -> list[list[AtomRecord]]:
        """Detect molecules using another structure's bonding graph as a template.

        This is intended for cases where a relaxed structure has drifted far
        enough that distance-only covalent detection fragments a molecule, while
        a corresponding reference structure (for example, the starting POSCAR)
        still has the correct topology. Atom count and element ordering must
        match between the current and template structures.
        """

        adjacency = build_connectivity_from_template(self, template_structure, covalent_scale)
        components = connected_components(adjacency)
        return [
            component_records(self, component, adjacency, unwrap=True)
            for component in components
        ]

    def deduplicate_molecules(
        self,
        covalent_scale: float = DEFAULT_COVALENT_SCALE,
    ) -> list[MoleculeGroup]:
        """Group chemically identical molecules after raw component detection."""

        adjacency = build_connectivity(self, covalent_scale)
        components = connected_components(adjacency)
        grouped: dict[str, dict[str, object]] = {}

        for component in components:
            ordered_original = sorted(component)
            mapping = {original: local for local, original in enumerate(ordered_original)}
            local_symbols = [self.atoms[index].element for index in ordered_original]
            graph = local_graph(component, adjacency, mapping)
            signature = chemical_signature(local_symbols, graph)
            molecule = component_records(self, component, adjacency, unwrap=True)

            if signature not in grouped:
                grouped[signature] = {
                    "representative_molecule": molecule,
                    "duplicate_molecules": [molecule],
                }
            else:
                grouped[signature]["duplicate_molecules"].append(molecule)

        return [
            MoleculeGroup(
                signature=signature,
                representative_molecule=group["representative_molecule"],
                duplicate_molecules=group["duplicate_molecules"],
            )
            for signature, group in grouped.items()
        ]

    def generate_gas_phase_vasp_structure(
        self,
        molecule: list[AtomRecord],
        box_length: float = DEFAULT_GAS_PHASE_BOX_LENGTH,
        *,
        name: str | None = None,
    ) -> "CrystalStructure":
        """Create a gas-phase P1 crystal structure centered in a cubic box."""

        if box_length <= 0.0:
            raise ValueError("box_length must be positive.")
        if not molecule:
            raise ValueError("molecule must contain at least one atom.")

        coordinates = np.array([atom.coordinates for atom in molecule], dtype=float)
        centroid = coordinates.mean(axis=0)
        target = np.array([box_length / 2.0, box_length / 2.0, box_length / 2.0])
        shift = target - centroid

        centered_atoms = [
            AtomRecord(
                label=atom.label,
                element=atom.element,
                coordinates=tuple(float(value) for value in (np.array(atom.coordinates) + shift)),
            )
            for atom in molecule
        ]

        return CrystalStructure(
            atoms=centered_atoms,
            cell_parameters=(box_length, box_length, box_length, 90.0, 90.0, 90.0),
            space_group="P 1",
            name=name or f"{self.name}_gas",
            explict_unit_cell=True,
        )

    # ----------------------------------------------------------- Z-matrices
    def generate_zmatrices(
        self,
        covalent_scale: float = DEFAULT_COVALENT_SCALE,
        linear_threshold: float = DEFAULT_LINEAR_THRESHOLD,
    ) -> list[ZMatrixRepresentation]:
        """Generate one Z-matrix per detected molecule (see :mod:`.zmatrix_builder`)."""

        return generate_zmatrices(self, covalent_scale, linear_threshold)

    def apply_zmatrix_topology(
        self,
        topology: ZMatrixRepresentation,
        covalent_scale: float = DEFAULT_COVALENT_SCALE,
        linear_threshold: float = DEFAULT_LINEAR_THRESHOLD,
    ) -> ZMatrixRepresentation:
        """Measure a fixed Z-matrix topology on this structure, matching atoms by label."""

        return apply_zmatrix_topology(self, topology, covalent_scale, linear_threshold)

    # --------------------------------------------------------------- images
    def write_unit_cell_molecule_image(
        self,
        destination: str | Path,
        covalent_scale: float = DEFAULT_COVALENT_SCALE,
        *,
        title: str | None = None,
        draw_box: bool = False,
    ) -> None:
        """Write one PNG image containing all detected molecules in the unit cell."""

        write_unit_cell_molecule_image(
            self, destination, covalent_scale, title=title, draw_box=draw_box
        )

    # ------------------------------------- kept for scripts written earlier
    def _to_ase_atoms(self) -> Atoms:
        return self.to_ase_atoms()

    def _build_connectivity(self, covalent_scale: float) -> dict[int, list[BondEdge]]:
        return build_connectivity(self, covalent_scale)

    @staticmethod
    def _connected_components(adjacency: dict[int, list[BondEdge]]) -> list[list[int]]:
        return connected_components(adjacency)
