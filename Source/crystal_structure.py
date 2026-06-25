"""Crystal structure container with lightweight multi-format IO.

This module provides a single in-memory representation for periodic crystal
structures used across the current CSP workflow. The class keeps:

- atom labels
- element types
- Cartesian coordinates in Angstrom
- cell parameters (a, b, c, alpha, beta, gamma)
- a Hermann-Mauguin space-group label

Notes:
- The ``.res`` writer preserves SHELX-style ``LATT``/``SYMM`` symmetry records
  when they are available on a reduced structure, so asymmetric-unit
  round-trips can reconstruct the original cell setting.
- The provided space-group label is preserved in a ``REM SPACE_GROUP``
  comment.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import math
import re
import shlex
from pathlib import Path
from typing import Iterable
import warnings

import networkx as nx
import numpy as np
from ase import Atoms
from ase.cell import Cell
from ase.data import atomic_numbers, covalent_radii
from ase.data.colors import jmol_colors
from ase.io import read as ase_read
from ase.neighborlist import neighbor_list
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from networkx.algorithms import isomorphism as nx_isomorphism
import spglib


CellParameters = tuple[float, float, float, float, float, float]
Vector3 = tuple[float, float, float]
DEFAULT_COVALENT_SCALE = 1.20
HYDROGEN_COVALENT_SCALE = 1.30
DEFAULT_GAS_PHASE_BOX_LENGTH = 20.0
DEFAULT_LINEAR_THRESHOLD = 15.0
DEFAULT_CIF_SITE_MERGE_TOLERANCE = 1e-2
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
class _ExpandedCifAtomSite:
    label: str
    element: str
    frac: tuple[float, float, float]
    occupancy: float | None


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


@dataclass(frozen=True)
class _ZMatrixTemplateRow:
    atom_index: int
    bond_to: int | None
    angle_to: int | None
    dihedral_to: int | None
    used_fallback_angle: bool
    used_fallback_dihedral: bool


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


@dataclass
class CrystalStructure:
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
                element=_normalize_element_symbol(atom.element or atom.label),
                coordinates=atom.coordinates,
            )
            for atom in self.atoms
        ]
        self.explict_unit_cell = bool(self.explict_unit_cell)
        self.symmetry_operations = tuple(str(operation).strip() for operation in self.symmetry_operations)

    @property
    def explicit_unit_cell(self) -> bool:
        """Compatibility alias for the canonical `explict_unit_cell` flag."""

        return self.explict_unit_cell

    @property
    def cell(self) -> Cell:
        if self.lattice_matrix is not None:
            return Cell(np.array(self.lattice_matrix, dtype=float))
        return Cell.fromcellpar(self.cell_parameters)

    @classmethod
    def from_file(cls, path: str | Path, fmt: str | None = None) -> "CrystalStructure":
        file_path = Path(path)
        file_format = _normalize_format(file_path, fmt)
        if file_format == "cif":
            return cls._from_cif(file_path)
        if file_format == "pdb":
            return cls._from_pdb(file_path)
        if file_format == "res":
            return cls._from_res(file_path)
        raise ValueError(f"Unsupported format: {file_format}")

    @classmethod
    def from_cif_unit_cell(
        cls,
        path: str | Path,
        *,
        site_merge_tolerance: float = DEFAULT_CIF_SITE_MERGE_TOLERANCE,
    ) -> "CrystalStructure":
        return cls.expand_cif_to_unit_cell(
            path,
            site_merge_tolerance=site_merge_tolerance,
        )

    def to_file(
        self,
        path: str | Path,
        fmt: str | None = None,
        *,
        rounding: bool = False,
    ) -> None:
        file_path = Path(path)
        file_format = _normalize_format(file_path, fmt)
        if file_format == "cif":
            self._write_cif(file_path)
            return
        if file_format == "pdb":
            self._write_pdb(file_path)
            return
        if file_format == "res":
            self._write_res(file_path, rounding=rounding)
            return
        raise ValueError(f"Unsupported format: {file_format}")

    def fractional_coordinates(self, cartesian: Vector3) -> Vector3:
        scaled = self.cell.scaled_positions(np.array([cartesian], dtype=float))[0]
        return tuple(float(value) for value in scaled)

    def cartesian_coordinates(self, fractional: Vector3) -> Vector3:
        cart = np.dot(np.array(fractional, dtype=float), self.cell.array)
        return tuple(float(value) for value in cart)

    def detect_space_group_symmetry(
        self,
        *,
        symprec: float = 0.05,
        angle_tolerance: float = -1.0,
    ) -> SpaceGroupDetection:
        dataset = self._symmetry_dataset(symprec=symprec, angle_tolerance=angle_tolerance)
        symbol = _normalize_space_group_symbol(str(dataset.international))
        shelx_latt_value, symmetry_operations = _shelx_symmetry_records_from_symmetry(
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
        detection = self.detect_space_group_symmetry(
            symprec=symprec,
            angle_tolerance=angle_tolerance,
        )
        representative_indices: list[int] = []
        seen_representatives: set[int] = set()
        for atom_index, representative in enumerate(detection.equivalent_atoms):
            if representative in seen_representatives:
                continue
            seen_representatives.add(representative)
            representative_indices.append(atom_index)

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
            latt_value, symmetry_operations = _shelx_symmetry_records(
                self.space_group,
                hall_number=self.hall_number,
            )

        expanded_atoms = _expand_shelx_atoms(
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
        if self.shelx_latt_value is not None:
            latt_value = self.shelx_latt_value
            symmetry_operations = tuple(self.symmetry_operations)
        else:
            latt_value, derived_operations = _shelx_symmetry_records(
                self.space_group,
                hall_number=self.hall_number,
            )
            symmetry_operations = tuple(derived_operations)

        parsed_operations: list[str] = []
        unsupported_operations: list[str] = []
        unsupported_reasons: list[str] = []

        for operation in symmetry_operations:
            try:
                parsed_operations.append(_normalize_csorm_symmetry_operation(operation))
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

        expanded_atoms = _expand_shelx_atoms(
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

        adjacency = self._build_connectivity(covalent_scale)
        components = self._connected_components(adjacency)
        return [
            self._component_records(component, adjacency, unwrap=True)
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

        adjacency = self._build_connectivity_from_template(template_structure, covalent_scale)
        components = self._connected_components(adjacency)
        return [
            self._component_records(component, adjacency, unwrap=True)
            for component in components
        ]

    def deduplicate_molecules(
        self,
        covalent_scale: float = DEFAULT_COVALENT_SCALE,
    ) -> list[MoleculeGroup]:
        """Group chemically identical molecules after raw component detection."""

        adjacency = self._build_connectivity(covalent_scale)
        components = self._connected_components(adjacency)
        grouped: dict[str, dict[str, object]] = {}

        for component in components:
            ordered_original = sorted(component)
            mapping = {original: local for local, original in enumerate(ordered_original)}
            local_symbols = [self.atoms[index].element for index in ordered_original]
            graph = self._local_graph(component, adjacency, mapping)
            signature = self._chemical_signature(local_symbols, graph)
            molecule = self._component_records(component, adjacency, unwrap=True)

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

    def generate_zmatrices(
        self,
        covalent_scale: float = DEFAULT_COVALENT_SCALE,
        linear_threshold: float = DEFAULT_LINEAR_THRESHOLD,
    ) -> list[ZMatrixRepresentation]:
        """Generate one Z-matrix representation for each detected molecule."""

        adjacency = self._build_connectivity(covalent_scale)
        components = self._connected_components(adjacency)
        component_data = []
        template_cache: dict[str, dict[str, object]] = {}

        for molecule_index, component in enumerate(components, start=1):
            ordered_original = sorted(component)
            mapping = {original: local for local, original in enumerate(ordered_original)}
            graph = self._local_graph(component, adjacency, mapping)
            molecule = self._component_records(component, adjacency, unwrap=True)
            symbols = [atom.element for atom in molecule]
            signature = self._chemical_signature(symbols, graph)
            component_data.append(
                {
                    "molecule_index": molecule_index,
                    "molecule": molecule,
                    "graph": graph,
                    "symbols": symbols,
                    "signature": signature,
                }
            )

            if signature not in template_cache:
                template, warnings = self._build_zmatrix_template(
                    molecule=molecule,
                    graph=graph,
                    linear_threshold=linear_threshold,
                )
                template_cache[signature] = {
                    "template": template,
                    "warnings": warnings,
                    "graph": graph,
                    "symbols": symbols,
                }

        zmatrices: list[ZMatrixRepresentation] = []
        for item in component_data:
            signature = item["signature"]
            cached = template_cache[signature]
            template = cached["template"]
            rep_graph = cached["graph"]
            rep_symbols = cached["symbols"]
            molecule = item["molecule"]
            graph = item["graph"]
            symbols = item["symbols"]

            mapping = self._graph_isomorphism_mapping(
                representative_symbols=rep_symbols,
                representative_graph=rep_graph,
                symbols=symbols,
                graph=graph,
            )
            remapped_template = self._remap_zmatrix_template(template, mapping)
            entries, ordered_labels, warnings = self._apply_zmatrix_template(
                molecule=molecule,
                template=remapped_template,
                linear_threshold=linear_threshold,
            )
            zmatrices.append(
                ZMatrixRepresentation(
                    molecule_index=item["molecule_index"],
                    entries=entries,
                    ordered_atom_labels=ordered_labels,
                    warnings=warnings,
                )
            )

        return zmatrices

    def write_unit_cell_molecule_image(
        self,
        destination: str | Path,
        covalent_scale: float = DEFAULT_COVALENT_SCALE,
        *,
        title: str | None = None,
        draw_box: bool = False,
    ) -> None:
        """Write one PNG image containing all detected molecules in the unit cell."""

        destination_path = Path(destination)
        destination_path.parent.mkdir(parents=True, exist_ok=True)

        adjacency = self._build_connectivity(covalent_scale)
        components = self._connected_components(adjacency)
        if not components:
            raise ValueError("No molecular components detected.")

        component_positions: dict[int, np.ndarray] = {}
        for component in components:
            unwrapped = self._unwrap_component(component, adjacency)
            wrapped_coords = np.array([self.atoms[index].coordinates for index in component], dtype=float)
            wrapped_centroid = wrapped_coords.mean(axis=0)
            ordered_indices = sorted(component)
            unwrapped_coords = np.vstack([unwrapped[index] for index in ordered_indices])
            translated_coords = unwrapped_coords - unwrapped_coords.mean(axis=0) + wrapped_centroid
            for index, coord in zip(ordered_indices, translated_coords):
                component_positions[index] = coord

        ordered_indices = list(range(len(self.atoms)))
        coords = np.vstack([component_positions[index] for index in ordered_indices])
        xy, projection_center, projection_basis = self._project_for_plot(coords)

        fig, ax = plt.subplots(figsize=(6, 6), dpi=200)

        if draw_box:
            box_corners = self._unit_cell_corners()
            box_xy = self._apply_projection(box_corners, projection_center, projection_basis)
            for left, right in self._unit_cell_edges():
                ax.plot(
                    [box_xy[left, 0], box_xy[right, 0]],
                    [box_xy[left, 1], box_xy[right, 1]],
                    color="dimgray",
                    linewidth=0.8,
                    linestyle="--",
                    zorder=0,
                )

        for atom_i, neighbors in adjacency.items():
            for edge in neighbors:
                atom_j = edge.neighbor
                if atom_i < atom_j:
                    ax.plot(
                        [xy[atom_i, 0], xy[atom_j, 0]],
                        [xy[atom_i, 1], xy[atom_j, 1]],
                        color="black",
                        linewidth=1.2,
                        zorder=1,
                    )

        numbers = [int(atomic_numbers[atom.element]) for atom in self.atoms]
        facecolors = [jmol_colors[number] for number in numbers]
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            s=220,
            facecolors=facecolors,
            edgecolors="black",
            linewidths=1.0,
            zorder=2,
        )

        span_x = float(np.ptp(xy[:, 0]))
        span_y = float(np.ptp(xy[:, 1]))
        span = max(span_x, span_y, 1.0)
        label_dx = 0.02 * span
        label_dy = 0.02 * span
        font_size = 6 if len(self.atoms) > 20 else 7

        for atom_index, atom in enumerate(self.atoms):
            ax.text(
                xy[atom_index, 0] + label_dx,
                xy[atom_index, 1] + label_dy,
                atom.label,
                fontsize=font_size,
                ha="left",
                va="bottom",
                color="black",
                bbox={
                    "boxstyle": "round,pad=0.12",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.85,
                },
                zorder=3,
            )

        margin = 0.18 * span
        ax.set_xlim(float(np.min(xy[:, 0]) - margin), float(np.max(xy[:, 0]) + margin))
        ax.set_ylim(float(np.min(xy[:, 1]) - margin), float(np.max(xy[:, 1]) + margin))
        ax.set_aspect("equal", adjustable="box")
        ax.set_axis_off()
        ax.set_title(title or f"{self.name} unit-cell molecules", fontsize=9, pad=6)
        fig.tight_layout(pad=0.1)
        fig.savefig(destination_path, bbox_inches="tight", pad_inches=0.05, transparent=False)
        plt.close(fig)

    @classmethod
    def _from_cif(cls, path: Path) -> "CrystalStructure":
        data, loops = _parse_cif_blocks(path.read_text())
        cell_parameters = (
            _parse_float(data["_cell_length_a"]),
            _parse_float(data["_cell_length_b"]),
            _parse_float(data["_cell_length_c"]),
            _parse_float(data["_cell_angle_alpha"]),
            _parse_float(data["_cell_angle_beta"]),
            _parse_float(data["_cell_angle_gamma"]),
        )
        space_group = (
            data.get("_symmetry_space_group_name_H-M")
            or data.get("_space_group_name_H-M_alt")
            or "P 1"
        )
        space_group = _strip_quotes(space_group)

        headers, rows = _find_loop(
            loops,
            {"_atom_site_label", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"},
        )
        label_idx = headers.index("_atom_site_label")
        type_idx = headers.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in headers else None
        frac_x_idx = headers.index("_atom_site_fract_x")
        frac_y_idx = headers.index("_atom_site_fract_y")
        frac_z_idx = headers.index("_atom_site_fract_z")

        cell = Cell.fromcellpar(cell_parameters)
        atoms: list[AtomRecord] = []
        for row in rows:
            label = _strip_quotes(row[label_idx])
            element = (
                _strip_quotes(row[type_idx]) if type_idx is not None else _guess_element_from_label(label)
            )
            frac = (
                _parse_float(row[frac_x_idx]),
                _parse_float(row[frac_y_idx]),
                _parse_float(row[frac_z_idx]),
            )
            cart = np.dot(np.array(frac, dtype=float), cell.array)
            atoms.append(
                AtomRecord(
                    label=label,
                    element=_normalize_element_symbol(element),
                    coordinates=tuple(float(value) for value in cart),
                )
            )

        return cls(
            atoms=atoms,
            cell_parameters=cell_parameters,
            space_group=space_group,
            name=path.stem or "CrystalStructure",
            explict_unit_cell=_parse_bool_tag(data.get("_csptoolbox_explict_unit_cell"), default=False),
        )

    @classmethod
    def expand_cif_to_unit_cell(
        cls,
        path: str | Path,
        *,
        site_merge_tolerance: float = DEFAULT_CIF_SITE_MERGE_TOLERANCE,
    ) -> "CrystalStructure":
        structure, _ = cls._expand_cif_unit_cell(
            path,
            site_merge_tolerance=site_merge_tolerance,
        )
        return structure

    @classmethod
    def inspect_cif_unit_cell_expansion(
        cls,
        path: str | Path,
        *,
        site_merge_tolerance: float = DEFAULT_CIF_SITE_MERGE_TOLERANCE,
    ) -> CifExpansionReport:
        _, report = cls._expand_cif_unit_cell(
            path,
            site_merge_tolerance=site_merge_tolerance,
        )
        return report

    @classmethod
    def _expand_cif_unit_cell(
        cls,
        path: str | Path,
        *,
        site_merge_tolerance: float = DEFAULT_CIF_SITE_MERGE_TOLERANCE,
    ) -> tuple["CrystalStructure", CifExpansionReport]:
        file_path = Path(path)
        site_merge_tolerance = float(site_merge_tolerance)
        if site_merge_tolerance < 0.0:
            raise ValueError("site_merge_tolerance must be non-negative.")

        data, loops = _parse_cif_blocks(file_path.read_text())
        cell_parameters = (
            _parse_float(data["_cell_length_a"]),
            _parse_float(data["_cell_length_b"]),
            _parse_float(data["_cell_length_c"]),
            _parse_float(data["_cell_angle_alpha"]),
            _parse_float(data["_cell_angle_beta"]),
            _parse_float(data["_cell_angle_gamma"]),
        )
        lattice = _lattice_matrix(*cell_parameters)
        space_group = (
            data.get("_symmetry_space_group_name_H-M")
            or data.get("_space_group_name_H-M_alt")
            or "P 1"
        )
        space_group = _strip_quotes(space_group)

        sym_headers, sym_rows = _find_loop_with_any_header(
            loops,
            ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"),
        )
        if "_space_group_symop_operation_xyz" in sym_headers:
            sym_idx = sym_headers.index("_space_group_symop_operation_xyz")
        else:
            sym_idx = sym_headers.index("_symmetry_equiv_pos_as_xyz")
        symmetry_operations = [row[sym_idx] for row in sym_rows]

        atom_headers, atom_rows = _find_loop(
            loops,
            {"_atom_site_label", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"},
        )
        label_idx = atom_headers.index("_atom_site_label")
        type_idx = atom_headers.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in atom_headers else None
        frac_x_idx = atom_headers.index("_atom_site_fract_x")
        frac_y_idx = atom_headers.index("_atom_site_fract_y")
        frac_z_idx = atom_headers.index("_atom_site_fract_z")
        occ_idx = atom_headers.index("_atom_site_occupancy") if "_atom_site_occupancy" in atom_headers else None

        expanded_sites: list[_ExpandedCifAtomSite] = []
        duplicate_sites_merged = 0
        partial_occupancies: list[tuple[str, float]] = []
        for row in atom_rows:
            label = _strip_quotes(row[label_idx])
            element = _normalize_element_symbol(
                _strip_quotes(row[type_idx]) if type_idx is not None else _guess_element_from_label(label)
            )
            occupancy = (
                _parse_float(row[occ_idx])
                if occ_idx is not None and row[occ_idx] not in {".", "?"}
                else None
            )
            if occupancy is not None and occupancy < 0.999:
                partial_occupancies.append((label, occupancy))
            frac = (
                _parse_float(row[frac_x_idx]),
                _parse_float(row[frac_y_idx]),
                _parse_float(row[frac_z_idx]),
            )

            for operation in symmetry_operations:
                transformed = _apply_symmetry_operation(operation, frac)
                canonical = _canonicalize_fractional(transformed)
                if _matches_existing_expanded_cif_site(
                    element=element,
                    frac=canonical,
                    expanded_sites=expanded_sites,
                    lattice=lattice,
                    site_merge_tolerance=site_merge_tolerance,
                ):
                    duplicate_sites_merged += 1
                    continue
                expanded_sites.append(
                    _ExpandedCifAtomSite(
                        label=label,
                        element=element,
                        frac=canonical,
                        occupancy=occupancy,
                    )
                )

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r".*This may result in wrong setting!.*",
            )
            ase_atoms = ase_read(str(file_path))
        ase_atom_count = len(ase_atoms)
        ase_matches_manual, ase_comparison_message = _compare_ase_and_manual_expansions(
            manual_atoms=expanded_sites,
            ase_atoms=ase_atoms,
        )
        if not ase_matches_manual:
            warnings.warn(
                (
                    f"ASE expansion sanity check differs from manual CIF expansion for {file_path.name}: "
                    f"{ase_comparison_message}"
                ),
                stacklevel=2,
            )

        label_counts: dict[str, int] = {}
        atoms: list[AtomRecord] = []
        for site in expanded_sites:
            label_counts[site.label] = label_counts.get(site.label, 0) + 1
            expanded_label = f"{site.label}_{label_counts[site.label]}"
            atoms.append(
                AtomRecord(
                    label=expanded_label,
                    element=site.element,
                    coordinates=_frac_to_cart(site.frac, lattice),
                )
            )

        structure = cls(
            atoms=atoms,
            cell_parameters=cell_parameters,
            lattice_matrix=tuple(tuple(float(value) for value in vector) for vector in lattice),
            space_group=space_group,
            name=file_path.stem or "CrystalStructure",
            explict_unit_cell=True,
        )
        report = CifExpansionReport(
            name=file_path.stem or "CrystalStructure",
            space_group=space_group,
            expanded_atom_count=len(expanded_sites),
            partial_occupancies=partial_occupancies,
            ase_expands_to_unit_cell=(ase_atom_count > len(atom_rows)),
            ase_atom_count=ase_atom_count,
            raw_atom_row_count=len(atom_rows),
            ase_matches_manual=ase_matches_manual,
            ase_comparison_message=ase_comparison_message,
            site_merge_tolerance=site_merge_tolerance,
            duplicate_sites_merged=duplicate_sites_merged,
        )
        return structure, report

    @classmethod
    def _from_pdb(cls, path: Path) -> "CrystalStructure":
        lines = path.read_text().splitlines()
        cryst1 = next((line for line in lines if line.startswith("CRYST1")), None)
        if cryst1 is None:
            raise ValueError(f"No CRYST1 record found in {path}")

        cell_parameters = (
            float(cryst1[6:15].strip()),
            float(cryst1[15:24].strip()),
            float(cryst1[24:33].strip()),
            float(cryst1[33:40].strip()),
            float(cryst1[40:47].strip()),
            float(cryst1[47:54].strip()),
        )
        space_group = cryst1[55:66].strip() or "P 1"

        atoms: list[AtomRecord] = []
        for line in lines:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            label = line[12:16].strip() or f"ATOM{len(atoms) + 1}"
            element = _normalize_element_symbol(line[76:78].strip() or _guess_element_from_label(label))
            coordinates = (
                float(line[30:38].strip()),
                float(line[38:46].strip()),
                float(line[46:54].strip()),
            )
            atoms.append(AtomRecord(label=label, element=element, coordinates=coordinates))

        return cls(
            atoms=atoms,
            cell_parameters=cell_parameters,
            space_group=space_group,
            name=path.stem or "CrystalStructure",
            explict_unit_cell=_parse_pdb_explict_unit_cell(lines),
        )

    @classmethod
    def _from_res(cls, path: Path) -> "CrystalStructure":
        lines = path.read_text().splitlines()
        sfac: list[str] = []
        space_group = "P 1"
        name = path.stem or "CrystalStructure"
        explict_unit_cell = False
        shelx_latt_value: int | None = None
        symmetry_operations: list[str] = []
        cell_parameters: CellParameters | None = None
        atoms: list[AtomRecord] = []

        for line in lines:
            stripped = line.strip()
            if not stripped:
                continue
            upper = stripped.upper()
            if upper.startswith("TITL"):
                parts = stripped.split(maxsplit=1)
                if len(parts) > 1:
                    name = parts[1].strip()
                continue
            if upper.startswith("REM SPACE_GROUP"):
                space_group = stripped.split("SPACE_GROUP", maxsplit=1)[1].strip()
                continue
            if upper.startswith("REM EXPLICT_UNIT_CELL"):
                explict_unit_cell = _parse_bool_tag(
                    stripped.split("EXPLICT_UNIT_CELL", maxsplit=1)[1].strip(),
                    default=False,
                )
                continue
            if upper.startswith("CELL"):
                parts = stripped.split()
                if len(parts) < 8:
                    raise ValueError(f"Invalid CELL line in {path}: {line}")
                cell_parameters = tuple(float(value) for value in parts[2:8])  # type: ignore[assignment]
                continue
            if upper.startswith("LATT"):
                parts = stripped.split()
                if len(parts) > 1:
                    shelx_latt_value = int(parts[1])
                continue
            if upper.startswith("SYMM"):
                symmetry_operations.append(stripped.split(maxsplit=1)[1].strip())
                continue
            if upper.startswith("SFAC"):
                sfac = stripped.split()[1:]
                continue
            if upper.startswith(("LATT", "SYMM", "ZERR", "UNIT", "END", "HKLF")):
                continue

            parts = stripped.split()
            if len(parts) < 5 or cell_parameters is None:
                continue
            label = parts[0]
            raw_species = parts[1]
            if raw_species.isdigit():
                species_index = int(raw_species) - 1
                if species_index < 0 or species_index >= len(sfac):
                    raise ValueError(f"SFAC index out of range in {path}: {line}")
                element = sfac[species_index]
                frac = tuple(float(value) for value in parts[2:5])
            else:
                element = raw_species
                frac = tuple(float(value) for value in parts[2:5])

            cell = Cell.fromcellpar(cell_parameters)
            cart = np.dot(np.array(frac, dtype=float), cell.array)
            atoms.append(
                AtomRecord(
                    label=label,
                    element=_normalize_element_symbol(element),
                    coordinates=tuple(float(value) for value in cart),
                )
            )

        if cell_parameters is None:
            raise ValueError(f"No CELL record found in {path}")

        return cls(
            atoms=atoms,
            cell_parameters=cell_parameters,
            space_group=space_group,
            name=name,
            explict_unit_cell=explict_unit_cell,
            shelx_latt_value=shelx_latt_value,
            symmetry_operations=tuple(symmetry_operations),
        )

    def _write_cif(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        a, b, c, alpha, beta, gamma = self.cell_parameters
        lines = [
            f"data_{_safe_data_name(self.name)}",
            "_audit_creation_method 'CSPToolbox CrystalStructure'",
            f"_csptoolbox_explict_unit_cell {'true' if self.explict_unit_cell else 'false'}",
            f"_cell_length_a {a:.10f}",
            f"_cell_length_b {b:.10f}",
            f"_cell_length_c {c:.10f}",
            f"_cell_angle_alpha {alpha:.10f}",
            f"_cell_angle_beta {beta:.10f}",
            f"_cell_angle_gamma {gamma:.10f}",
            f"_symmetry_space_group_name_H-M '{self.space_group}'",
            "loop_",
            "_atom_site_label",
            "_atom_site_type_symbol",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
            "_atom_site_occupancy",
        ]
        for atom in self.atoms:
            frac = self.fractional_coordinates(atom.coordinates)
            lines.append(
                f"{atom.label} {atom.element} "
                f"{frac[0]:.10f} {frac[1]:.10f} {frac[2]:.10f} 1.0"
            )
        path.write_text("\n".join(lines) + "\n")

    def _write_pdb(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        a, b, c, alpha, beta, gamma = self.cell_parameters
        lines = [
            f"HEADER    {self.name}",
            f"REMARK   1 CSPTOOLBOX_EXPLICT_UNIT_CELL {'TRUE' if self.explict_unit_cell else 'FALSE'}",
            f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} "
            f"{self.space_group:<11}{1:>4d}",
        ]
        for index, atom in enumerate(self.atoms, start=1):
            x, y, z = atom.coordinates
            atom_name = atom.label[:4]
            lines.append(
                f"HETATM{index:5d} {atom_name:<4s} MOL A   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          "
                f"{atom.element:>2s}"
            )
        lines.extend(["END", ""])
        path.write_text("\n".join(lines))

    def _write_res(self, path: Path, *, rounding: bool = False) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        ordered_elements = _element_order(atom.element for atom in self.atoms)
        if self.shelx_latt_value is not None:
            latt_value = self.shelx_latt_value
            symmetry_operations = _deduplicate_shelx_symmetry_operations(
                self.symmetry_operations,
                latt_value=latt_value,
                rounding=rounding,
            )
        else:
            latt_value, symmetry_operations = _shelx_symmetry_records(
                self.space_group,
                hall_number=self.hall_number,
                rounding=rounding,
            )
        lines = [
            f"TITL {self.name}",
            f"REM SPACE_GROUP {self.space_group}",
            f"REM EXPLICT_UNIT_CELL {'TRUE' if self.explict_unit_cell else 'FALSE'}",
            "CELL 1.54184 "
            f"{self.cell_parameters[0]:.10f} {self.cell_parameters[1]:.10f} {self.cell_parameters[2]:.10f} "
            f"{self.cell_parameters[3]:.10f} {self.cell_parameters[4]:.10f} {self.cell_parameters[5]:.10f}",
            f"LATT {latt_value}",
        ]
        lines.extend(f"SYMM {operation}" for operation in symmetry_operations)
        lines.append("SFAC " + " ".join(ordered_elements))
        for atom in self.atoms:
            frac = self.fractional_coordinates(atom.coordinates)
            sfac_index = ordered_elements.index(atom.element) + 1
            lines.append(
                f"{atom.label:<8s} {sfac_index:2d} "
                f"{frac[0]:.10f} {frac[1]:.10f} {frac[2]:.10f} 11.00000 0.05000"
            )
        lines.extend(["END", ""])
        path.write_text("\n".join(lines))

    def _to_ase_atoms(self) -> Atoms:
        return Atoms(
            symbols=[atom.element for atom in self.atoms],
            positions=np.array([atom.coordinates for atom in self.atoms], dtype=float),
            cell=self.cell.array,
            pbc=True,
        )

    def _spglib_cell(self) -> tuple[np.ndarray, np.ndarray, list[int]]:
        atoms = self._to_ase_atoms()
        return (
            np.asarray(atoms.cell, dtype=float),
            np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float),
            [int(atomic_numbers[atom.element]) for atom in self.atoms],
        )

    def _symmetry_dataset(
        self,
        *,
        symprec: float,
        angle_tolerance: float,
    ) -> spglib.SpglibDataset:
        dataset = spglib.get_symmetry_dataset(
            self._spglib_cell(),
            symprec=symprec,
            angle_tolerance=angle_tolerance,
        )
        if dataset is None:
            raise ValueError("Could not determine space-group symmetry for this structure.")
        return dataset

    def _build_connectivity(self, covalent_scale: float) -> dict[int, list[BondEdge]]:
        atoms = self._to_ase_atoms()
        numbers = [int(atomic_numbers[atom.element]) for atom in self.atoms]

        radii: list[float] = []
        for number in numbers:
            radius = float(covalent_radii[number])
            if np.isnan(radius) or radius <= 0.0:
                raise ValueError(f"No covalent radius available for atomic number {number}.")
            radii.append(radius)

        max_scale = max(float(covalent_scale), HYDROGEN_COVALENT_SCALE)
        cutoffs = np.asarray(radii) * max_scale
        i_list, j_list, shifts = neighbor_list(
            "ijS",
            atoms,
            cutoff=cutoffs,
            self_interaction=False,
        )

        positions = atoms.get_positions()
        cell = np.asarray(atoms.cell)
        best_edges: dict[tuple[int, int], tuple[tuple[int, int, int], float, float]] = {}

        for left, right, shift in zip(i_list, j_list, shifts):
            if left == right:
                continue
            shift_tuple = tuple(int(value) for value in shift)
            if left < right:
                key = (int(left), int(right))
                stored_shift = shift_tuple
            else:
                key = (int(right), int(left))
                stored_shift = tuple(-value for value in shift_tuple)

            displacement = positions[key[1]] + np.dot(stored_shift, cell) - positions[key[0]]
            distance = float(np.linalg.norm(displacement))
            denom = radii[key[0]] + radii[key[1]]
            cutoff_distance = denom * self._pair_covalent_scale(
                numbers[key[0]],
                numbers[key[1]],
                covalent_scale,
            )
            if distance > cutoff_distance:
                continue
            ratio = distance / denom if denom > 0.0 else float("inf")

            current = best_edges.get(key)
            if current is None or distance < current[1]:
                best_edges[key] = (stored_shift, distance, ratio)

        allowed_keys = set(best_edges)
        hydrogen_indices = [index for index, number in enumerate(numbers) if number == 1]
        for hydrogen_index in hydrogen_indices:
            incident = [
                (key, best_edges[key])
                for key in allowed_keys
                if hydrogen_index in key
            ]
            if len(incident) <= 1:
                continue
            best_key, _ = min(
                incident,
                key=lambda item: (item[1][1], item[1][2], item[0]),
            )
            for key, _ in incident:
                if key != best_key:
                    allowed_keys.discard(key)

        adjacency = {index: [] for index in range(len(self.atoms))}
        for left, right in sorted(allowed_keys):
            shift, distance, ratio = best_edges[(left, right)]
            adjacency[left].append(BondEdge(right, shift, distance, ratio))
            adjacency[right].append(
                BondEdge(left, tuple(-value for value in shift), distance, ratio)
            )
        return adjacency

    def _build_connectivity_from_template(
        self,
        template_structure: "CrystalStructure",
        covalent_scale: float,
    ) -> dict[int, list[BondEdge]]:
        if len(self.atoms) != len(template_structure.atoms):
            raise ValueError(
                "Template-based molecule detection requires the same atom count in "
                "the current and template structures."
            )

        current_elements = [atom.element for atom in self.atoms]
        template_elements = [atom.element for atom in template_structure.atoms]
        if current_elements != template_elements:
            raise ValueError(
                "Template-based molecule detection requires the same atom ordering "
                "and element sequence in the current and template structures."
            )

        template_adjacency = template_structure._build_connectivity(covalent_scale)
        positions = np.array([atom.coordinates for atom in self.atoms], dtype=float)
        cell = np.asarray(self.cell.array, dtype=float)
        fractional = self.cell.scaled_positions(positions)
        numbers = [int(atomic_numbers[atom.element]) for atom in self.atoms]
        radii = [float(covalent_radii[number]) for number in numbers]

        adjacency = {index: [] for index in range(len(self.atoms))}
        added_keys: set[tuple[int, int]] = set()

        for left, edges in template_adjacency.items():
            for edge in edges:
                right = edge.neighbor
                key = (min(left, right), max(left, right))
                if key in added_keys:
                    continue
                added_keys.add(key)

                delta = fractional[key[1]] - fractional[key[0]]
                shift = tuple(int(value) for value in -np.rint(delta))
                displacement = positions[key[1]] + np.dot(shift, cell) - positions[key[0]]
                distance = float(np.linalg.norm(displacement))
                denom = radii[key[0]] + radii[key[1]]
                ratio = distance / denom if denom > 0.0 else float("inf")

                adjacency[key[0]].append(BondEdge(key[1], shift, distance, ratio))
                adjacency[key[1]].append(
                    BondEdge(key[0], tuple(-value for value in shift), distance, ratio)
                )

        return adjacency

    @staticmethod
    def _connected_components(adjacency: dict[int, list[BondEdge]]) -> list[list[int]]:
        remaining = set(adjacency)
        components: list[list[int]] = []
        while remaining:
            root = min(remaining)
            stack = [root]
            component: list[int] = []
            remaining.remove(root)
            while stack:
                node = stack.pop()
                component.append(node)
                for edge in adjacency[node]:
                    if edge.neighbor in remaining:
                        remaining.remove(edge.neighbor)
                        stack.append(edge.neighbor)
            components.append(sorted(component))
        return components

    @staticmethod
    def _pair_covalent_scale(
        left_number: int,
        right_number: int,
        covalent_scale: float,
    ) -> float:
        if left_number == 1 or right_number == 1:
            return max(float(covalent_scale), HYDROGEN_COVALENT_SCALE)
        return float(covalent_scale)

    def _component_records(
        self,
        component: list[int],
        adjacency: dict[int, list[BondEdge]],
        *,
        unwrap: bool,
    ) -> list[AtomRecord]:
        if not unwrap:
            return [self.atoms[index] for index in component]

        unwrapped = self._unwrap_component(component, adjacency)
        ordered_indices = sorted(component)
        return [
            AtomRecord(
                label=self.atoms[index].label,
                element=self.atoms[index].element,
                coordinates=tuple(float(value) for value in unwrapped[index]),
            )
            for index in ordered_indices
        ]

    def _unwrap_component(
        self,
        component: list[int],
        adjacency: dict[int, list[BondEdge]],
    ) -> dict[int, np.ndarray]:
        component_set = set(component)
        root = component[0]
        cell = np.asarray(self.cell.array)
        wrapped_positions = np.array([atom.coordinates for atom in self.atoms], dtype=float)

        unwrapped: dict[int, np.ndarray] = {root: wrapped_positions[root].copy()}
        stack = [root]
        while stack:
            node = stack.pop()
            base = unwrapped[node]
            for edge in adjacency[node]:
                if edge.neighbor not in component_set or edge.neighbor in unwrapped:
                    continue
                displacement = (
                    wrapped_positions[edge.neighbor]
                    + np.dot(edge.shift, cell)
                    - wrapped_positions[node]
                )
                unwrapped[edge.neighbor] = base + displacement
                stack.append(edge.neighbor)
        return unwrapped

    @staticmethod
    def _local_graph(
        component: list[int],
        adjacency: dict[int, list[BondEdge]],
        mapping: dict[int, int],
    ) -> dict[int, set[int]]:
        component_set = set(component)
        graph = {mapping[index]: set() for index in component}
        for original in component:
            local = mapping[original]
            for edge in adjacency[original]:
                if edge.neighbor in component_set:
                    graph[local].add(mapping[edge.neighbor])
        return graph

    @staticmethod
    def _chemical_signature(symbols: list[str], graph: dict[int, set[int]]) -> str:
        labels: dict[int, object] = {
            atom_index: (symbols[atom_index], len(graph[atom_index]))
            for atom_index in graph
        }
        for _ in range(len(graph)):
            updated: dict[int, object] = {}
            pattern_to_id: dict[tuple[object, ...], int] = {}
            for atom_index in graph:
                pattern = (
                    labels[atom_index],
                    tuple(sorted(labels[neighbor] for neighbor in graph[atom_index])),
                )
                if pattern not in pattern_to_id:
                    pattern_to_id[pattern] = len(pattern_to_id)
                updated[atom_index] = pattern_to_id[pattern]
            if updated == labels:
                break
            labels = updated

        edge_labels: list[tuple[object, object]] = []
        for atom_index, neighbors in graph.items():
            for neighbor in neighbors:
                if atom_index < neighbor:
                    edge_labels.append(tuple(sorted((labels[atom_index], labels[neighbor]))))

        atom_terms = sorted(str(value) for value in labels.values())
        edge_terms = sorted(f"{left}--{right}" for left, right in edge_labels)
        return "atoms:" + ";".join(atom_terms) + "|edges:" + ";".join(edge_terms)

    def _build_zmatrix_template(
        self,
        molecule: list[AtomRecord],
        graph: dict[int, set[int]],
        linear_threshold: float,
    ) -> tuple[list[_ZMatrixTemplateRow], list[str]]:
        coords = np.array([atom.coordinates for atom in molecule], dtype=float)
        symbols = [atom.element for atom in molecule]
        order, parent = self._atom_order(symbols, coords, graph)
        order, branch_kept_child = self._branch_improper_atom_order(
            order,
            parent,
            symbols,
        )
        defined: set[int] = set()
        template: list[_ZMatrixTemplateRow] = []
        warnings: list[str] = []

        for step, atom_i in enumerate(order):
            defined.add(atom_i)
            if step == 0:
                template.append(_ZMatrixTemplateRow(atom_i, None, None, None, False, False))
                continue

            atom_b = parent[atom_i]
            if atom_b is None or atom_b not in graph[atom_i]:
                raise ValueError(f"Atom index {atom_i} is missing a bonded parent reference.")

            if step == 1:
                template.append(_ZMatrixTemplateRow(atom_i, atom_b, None, None, False, False))
                continue

            defined_refs = defined - {atom_i}
            branch_improper = None
            if step >= 3:
                branch_improper = self._choose_branch_improper_reference(
                    atom_i,
                    atom_b,
                    defined_refs,
                    parent,
                    branch_kept_child,
                    symbols,
                    coords,
                    graph,
                )
            used_fallback_angle = False
            used_fallback_dihedral = False
            if branch_improper is not None:
                atom_a, atom_d = branch_improper
                angle_deg = self._angle_value(coords, atom_i, atom_b, atom_a)
                dihedral_deg = self._dihedral_value(coords, atom_i, atom_b, atom_a, atom_d)
            else:
                atom_a, angle_deg, angle_is_bonded_path = self._choose_angle_reference(
                    atom_i,
                    atom_b,
                    defined_refs,
                    parent,
                    symbols,
                    coords,
                    graph,
                    linear_threshold,
                )
                if not angle_is_bonded_path:
                    used_fallback_angle = True

            angle_margin = self._linear_margin(angle_deg)
            if angle_margin < linear_threshold:
                warnings.append(f"Atom index {atom_i} uses a near-linear angle reference.")

            if step == 2:
                template.append(_ZMatrixTemplateRow(atom_i, atom_b, atom_a, None, used_fallback_angle, False))
                continue

            if branch_improper is None:
                atom_d, dihedral_deg, dihedral_is_bonded_path = self._choose_dihedral_reference(
                    atom_i,
                    atom_b,
                    atom_a,
                    defined_refs,
                    parent,
                    symbols,
                    coords,
                    graph,
                    linear_threshold,
                )
                if not dihedral_is_bonded_path:
                    used_fallback_dihedral = True

            dihedral_margin = self._linear_margin(self._angle_value(coords, atom_b, atom_a, atom_d))
            if dihedral_margin < linear_threshold:
                warnings.append(f"Atom index {atom_i} uses a near-linear dihedral anchor.")

            template.append(
                _ZMatrixTemplateRow(
                    atom_i,
                    atom_b,
                    atom_a,
                    atom_d,
                    used_fallback_angle,
                    used_fallback_dihedral,
                )
            )

        return template, warnings

    def _apply_zmatrix_template(
        self,
        molecule: list[AtomRecord],
        template: list[_ZMatrixTemplateRow],
        linear_threshold: float,
    ) -> tuple[list[ZMatrixEntry], list[str], list[str]]:
        coords = np.array([atom.coordinates for atom in molecule], dtype=float)
        labels = [atom.label for atom in molecule]
        symbols = [atom.element for atom in molecule]
        index_in_order = {row.atom_index: position + 1 for position, row in enumerate(template)}

        entries: list[ZMatrixEntry] = []
        warnings: list[str] = []

        for row in template:
            atom_i = row.atom_index
            label = labels[atom_i]
            element = symbols[atom_i]

            bond_length = None
            angle_degrees = None
            dihedral_degrees = None

            if row.bond_to is not None:
                bond_length = float(np.linalg.norm(coords[atom_i] - coords[row.bond_to]))
            if row.angle_to is not None and row.bond_to is not None:
                angle_degrees = self._angle_value(coords, atom_i, row.bond_to, row.angle_to)
            if row.dihedral_to is not None and row.angle_to is not None and row.bond_to is not None:
                dihedral_degrees = self._dihedral_value(
                    coords,
                    atom_i,
                    row.bond_to,
                    row.angle_to,
                    row.dihedral_to,
                )

            if row.used_fallback_angle:
                warnings.append(f"Atom {label} used a fallback non-bonded angle reference.")
            if row.used_fallback_dihedral:
                warnings.append(f"Atom {label} used a fallback non-bonded dihedral reference.")
            if angle_degrees is not None and self._linear_margin(angle_degrees) < linear_threshold:
                warnings.append(f"Atom {label} uses a near-linear angle reference ({angle_degrees:.2f} deg).")
            if (
                row.dihedral_to is not None
                and row.angle_to is not None
                and row.bond_to is not None
                and self._linear_margin(self._angle_value(coords, row.bond_to, row.angle_to, row.dihedral_to))
                < linear_threshold
            ):
                warnings.append(f"Atom {label} uses a near-linear dihedral anchor.")

            entries.append(
                ZMatrixEntry(
                    label=label,
                    element=element,
                    bond_to=index_in_order[row.bond_to] if row.bond_to is not None else None,
                    bond_length=bond_length,
                    angle_to=index_in_order[row.angle_to] if row.angle_to is not None else None,
                    angle_degrees=angle_degrees,
                    dihedral_to=index_in_order[row.dihedral_to] if row.dihedral_to is not None else None,
                    dihedral_degrees=dihedral_degrees,
                )
            )

        ordered_labels = [labels[row.atom_index] for row in template]
        return entries, ordered_labels, warnings

    @staticmethod
    def _remap_zmatrix_template(
        template: list[_ZMatrixTemplateRow],
        mapping: dict[int, int],
    ) -> list[_ZMatrixTemplateRow]:
        return [
            _ZMatrixTemplateRow(
                atom_index=mapping[row.atom_index],
                bond_to=mapping[row.bond_to] if row.bond_to is not None else None,
                angle_to=mapping[row.angle_to] if row.angle_to is not None else None,
                dihedral_to=mapping[row.dihedral_to] if row.dihedral_to is not None else None,
                used_fallback_angle=row.used_fallback_angle,
                used_fallback_dihedral=row.used_fallback_dihedral,
            )
            for row in template
        ]

    def _graph_isomorphism_mapping(
        self,
        representative_symbols: list[str],
        representative_graph: dict[int, set[int]],
        symbols: list[str],
        graph: dict[int, set[int]],
    ) -> dict[int, int]:
        if representative_symbols == symbols and representative_graph == graph:
            return {index: index for index in range(len(symbols))}

        left = self._graph_to_networkx(representative_symbols, representative_graph)
        right = self._graph_to_networkx(symbols, graph)
        matcher = nx_isomorphism.GraphMatcher(
            left,
            right,
            node_match=lambda a, b: a["element"] == b["element"],
        )
        try:
            return next(matcher.isomorphisms_iter())
        except StopIteration as error:
            raise ValueError("Could not find an isomorphism between chemically identical molecules.") from error

    @staticmethod
    def _graph_to_networkx(symbols: list[str], graph: dict[int, set[int]]) -> nx.Graph:
        nx_graph = nx.Graph()
        for atom_index, symbol in enumerate(symbols):
            nx_graph.add_node(atom_index, element=symbol)
        for atom_index, neighbors in graph.items():
            for neighbor in neighbors:
                if atom_index < neighbor:
                    nx_graph.add_edge(atom_index, neighbor)
        return nx_graph

    @staticmethod
    def _project_for_plot(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if len(coords) == 0:
            return np.zeros((0, 2)), np.zeros(3), np.eye(3)
        if len(coords) <= 2:
            centered = coords - coords.mean(axis=0)
            basis = np.eye(3)
            return centered[:, :2], coords.mean(axis=0), basis

        centered = coords - coords.mean(axis=0)
        _, _, right_vectors = np.linalg.svd(centered, full_matrices=False)
        rotated = centered @ right_vectors.T
        return rotated[:, :2], coords.mean(axis=0), right_vectors

    @staticmethod
    def _apply_projection(
        coords: np.ndarray,
        center: np.ndarray,
        basis: np.ndarray,
    ) -> np.ndarray:
        centered = coords - center
        rotated = centered @ basis.T
        return rotated[:, :2]

    def _unit_cell_corners(self) -> np.ndarray:
        a_vec, b_vec, c_vec = np.asarray(self.cell.array)
        origin = np.zeros(3, dtype=float)
        return np.array(
            [
                origin,
                a_vec,
                b_vec,
                c_vec,
                a_vec + b_vec,
                a_vec + c_vec,
                b_vec + c_vec,
                a_vec + b_vec + c_vec,
            ],
            dtype=float,
        )

    @staticmethod
    def _unit_cell_edges() -> list[tuple[int, int]]:
        return [
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 4),
            (1, 5),
            (2, 4),
            (2, 6),
            (3, 5),
            (3, 6),
            (4, 7),
            (5, 7),
            (6, 7),
        ]

    @staticmethod
    def _choose_root(symbols: list[str], coords: np.ndarray) -> int:
        centroid = coords.mean(axis=0)
        heavy = [index for index, symbol in enumerate(symbols) if symbol != "H"]
        pool = heavy or list(range(len(symbols)))
        return min(
            pool,
            key=lambda index: (
                ELEMENT_PRIORITY.get(symbols[index].upper(), 7),
                float(np.linalg.norm(coords[index] - centroid)),
                index,
            ),
        )

    @staticmethod
    def _traversal_priority(
        atom_index: int,
        symbols: list[str],
        graph: dict[int, set[int]],
        coords: np.ndarray,
    ) -> tuple[int, int, float, int]:
        return (
            ELEMENT_PRIORITY.get(symbols[atom_index].upper(), 7),
            -len(graph[atom_index]),
            float(np.linalg.norm(coords[atom_index] - coords.mean(axis=0))),
            atom_index,
        )

    def _atom_order(
        self,
        symbols: list[str],
        coords: np.ndarray,
        graph: dict[int, set[int]],
    ) -> tuple[list[int], dict[int, int | None]]:
        heavy_atoms = [index for index, symbol in enumerate(symbols) if symbol != "H"]
        if not heavy_atoms:
            root = self._choose_root(symbols, coords)
            visited = {root}
            parent: dict[int, int | None] = {root: None}
            order = [root]
            frontier = {neighbor: root for neighbor in graph[root]}

            while frontier:
                next_atom = min(
                    frontier,
                    key=lambda index: self._traversal_priority(index, symbols, graph, coords),
                )
                next_parent = frontier.pop(next_atom)
                if next_atom in visited:
                    continue
                visited.add(next_atom)
                parent[next_atom] = next_parent
                order.append(next_atom)
                for neighbor in graph[next_atom]:
                    if neighbor not in visited and neighbor not in frontier:
                        frontier[neighbor] = next_atom
            return order, parent

        heavy_graph = {
            atom_index: {neighbor for neighbor in graph[atom_index] if symbols[neighbor] != "H"}
            for atom_index in heavy_atoms
        }
        root = self._choose_root(symbols, coords)
        visited_heavy = {root}
        parent: dict[int, int | None] = {root: None}
        heavy_order = [root]
        frontier = {neighbor: root for neighbor in heavy_graph[root]}

        while frontier:
            next_atom = min(
                frontier,
                key=lambda index: self._traversal_priority(index, symbols, graph, coords),
            )
            next_parent = frontier.pop(next_atom)
            if next_atom in visited_heavy:
                continue
            visited_heavy.add(next_atom)
            parent[next_atom] = next_parent
            heavy_order.append(next_atom)
            for neighbor in heavy_graph[next_atom]:
                if neighbor not in visited_heavy and neighbor not in frontier:
                    frontier[neighbor] = next_atom

        order = list(heavy_order)
        for atom_index in heavy_order:
            attached_h = sorted(
                (neighbor for neighbor in graph[atom_index] if symbols[neighbor] == "H"),
                key=lambda index: index,
            )
            for hydrogen_index in attached_h:
                if hydrogen_index not in parent:
                    parent[hydrogen_index] = atom_index
                    order.append(hydrogen_index)

        return order, parent

    def _branch_improper_atom_order(
        self,
        order: list[int],
        parent: dict[int, int | None],
        symbols: list[str],
    ) -> tuple[list[int], dict[int, int]]:
        """Reorder a parent tree so the retained branch child is emitted first."""

        if not order:
            return order, {}

        children_by_parent: dict[int, list[int]] = {atom_index: [] for atom_index in order}
        for atom_index in order:
            parent_index = parent.get(atom_index)
            if parent_index is not None:
                children_by_parent.setdefault(parent_index, []).append(atom_index)

        original_position = {atom_index: position for position, atom_index in enumerate(order)}
        subtree_size_cache: dict[int, int] = {}

        def subtree_size(atom_index: int) -> int:
            if atom_index in subtree_size_cache:
                return subtree_size_cache[atom_index]
            size = 1 + sum(subtree_size(child) for child in children_by_parent.get(atom_index, []))
            subtree_size_cache[atom_index] = size
            return size

        kept_child: dict[int, int] = {}
        for center, children in children_by_parent.items():
            if len(children) <= 1:
                continue
            kept_child[center] = max(
                children,
                key=lambda child: (
                    symbols[child] != "H",
                    subtree_size(child),
                    -original_position[child],
                ),
            )

        def child_sort_key(center: int, child: int) -> tuple[int, int, int, int]:
            keep = kept_child.get(center)
            return (
                0 if child == keep else 1,
                0 if symbols[child] == "H" else 1,
                subtree_size(child),
                original_position[child],
            )

        reordered: list[int] = []
        visited: set[int] = set()

        def visit_once(atom_index: int) -> None:
            if atom_index in visited:
                return
            visited.add(atom_index)
            reordered.append(atom_index)
            for child in sorted(
                children_by_parent.get(atom_index, []),
                key=lambda child: child_sort_key(atom_index, child),
            ):
                visit_once(child)

        visit_once(order[0])
        missing = [atom_index for atom_index in order if atom_index not in visited]
        reordered.extend(missing)
        return reordered, kept_child

    @staticmethod
    def _angle_value(coords: np.ndarray, atom_i: int, atom_j: int, atom_k: int) -> float:
        vec_ji = coords[atom_i] - coords[atom_j]
        vec_jk = coords[atom_k] - coords[atom_j]
        denom = np.linalg.norm(vec_ji) * np.linalg.norm(vec_jk)
        if denom < 1e-12:
            return float("nan")
        cosine = float(np.clip(np.dot(vec_ji, vec_jk) / denom, -1.0, 1.0))
        return float(np.degrees(np.arccos(cosine)))

    @staticmethod
    def _dihedral_value(
        coords: np.ndarray,
        atom_i: int,
        atom_j: int,
        atom_k: int,
        atom_l: int,
    ) -> float:
        p0, p1, p2, p3 = coords[[atom_i, atom_j, atom_k, atom_l]]
        b0 = p0 - p1
        b1 = p2 - p1
        b2 = p3 - p2

        b1_norm = np.linalg.norm(b1)
        if b1_norm < 1e-12:
            return float("nan")
        b1_unit = b1 / b1_norm

        v = b0 - np.dot(b0, b1_unit) * b1_unit
        w = b2 - np.dot(b2, b1_unit) * b1_unit
        v_norm = np.linalg.norm(v)
        w_norm = np.linalg.norm(w)
        if v_norm < 1e-12 or w_norm < 1e-12:
            return float("nan")

        x_value = np.dot(v, w)
        y_value = np.dot(np.cross(b1_unit, v), w)
        return float(np.degrees(np.arctan2(y_value, x_value)))

    @staticmethod
    def _linear_margin(angle_degrees: float) -> float:
        if np.isnan(angle_degrees):
            return -1.0
        return float(min(abs(angle_degrees), abs(180.0 - angle_degrees)))

    @staticmethod
    def _unique(sequence: Iterable[int]) -> list[int]:
        seen: set[int] = set()
        result: list[int] = []
        for value in sequence:
            if value not in seen:
                seen.add(value)
                result.append(value)
        return result

    @staticmethod
    def _candidate_priority(
        candidate: int,
        anchor: int,
        symbols: list[str],
        graph: dict[int, set[int]],
        coords: np.ndarray,
    ) -> tuple[bool, int, float, int]:
        return (
            symbols[candidate] != "H",
            len(graph[candidate]),
            -float(np.linalg.norm(coords[candidate] - coords[anchor])),
            -candidate,
        )

    def _choose_angle_reference(
        self,
        atom_i: int,
        atom_b: int,
        defined: set[int],
        parent: dict[int, int | None],
        symbols: list[str],
        coords: np.ndarray,
        graph: dict[int, set[int]],
        threshold: float,
    ) -> tuple[int, float, bool]:
        bonded_candidates: list[int] = []
        if parent.get(atom_b) is not None and parent[atom_b] in defined:
            bonded_candidates.append(parent[atom_b])  # type: ignore[arg-type]
        bonded_candidates.extend(
            sorted(
                (graph[atom_b] & defined) - {atom_b, atom_i},
                key=lambda index: self._candidate_priority(index, atom_b, symbols, graph, coords),
                reverse=True,
            )
        )
        fallback_candidates = sorted(
            defined - {atom_b, atom_i} - set(bonded_candidates),
            key=lambda index: self._candidate_priority(index, atom_b, symbols, graph, coords),
            reverse=True,
        )

        best_candidate = None
        best_margin = -1.0
        for candidate in self._unique(bonded_candidates):
            angle_deg = self._angle_value(coords, atom_i, atom_b, candidate)
            margin = self._linear_margin(angle_deg)
            if margin > best_margin:
                best_margin = margin
                best_candidate = candidate
        if best_candidate is not None:
            return best_candidate, self._angle_value(coords, atom_i, atom_b, best_candidate), True

        best_candidate = None
        best_margin = -1.0
        for candidate in self._unique(fallback_candidates):
            angle_deg = self._angle_value(coords, atom_i, atom_b, candidate)
            margin = self._linear_margin(angle_deg)
            if margin > best_margin:
                best_margin = margin
                best_candidate = candidate
            if margin >= threshold:
                return candidate, angle_deg, False

        if best_candidate is None:
            raise ValueError(f"Failed to choose an angle reference for atom {atom_i}.")
        return best_candidate, self._angle_value(coords, atom_i, atom_b, best_candidate), False

    def _choose_branch_improper_reference(
        self,
        atom_i: int,
        atom_b: int,
        defined: set[int],
        parent: dict[int, int | None],
        branch_kept_child: dict[int, int],
        symbols: list[str],
        coords: np.ndarray,
        graph: dict[int, set[int]],
    ) -> tuple[int, int] | None:
        kept_child = branch_kept_child.get(atom_b)
        if kept_child is None or kept_child == atom_i or kept_child not in defined:
            return None

        candidate_refs: list[int] = []
        parent_b = parent.get(atom_b)
        if parent_b is not None and parent_b in graph[atom_b] and parent_b in defined:
            candidate_refs.append(parent_b)
        candidate_refs.append(kept_child)
        candidate_refs.extend(
            sorted(
                (graph[atom_b] & defined) - {atom_i, atom_b, parent_b, kept_child},
                key=lambda index: self._branch_improper_reference_priority(
                    index,
                    atom_b,
                    kept_child,
                    symbols,
                    coords,
                ),
            )
        )

        refs = self._unique(candidate_refs)
        if len(refs) < 2:
            return None

        best_pair: tuple[int, int] | None = None
        best_score: tuple[bool, float, float, int, int] | None = None
        for left_index, left in enumerate(refs):
            for right_index, right in enumerate(refs[left_index + 1:], start=left_index + 1):
                angle_deg = self._angle_value(coords, atom_i, atom_b, left)
                anchor_angle = self._angle_value(coords, atom_b, left, right)
                score = (
                    right not in graph[left],
                    self._linear_margin(angle_deg),
                    self._linear_margin(anchor_angle),
                    -left_index,
                    -right_index,
                )
                if best_score is None or score > best_score:
                    best_score = score
                    best_pair = (left, right)

        return best_pair

    @staticmethod
    def _branch_improper_reference_priority(
        atom_index: int,
        center: int,
        kept_child: int,
        symbols: list[str],
        coords: np.ndarray,
    ) -> tuple[int, float, int]:
        if atom_index == kept_child:
            group = 0
        elif symbols[atom_index] != "H":
            group = 1
        else:
            group = 2
        return (
            group,
            float(np.linalg.norm(coords[atom_index] - coords[center])),
            atom_index,
        )

    def _choose_dihedral_reference(
        self,
        atom_i: int,
        atom_b: int,
        atom_a: int,
        defined: set[int],
        parent: dict[int, int | None],
        symbols: list[str],
        coords: np.ndarray,
        graph: dict[int, set[int]],
        threshold: float,
    ) -> tuple[int, float, bool]:
        bonded_candidates: list[int] = []
        parent_a = parent.get(atom_a)
        if parent_a is not None and parent_a in defined and parent_a not in {atom_i, atom_b}:
            bonded_candidates.append(parent_a)
        bonded_candidates.extend(
            sorted(
                (graph[atom_a] & defined) - {atom_i, atom_b, atom_a},
                key=lambda index: self._candidate_priority(index, atom_a, symbols, graph, coords),
                reverse=True,
            )
        )
        center_bonded_exclusions = (graph[atom_b] & defined) - {atom_i, atom_b, atom_a}
        fallback_candidates: list[int] = []
        fallback_candidates.extend(
            sorted(
                defined
                - {atom_i, atom_b, atom_a}
                - set(bonded_candidates)
                - center_bonded_exclusions,
                key=lambda index: self._candidate_priority(index, atom_a, symbols, graph, coords),
                reverse=True,
            )
        )

        best_candidate = None
        best_margin = -1.0
        for candidate in self._unique(bonded_candidates):
            anchor_angle = self._angle_value(coords, atom_b, atom_a, candidate)
            margin = self._linear_margin(anchor_angle)
            if margin > best_margin:
                best_margin = margin
                best_candidate = candidate
        if best_candidate is not None:
            return best_candidate, self._dihedral_value(coords, atom_i, atom_b, atom_a, best_candidate), True

        best_candidate = None
        best_margin = -1.0
        for candidate in self._unique(fallback_candidates):
            anchor_angle = self._angle_value(coords, atom_b, atom_a, candidate)
            margin = self._linear_margin(anchor_angle)
            if margin > best_margin:
                best_margin = margin
                best_candidate = candidate
            if margin >= threshold:
                return candidate, self._dihedral_value(coords, atom_i, atom_b, atom_a, candidate), False

        if best_candidate is None:
            raise ValueError(f"Failed to choose a dihedral reference for atom {atom_i}.")
        return best_candidate, self._dihedral_value(coords, atom_i, atom_b, atom_a, best_candidate), False


def _normalize_format(path: Path, fmt: str | None) -> str:
    if fmt is not None:
        return fmt.lower()
    suffix = path.suffix.lower()
    if suffix == ".cif":
        return "cif"
    if suffix == ".pdb":
        return "pdb"
    if suffix == ".res":
        return "res"
    raise ValueError(f"Could not infer format from path: {path}")


def _strip_quotes(value: str) -> str:
    return value.strip().strip("'").strip('"')


def _strip_uncertainty(value: str) -> str:
    return re.sub(r"\([^)]+\)", "", value).strip()


def _parse_float(value: str) -> float:
    cleaned = _strip_uncertainty(_strip_quotes(value))
    if "/" in cleaned and re.fullmatch(r"[+-]?\d+/\d+", cleaned):
        return float(Fraction(cleaned))
    return float(cleaned)


def _parse_bool_tag(value: str | None, *, default: bool) -> bool:
    if value is None:
        return default

    normalized = _strip_quotes(value).strip().lower()
    if normalized in {"true", "t", "yes", "y", "1"}:
        return True
    if normalized in {"false", "f", "no", "n", "0"}:
        return False
    return default


def _parse_pdb_explict_unit_cell(lines: list[str]) -> bool:
    prefix = "REMARK   1 CSPTOOLBOX_EXPLICT_UNIT_CELL "
    for line in lines:
        if line.startswith(prefix):
            return _parse_bool_tag(line[len(prefix):], default=False)
    return False


def _lattice_matrix(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> list[tuple[float, float, float]]:
    alpha_r = math.radians(alpha)
    beta_r = math.radians(beta)
    gamma_r = math.radians(gamma)

    lattice = [
        (a, 0.0, 0.0),
        (b * math.cos(gamma_r), b * math.sin(gamma_r), 0.0),
        (0.0, 0.0, 0.0),
    ]
    c_x = c * math.cos(beta_r)
    c_y = c * (math.cos(alpha_r) - math.cos(beta_r) * math.cos(gamma_r)) / math.sin(gamma_r)
    c_z_sq = c**2 - c_x**2 - c_y**2
    lattice[2] = (c_x, c_y, math.sqrt(max(c_z_sq, 0.0)))
    return lattice


def _canonicalize_fractional(
    frac_coords: tuple[float, float, float],
    decimals: int = 8,
) -> tuple[float, float, float]:
    normalized: list[float] = []
    tolerance = 10 ** (-decimals)
    for value in frac_coords:
        wrapped = value % 1.0
        if math.isclose(wrapped, 1.0, abs_tol=tolerance):
            wrapped = 0.0
        normalized.append(round(wrapped, decimals))
    return tuple(normalized)  # type: ignore[return-value]


def _fractional_minimum_image_distance(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    lattice: list[tuple[float, float, float]],
) -> float:
    delta = np.array(right, dtype=float) - np.array(left, dtype=float)
    delta -= np.rint(delta)
    cartesian = np.dot(delta, np.array(lattice, dtype=float))
    return float(np.linalg.norm(cartesian))


def _matches_existing_expanded_cif_site(
    *,
    element: str,
    frac: tuple[float, float, float],
    expanded_sites: list[_ExpandedCifAtomSite],
    lattice: list[tuple[float, float, float]],
    site_merge_tolerance: float,
) -> bool:
    for site in expanded_sites:
        if site.element != element:
            continue
        if site.frac == frac:
            return True
        if (
            site_merge_tolerance > 0.0
            and _fractional_minimum_image_distance(site.frac, frac, lattice) <= site_merge_tolerance
        ):
            return True
    return False


def _apply_symmetry_operation(
    operation: str,
    frac_coords: tuple[float, float, float],
) -> tuple[float, float, float]:
    x, y, z = frac_coords
    components = [part.strip().replace(" ", "").lower() for part in operation.split(",")]
    values = [eval(component, {"__builtins__": {}}, {"x": x, "y": y, "z": z}) for component in components]
    return tuple(float(value) for value in values)  # type: ignore[return-value]


def _expand_shelx_atoms(
    *,
    atoms: list[AtomRecord],
    cell: np.ndarray,
    latt_value: int,
    symmetry_operations: Iterable[str],
) -> list[AtomRecord]:
    generated_atoms: list[AtomRecord] = []
    centrosymmetric = latt_value > 0
    centering_translations = _shelx_centering_translations(latt_value)
    normalized_operations = [str(operation).strip() for operation in symmetry_operations]

    for atom in atoms:
        frac = tuple(float(value) for value in Cell(cell).scaled_positions(np.array([atom.coordinates], dtype=float))[0])
        generated_fractionals = [frac]
        generated_fractionals.extend(
            _apply_symmetry_operation(operation, frac)
            for operation in normalized_operations
        )
        if centrosymmetric:
            generated_fractionals.extend(
                tuple(-value for value in coords)
                for coords in list(generated_fractionals)
            )

        unique_positions: list[tuple[float, float, float]] = []
        seen_positions: set[tuple[float, float, float]] = set()
        for coords in generated_fractionals:
            for translation in centering_translations:
                canonical = _canonicalize_fractional(
                    tuple(float(left + right) for left, right in zip(coords, translation))
                )
                if canonical in seen_positions:
                    continue
                seen_positions.add(canonical)
                unique_positions.append(canonical)

        use_suffix = len(unique_positions) > 1
        for index, canonical in enumerate(unique_positions, start=1):
            cart = np.dot(np.array(canonical, dtype=float), cell)
            generated_atoms.append(
                AtomRecord(
                    label=f"{atom.label}_{index}" if use_suffix else atom.label,
                    element=atom.element,
                    coordinates=tuple(float(value) for value in cart),
                )
            )

    return generated_atoms


def _frac_to_cart(
    frac_coords: tuple[float, float, float],
    lattice: list[tuple[float, float, float]],
) -> tuple[float, float, float]:
    fx, fy, fz = frac_coords
    a_vec, b_vec, c_vec = lattice
    return (
        fx * a_vec[0] + fy * b_vec[0] + fz * c_vec[0],
        fx * a_vec[1] + fy * b_vec[1] + fz * c_vec[1],
        fx * a_vec[2] + fy * b_vec[2] + fz * c_vec[2],
    )


def _parse_cif_blocks(text: str) -> tuple[dict[str, str], list[tuple[list[str], list[list[str]]]]]:
    data: dict[str, str] = {}
    loops: list[tuple[list[str], list[list[str]]]] = []
    lines = text.splitlines()
    index = 0

    while index < len(lines):
        line = lines[index].strip()
        if not line or line.startswith("#"):
            index += 1
            continue

        if line.startswith("loop_"):
            index += 1
            headers: list[str] = []
            while index < len(lines) and lines[index].lstrip().startswith("_"):
                headers.append(lines[index].strip())
                index += 1

            rows: list[list[str]] = []
            pending: list[str] = []
            while index < len(lines):
                row_line = lines[index].strip()
                if not row_line or row_line.startswith("#"):
                    index += 1
                    continue
                if row_line.startswith(("loop_", "_", "data_")):
                    break
                pending.extend(_split_cif_tokens(row_line))
                while len(pending) >= len(headers):
                    rows.append(pending[: len(headers)])
                    pending = pending[len(headers) :]
                index += 1

            loops.append((headers, rows))
            continue

        if line.startswith("_"):
            tokens = _split_cif_tokens(line)
            if len(tokens) > 1:
                data[tokens[0]] = " ".join(tokens[1:])
            index += 1
            continue

        index += 1

    return data, loops


def _split_cif_tokens(line: str) -> list[str]:
    try:
        return shlex.split(line, posix=True)
    except ValueError:
        return line.split()


def _find_loop(
    loops: list[tuple[list[str], list[list[str]]]],
    required_headers: set[str],
) -> tuple[list[str], list[list[str]]]:
    for headers, rows in loops:
        if required_headers.issubset(set(headers)):
            return headers, rows
    required = ", ".join(sorted(required_headers))
    raise ValueError(f"Could not find CIF atom loop with headers: {required}")


def _find_loop_with_any_header(
    loops: list[tuple[list[str], list[list[str]]]],
    candidate_headers: tuple[str, ...],
) -> tuple[list[str], list[list[str]]]:
    for headers, rows in loops:
        if any(header in headers for header in candidate_headers):
            return headers, rows
    required = ", ".join(candidate_headers)
    raise ValueError(f"Could not find CIF loop containing one of: {required}")


def _guess_element_from_label(label: str) -> str:
    match = re.match(r"([A-Za-z]+)", label)
    if not match:
        raise ValueError(f"Could not infer element from label: {label}")
    return _normalize_element_symbol(match.group(1))


def _normalize_element_symbol(symbol: str) -> str:
    raw = re.sub(r"[^A-Za-z]", "", str(symbol).strip())
    if not raw:
        raise ValueError(f"Could not normalize element symbol: {symbol!r}")
    if raw.upper() in {"D", "T"}:
        return "H"

    candidates = [raw]
    title_case = raw[0].upper() + raw[1:].lower()
    upper_case = raw.upper()
    lower_case = raw.lower()
    first_letter = raw[0].upper()
    for candidate in (title_case, upper_case, lower_case, first_letter):
        if candidate not in candidates:
            candidates.append(candidate)

    for candidate in candidates:
        if candidate in atomic_numbers:
            return candidate

    raise ValueError(f"Unsupported element symbol: {symbol!r}")


def _safe_data_name(name: str) -> str:
    return re.sub(r"\s+", "_", name.strip()) or "CrystalStructure"


def _normalize_space_group_symbol(symbol: str) -> str:
    normalized = " ".join(str(symbol).strip().split())
    match = re.match(r"^([A-Z])([^\s].*)$", normalized)
    if match:
        return f"{match.group(1)} {match.group(2)}"
    return normalized


def _normalize_space_group_lookup_key(symbol: str) -> str:
    return _normalize_space_group_symbol(symbol).replace(" ", "").replace("_", "").upper()


def _space_group_hall_number(space_group_symbol: str) -> int:
    target = _normalize_space_group_lookup_key(space_group_symbol)
    for hall_number in range(1, 531):
        space_group_type = spglib.get_spacegroup_type(hall_number)
        candidates: set[str] = set()
        for raw in (
            space_group_type.international_short,
            space_group_type.international_full,
            space_group_type.international,
        ):
            for piece in str(raw).split("="):
                candidates.add(_normalize_space_group_lookup_key(piece))
        if target in candidates:
            return hall_number
    raise ValueError(f"Unsupported or unknown space group for SHELX output: {space_group_symbol!r}")


def _shelx_symmetry_records(
    space_group_symbol: str,
    *,
    hall_number: int | None = None,
    rounding: bool = False,
) -> tuple[int, list[str]]:
    if hall_number is None:
        hall_number = _space_group_hall_number(space_group_symbol)
    symmetry = spglib.get_symmetry_from_database(hall_number)
    return _shelx_symmetry_records_from_symmetry(
        space_group_symbol=space_group_symbol,
        rotations=np.asarray(symmetry["rotations"], dtype=int),
        translations=np.asarray(symmetry["translations"], dtype=float),
        rounding=rounding,
    )


def _shelx_symmetry_records_from_symmetry(
    *,
    space_group_symbol: str,
    rotations: np.ndarray,
    translations: np.ndarray,
    rounding: bool = False,
) -> tuple[int, list[str]]:
    rotations = np.asarray(rotations, dtype=int)
    translations = np.asarray(translations, dtype=float)

    inversion_index = None
    for index, rotation in enumerate(rotations):
        if np.array_equal(rotation, -np.eye(3, dtype=int)):
            inversion_index = index
            break

    centrosymmetric = inversion_index is not None
    latt_value = _shelx_latt_value(
        space_group_symbol=space_group_symbol,
        centrosymmetric=centrosymmetric,
    )
    centering_translations = _shelx_centering_translations(latt_value)

    operations: list[tuple[np.ndarray, np.ndarray]] = []
    seen_keys: set[tuple[tuple[int, ...], tuple[int, int, int]]] = set()
    inversion_translation = (
        translations[inversion_index] if inversion_index is not None else None
    )

    for rotation, translation in zip(rotations, translations):
        if np.array_equal(rotation, np.eye(3, dtype=int)):
            continue
        if centrosymmetric and np.array_equal(rotation, -np.eye(3, dtype=int)):
            continue

        key = _centering_reduced_symmetry_operation_key(
            rotation,
            translation,
            centering_translations=centering_translations,
        )
        if key in seen_keys:
            continue
        seen_keys.add(key)

        if centrosymmetric and inversion_translation is not None:
            inversion_mate_key = _centering_reduced_symmetry_operation_key(
                -rotation,
                translation + inversion_translation,
                centering_translations=centering_translations,
            )
            seen_keys.add(inversion_mate_key)

        operations.append((rotation, translation))

    operation_strings = [
        _format_shelx_symmetry_operation(rotation, translation, rounding=rounding)
        for rotation, translation in operations
    ]
    return latt_value, operation_strings


def _shelx_centering_translations(latt_value: int) -> tuple[tuple[float, float, float], ...]:
    centering = abs(latt_value)
    centering_map: dict[int, tuple[tuple[float, float, float], ...]] = {
        1: ((0.0, 0.0, 0.0),),
        2: ((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)),
        3: ((0.0, 0.0, 0.0), (2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0), (1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)),
        4: ((0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)),
        5: ((0.0, 0.0, 0.0), (0.0, 0.5, 0.5)),
        6: ((0.0, 0.0, 0.0), (0.5, 0.0, 0.5)),
        7: ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0)),
    }
    if centering not in centering_map:
        raise ValueError(f"Unsupported SHELX LATT value: {latt_value}")
    return centering_map[centering]


def _shelx_latt_value(space_group_symbol: str, *, centrosymmetric: bool) -> int:
    leading = _normalize_space_group_symbol(space_group_symbol).split()[0].upper()
    centering_map = {
        "P": 1,
        "I": 2,
        "R": 3,
        "F": 4,
        "A": 5,
        "B": 6,
        "C": 7,
    }
    if leading not in centering_map:
        raise ValueError(f"Unsupported lattice centering for SHELX output: {space_group_symbol!r}")
    value = centering_map[leading]
    return value if centrosymmetric else -value


def _symmetry_operation_key(
    rotation: np.ndarray,
    translation: np.ndarray | Iterable[float],
) -> tuple[tuple[int, ...], tuple[int, int, int]]:
    normalized_translation = tuple(
        int(round(value * 24)) % 24
        for value in np.mod(np.asarray(list(translation), dtype=float), 1.0)
    )
    return tuple(rotation.reshape(-1).tolist()), normalized_translation


def _centering_reduced_symmetry_operation_key(
    rotation: np.ndarray,
    translation: np.ndarray | Iterable[float],
    *,
    centering_translations: Iterable[tuple[float, float, float]],
) -> tuple[tuple[int, ...], tuple[int, int, int]]:
    translation_array = np.asarray(list(translation), dtype=float)
    candidate_keys = [
        _symmetry_operation_key(rotation, translation_array - np.asarray(centering, dtype=float))
        for centering in centering_translations
    ]
    return min(candidate_keys)


def _deduplicate_shelx_symmetry_operations(
    symmetry_operations: Iterable[str],
    *,
    latt_value: int,
    rounding: bool = False,
) -> list[str]:
    centering_translations = _shelx_centering_translations(latt_value)
    deduplicated: list[str] = []
    seen_keys: set[tuple[tuple[int, ...], tuple[int, int, int]]] = set()

    for operation in symmetry_operations:
        normalized = _normalize_shelx_symmetry_operation(operation, rounding=rounding)
        rotation_rows: list[tuple[int, int, int]] = []
        translation_terms: list[float] = []
        for component in normalized.split(","):
            rotation_row, translation = _parse_general_symmetry_component(component)
            rotation_rows.append(rotation_row)
            translation_terms.append(translation)
        key = _centering_reduced_symmetry_operation_key(
            np.asarray(rotation_rows, dtype=int),
            np.asarray(translation_terms, dtype=float),
            centering_translations=centering_translations,
        )
        if key in seen_keys:
            continue
        seen_keys.add(key)
        deduplicated.append(normalized if rounding else str(operation).strip())

    return deduplicated


_CSORM_SUPPORTED_TRANSLATIONS = frozenset(
    {
        Fraction(0, 1),
        Fraction(1, 4),
        Fraction(1, 3),
        Fraction(1, 2),
        Fraction(2, 3),
        Fraction(3, 4),
    }
)
_CSORM_DECIMAL_TRANSLATION_MAP = {
    "2": Fraction(1, 4),
    "3": Fraction(1, 3),
    "5": Fraction(1, 2),
    "6": Fraction(2, 3),
    "7": Fraction(3, 4),
}
_RES_TRANSLATION_ROUNDING_TOLERANCE = 0.01
_RES_TRANSLATION_TARGETS = (
    0.0,
    1.0 / 4.0,
    1.0 / 3.0,
    1.0 / 2.0,
    2.0 / 3.0,
    3.0 / 4.0,
    1.0,
)


def _normalize_csorm_symmetry_operation(operation: str) -> str:
    components = [component.strip() for component in operation.split(",")]
    if len(components) != 3:
        raise ValueError(f"Expected 3 symmetry components, got {len(components)} in {operation!r}")

    rotations: list[tuple[int, int, int]] = []
    translations: list[float] = []
    for component in components:
        rotation_row, translation = _parse_csorm_symmetry_component(component)
        rotations.append(rotation_row)
        translations.append(float(translation))

    return _format_shelx_symmetry_operation(
        np.asarray(rotations, dtype=int),
        np.asarray(translations, dtype=float),
        rounding=True,
    )


def _normalize_shelx_symmetry_operation(operation: str, *, rounding: bool = False) -> str:
    components = [component.strip() for component in operation.split(",")]
    if len(components) != 3:
        raise ValueError(f"Expected 3 symmetry components, got {len(components)} in {operation!r}")

    rotations: list[tuple[int, int, int]] = []
    translations: list[float] = []
    for component in components:
        rotation_row, translation = _parse_general_symmetry_component(component)
        rotations.append(rotation_row)
        translations.append(float(translation))

    return _format_shelx_symmetry_operation(
        np.asarray(rotations, dtype=int),
        np.asarray(translations, dtype=float),
        rounding=rounding,
    )


def _parse_general_symmetry_component(component: str) -> tuple[tuple[int, int, int], float]:
    expr = component.replace(" ", "").upper()
    if not expr:
        raise ValueError("Empty symmetry component.")

    coeffs = {"X": 0, "Y": 0, "Z": 0}
    translation = 0.0
    expr = expr.replace("-", "+-")
    if expr.startswith("+-"):
        expr = "-" + expr[2:]
    terms = [term for term in expr.split("+") if term]
    if not terms:
        raise ValueError(f"Could not parse symmetry component {component!r}")

    for term in terms:
        if term in ("X", "Y", "Z"):
            coeffs[term] += 1
            continue
        if term in ("-X", "-Y", "-Z"):
            coeffs[term[1:]] -= 1
            continue
        translation += _parse_general_translation_token(term)

    rotation = (coeffs["X"], coeffs["Y"], coeffs["Z"])
    if any(abs(value) > 1 for value in rotation):
        raise ValueError(f"Unsupported symmetry rotation term in {component!r}")

    return rotation, translation


def _parse_general_translation_token(token: str) -> float:
    cleaned = token.strip().upper()
    sign = 1.0
    if cleaned.startswith("-"):
        sign = -1.0
        cleaned = cleaned[1:]
    elif cleaned.startswith("+"):
        cleaned = cleaned[1:]

    if not cleaned:
        raise ValueError("Empty translation token.")

    if "/" in cleaned:
        return sign * float(Fraction(cleaned))
    return sign * float(cleaned)


def _parse_csorm_symmetry_component(component: str) -> tuple[tuple[int, int, int], Fraction]:
    expr = component.replace(" ", "").upper()
    if not expr:
        raise ValueError("Empty symmetry component.")

    coeffs = {"X": 0, "Y": 0, "Z": 0}
    translation = Fraction(0, 1)
    expr = expr.replace("-", "+-")
    if expr.startswith("+-"):
        expr = "-" + expr[2:]
    terms = [term for term in expr.split("+") if term]
    if not terms:
        raise ValueError(f"Could not parse symmetry component {component!r}")

    for term in terms:
        if term in ("X", "Y", "Z"):
            coeffs[term] += 1
            continue
        if term in ("-X", "-Y", "-Z"):
            coeffs[term[1:]] -= 1
            continue
        translation += _parse_csorm_translation_token(term)

    normalized_translation = translation % 1
    if normalized_translation not in _CSORM_SUPPORTED_TRANSLATIONS:
        raise ValueError(
            f"Unsupported CSORM translation {component!r}: final offset {normalized_translation}"
        )

    rotation = (coeffs["X"], coeffs["Y"], coeffs["Z"])
    if any(abs(value) > 1 for value in rotation):
        raise ValueError(f"Unsupported CSORM rotation term in {component!r}")

    return rotation, normalized_translation


def _parse_csorm_translation_token(token: str) -> Fraction:
    cleaned = token.strip().upper()
    sign = 1
    if cleaned.startswith("-"):
        sign = -1
        cleaned = cleaned[1:]
    elif cleaned.startswith("+"):
        cleaned = cleaned[1:]

    if not cleaned:
        raise ValueError("Empty translation token.")

    if "/" in cleaned:
        fraction = Fraction(cleaned)
        normalized = fraction % 1
        if normalized not in _CSORM_SUPPORTED_TRANSLATIONS:
            raise ValueError(f"Unsupported CSORM rational translation {token!r}")
        return sign * fraction

    if "." in cleaned:
        dot_index = cleaned.find(".")
        if dot_index == len(cleaned) - 1:
            raise ValueError(f"Invalid CSORM decimal translation {token!r}")
        decimal_digit = cleaned[dot_index + 1]
        if decimal_digit not in _CSORM_DECIMAL_TRANSLATION_MAP:
            raise ValueError(f"Unsupported CSORM decimal translation {token!r}")
        integer_part = cleaned[:dot_index]
        whole = int(integer_part) if integer_part else 0
        return sign * (Fraction(whole, 1) + _CSORM_DECIMAL_TRANSLATION_MAP[decimal_digit])

    return sign * Fraction(cleaned)


def _format_shelx_symmetry_operation(
    rotation: np.ndarray,
    translation: np.ndarray,
    *,
    rounding: bool = False,
) -> str:
    axes = ("X", "Y", "Z")
    components: list[str] = []
    for row_index in range(3):
        terms: list[str] = []
        offset = _format_fractional_offset(float(translation[row_index]), rounding=rounding)
        if offset:
            terms.append(offset)

        for coeff, axis in zip(rotation[row_index], axes):
            if coeff == 0:
                continue
            if coeff == 1:
                if terms:
                    terms.append(f"+{axis}")
                else:
                    terms.append(axis)
            elif coeff == -1:
                if terms:
                    terms.append(f"-{axis}")
                else:
                    terms.append(f"-{axis}")
            else:
                raise ValueError("Only unit rotation coefficients are supported for SHELX symmetry output.")

        component = "".join(terms) if terms else "0"
        components.append(component)
    return ",".join(components)


def _format_fractional_offset(value: float, *, rounding: bool = False) -> str:
    normalized = float(np.mod(value, 1.0))
    if rounding:
        normalized = _round_res_fractional_translation(normalized)
    if abs(normalized) < 1e-8 or abs(normalized - 1.0) < 1e-8:
        return ""
    fraction = Fraction(normalized).limit_denominator(24)
    if abs(float(fraction) - normalized) > 1e-6:
        return f"{normalized:.6f}".rstrip("0").rstrip(".")
    if fraction.denominator == 1:
        return str(fraction.numerator)
    return f"{fraction.numerator}/{fraction.denominator}"


def _round_res_fractional_translation(value: float) -> float:
    normalized = float(np.mod(value, 1.0))
    for target in _RES_TRANSLATION_TARGETS:
        if abs(normalized - target) <= _RES_TRANSLATION_ROUNDING_TOLERANCE:
            return 0.0 if math.isclose(target, 1.0, abs_tol=1e-12) else target
    return normalized


def _element_order(elements: Iterable[str]) -> list[str]:
    ordered: list[str] = []
    for element in elements:
        if element not in ordered:
            ordered.append(element)
    return ordered


def _compare_ase_and_manual_expansions(
    manual_atoms: list[_ExpandedCifAtomSite],
    ase_atoms: Atoms,
) -> tuple[bool, str | None]:
    manual_counter = Counter(
        (atom.element, _canonicalize_fractional(atom.frac))
        for atom in manual_atoms
    )
    ase_counter = Counter(
        (
            _normalize_element_symbol(symbol),
            _canonicalize_fractional(
                tuple(float(value) for value in frac)
            ),
        )
        for symbol, frac in zip(
            ase_atoms.get_chemical_symbols(),
            ase_atoms.get_scaled_positions(wrap=True),
        )
    )

    if manual_counter == ase_counter:
        return True, None

    manual_only = manual_counter - ase_counter
    ase_only = ase_counter - manual_counter
    messages: list[str] = []
    if manual_only:
        messages.append(
            "manual-only sites: "
            + ", ".join(
                f"{element}@{frac} x{count}"
                for (element, frac), count in sorted(manual_only.items())
            )
        )
    if ase_only:
        messages.append(
            "ase-only sites: "
            + ", ".join(
                f"{element}@{frac} x{count}"
                for (element, frac), count in sorted(ase_only.items())
            )
        )
    return False, "; ".join(messages) or "expanded atom sets differ"
