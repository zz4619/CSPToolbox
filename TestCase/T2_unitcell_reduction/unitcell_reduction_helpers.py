"""Helpers for the T2 unit-cell reduction regression tests."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import CrystalStructure
from Source.pdd_descriptor import pdd_distance


CASE_ROOT = Path(__file__).resolve().parent
ASYMMETRIC_DIR = CASE_ROOT / "Experimental"
FULL_UNIT_CELL_DIR = CASE_ROOT / "Experimental_FullUnitCell"

DEFAULT_SMOKE_REFCODES = (
    "UMIQEO",
    "ABALAS",
    "DMSULO04",
    "CYCYPR",
    "TRIZIN01",
)

DEFAULT_SYMPREC = 0.05
DEFAULT_PDD_K = 100
DEFAULT_PDD_TOLERANCE = 1e-6
DEFAULT_SITE_TOLERANCE_ANGSTROM = 1e-6
DEFAULT_LATTICE_TOLERANCE = 1e-6


@dataclass(frozen=True)
class ReductionComparison:
    """Summary of one full-cell -> asymmetric-unit -> full-cell comparison."""

    refcode: str
    full_atom_count: int
    reference_asymmetric_atom_count: int
    reduced_atom_count: int
    reexpanded_atom_count: int
    detected_space_group: str
    detected_hall_number: int | None
    detected_symmetry_operation_count: int
    pdd_to_full: float
    pdd_to_t1_expansion: float
    max_site_displacement_to_full: float
    max_site_displacement_to_t1_expansion: float


def paired_refcodes(
    *,
    asymmetric_dir: Path = ASYMMETRIC_DIR,
    full_unit_cell_dir: Path = FULL_UNIT_CELL_DIR,
) -> list[str]:
    """Return sorted refcodes that exist in both T2 CIF folders."""

    asymmetric = {path.stem for path in asymmetric_dir.glob("*.cif")}
    full_cell = {path.stem for path in full_unit_cell_dir.glob("*.cif")}
    missing_asymmetric = sorted(full_cell - asymmetric)
    missing_full_cell = sorted(asymmetric - full_cell)
    if missing_asymmetric or missing_full_cell:
        raise AssertionError(
            "T2 CIF folders are not paired: "
            f"missing asymmetric={missing_asymmetric[:10]}, "
            f"missing full_cell={missing_full_cell[:10]}"
        )
    return sorted(asymmetric & full_cell)


def paths_for_refcode(
    refcode: str,
    *,
    asymmetric_dir: Path = ASYMMETRIC_DIR,
    full_unit_cell_dir: Path = FULL_UNIT_CELL_DIR,
) -> tuple[Path, Path]:
    """Return `(asymmetric_path, full_unit_cell_path)` for one refcode."""

    asymmetric_path = asymmetric_dir / f"{refcode}.cif"
    full_unit_cell_path = full_unit_cell_dir / f"{refcode}.cif"
    if not asymmetric_path.is_file():
        raise FileNotFoundError(f"Missing T2 asymmetric CIF: {asymmetric_path}")
    if not full_unit_cell_path.is_file():
        raise FileNotFoundError(f"Missing T2 full-unit-cell CIF: {full_unit_cell_path}")
    return asymmetric_path, full_unit_cell_path


def compare_reduction_round_trip(
    refcode: str,
    *,
    symprec: float = DEFAULT_SYMPREC,
    pdd_k: int = DEFAULT_PDD_K,
) -> ReductionComparison:
    """Reduce a full-cell CIF, re-expand it, and compare against references."""

    asymmetric_path, full_unit_cell_path = paths_for_refcode(refcode)

    full_structure = CrystalStructure.from_file(full_unit_cell_path)
    reference_asymmetric = CrystalStructure.from_file(asymmetric_path)
    t1_expanded = CrystalStructure.expand_cif_to_unit_cell(asymmetric_path)

    reduced = full_structure.reduce_to_asymmetric_unit(symprec=symprec)
    reexpanded = reduced.expand_to_explicit_unit_cell()

    return ReductionComparison(
        refcode=refcode,
        full_atom_count=len(full_structure.atoms),
        reference_asymmetric_atom_count=len(reference_asymmetric.atoms),
        reduced_atom_count=len(reduced.atoms),
        reexpanded_atom_count=len(reexpanded.atoms),
        detected_space_group=reduced.space_group,
        detected_hall_number=reduced.hall_number,
        detected_symmetry_operation_count=len(reduced.symmetry_operations),
        pdd_to_full=pdd_distance(
            reexpanded,
            full_structure,
            k=pdd_k,
            metric="chebyshev",
        ),
        pdd_to_t1_expansion=pdd_distance(
            reexpanded,
            t1_expanded,
            k=pdd_k,
            metric="chebyshev",
        ),
        max_site_displacement_to_full=max_same_element_site_displacement(
            reexpanded,
            full_structure,
        ),
        max_site_displacement_to_t1_expansion=max_same_element_site_displacement(
            reexpanded,
            t1_expanded,
        ),
    )


def assert_reduction_round_trip(
    refcode: str,
    *,
    symprec: float = DEFAULT_SYMPREC,
    pdd_tolerance: float = DEFAULT_PDD_TOLERANCE,
    site_tolerance: float = DEFAULT_SITE_TOLERANCE_ANGSTROM,
) -> ReductionComparison:
    """Assert the T2 reduction invariant for one refcode."""

    comparison = compare_reduction_round_trip(refcode, symprec=symprec)
    errors: list[str] = []

    if comparison.reduced_atom_count < 1:
        errors.append("reduced structure has no atoms")
    if comparison.reduced_atom_count > comparison.full_atom_count:
        errors.append(
            "reduced structure has more atoms than the full-cell input "
            f"({comparison.reduced_atom_count} > {comparison.full_atom_count})"
        )
    if comparison.reexpanded_atom_count != comparison.full_atom_count:
        errors.append(
            "re-expansion did not restore the full-cell atom count "
            f"({comparison.reexpanded_atom_count} != {comparison.full_atom_count})"
        )
    if comparison.detected_hall_number is None:
        errors.append("reduced structure has no detected Hall number")
    if comparison.detected_space_group == "P 1" and comparison.reduced_atom_count != comparison.full_atom_count:
        errors.append("P 1 detection reduced the atom count unexpectedly")
    if comparison.pdd_to_full > pdd_tolerance:
        errors.append(
            f"PDD to full-cell reference is {comparison.pdd_to_full:.3e}, "
            f"above {pdd_tolerance:.3e}"
        )
    if comparison.pdd_to_t1_expansion > pdd_tolerance:
        errors.append(
            f"PDD to T1 expansion is {comparison.pdd_to_t1_expansion:.3e}, "
            f"above {pdd_tolerance:.3e}"
        )
    if comparison.max_site_displacement_to_full > site_tolerance:
        errors.append(
            "site displacement to full-cell reference is "
            f"{comparison.max_site_displacement_to_full:.3e} A, above {site_tolerance:.3e} A"
        )
    if comparison.max_site_displacement_to_t1_expansion > site_tolerance:
        errors.append(
            "site displacement to T1 expansion is "
            f"{comparison.max_site_displacement_to_t1_expansion:.3e} A, above {site_tolerance:.3e} A"
        )

    if errors:
        detail = "\n".join(f"- {message}" for message in errors)
        raise AssertionError(f"T2 reduction round-trip failed for {refcode}:\n{detail}")

    return comparison


def max_same_element_site_displacement(
    left: CrystalStructure,
    right: CrystalStructure,
    *,
    lattice_tolerance: float = DEFAULT_LATTICE_TOLERANCE,
) -> float:
    """Return max same-element minimum-image displacement between two structures."""

    if not np.allclose(left.cell.array, right.cell.array, atol=lattice_tolerance, rtol=0.0):
        raise AssertionError(
            f"Lattice mismatch between {left.name!r} and {right.name!r}."
        )

    left_counts = element_counts(left)
    right_counts = element_counts(right)
    if left_counts != right_counts:
        raise AssertionError(
            f"Element-count mismatch between {left.name!r} and {right.name!r}: "
            f"{left_counts} != {right_counts}"
        )

    lattice = np.asarray(left.cell.array, dtype=float)
    max_displacement = 0.0
    for element in sorted(left_counts):
        left_frac = _fractional_positions_for_element(left, element)
        right_frac = _fractional_positions_for_element(right, element)
        if len(left_frac) == 0:
            continue

        distances = _minimum_image_distance_matrix(left_frac, right_frac, lattice)
        row_indices, column_indices = linear_sum_assignment(distances)
        matched = distances[row_indices, column_indices]
        if len(matched):
            max_displacement = max(max_displacement, float(np.max(matched)))

    return max_displacement


def element_counts(structure: CrystalStructure) -> Counter[str]:
    """Return normalized element counts for a structure."""

    return Counter(atom.element for atom in structure.atoms)


def _fractional_positions_for_element(
    structure: CrystalStructure,
    element: str,
) -> np.ndarray:
    positions = [
        structure.fractional_coordinates(atom.coordinates)
        for atom in structure.atoms
        if atom.element == element
    ]
    return np.asarray(positions, dtype=float)


def _minimum_image_distance_matrix(
    left_frac: np.ndarray,
    right_frac: np.ndarray,
    lattice: np.ndarray,
) -> np.ndarray:
    deltas = right_frac[None, :, :] - left_frac[:, None, :]
    deltas -= np.rint(deltas)
    cartesian = deltas @ lattice
    return np.linalg.norm(cartesian, axis=2)

