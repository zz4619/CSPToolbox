"""Compare two explicit unit cells site by site."""

from __future__ import annotations

from collections import Counter

import numpy as np
from scipy.optimize import linear_sum_assignment

from Source.crystal_structure import CrystalStructure


def element_counts(structure: CrystalStructure) -> Counter[str]:
    """Return element counts for a structure."""

    return Counter(atom.element for atom in structure.atoms)


def max_same_element_site_displacement(
    left: CrystalStructure,
    right: CrystalStructure,
    *,
    lattice_tolerance: float = 1e-6,
) -> float:
    """Return the largest same-element, minimum-image displacement in Å.

    Sites are paired by an optimal assignment per element, so atom order and
    labels do not matter. Both structures must share the same lattice and
    element counts; otherwise an AssertionError explains the mismatch.
    """

    if not np.allclose(left.cell.array, right.cell.array, atol=lattice_tolerance, rtol=0.0):
        raise AssertionError(f"Lattice mismatch between {left.name!r} and {right.name!r}.")
    left_counts = element_counts(left)
    right_counts = element_counts(right)
    if left_counts != right_counts:
        raise AssertionError(
            f"Element-count mismatch between {left.name!r} and {right.name!r}: "
            f"{dict(left_counts)} != {dict(right_counts)}"
        )

    lattice = np.asarray(left.cell.array, dtype=float)
    largest = 0.0
    for element in sorted(left_counts):
        left_frac = _fractional_positions(left, element)
        right_frac = _fractional_positions(right, element)
        deltas = right_frac[None, :, :] - left_frac[:, None, :]
        deltas -= np.rint(deltas)
        distances = np.linalg.norm(deltas @ lattice, axis=2)
        rows, columns = linear_sum_assignment(distances)
        largest = max(largest, float(np.max(distances[rows, columns])))
    return largest


def _fractional_positions(structure: CrystalStructure, element: str) -> np.ndarray:
    return np.asarray(
        [
            structure.fractional_coordinates(atom.coordinates)
            for atom in structure.atoms
            if atom.element == element
        ],
        dtype=float,
    )
