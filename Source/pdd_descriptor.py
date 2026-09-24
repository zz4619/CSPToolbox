"""PDD descriptor and distance utilities for explicit unit-cell structures."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache, reduce
from itertools import product
from math import gcd
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import linprog, linear_sum_assignment
from scipy.sparse import coo_matrix
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist, pdist, squareform

if TYPE_CHECKING:
    from .crystal_structure import CrystalStructure

__all__ = ['PDDDescriptor', 'calculate_pdd', 'calculate_pdd_from_arrays',
           'pdd_distance', 'pdd_distance_breakdown', 'pdd_descriptor_distance', 'pdd_to_amd']


@dataclass(frozen=True)
class PDDDescriptor:
    source_name: str
    k: int
    weights: np.ndarray
    distances: np.ndarray
    center_elements: tuple[str, ...]
    typed: bool
    collapse: bool
    collapse_tol: float

    @property
    def matrix(self) -> np.ndarray:
        return np.hstack((self.weights[:, None], self.distances))

    def rows_for_element(self, element: str) -> tuple[np.ndarray, np.ndarray]:
        mask = np.array([center == element for center in self.center_elements], dtype=bool)
        return self.weights[mask], self.distances[mask]


def calculate_pdd(
    structure: CrystalStructure,
    *,
    k: int = 100,
    typed: bool = True,
    lexsort: bool = True,
    collapse: bool = True,
    collapse_tol: float = 1e-4,
    workers: int = 1,
) -> PDDDescriptor:
    """Calculate the PDD descriptor for an explicit unit-cell structure."""

    _require_explict_unit_cell(
        structure,
        context="PDD calculation requires an explict unit cell structure.",
    )
    if k < 1:
        raise ValueError("k must be at least 1.")

    motif, cell, center_numbers = _explicit_structure_to_pdd_input(structure)
    return calculate_pdd_from_arrays(
        motif, cell, [atom.element for atom in structure.atoms],
        source_name=structure.name, k=k, typed=typed, lexsort=lexsort,
        collapse=collapse, collapse_tol=collapse_tol, workers=workers,
    )


def calculate_pdd_from_arrays(
    coordinates, cell, elements, *, source_name="CrystalStructure", k=100,
    typed=True, lexsort=True, collapse=True, collapse_tol=1e-4, workers=1,
) -> PDDDescriptor:
    """Describe a FULL periodic cell; Cartesian coordinates and cell rows are in Å.

    This entry point deliberately does not infer or expand space-group symmetry.
    """
    from ase.data import atomic_numbers, chemical_symbols
    motif = np.asarray(coordinates, dtype=float)
    cell = np.asarray(cell, dtype=float)
    elements = tuple(elements)
    if motif.ndim != 2 or motif.shape[1:] != (3,) or not len(motif):
        raise ValueError("Coordinates must be a nonempty (N, 3) array.")
    if cell.shape != (3, 3) or not np.all(np.isfinite(cell)) or abs(np.linalg.det(cell)) < 1e-12:
        raise ValueError("Cell must be finite, nonsingular and 3 by 3.")
    if not np.all(np.isfinite(motif)) or len(elements) != len(motif):
        raise ValueError("Coordinates must be finite and match the element count.")
    if isinstance(k, bool) or int(k) != k or k < 1:
        raise ValueError("k must be a positive integer.")
    if not np.isfinite(collapse_tol) or collapse_tol < 0:
        raise ValueError("collapse_tol must be finite and nonnegative.")
    center_numbers = np.array([atomic_numbers[e] for e in elements], dtype=int)
    distances = _nearest_neighbours(motif, cell, int(k), workers=workers)
    weights = np.full((distances.shape[0],), 1.0 / distances.shape[0], dtype=float)
    center_elements = [chemical_symbols[number] for number in center_numbers]

    if collapse:
        if typed:
            weights, distances, center_elements = _collapse_typed_pdd_rows(
                weights,
                distances,
                center_numbers,
                collapse_tol=collapse_tol,
            )
        else:
            weights, distances, center_elements = _collapse_untyped_pdd_rows(
                weights,
                distances,
                center_elements,
                collapse_tol=collapse_tol,
            )

    if lexsort:
        if typed:
            ordering = sorted(
                range(len(weights)),
                key=lambda index: (center_elements[index], *distances[index].tolist()),
            )
        else:
            ordering = sorted(
                range(len(weights)),
                key=lambda index: tuple(distances[index].tolist()),
            )
        weights = weights[ordering]
        distances = distances[ordering]
        center_elements = [center_elements[index] for index in ordering]

    return PDDDescriptor(
        source_name=source_name,
        k=k,
        weights=weights,
        distances=distances,
        center_elements=tuple(center_elements),
        typed=typed,
        collapse=collapse,
        collapse_tol=collapse_tol,
    )


def _validate_descriptor(descriptor: PDDDescriptor) -> None:
    w, d = descriptor.weights, descriptor.distances
    if (w.ndim != 1 or not len(w) or d.shape != (len(w), descriptor.k)
            or len(descriptor.center_elements) != len(w) or descriptor.k < 1):
        raise ValueError("Invalid PDD array dimensions.")
    if (not np.all(np.isfinite(w)) or not np.all(np.isfinite(d))
            or np.any(w <= 0) or np.any(d < 0)
            or not np.isclose(w.sum(), 1.0, atol=1e-12, rtol=0)):
        raise ValueError("PDD needs finite nonnegative distances and positive normalized weights.")
    if np.any(np.diff(d, axis=1) < -1e-12):
        raise ValueError("PDD neighbour distances must be sorted within each row.")


def pdd_to_amd(descriptor: PDDDescriptor) -> np.ndarray:
    """Return the weighted column means, including all central elements."""
    _validate_descriptor(descriptor)
    return descriptor.weights @ descriptor.distances


def pdd_descriptor_distance(a: PDDDescriptor, b: PDDDescriptor, *, metric="chebyshev") -> float:
    """Compare cached descriptors, allowing equivalent cells of different sizes.

    Typed comparisons require equal ELEMENT FRACTIONS, not equal cell atom counts.
    Default EMD and AMD distances are in Å. Equal descriptors need not uniquely
    identify a crystal at finite k; this is a geometric descriptor metric.
    """
    _validate_descriptor(a)
    _validate_descriptor(b)
    if a.k != b.k or a.typed != b.typed:
        raise ValueError("PDDs must use the same k and typing mode.")
    composition = None
    if a.typed:
        elements = sorted(set(a.center_elements) | set(b.center_elements))
        for e in elements:
            if not np.isclose(a.rows_for_element(e)[0].sum(), b.rows_for_element(e)[0].sum(), atol=1e-12, rtol=0):
                raise ValueError("Typed PDDs require matching element fractions.")
        composition = tuple((e, 1) for e in elements)
    return _pdd_distance_breakdown(a, b, metric=metric, expected_composition=composition)[0]


def pdd_distance(
    structure_a: CrystalStructure,
    structure_b: CrystalStructure,
    *,
    k: int = 100,
    typed: bool = True,
    metric: str = "chebyshev",
    collapse: bool = True,
    collapse_tol: float = 1e-4,
    **metric_kwargs,
) -> float:
    """Compare two explicit unit-cell structures of the same chemistry."""

    descriptor_a = calculate_pdd(
        structure_a,
        k=k,
        typed=typed,
        collapse=collapse,
        collapse_tol=collapse_tol,
    )
    descriptor_b = calculate_pdd(
        structure_b,
        k=k,
        typed=typed,
        collapse=collapse,
        collapse_tol=collapse_tol,
    )
    total, _ = _pdd_distance_breakdown(
        descriptor_a,
        descriptor_b,
        metric=metric,
        expected_composition=_composition_signature(structure_a, structure_b) if typed else None,
        **metric_kwargs,
    )
    return total


def pdd_distance_breakdown(
    structure_a: CrystalStructure,
    structure_b: CrystalStructure,
    *,
    k: int = 100,
    typed: bool = True,
    metric: str = "chebyshev",
    collapse: bool = True,
    collapse_tol: float = 1e-4,
    **metric_kwargs,
) -> tuple[float, dict[str, float]]:
    """Return the total PDD distance and the per-element contributions."""

    descriptor_a = calculate_pdd(
        structure_a,
        k=k,
        typed=typed,
        collapse=collapse,
        collapse_tol=collapse_tol,
    )
    descriptor_b = calculate_pdd(
        structure_b,
        k=k,
        typed=typed,
        collapse=collapse,
        collapse_tol=collapse_tol,
    )
    return _pdd_distance_breakdown(
        descriptor_a,
        descriptor_b,
        metric=metric,
        expected_composition=_composition_signature(structure_a, structure_b) if typed else None,
        **metric_kwargs,
    )


def _pdd_distance_breakdown(
    descriptor_a: PDDDescriptor,
    descriptor_b: PDDDescriptor,
    *,
    metric: str,
    expected_composition: tuple[tuple[str, int], ...] | None,
    **metric_kwargs,
) -> tuple[float, dict[str, float]]:
    if descriptor_a.typed != descriptor_b.typed:
        raise ValueError("PDD descriptors must use the same typed/geometry-only mode.")

    if not descriptor_a.typed:
        cost_matrix = cdist(descriptor_a.distances, descriptor_b.distances, metric=metric, **metric_kwargs)
        total = _earth_movers_distance(descriptor_a.weights, descriptor_b.weights, cost_matrix)
        return total, {"ALL": total}

    if expected_composition is None:
        raise ValueError("Typed PDD comparison requires an expected composition.")

    expected_elements = {element for element, _ in expected_composition}
    if set(descriptor_a.center_elements) != expected_elements or set(descriptor_b.center_elements) != expected_elements:
        raise ValueError("PDD descriptors do not cover the expected element types.")

    by_element: dict[str, float] = {}
    total = 0.0
    for element, _ in expected_composition:
        weights_a, distances_a = descriptor_a.rows_for_element(element)
        weights_b, distances_b = descriptor_b.rows_for_element(element)
        cost_matrix = cdist(distances_a, distances_b, metric=metric, **metric_kwargs)
        element_distance = _earth_movers_distance(weights_a, weights_b, cost_matrix)
        by_element[element] = element_distance
        total += element_distance

    return total, by_element


def _explicit_structure_to_pdd_input(
    structure: CrystalStructure,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from ase.data import atomic_numbers
    coordinates = np.array([atom.coordinates for atom in structure.atoms], dtype=float)
    fractional = np.mod(structure.cell.scaled_positions(coordinates), 1.0)
    cell = np.asarray(structure.cell.array, dtype=float)
    motif = fractional @ cell
    center_numbers = np.array([atomic_numbers[atom.element] for atom in structure.atoms], dtype=int)
    return motif, cell, center_numbers


def _collapse_typed_pdd_rows(
    weights: np.ndarray,
    distances: np.ndarray,
    center_numbers: np.ndarray,
    *,
    collapse_tol: float,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    from ase.data import chemical_symbols
    collapsed_weights: list[float] = []
    collapsed_distances: list[np.ndarray] = []
    collapsed_centers: list[str] = []

    for atomic_number in sorted(set(center_numbers.tolist())):
        element_indexes = np.flatnonzero(center_numbers == atomic_number)
        element_weights = weights[element_indexes]
        element_distances = distances[element_indexes]
        groups = [[index] for index in range(len(element_indexes))]

        if len(element_indexes) > 1:
            overlapping = pdist(element_distances, metric="chebyshev") < collapse_tol
            if overlapping.any():
                groups = _collapse_into_groups(overlapping)

        for group in groups:
            group_weights = element_weights[group]
            group_distances = element_distances[group]
            collapsed_weights.append(float(np.sum(group_weights)))
            collapsed_distances.append(
                np.average(group_distances, axis=0, weights=group_weights)
            )
            collapsed_centers.append(chemical_symbols[atomic_number])

    return (
        np.array(collapsed_weights, dtype=float),
        np.array(collapsed_distances, dtype=float),
        collapsed_centers,
    )


def _collapse_untyped_pdd_rows(
    weights: np.ndarray,
    distances: np.ndarray,
    center_elements: list[str],
    *,
    collapse_tol: float,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    groups = [[index] for index in range(len(weights))]

    if len(weights) > 1:
        overlapping = pdist(distances, metric="chebyshev") < collapse_tol
        if overlapping.any():
            groups = _collapse_into_groups(overlapping)

    collapsed_weights: list[float] = []
    collapsed_distances: list[np.ndarray] = []
    collapsed_centers: list[str] = []

    for group in groups:
        group_weights = weights[group]
        group_distances = distances[group]
        group_elements = {center_elements[index] for index in group}
        collapsed_weights.append(float(np.sum(group_weights)))
        collapsed_distances.append(
            np.average(group_distances, axis=0, weights=group_weights)
        )
        if len(group_elements) == 1:
            collapsed_centers.append(next(iter(group_elements)))
        else:
            collapsed_centers.append("MIXED")

    return (
        np.array(collapsed_weights, dtype=float),
        np.array(collapsed_distances, dtype=float),
        collapsed_centers,
    )


def _composition_signature(
    structure_a: CrystalStructure,
    structure_b: CrystalStructure,
) -> tuple[tuple[str, int], ...]:
    composition_a = Counter(atom.element for atom in structure_a.atoms)
    composition_b = Counter(atom.element for atom in structure_b.atoms)
    divisor_a = reduce(gcd, composition_a.values())
    divisor_b = reduce(gcd, composition_b.values())
    signature_a = tuple(sorted((e, n // divisor_a) for e, n in composition_a.items()))
    signature_b = tuple(sorted((e, n // divisor_b) for e, n in composition_b.items()))
    if signature_a != signature_b:
        raise ValueError("PDD comparison requires matching element ratios in the explicit unit cell.")
    return signature_a


def _collapse_into_groups(overlapping: np.ndarray) -> list[list[int]]:
    adjacency = squareform(overlapping)
    visited: set[int] = set()
    groups: list[list[int]] = []

    for start in range(adjacency.shape[0]):
        if start in visited:
            continue
        queue = [start]
        component: list[int] = []
        visited.add(start)
        while queue:
            node = queue.pop()
            component.append(node)
            neighbors = np.flatnonzero(adjacency[node])
            for neighbor in neighbors.tolist():
                if neighbor in visited:
                    continue
                visited.add(neighbor)
                queue.append(neighbor)
        groups.append(sorted(component))

    return groups


def _nearest_neighbours(motif: np.ndarray, cell: np.ndarray, k: int, *, workers=1) -> np.ndarray:
    """Certify image coverage at the returned kth-neighbour radius.

    Previously: append integer shells until distances stop changing. A short
    lattice vector with large coefficients in a skew basis can occur much later,
    so a temporary plateau did not establish completeness.

    With wrapped fractional positions, any pair within R has image coefficients
    |n_j| < 1 + R*||inverse(cell)[:,j]||. The integer bounds below therefore cover
    every such pair. Reduction changes the basis only, not the periodic point set.
    """
    from ase.geometry import minkowski_reduce
    cell, _ = minkowski_reduce(np.asarray(cell, dtype=float))
    inverse = np.linalg.inv(cell)
    motif = np.mod(np.asarray(motif) @ inverse, 1.0) @ cell
    reciprocal_lengths = np.linalg.norm(inverse, axis=0)
    volume = abs(np.linalg.det(cell))
    radius = (3.0 * (k + 1) * volume / (4.0 * np.pi * len(motif))) ** (1.0 / 3.0)
    radius = max(radius * 1.2, 1e-6)
    while True:
        bounds = np.ceil(radius * reciprocal_lengths + 1e-12).astype(int)
        shifts = np.array(list(product(*(range(-b, b + 1) for b in bounds)))) @ cell
        cloud = (shifts[:, None, :] + motif[None, :, :]).reshape(-1, 3)
        distances, _ = KDTree(cloud).query(motif, k=k + 1, workers=workers)
        kth = float(np.max(distances[:, -1]))
        if np.isfinite(kth) and kth <= radius:
            return distances[:, 1:]
        radius = max(radius * 1.5, kth * (1.0 + 1e-12)) if np.isfinite(kth) else radius * 2


@lru_cache(maxsize=32)
def _transport_constraints(n_sources, n_sinks):
    columns = np.arange(n_sources * n_sinks)
    rows = np.concatenate((np.repeat(np.arange(n_sources), n_sinks),
                           n_sources + np.tile(np.arange(n_sinks), n_sources)))
    return coo_matrix((np.ones(2 * len(columns)), (rows, np.tile(columns, 2))),
                      shape=(n_sources + n_sinks, len(columns))).tocsr()


def _earth_movers_distance(
    source_weights: np.ndarray,
    sink_weights: np.ndarray,
    cost_matrix: np.ndarray,
) -> float:
    if source_weights.ndim != 1 or sink_weights.ndim != 1:
        raise ValueError("EMD weights must be one-dimensional.")
    if cost_matrix.shape != (source_weights.shape[0], sink_weights.shape[0]):
        raise ValueError("Cost matrix shape does not match the source/sink weights.")

    source_total = float(np.sum(source_weights))
    sink_total = float(np.sum(sink_weights))
    if not np.isclose(source_total, sink_total, atol=1e-10, rtol=0):
        raise ValueError("EMD requires matching total weights.")

    n_sources, n_sinks = cost_matrix.shape
    if n_sources == 1:
        return float(cost_matrix[0] @ sink_weights)
    if n_sinks == 1:
        return float(cost_matrix[:, 0] @ source_weights)
    # Equal masses reduce exactly to bipartite assignment, avoiding a general LP.
    if (n_sources == n_sinks and np.all(source_weights == source_weights[0])
            and np.array_equal(source_weights, sink_weights)):
        rows, cols = linear_sum_assignment(cost_matrix)
        return float(source_weights[0] * cost_matrix[rows, cols].sum())
    c = cost_matrix.ravel()
    a_eq = _transport_constraints(n_sources, n_sinks)

    b_eq = np.concatenate((source_weights, sink_weights))
    attempts: list[tuple[str, object]] = []
    for method in ("highs", "highs-ds", "highs-ipm"):
        result = linprog(
            c,
            A_eq=a_eq.tocsr(),
            b_eq=b_eq,
            bounds=(0.0, None),
            method=method,
        )
        attempts.append((method, result))
        if result.success:
            return max(0.0, float(result.fun))

    details = "; ".join(
        f"{method}: {result.message}"
        for method, result in attempts
    )
    raise ValueError(f"Could not solve the PDD transport problem: {details}")


def _require_explict_unit_cell(structure: CrystalStructure, *, context: str) -> None:
    if not structure.explict_unit_cell:
        raise ValueError(context)
