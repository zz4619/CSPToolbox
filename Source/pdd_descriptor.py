"""PDD descriptor and distance utilities for explicit unit-cell structures."""

from __future__ import annotations

import collections
from collections import Counter
from dataclasses import dataclass
from itertools import combinations, product

import numpy as np
from ase.data import atomic_numbers, chemical_symbols
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist, pdist, squareform

from .crystal_structure import CrystalStructure


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
) -> PDDDescriptor:
    """Calculate the PDD descriptor for an explicit unit-cell structure."""

    _require_explict_unit_cell(
        structure,
        context="PDD calculation requires an explict unit cell structure.",
    )
    if k < 1:
        raise ValueError("k must be at least 1.")

    motif, cell, center_numbers = _explicit_structure_to_pdd_input(structure)
    distances = _nearest_neighbours(motif, cell, k)
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
        source_name=structure.name,
        k=k,
        weights=weights,
        distances=distances,
        center_elements=tuple(center_elements),
        typed=typed,
        collapse=collapse,
        collapse_tol=collapse_tol,
    )


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
    signature_a = tuple(sorted(composition_a.items()))
    signature_b = tuple(sorted(composition_b.items()))
    if signature_a != signature_b:
        raise ValueError("PDD comparison requires matching element counts in the explicit unit cell.")
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


def _distance_squared(xy: tuple[int, ...], z: int) -> int:
    return int(z**2 + sum(value**2 for value in xy))


def _distkey(point: tuple[int, ...] | list[int]) -> int:
    return int(sum(value**2 for value in point))


def _generate_integer_lattice(dims: int):
    ymax = collections.defaultdict(int)
    d = 0

    if dims == 1:
        yield np.array([[0]], dtype=int)
        while True:
            d += 1
            yield np.array([[-d], [d]], dtype=int)

    while True:
        positive_int_lattice: list[tuple[int, ...]] = []
        while True:
            batch: list[tuple[int, ...]] = []
            for xy in product(range(d + 1), repeat=dims - 1):
                if _distance_squared(xy, ymax[xy]) <= d**2:
                    batch.append((*xy, ymax[xy]))
                    ymax[xy] += 1
            if not batch:
                break
            positive_int_lattice += batch
        positive_int_lattice.sort(key=_distkey)

        int_lattice: list[tuple[int, ...] | list[int]] = []
        for point in positive_int_lattice:
            int_lattice.append(point)
            for n_reflections in range(1, dims + 1):
                for indexes in combinations(range(dims), n_reflections):
                    if all(point[i] for i in indexes):
                        reflected = list(point)
                        for i in indexes:
                            reflected[i] *= -1
                        int_lattice.append(reflected)

        yield np.array(int_lattice, dtype=int)
        d += 1


def _generate_concentric_cloud(motif: np.ndarray, cell: np.ndarray):
    int_lattice_generator = _generate_integer_lattice(cell.shape[0])

    while True:
        int_lattice = next(int_lattice_generator) @ cell
        yield np.concatenate([motif + translation for translation in int_lattice])


def _nearest_neighbours(motif: np.ndarray, cell: np.ndarray, k: int) -> np.ndarray:
    cloud_generator = _generate_concentric_cloud(motif, cell)
    n_points = 0
    cloud_batches: list[np.ndarray] = []
    while n_points <= k:
        batch = next(cloud_generator)
        n_points += batch.shape[0]
        cloud_batches.append(batch)
    cloud_batches.append(next(cloud_generator))
    cloud = np.concatenate(cloud_batches)

    tree = KDTree(cloud, compact_nodes=False, balanced_tree=False)
    candidate_distances, _ = tree.query(motif, k=k + 1, workers=-1)
    previous = np.zeros_like(candidate_distances)

    while not np.allclose(previous, candidate_distances, atol=1e-12, rtol=0):
        previous = candidate_distances
        cloud = np.vstack((cloud, next(cloud_generator), next(cloud_generator)))
        tree = KDTree(cloud, compact_nodes=False, balanced_tree=False)
        candidate_distances, _ = tree.query(motif, k=k + 1, workers=-1)

    return candidate_distances[:, 1:]


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
    n_variables = n_sources * n_sinks
    c = cost_matrix.ravel()
    a_eq = lil_matrix((n_sources + n_sinks, n_variables), dtype=float)

    for source_index in range(n_sources):
        a_eq[source_index, source_index * n_sinks : (source_index + 1) * n_sinks] = 1.0
    for sink_index in range(n_sinks):
        a_eq[n_sources + sink_index, sink_index::n_sinks] = 1.0

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
