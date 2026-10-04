"""Covalent bond graphs, molecular components and periodic unwrapping."""

from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING, Iterable

import networkx as nx
import numpy as np
from ase.data import atomic_numbers, covalent_radii
from ase.neighborlist import neighbor_list
from networkx.algorithms import isomorphism as nx_isomorphism

from .records import DEFAULT_COVALENT_SCALE, HYDROGEN_COVALENT_SCALE, AtomRecord, BondEdge

if TYPE_CHECKING:
    from .structure import CrystalStructure


def molecule_aware_representative_indices(
    structure: CrystalStructure,
    equivalent_atoms: Iterable[int],
) -> list[int]:
    """Choose one site from each symmetry orbit while preserving molecules.

    spglib reports atom orbits independently. Picking the first atom from
    every orbit can combine atoms from different symmetry-equivalent
    molecular copies, which leaves otherwise intact molecules fragmented in
    the reduced structure. Prefer representatives from the same connected
    component where possible; this still selects exactly one atom per orbit.
    """

    orbit_by_index = [int(value) for value in equivalent_atoms]
    orbit_to_indices: dict[int, list[int]] = defaultdict(list)
    for atom_index, orbit in enumerate(orbit_by_index):
        orbit_to_indices[orbit].append(atom_index)

    remaining_orbits = set(orbit_to_indices)
    if not remaining_orbits:
        return []

    adjacency = build_connectivity(structure, DEFAULT_COVALENT_SCALE)
    components = connected_components(adjacency)
    selected_indices: list[int] = []

    while remaining_orbits:
        best_component: list[int] | None = None
        best_cover: set[int] = set()
        for component in components:
            cover = {orbit_by_index[index] for index in component} & remaining_orbits
            if len(cover) > len(best_cover):
                best_component = component
                best_cover = cover

        if best_component is None or not best_cover:
            orbit = min(remaining_orbits)
            selected_indices.append(orbit_to_indices[orbit][0])
            remaining_orbits.remove(orbit)
            continue

        for atom_index in sorted(best_component):
            orbit = orbit_by_index[atom_index]
            if orbit not in remaining_orbits:
                continue
            selected_indices.append(atom_index)
            remaining_orbits.remove(orbit)

    return selected_indices


def build_connectivity(structure: CrystalStructure, covalent_scale: float) -> dict[int, list[BondEdge]]:
    """Periodic covalent bond graph: atom index -> bonds (neighbour, image shift, distance)."""

    atoms = structure.to_ase_atoms()
    numbers = [int(atomic_numbers[atom.element]) for atom in structure.atoms]

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
        cutoff_distance = denom * _pair_covalent_scale(
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

    adjacency = {index: [] for index in range(len(structure.atoms))}
    for left, right in sorted(allowed_keys):
        shift, distance, ratio = best_edges[(left, right)]
        adjacency[left].append(BondEdge(right, shift, distance, ratio))
        adjacency[right].append(
            BondEdge(left, tuple(-value for value in shift), distance, ratio)
        )
    return adjacency


def build_connectivity_from_template(
    structure: CrystalStructure,
    template_structure: "CrystalStructure",
    covalent_scale: float,
) -> dict[int, list[BondEdge]]:
    """Bond graph taken from a template structure with the same atom order."""

    if len(structure.atoms) != len(template_structure.atoms):
        raise ValueError(
            "Template-based molecule detection requires the same atom count in "
            "the current and template structures."
        )

    current_elements = [atom.element for atom in structure.atoms]
    template_elements = [atom.element for atom in template_structure.atoms]
    if current_elements != template_elements:
        raise ValueError(
            "Template-based molecule detection requires the same atom ordering "
            "and element sequence in the current and template structures."
        )

    template_adjacency = build_connectivity(template_structure, covalent_scale)
    positions = np.array([atom.coordinates for atom in structure.atoms], dtype=float)
    cell = np.asarray(structure.cell.array, dtype=float)
    fractional = structure.cell.scaled_positions(positions)
    numbers = [int(atomic_numbers[atom.element]) for atom in structure.atoms]
    radii = [float(covalent_radii[number]) for number in numbers]

    adjacency = {index: [] for index in range(len(structure.atoms))}
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


def connected_components(adjacency: dict[int, list[BondEdge]]) -> list[list[int]]:
    """Molecules as sorted lists of atom indices, ordered by their lowest index."""

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


def component_records(
    structure: CrystalStructure,
    component: list[int],
    adjacency: dict[int, list[BondEdge]],
    *,
    unwrap: bool,
) -> list[AtomRecord]:
    """AtomRecords of one molecule, optionally unwrapped across cell boundaries."""

    if not unwrap:
        return [structure.atoms[index] for index in component]

    unwrapped = unwrap_component(structure, component, adjacency)
    ordered_indices = sorted(component)
    return [
        AtomRecord(
            label=structure.atoms[index].label,
            element=structure.atoms[index].element,
            coordinates=tuple(float(value) for value in unwrapped[index]),
        )
        for index in ordered_indices
    ]


def unwrap_component(
    structure: CrystalStructure,
    component: list[int],
    adjacency: dict[int, list[BondEdge]],
) -> dict[int, np.ndarray]:
    """Cartesian positions of one molecule made contiguous across periodic images."""

    component_set = set(component)
    root = component[0]
    cell = np.asarray(structure.cell.array)
    wrapped_positions = np.array([atom.coordinates for atom in structure.atoms], dtype=float)

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


def local_graph(
    component: list[int],
    adjacency: dict[int, list[BondEdge]],
    mapping: dict[int, int],
) -> dict[int, set[int]]:
    """Bond graph of one molecule with atoms renumbered by ``mapping``."""

    component_set = set(component)
    graph = {mapping[index]: set() for index in component}
    for original in component:
        local = mapping[original]
        for edge in adjacency[original]:
            if edge.neighbor in component_set:
                graph[local].add(mapping[edge.neighbor])
    return graph


def chemical_signature(symbols: list[str], graph: dict[int, set[int]]) -> str:
    """Label-independent signature of a molecular graph (Weisfeiler-Lehman refinement of element and degree)."""

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


def graph_isomorphism_mapping(
    representative_symbols: list[str],
    representative_graph: dict[int, set[int]],
    symbols: list[str],
    graph: dict[int, set[int]],
) -> dict[int, int]:
    """Map representative-molecule indices onto an isomorphic molecule (first match)."""

    if representative_symbols == symbols and representative_graph == graph:
        return {index: index for index in range(len(symbols))}

    left = _graph_to_networkx(representative_symbols, representative_graph)
    right = _graph_to_networkx(symbols, graph)
    matcher = nx_isomorphism.GraphMatcher(
        left,
        right,
        node_match=lambda a, b: a["element"] == b["element"],
    )
    try:
        return next(matcher.isomorphisms_iter())
    except StopIteration as error:
        raise ValueError("Could not find an isomorphism between chemically identical molecules.") from error


def _pair_covalent_scale(
    left_number: int,
    right_number: int,
    covalent_scale: float,
) -> float:
    if left_number == 1 or right_number == 1:
        return max(float(covalent_scale), HYDROGEN_COVALENT_SCALE)
    return float(covalent_scale)


def _graph_to_networkx(symbols: list[str], graph: dict[int, set[int]]) -> nx.Graph:
    nx_graph = nx.Graph()
    for atom_index, symbol in enumerate(symbols):
        nx_graph.add_node(atom_index, element=symbol)
    for atom_index, neighbors in graph.items():
        for neighbor in neighbors:
            if atom_index < neighbor:
                nx_graph.add_edge(atom_index, neighbor)
    return nx_graph
