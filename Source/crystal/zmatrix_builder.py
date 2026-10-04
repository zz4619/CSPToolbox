"""Z-matrix construction from Cartesian molecules, and evaluation of a shared
Z-matrix topology on other conformers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable

import numpy as np

from .connectivity import (
    build_connectivity,
    chemical_signature,
    component_records,
    connected_components,
    graph_isomorphism_mapping,
    local_graph,
)
from .records import (
    DEFAULT_COVALENT_SCALE,
    DEFAULT_LINEAR_THRESHOLD,
    ELEMENT_PRIORITY,
    AtomRecord,
    ZMatrixEntry,
    ZMatrixRepresentation,
    zmatrix_topology_signature,
)

if TYPE_CHECKING:
    from .structure import CrystalStructure


@dataclass(frozen=True)
class _ZMatrixTemplateRow:
    atom_index: int
    bond_to: int | None
    angle_to: int | None
    dihedral_to: int | None
    used_fallback_angle: bool
    used_fallback_dihedral: bool


def generate_zmatrices(
    structure: CrystalStructure,
    covalent_scale: float = DEFAULT_COVALENT_SCALE,
    linear_threshold: float = DEFAULT_LINEAR_THRESHOLD,
) -> list[ZMatrixRepresentation]:
    """Generate one Z-matrix representation for each detected molecule."""

    adjacency = build_connectivity(structure, covalent_scale)
    components = connected_components(adjacency)
    component_data = []
    template_cache: dict[str, dict[str, object]] = {}

    for molecule_index, component in enumerate(components, start=1):
        ordered_original = sorted(component)
        mapping = {original: local for local, original in enumerate(ordered_original)}
        graph = local_graph(component, adjacency, mapping)
        molecule = component_records(structure, component, adjacency, unwrap=True)
        symbols = [atom.element for atom in molecule]
        signature = chemical_signature(symbols, graph)
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
            template, warnings = _build_zmatrix_template(
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

        mapping = graph_isomorphism_mapping(
            representative_symbols=rep_symbols,
            representative_graph=rep_graph,
            symbols=symbols,
            graph=graph,
        )
        remapped_template = _remap_zmatrix_template(template, mapping)
        entries, ordered_labels, warnings = _apply_zmatrix_template(
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


def apply_zmatrix_topology(
    structure: CrystalStructure,
    topology: ZMatrixRepresentation,
    covalent_scale: float = DEFAULT_COVALENT_SCALE,
    linear_threshold: float = DEFAULT_LINEAR_THRESHOLD,
) -> ZMatrixRepresentation:
    """Evaluate one fixed Z-matrix topology on this structure.

    Atom order and all reference indices are copied unchanged from
    ``topology``. Only bond lengths, bond angles, and dihedral angles are
    evaluated from this structure's Cartesian coordinates. Atom labels are
    therefore the required identity map between conformers.
    """

    signature = zmatrix_topology_signature(topology)
    topology_labels = [row[0] for row in signature]
    if len(topology_labels) != len(set(topology_labels)):
        raise ValueError("Z-matrix topology contains duplicate atom labels.")

    atoms_by_label = {atom.label: atom for atom in structure.atoms}
    if len(atoms_by_label) != len(structure.atoms):
        raise ValueError("Crystal structure contains duplicate atom labels.")

    expected_labels = set(topology_labels)
    observed_labels = set(atoms_by_label)
    if observed_labels != expected_labels:
        missing = sorted(expected_labels - observed_labels)
        extra = sorted(observed_labels - expected_labels)
        raise ValueError(
            "Crystal structure atom labels do not match the Z-matrix "
            f"topology; missing={missing}, extra={extra}."
        )

    ordered_atoms = [atoms_by_label[label] for label in topology_labels]
    for position, (atom, row) in enumerate(zip(ordered_atoms, signature), start=1):
        _label, expected_element, _bond_to, _angle_to, _dihedral_to = row
        if atom.element != expected_element:
            raise ValueError(
                f"Z-matrix row {position} label {atom.label!r} expects "
                f"element {expected_element}, found {atom.element}."
            )

    atom_index_by_label = {
        atom.label: index for index, atom in enumerate(structure.atoms)
    }
    adjacency = build_connectivity(structure, covalent_scale)
    ordered_coordinates = np.array(
        [atom.coordinates for atom in ordered_atoms],
        dtype=float,
    )
    entries: list[ZMatrixEntry] = []
    warnings: list[str] = []

    for position, (label, element, bond_to, angle_to, dihedral_to) in enumerate(
        signature,
        start=1,
    ):
        expected_reference_count = (
            0 if position == 1 else 1 if position == 2 else 2 if position == 3 else 3
        )
        observed_reference_count = sum(
            reference is not None
            for reference in (bond_to, angle_to, dihedral_to)
        )
        if observed_reference_count != expected_reference_count:
            raise ValueError(
                f"Z-matrix row {position} has {observed_reference_count} "
                f"references; expected {expected_reference_count}."
            )

        for reference_name, reference in (
            ("bond", bond_to),
            ("angle", angle_to),
            ("dihedral", dihedral_to),
        ):
            if reference is not None and not 1 <= reference < position:
                raise ValueError(
                    f"Z-matrix row {position} {reference_name} reference "
                    f"{reference} does not point to an earlier row."
                )

        bond_length = None
        angle_degrees = None
        dihedral_degrees = None
        atom_position = position - 1

        if bond_to is not None:
            bond_position = bond_to - 1
            atom_index = atom_index_by_label[label]
            bond_label = topology_labels[bond_position]
            bond_atom_index = atom_index_by_label[bond_label]
            bonded_neighbors = {
                edge.neighbor for edge in adjacency[atom_index]
            }
            if bond_atom_index not in bonded_neighbors:
                raise ValueError(
                    f"Z-matrix row {position} bond reference {bond_to} "
                    f"({label}-{bond_label}) is not bonded in this structure."
                )
            bond_length = float(
                np.linalg.norm(
                    ordered_coordinates[atom_position]
                    - ordered_coordinates[bond_position]
                )
            )

        if angle_to is not None and bond_to is not None:
            angle_position = angle_to - 1
            angle_degrees = _angle_value(
                ordered_coordinates,
                atom_position,
                bond_to - 1,
                angle_position,
            )
            if _linear_margin(angle_degrees) < linear_threshold:
                warnings.append(
                    f"Atom {label} uses a near-linear angle reference "
                    f"({angle_degrees:.2f} deg)."
                )

        if (
            dihedral_to is not None
            and angle_to is not None
            and bond_to is not None
        ):
            dihedral_position = dihedral_to - 1
            dihedral_degrees = _dihedral_value(
                ordered_coordinates,
                atom_position,
                bond_to - 1,
                angle_to - 1,
                dihedral_position,
            )
            anchor_angle = _angle_value(
                ordered_coordinates,
                bond_to - 1,
                angle_to - 1,
                dihedral_position,
            )
            if _linear_margin(anchor_angle) < linear_threshold:
                warnings.append(
                    f"Atom {label} uses a near-linear dihedral anchor "
                    f"({anchor_angle:.2f} deg)."
                )

        entries.append(
            ZMatrixEntry(
                label=label,
                element=element,
                bond_to=bond_to,
                bond_length=bond_length,
                angle_to=angle_to,
                angle_degrees=angle_degrees,
                dihedral_to=dihedral_to,
                dihedral_degrees=dihedral_degrees,
            )
        )

    return ZMatrixRepresentation(
        molecule_index=topology.molecule_index,
        entries=entries,
        ordered_atom_labels=list(topology_labels),
        warnings=warnings,
    )


def _build_zmatrix_template(
    molecule: list[AtomRecord],
    graph: dict[int, set[int]],
    linear_threshold: float,
) -> tuple[list[_ZMatrixTemplateRow], list[str]]:
    coords = np.array([atom.coordinates for atom in molecule], dtype=float)
    symbols = [atom.element for atom in molecule]
    order, parent = _atom_order(symbols, coords, graph)
    order, branch_kept_child = _branch_improper_atom_order(
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
            branch_improper = _choose_branch_improper_reference(
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
            angle_deg = _angle_value(coords, atom_i, atom_b, atom_a)
            dihedral_deg = _dihedral_value(coords, atom_i, atom_b, atom_a, atom_d)
        else:
            atom_a, angle_deg, angle_is_bonded_path = _choose_angle_reference(
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

        angle_margin = _linear_margin(angle_deg)
        if angle_margin < linear_threshold:
            warnings.append(f"Atom index {atom_i} uses a near-linear angle reference.")

        if step == 2:
            template.append(_ZMatrixTemplateRow(atom_i, atom_b, atom_a, None, used_fallback_angle, False))
            continue

        if branch_improper is None:
            atom_d, dihedral_deg, dihedral_is_bonded_path = _choose_dihedral_reference(
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

        dihedral_margin = _linear_margin(_angle_value(coords, atom_b, atom_a, atom_d))
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
            angle_degrees = _angle_value(coords, atom_i, row.bond_to, row.angle_to)
        if row.dihedral_to is not None and row.angle_to is not None and row.bond_to is not None:
            dihedral_degrees = _dihedral_value(
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
        if angle_degrees is not None and _linear_margin(angle_degrees) < linear_threshold:
            warnings.append(f"Atom {label} uses a near-linear angle reference ({angle_degrees:.2f} deg).")
        if (
            row.dihedral_to is not None
            and row.angle_to is not None
            and row.bond_to is not None
            and _linear_margin(_angle_value(coords, row.bond_to, row.angle_to, row.dihedral_to))
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
    symbols: list[str],
    coords: np.ndarray,
    graph: dict[int, set[int]],
) -> tuple[list[int], dict[int, int | None]]:
    heavy_atoms = [index for index, symbol in enumerate(symbols) if symbol != "H"]
    if not heavy_atoms:
        root = _choose_root(symbols, coords)
        visited = {root}
        parent: dict[int, int | None] = {root: None}
        order = [root]
        frontier = {neighbor: root for neighbor in graph[root]}

        while frontier:
            next_atom = min(
                frontier,
                key=lambda index: _traversal_priority(index, symbols, graph, coords),
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
    root = _choose_root(symbols, coords)
    visited_heavy = {root}
    parent: dict[int, int | None] = {root: None}
    heavy_order = [root]
    frontier = {neighbor: root for neighbor in heavy_graph[root]}

    while frontier:
        next_atom = min(
            frontier,
            key=lambda index: _traversal_priority(index, symbols, graph, coords),
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


def _angle_value(coords: np.ndarray, atom_i: int, atom_j: int, atom_k: int) -> float:
    vec_ji = coords[atom_i] - coords[atom_j]
    vec_jk = coords[atom_k] - coords[atom_j]
    denom = np.linalg.norm(vec_ji) * np.linalg.norm(vec_jk)
    if denom < 1e-12:
        return float("nan")
    cosine = float(np.clip(np.dot(vec_ji, vec_jk) / denom, -1.0, 1.0))
    return float(np.degrees(np.arccos(cosine)))


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


def _linear_margin(angle_degrees: float) -> float:
    if np.isnan(angle_degrees):
        return -1.0
    return float(min(abs(angle_degrees), abs(180.0 - angle_degrees)))


def _unique(sequence: Iterable[int]) -> list[int]:
    seen: set[int] = set()
    result: list[int] = []
    for value in sequence:
        if value not in seen:
            seen.add(value)
            result.append(value)
    return result


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
            key=lambda index: _candidate_priority(index, atom_b, symbols, graph, coords),
            reverse=True,
        )
    )
    fallback_candidates = sorted(
        defined - {atom_b, atom_i} - set(bonded_candidates),
        key=lambda index: _candidate_priority(index, atom_b, symbols, graph, coords),
        reverse=True,
    )

    best_candidate = None
    best_margin = -1.0
    for candidate in _unique(bonded_candidates):
        angle_deg = _angle_value(coords, atom_i, atom_b, candidate)
        margin = _linear_margin(angle_deg)
        if margin > best_margin:
            best_margin = margin
            best_candidate = candidate
    if best_candidate is not None:
        return best_candidate, _angle_value(coords, atom_i, atom_b, best_candidate), True

    best_candidate = None
    best_margin = -1.0
    for candidate in _unique(fallback_candidates):
        angle_deg = _angle_value(coords, atom_i, atom_b, candidate)
        margin = _linear_margin(angle_deg)
        if margin > best_margin:
            best_margin = margin
            best_candidate = candidate
        if margin >= threshold:
            return candidate, angle_deg, False

    if best_candidate is None:
        raise ValueError(f"Failed to choose an angle reference for atom {atom_i}.")
    return best_candidate, _angle_value(coords, atom_i, atom_b, best_candidate), False


def _choose_branch_improper_reference(
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
            key=lambda index: _branch_improper_reference_priority(
                index,
                atom_b,
                kept_child,
                symbols,
                coords,
            ),
        )
    )

    refs = _unique(candidate_refs)
    if len(refs) < 2:
        return None

    best_pair: tuple[int, int] | None = None
    best_score: tuple[bool, float, float, int, int] | None = None
    for left_index, left in enumerate(refs):
        for right_index, right in enumerate(refs[left_index + 1:], start=left_index + 1):
            angle_deg = _angle_value(coords, atom_i, atom_b, left)
            anchor_angle = _angle_value(coords, atom_b, left, right)
            score = (
                right not in graph[left],
                _linear_margin(angle_deg),
                _linear_margin(anchor_angle),
                -left_index,
                -right_index,
            )
            if best_score is None or score > best_score:
                best_score = score
                best_pair = (left, right)

    return best_pair


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
            key=lambda index: _candidate_priority(index, atom_a, symbols, graph, coords),
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
            key=lambda index: _candidate_priority(index, atom_a, symbols, graph, coords),
            reverse=True,
        )
    )

    best_candidate = None
    best_margin = -1.0
    for candidate in _unique(bonded_candidates):
        anchor_angle = _angle_value(coords, atom_b, atom_a, candidate)
        margin = _linear_margin(anchor_angle)
        if margin > best_margin:
            best_margin = margin
            best_candidate = candidate
    if best_candidate is not None:
        return best_candidate, _dihedral_value(coords, atom_i, atom_b, atom_a, best_candidate), True

    best_candidate = None
    best_margin = -1.0
    for candidate in _unique(fallback_candidates):
        anchor_angle = _angle_value(coords, atom_b, atom_a, candidate)
        margin = _linear_margin(anchor_angle)
        if margin > best_margin:
            best_margin = margin
            best_candidate = candidate
        if margin >= threshold:
            return candidate, _dihedral_value(coords, atom_i, atom_b, atom_a, candidate), False

    if best_candidate is None:
        raise ValueError(f"Failed to choose a dihedral reference for atom {atom_i}.")
    return best_candidate, _dihedral_value(coords, atom_i, atom_b, atom_a, best_candidate), False
