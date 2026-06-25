"""Cartesian reconstruction and viewer data preparation for Z-matrices."""

from __future__ import annotations

import math
from typing import Iterable

from .model import (
    ViewerAtom,
    ViewerBond,
    ViewerDihedral,
    ViewerMolecule,
    ZMatrixAtom,
    ZMatrixDocument,
)


ELEMENT_COLORS = {
    "H": "#f7f7f2",
    "C": "#30343b",
    "N": "#2f66d0",
    "O": "#d83b36",
    "F": "#3dbf78",
    "P": "#d46b18",
    "S": "#d8b52f",
    "Cl": "#43a047",
    "Br": "#8b3a2f",
    "I": "#6b4bb3",
}

DISPLAY_RADII = {
    "H": 0.13,
    "C": 0.20,
    "N": 0.19,
    "O": 0.18,
    "F": 0.17,
    "P": 0.24,
    "S": 0.24,
    "Cl": 0.23,
    "Br": 0.25,
    "I": 0.27,
}

COVALENT_RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
}

DEFAULT_COLOR = "#7b8794"
DEFAULT_DISPLAY_RADIUS = 0.20
DEFAULT_COVALENT_RADIUS = 0.77

Vec3 = tuple[float, float, float]


def build_viewer_molecule(
    document: ZMatrixDocument,
    *,
    infer_bonds: bool = True,
    covalent_scale: float = 1.25,
) -> ViewerMolecule:
    """Reconstruct a parsed Z-matrix and prepare the browser payload."""

    coordinates = reconstruct_coordinates(document.atoms)
    viewer_atoms = tuple(
        ViewerAtom(
            index=atom.row_index,
            label=atom.label,
            element=atom.element,
            coordinates=coordinates[atom.row_index - 1],
            color=ELEMENT_COLORS.get(atom.element, DEFAULT_COLOR),
            display_radius=DISPLAY_RADII.get(atom.element, DEFAULT_DISPLAY_RADIUS),
        )
        for atom in document.atoms
    )

    construction_bonds = {
        _normal_pair(atom.row_index, atom.bond_to)
        for atom in document.atoms
        if atom.bond_to is not None
    }
    explicit_bonds = set(document.explicit_bonds)
    inferred_bonds = (
        _infer_bonds(document.atoms, coordinates, covalent_scale=covalent_scale)
        if infer_bonds
        else set()
    )
    all_bonds = construction_bonds | explicit_bonds | inferred_bonds
    viewer_bonds = tuple(
        ViewerBond(left=left, right=right, kind=_bond_kind((left, right), construction_bonds, explicit_bonds))
        for left, right in sorted(all_bonds)
    )
    bond_set = frozenset(all_bonds)
    dihedrals = tuple(_build_dihedral(atom, document.atoms, bond_set) for atom in document.atoms if atom.has_dihedral)

    return ViewerMolecule(
        title=document.title,
        atoms=viewer_atoms,
        bonds=viewer_bonds,
        dihedrals=dihedrals,
        warnings=tuple(_geometry_warnings(document.atoms, coordinates)),
        source_name=document.source_name,
    )


def reconstruct_coordinates(atoms: Iterable[ZMatrixAtom]) -> tuple[Vec3, ...]:
    """Reconstruct Cartesian coordinates from 1-based Z-matrix rows."""

    rows = tuple(atoms)
    if not rows:
        raise ValueError("Cannot reconstruct an empty Z-matrix.")

    coordinates: list[Vec3] = []
    for atom in rows:
        row = atom.row_index
        if row == 1:
            coordinates.append((0.0, 0.0, 0.0))
            continue

        if atom.bond_to is None or atom.bond_length is None:
            raise ValueError(f"Row {row} is missing a bond reference.")

        bond_origin = coordinates[atom.bond_to - 1]
        if row == 2:
            coordinates.append(_add(bond_origin, (atom.bond_length, 0.0, 0.0)))
            continue

        if atom.angle_to is None or atom.angle_degrees is None:
            raise ValueError(f"Row {row} is missing an angle reference.")

        angle_origin = coordinates[atom.angle_to - 1]
        bond_axis = _unit(_sub(angle_origin, bond_origin))
        theta = math.radians(atom.angle_degrees)

        if row == 3:
            perpendicular = _perpendicular_unit(bond_axis)
            direction = _add(_scale(bond_axis, math.cos(theta)), _scale(perpendicular, math.sin(theta)))
            coordinates.append(_add(bond_origin, _scale(direction, atom.bond_length)))
            continue

        if atom.dihedral_to is None or atom.dihedral_degrees is None:
            raise ValueError(f"Row {row} is missing a dihedral reference.")

        dihedral_origin = coordinates[atom.dihedral_to - 1]
        phi = math.radians(atom.dihedral_degrees)
        direction = _direction_from_internal(
            bond_origin,
            angle_origin,
            dihedral_origin,
            theta=theta,
            phi=phi,
        )
        coordinates.append(_add(bond_origin, _scale(direction, atom.bond_length)))

    return tuple(coordinates)


def measure_angle_degrees(coords: tuple[Vec3, ...], atom_i: int, atom_j: int, atom_k: int) -> float:
    """Measure angle i-j-k from 1-based indices."""

    vec_ji = _sub(coords[atom_i - 1], coords[atom_j - 1])
    vec_jk = _sub(coords[atom_k - 1], coords[atom_j - 1])
    denom = _norm(vec_ji) * _norm(vec_jk)
    if denom < 1e-12:
        return float("nan")
    cosine = max(-1.0, min(1.0, _dot(vec_ji, vec_jk) / denom))
    return math.degrees(math.acos(cosine))


def measure_dihedral_degrees(
    coords: tuple[Vec3, ...],
    atom_i: int,
    atom_j: int,
    atom_k: int,
    atom_l: int,
) -> float:
    """Measure the CSPToolbox dihedral convention from 1-based indices."""

    p0, p1, p2, p3 = (
        coords[atom_i - 1],
        coords[atom_j - 1],
        coords[atom_k - 1],
        coords[atom_l - 1],
    )
    b0 = _sub(p0, p1)
    b1 = _sub(p2, p1)
    b2 = _sub(p3, p2)

    b1_norm = _norm(b1)
    if b1_norm < 1e-12:
        return float("nan")
    b1_unit = _scale(b1, 1.0 / b1_norm)

    v = _sub(b0, _scale(b1_unit, _dot(b0, b1_unit)))
    w = _sub(b2, _scale(b1_unit, _dot(b2, b1_unit)))
    if _norm(v) < 1e-12 or _norm(w) < 1e-12:
        return float("nan")

    x_value = _dot(v, w)
    y_value = _dot(_cross(b1_unit, v), w)
    return math.degrees(math.atan2(y_value, x_value))


def _direction_from_internal(
    bond_origin: Vec3,
    angle_origin: Vec3,
    dihedral_origin: Vec3,
    *,
    theta: float,
    phi: float,
) -> Vec3:
    bond_axis = _unit(_sub(angle_origin, bond_origin))
    reference = _sub(dihedral_origin, angle_origin)
    projected = _sub(reference, _scale(bond_axis, _dot(reference, bond_axis)))
    if _norm(projected) < 1e-10:
        projected = _perpendicular_unit(bond_axis)
    else:
        projected = _unit(projected)

    handed = _unit(_cross(bond_axis, projected))
    dihedral_plane_direction = _sub(
        _scale(projected, math.cos(phi)),
        _scale(handed, math.sin(phi)),
    )
    direction = _add(
        _scale(bond_axis, math.cos(theta)),
        _scale(dihedral_plane_direction, math.sin(theta)),
    )
    return _unit(direction)


def _infer_bonds(
    atoms: tuple[ZMatrixAtom, ...],
    coords: tuple[Vec3, ...],
    *,
    covalent_scale: float,
) -> set[tuple[int, int]]:
    bonds: set[tuple[int, int]] = set()
    for left_index, left_atom in enumerate(atoms, start=1):
        left_radius = COVALENT_RADII.get(left_atom.element, DEFAULT_COVALENT_RADIUS)
        for right_index in range(left_index + 1, len(atoms) + 1):
            right_atom = atoms[right_index - 1]
            right_radius = COVALENT_RADII.get(right_atom.element, DEFAULT_COVALENT_RADIUS)
            cutoff = covalent_scale * (left_radius + right_radius)
            if _distance(coords[left_index - 1], coords[right_index - 1]) <= cutoff:
                bonds.add((left_index, right_index))
    return bonds


def _build_dihedral(
    atom: ZMatrixAtom,
    atoms: tuple[ZMatrixAtom, ...],
    bond_set: frozenset[tuple[int, int]],
) -> ViewerDihedral:
    if atom.bond_to is None or atom.angle_to is None or atom.dihedral_to is None:
        raise ValueError(f"Row {atom.row_index} is missing dihedral references.")

    atom_indices = (atom.row_index, atom.bond_to, atom.angle_to, atom.dihedral_to)
    labels = tuple(atoms[index - 1].label for index in atom_indices)
    kind = _classify_dihedral(atom_indices, bond_set)
    return ViewerDihedral(
        id=f"row-{atom.row_index}",
        row_index=atom.row_index,
        atom_indices=atom_indices,
        atom_labels=labels,
        value_degrees=float(atom.dihedral_degrees),
        kind=kind,
        links=(_normal_pair(atom.row_index, atom.bond_to), _normal_pair(atom.bond_to, atom.angle_to), _normal_pair(atom.angle_to, atom.dihedral_to)),
    )


def _classify_dihedral(
    atom_indices: tuple[int, int, int, int],
    bond_set: frozenset[tuple[int, int]],
) -> str:
    atom_i, atom_j, atom_k, atom_l = atom_indices
    has_ij = _normal_pair(atom_i, atom_j) in bond_set
    has_jk = _normal_pair(atom_j, atom_k) in bond_set
    proper = has_ij and has_jk and _normal_pair(atom_k, atom_l) in bond_set
    improper = has_ij and has_jk and _normal_pair(atom_j, atom_l) in bond_set
    if proper and improper:
        return "ambiguous"
    if proper:
        return "proper"
    if improper:
        return "improper"
    return "fallback"


def _geometry_warnings(
    atoms: tuple[ZMatrixAtom, ...],
    coords: tuple[Vec3, ...],
) -> list[str]:
    warnings: list[str] = []
    for atom in atoms:
        if atom.angle_to is None or atom.bond_to is None or atom.angle_degrees is None:
            continue
        angle = measure_angle_degrees(coords, atom.row_index, atom.bond_to, atom.angle_to)
        if abs(angle - atom.angle_degrees) > 1e-5:
            warnings.append(f"Row {atom.row_index} reconstructed angle differs from source.")
        if atom.has_dihedral:
            measured = measure_dihedral_degrees(
                coords,
                atom.row_index,
                atom.bond_to,
                atom.angle_to,
                atom.dihedral_to,
            )
            if abs(_angle_delta(measured, atom.dihedral_degrees)) > 1e-5:
                warnings.append(f"Row {atom.row_index} reconstructed dihedral differs from source.")
    return warnings


def _bond_kind(
    pair: tuple[int, int],
    construction_bonds: set[tuple[int, int]],
    explicit_bonds: set[tuple[int, int]],
) -> str:
    if pair in construction_bonds:
        return "construction"
    if pair in explicit_bonds:
        return "explicit"
    return "inferred"


def _angle_delta(left: float, right: float) -> float:
    return ((left - right + 180.0) % 360.0) - 180.0


def _normal_pair(left: int, right: int | None) -> tuple[int, int]:
    if right is None:
        raise ValueError("Cannot build a bond pair with a missing reference.")
    if left == right:
        raise ValueError("A bond cannot connect an atom to itself.")
    return (left, right) if left < right else (right, left)


def _perpendicular_unit(vec: Vec3) -> Vec3:
    axis = (1.0, 0.0, 0.0) if abs(vec[0]) < 0.8 else (0.0, 1.0, 0.0)
    return _unit(_cross(vec, axis))


def _add(left: Vec3, right: Vec3) -> Vec3:
    return (left[0] + right[0], left[1] + right[1], left[2] + right[2])


def _sub(left: Vec3, right: Vec3) -> Vec3:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def _scale(vec: Vec3, factor: float) -> Vec3:
    return (vec[0] * factor, vec[1] * factor, vec[2] * factor)


def _dot(left: Vec3, right: Vec3) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _cross(left: Vec3, right: Vec3) -> Vec3:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _norm(vec: Vec3) -> float:
    return math.sqrt(_dot(vec, vec))


def _unit(vec: Vec3) -> Vec3:
    norm = _norm(vec)
    if norm < 1e-12:
        raise ValueError("Cannot normalize a zero-length vector during Z-matrix reconstruction.")
    return _scale(vec, 1.0 / norm)


def _distance(left: Vec3, right: Vec3) -> float:
    return _norm(_sub(left, right))

