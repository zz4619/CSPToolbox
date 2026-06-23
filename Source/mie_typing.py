"""Mie/FIT atom typing compatible with the CSO-RM labelling rules."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
from ase.data import atomic_numbers

from .crystal_structure import AtomRecord, CrystalStructure


SUPPORTED_LABEL_TYPES = {"fit", "fit_hp"}

CSORM_BOND_DISTANCE: dict[frozenset[int], float] = {
    frozenset((1, 6)): 1.2,
    frozenset((1, 7)): 1.2,
    frozenset((1, 8)): 1.3,
    frozenset((1, 16)): 1.4,
    frozenset((6, 6)): 1.7,
    frozenset((6, 7)): 1.7,
    frozenset((6, 8)): 1.5,
    frozenset((6, 9)): 1.5,
    frozenset((6, 14)): 2.0,
    frozenset((6, 16)): 2.0,
    frozenset((6, 17)): 1.9,
    frozenset((6, 53)): 2.3,
    frozenset((6, 55)): 2.0,
    frozenset((7, 7)): 1.7,
    frozenset((7, 8)): 1.7,
    frozenset((7, 16)): 1.75,
    frozenset((7, 17)): 2.0,
    frozenset((8, 16)): 1.7,
    frozenset((16, 16)): 2.2,
}


@dataclass(frozen=True)
class InterSpec:
    potential_form: str
    label_type: str
    allowed_site_types: frozenset[str]


@dataclass(frozen=True)
class FITAtomType:
    atom_label: str
    element: str
    fit_type: str
    molecule_index: int
    bonded_parent_label: str | None = None
    bonded_parent_element: str | None = None


@dataclass(frozen=True)
class MieAtomType:
    atom_label: str
    element: str
    mie_type: str
    molecule_index: int
    bonded_parent_label: str | None = None
    bonded_parent_element: str | None = None


def read_inter_spec(path: str | Path) -> InterSpec:
    """Read potential form, label type, and site labels from a CSO-RM .inter file."""

    inter_path = Path(path)
    potential_form = ""
    label_type = ""
    allowed_site_types: set[str] = set()
    section: str | None = None

    for raw_line in inter_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("!", maxsplit=1)[0].strip()
        if not line:
            continue
        upper = line.upper()
        if upper.startswith("START "):
            section = line.split(maxsplit=1)[1].strip().lower()
            continue
        if upper == "END":
            section = None
            continue

        if section == "spec":
            if ":" not in line:
                continue
            key, value = line.split(":", maxsplit=1)
            key = key.strip().lower()
            value = value.strip()
            if key == "potential form":
                potential_form = value.split()[0].lower()
            elif key == "label type":
                label_type = value.split()[0].lower()
        elif section == "params":
            parts = line.split()
            if len(parts) >= 2:
                allowed_site_types.add(parts[0])
                allowed_site_types.add(parts[1])

    if not potential_form:
        raise ValueError(f"No potential form found in {inter_path}")
    if not label_type:
        raise ValueError(f"No label type found in {inter_path}")
    if not allowed_site_types:
        raise ValueError(f"No potential site types found in {inter_path}")

    return InterSpec(
        potential_form=potential_form,
        label_type=label_type,
        allowed_site_types=frozenset(allowed_site_types),
    )


def detect_fit_atom_types(
    structure: CrystalStructure,
    *,
    label_type: str = "fit",
) -> list[FITAtomType]:
    """Assign CSO-RM FIT/FIT_Hp atom labels to each atom in a structure."""

    normalized_label_type = label_type.lower()
    if normalized_label_type not in SUPPORTED_LABEL_TYPES:
        raise ValueError(f"Unsupported FIT label type: {label_type}")

    typed_atoms: list[FITAtomType] = []
    adjacency = _build_csorm_connectivity(structure)
    components = _connected_components(adjacency)
    for molecule_index, component in enumerate(components, start=1):
        for atom_index in component:
            atom = structure.atoms[atom_index]
            if atom.element == "H":
                parent = _single_hydrogen_parent(atom_index, structure, adjacency)
                fit_type = _hydrogen_fit_type(
                    parent.element,
                    label_type=normalized_label_type,
                )
                typed_atoms.append(
                    FITAtomType(
                        atom_label=atom.label,
                        element=atom.element,
                        fit_type=fit_type,
                        molecule_index=molecule_index,
                        bonded_parent_label=parent.label,
                        bonded_parent_element=parent.element,
                    )
                )
            else:
                typed_atoms.append(
                    FITAtomType(
                        atom_label=atom.label,
                        element=atom.element,
                        fit_type=_heavy_atom_fit_type(atom.element),
                        molecule_index=molecule_index,
                    )
                )
    return typed_atoms


def detect_mie_atom_types(
    structure: CrystalStructure,
    *,
    label_type: str = "fit",
) -> list[MieAtomType]:
    """Assign Mie atom types using the requested FIT-style label scheme."""

    return [
        MieAtomType(
            atom_label=atom.atom_label,
            element=atom.element,
            mie_type=atom.fit_type,
            molecule_index=atom.molecule_index,
            bonded_parent_label=atom.bonded_parent_label,
            bonded_parent_element=atom.bonded_parent_element,
        )
        for atom in detect_fit_atom_types(structure, label_type=label_type)
    ]


def validate_mie_atom_types(
    structure: CrystalStructure,
    inter_path: str | Path,
) -> list[MieAtomType]:
    """Validate that detected Mie atom types are present in a .inter potential."""

    spec = read_inter_spec(inter_path)
    if spec.potential_form != "mie":
        raise ValueError(f"Expected Mie potential, found {spec.potential_form!r}")

    typed_atoms = detect_mie_atom_types(structure, label_type=spec.label_type)
    _validate_detected_types(
        {atom.mie_type for atom in typed_atoms},
        spec.allowed_site_types,
        context="Mie",
    )
    return typed_atoms


def validate_fit_atom_types(
    structure: CrystalStructure,
    inter_path: str | Path,
) -> list[FITAtomType]:
    """Validate FIT/FIT_Hp atom types against any CSO-RM .inter potential."""

    spec = read_inter_spec(inter_path)
    typed_atoms = detect_fit_atom_types(structure, label_type=spec.label_type)
    _validate_detected_types(
        {atom.fit_type for atom in typed_atoms},
        spec.allowed_site_types,
        context="FIT",
    )
    return typed_atoms


def validate_inter_atom_types(
    structure: CrystalStructure,
    inter_path: str | Path,
) -> list[FITAtomType] | list[MieAtomType]:
    """Validate atom labels required by a .inter file.

    Mie potentials return MieAtomType records. Other potentials using FIT-style
    labels return FITAtomType records.
    """

    spec = read_inter_spec(inter_path)
    if spec.potential_form == "mie":
        typed_atoms = detect_mie_atom_types(structure, label_type=spec.label_type)
        _validate_detected_types(
            {atom.mie_type for atom in typed_atoms},
            spec.allowed_site_types,
            context="Mie",
        )
        return typed_atoms

    typed_atoms = detect_fit_atom_types(structure, label_type=spec.label_type)
    _validate_detected_types(
        {atom.fit_type for atom in typed_atoms},
        spec.allowed_site_types,
        context="FIT",
    )
    return typed_atoms


def _validate_detected_types(
    detected_types: set[str],
    allowed_types: frozenset[str],
    *,
    context: str,
) -> None:
    missing_types = sorted(detected_types - allowed_types)
    if missing_types:
        raise ValueError(
            f"{context} atom types are not present in the potential: "
            + ", ".join(missing_types)
        )


def count_atom_types(typed_atoms: list[FITAtomType] | list[MieAtomType]) -> Counter[str]:
    """Count detected FIT or Mie labels."""

    counts: Counter[str] = Counter()
    for atom in typed_atoms:
        label = atom.fit_type if isinstance(atom, FITAtomType) else atom.mie_type
        counts[label] += 1
    return counts


def _build_csorm_connectivity(structure: CrystalStructure) -> dict[int, set[int]]:
    positions = np.array([atom.coordinates for atom in structure.atoms], dtype=float)
    fractional = structure.cell.scaled_positions(positions)
    cell = np.asarray(structure.cell.array, dtype=float)
    numbers = [atomic_numbers[atom.element] for atom in structure.atoms]
    adjacency: dict[int, set[int]] = {index: set() for index in range(len(structure.atoms))}

    for left in range(len(structure.atoms)):
        for right in range(left + 1, len(structure.atoms)):
            cutoff = CSORM_BOND_DISTANCE.get(frozenset((numbers[left], numbers[right])), 0.0)
            if cutoff <= 0.0:
                continue
            delta = fractional[right] - fractional[left]
            delta -= np.rint(delta)
            distance = float(np.linalg.norm(np.dot(delta, cell)))
            if distance < cutoff:
                adjacency[left].add(right)
                adjacency[right].add(left)
    return adjacency


def _connected_components(adjacency: dict[int, set[int]]) -> list[list[int]]:
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
            for neighbor in sorted(adjacency[node]):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    return components


def _single_hydrogen_parent(
    hydrogen_index: int,
    structure: CrystalStructure,
    adjacency: dict[int, set[int]],
) -> AtomRecord:
    hydrogen = structure.atoms[hydrogen_index]
    parents: list[AtomRecord] = []
    for neighbor in sorted(adjacency[hydrogen_index]):
        candidate = structure.atoms[neighbor]
        if candidate.element != "H":
            parents.append(candidate)

    if not parents:
        raise ValueError(f"Hydrogen {hydrogen.label} is not connected to any supported atom.")
    if len(parents) > 1:
        labels = ", ".join(parent.label for parent in parents)
        raise ValueError(f"Hydrogen {hydrogen.label} is connected to multiple atoms: {labels}")
    return parents[0]


def _hydrogen_fit_type(parent_element: str, *, label_type: str) -> str:
    if label_type == "fit":
        if parent_element in {"C", "S"}:
            return "H_c"
        if parent_element == "N":
            return "H_n"
        if parent_element == "O":
            return "H_o"
    elif label_type == "fit_hp":
        if parent_element in {"C", "S"}:
            return "H_c"
        if parent_element in {"N", "O"}:
            return "H_p"
    raise ValueError(f"Hydrogen bonded to unsupported parent element: {parent_element}")


def _heavy_atom_fit_type(element: str) -> str:
    symbol = _fit_element_symbol(element)
    if len(symbol) == 1:
        return f"{symbol}__"
    if len(symbol) == 2:
        return f"{symbol}_"
    raise ValueError(f"Unsupported FIT atom element: {element}")


def _fit_element_symbol(element: str) -> str:
    if not re.fullmatch(r"[A-Z][a-z]?", element):
        raise ValueError(f"Unsupported element symbol: {element}")
    return element
