"""CIF reading, symmetry expansion to an explicit unit cell, and writing."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import re
import shlex
from pathlib import Path
from typing import TYPE_CHECKING
import warnings

import numpy as np
from ase import Atoms
from ase.cell import Cell

from .elements import guess_element_from_label, normalize_element_symbol
from .fileio import parse_bool_tag, strip_quotes
from .lattice import canonicalize_fractional, frac_to_cart, fractional_minimum_image_distance, lattice_matrix
from .records import DEFAULT_CIF_SITE_MERGE_TOLERANCE, AtomRecord, CifExpansionReport
from .symmetry import apply_symmetry_operation

if TYPE_CHECKING:
    from .structure import CrystalStructure


@dataclass(frozen=True)
class _ExpandedCifAtomSite:
    label: str
    element: str
    frac: tuple[float, float, float]
    occupancy: float | None


def read_cif(cls: type[CrystalStructure], path: Path) -> "CrystalStructure":
    """Read a CIF's atom-site loop as written (no symmetry expansion)."""

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
    space_group = strip_quotes(space_group)

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
        label = strip_quotes(row[label_idx])
        element = (
            strip_quotes(row[type_idx]) if type_idx is not None else guess_element_from_label(label)
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
                element=normalize_element_symbol(element),
                coordinates=tuple(float(value) for value in cart),
            )
        )

    return cls(
        atoms=atoms,
        cell_parameters=cell_parameters,
        space_group=space_group,
        name=path.stem or "CrystalStructure",
        explict_unit_cell=parse_bool_tag(data.get("_csptoolbox_explict_unit_cell"), default=False),
    )


def expand_cif_unit_cell(
    cls: type[CrystalStructure],
    path: str | Path,
    *,
    site_merge_tolerance: float = DEFAULT_CIF_SITE_MERGE_TOLERANCE,
) -> tuple["CrystalStructure", CifExpansionReport]:
    """Expand a CIF's sites with its symmetry operations; return the cell and a report."""

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
    lattice = lattice_matrix(*cell_parameters)
    space_group = (
        data.get("_symmetry_space_group_name_H-M")
        or data.get("_space_group_name_H-M_alt")
        or "P 1"
    )
    space_group = strip_quotes(space_group)

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
        label = strip_quotes(row[label_idx])
        element = normalize_element_symbol(
            strip_quotes(row[type_idx]) if type_idx is not None else guess_element_from_label(label)
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
            transformed = apply_symmetry_operation(operation, frac)
            canonical = canonicalize_fractional(transformed)
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
        from ase.io import read as ase_read  # imported here: ase.io is slow to import

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
                coordinates=frac_to_cart(site.frac, lattice),
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


def write_cif(structure: CrystalStructure, path: Path) -> None:
    """Write a P1-style CIF of the atoms as stored, with the space-group label."""

    path.parent.mkdir(parents=True, exist_ok=True)
    a, b, c, alpha, beta, gamma = structure.cell_parameters
    lines = [
        f"data_{_safe_data_name(structure.name)}",
        "_audit_creation_method 'CSPToolbox CrystalStructure'",
        f"_csptoolbox_explict_unit_cell {'true' if structure.explict_unit_cell else 'false'}",
        f"_cell_length_a {a:.10f}",
        f"_cell_length_b {b:.10f}",
        f"_cell_length_c {c:.10f}",
        f"_cell_angle_alpha {alpha:.10f}",
        f"_cell_angle_beta {beta:.10f}",
        f"_cell_angle_gamma {gamma:.10f}",
        f"_symmetry_space_group_name_H-M '{structure.space_group}'",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
    ]
    for atom in structure.atoms:
        frac = structure.fractional_coordinates(atom.coordinates)
        lines.append(
            f"{atom.label} {atom.element} "
            f"{frac[0]:.10f} {frac[1]:.10f} {frac[2]:.10f} 1.0"
        )
    path.write_text("\n".join(lines) + "\n")


def _strip_uncertainty(value: str) -> str:
    return re.sub(r"\([^)]+\)", "", value).strip()


def _parse_float(value: str) -> float:
    cleaned = _strip_uncertainty(strip_quotes(value))
    if "/" in cleaned and re.fullmatch(r"[+-]?\d+/\d+", cleaned):
        return float(Fraction(cleaned))
    return float(cleaned)


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
            and fractional_minimum_image_distance(site.frac, frac, lattice) <= site_merge_tolerance
        ):
            return True
    return False


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


def _safe_data_name(name: str) -> str:
    return re.sub(r"\s+", "_", name.strip()) or "CrystalStructure"


def _compare_ase_and_manual_expansions(
    manual_atoms: list[_ExpandedCifAtomSite],
    ase_atoms: Atoms,
) -> tuple[bool, str | None]:
    manual_counter = Counter(
        (atom.element, canonicalize_fractional(atom.frac))
        for atom in manual_atoms
    )
    ase_counter = Counter(
        (
            normalize_element_symbol(symbol),
            canonicalize_fractional(
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
