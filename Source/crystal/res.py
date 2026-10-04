"""SHELX .res reading and writing; LATT/SYMM records are preserved."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.cell import Cell

from .elements import element_order, normalize_element_symbol
from .fileio import parse_bool_tag
from .records import AtomRecord, CellParameters
from .symmetry import deduplicate_shelx_symmetry_operations, shelx_symmetry_records

if TYPE_CHECKING:
    from .structure import CrystalStructure


def read_res(cls: type[CrystalStructure], path: Path) -> "CrystalStructure":
    """Read a SHELX .res file, keeping LATT/SYMM records for later expansion."""

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
            explict_unit_cell = parse_bool_tag(
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
                element=normalize_element_symbol(element),
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


def write_res(structure: CrystalStructure, path: Path, *, rounding: bool = False) -> None:
    """Write a SHELX .res file with LATT/SYMM records and fractional coordinates."""

    path.parent.mkdir(parents=True, exist_ok=True)
    ordered_elements = element_order(atom.element for atom in structure.atoms)
    if structure.shelx_latt_value is not None:
        latt_value = structure.shelx_latt_value
        symmetry_operations = deduplicate_shelx_symmetry_operations(
            structure.symmetry_operations,
            latt_value=latt_value,
            rounding=rounding,
        )
    else:
        latt_value, symmetry_operations = shelx_symmetry_records(
            structure.space_group,
            hall_number=structure.hall_number,
            rounding=rounding,
        )
    lines = [
        f"TITL {structure.name}",
        f"REM SPACE_GROUP {structure.space_group}",
        f"REM EXPLICT_UNIT_CELL {'TRUE' if structure.explict_unit_cell else 'FALSE'}",
        "CELL 1.54184 "
        f"{structure.cell_parameters[0]:.10f} {structure.cell_parameters[1]:.10f} {structure.cell_parameters[2]:.10f} "
        f"{structure.cell_parameters[3]:.10f} {structure.cell_parameters[4]:.10f} {structure.cell_parameters[5]:.10f}",
        f"LATT {latt_value}",
    ]
    lines.extend(f"SYMM {operation}" for operation in symmetry_operations)
    lines.append("SFAC " + " ".join(ordered_elements))
    for atom in structure.atoms:
        frac = structure.fractional_coordinates(atom.coordinates)
        sfac_index = ordered_elements.index(atom.element) + 1
        lines.append(
            f"{atom.label:<8s} {sfac_index:2d} "
            f"{frac[0]:.10f} {frac[1]:.10f} {frac[2]:.10f} 11.00000 0.05000"
        )
    lines.extend(["END", ""])
    path.write_text("\n".join(lines))
