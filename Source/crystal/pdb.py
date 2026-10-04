"""PDB reading and writing (CRYST1 cell, ATOM/HETATM records)."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from .elements import guess_element_from_label, normalize_element_symbol
from .fileio import parse_bool_tag
from .records import AtomRecord

if TYPE_CHECKING:
    from .structure import CrystalStructure


def read_pdb(cls: type[CrystalStructure], path: Path) -> "CrystalStructure":
    """Read CRYST1 and ATOM/HETATM records (Cartesian Å)."""

    lines = path.read_text().splitlines()
    cryst1 = next((line for line in lines if line.startswith("CRYST1")), None)
    if cryst1 is None:
        raise ValueError(f"No CRYST1 record found in {path}")

    cell_parameters = (
        float(cryst1[6:15].strip()),
        float(cryst1[15:24].strip()),
        float(cryst1[24:33].strip()),
        float(cryst1[33:40].strip()),
        float(cryst1[40:47].strip()),
        float(cryst1[47:54].strip()),
    )
    space_group = cryst1[55:66].strip() or "P 1"

    atoms: list[AtomRecord] = []
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        label = line[12:16].strip() or f"ATOM{len(atoms) + 1}"
        element = normalize_element_symbol(line[76:78].strip() or guess_element_from_label(label))
        coordinates = (
            float(line[30:38].strip()),
            float(line[38:46].strip()),
            float(line[46:54].strip()),
        )
        atoms.append(AtomRecord(label=label, element=element, coordinates=coordinates))

    return cls(
        atoms=atoms,
        cell_parameters=cell_parameters,
        space_group=space_group,
        name=path.stem or "CrystalStructure",
        explict_unit_cell=_parse_pdb_explict_unit_cell(lines),
    )


def write_pdb(structure: CrystalStructure, path: Path) -> None:
    """Write CRYST1 and HETATM records (labels truncated to 4 characters)."""

    path.parent.mkdir(parents=True, exist_ok=True)
    a, b, c, alpha, beta, gamma = structure.cell_parameters
    lines = [
        f"HEADER    {structure.name}",
        f"REMARK   1 CSPTOOLBOX_EXPLICT_UNIT_CELL {'TRUE' if structure.explict_unit_cell else 'FALSE'}",
        f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} "
        f"{structure.space_group:<11}{1:>4d}",
    ]
    for index, atom in enumerate(structure.atoms, start=1):
        x, y, z = atom.coordinates
        atom_name = atom.label[:4]
        lines.append(
            f"HETATM{index:5d} {atom_name:<4s} MOL A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          "
            f"{atom.element:>2s}"
        )
    lines.extend(["END", ""])
    path.write_text("\n".join(lines))


def _parse_pdb_explict_unit_cell(lines: list[str]) -> bool:
    prefix = "REMARK   1 CSPTOOLBOX_EXPLICT_UNIT_CELL "
    for line in lines:
        if line.startswith(prefix):
            return parse_bool_tag(line[len(prefix):], default=False)
    return False
