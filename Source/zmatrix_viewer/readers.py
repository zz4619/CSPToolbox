"""File adapters for the viewer. Geometry stays in the input frame and in angstrom."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
from pathlib import Path

from .geometry import COVALENT_RADII, DISPLAY_RADII, ELEMENT_COLORS
from .model import ViewerAtom, ViewerBond, ViewerMolecule


@dataclass(frozen=True)
class StructureInput:
    molecule: ViewerMolecule
    labels_explicit: bool = True
    crystal: object | None = None
    notices: tuple[str, ...] = ()


def molecule_from_coordinates(labels, elements, coordinates, *, title="Molecule", source=None,
                              bonds=None, infer_bonds=True, covalent_scale=1.25):
    labels, elements = list(labels), list(elements)
    coordinates = [tuple(map(float, point)) for point in coordinates]
    if not labels or not (len(labels) == len(elements) == len(coordinates)):
        raise ValueError("Atom labels, elements and coordinates must have the same nonzero length.")
    if len(set(labels)) != len(labels):
        raise ValueError("Atom labels are not unique; use qualified labels for repeated molecules.")
    if any(len(point) != 3 or not all(math.isfinite(x) for x in point) for point in coordinates):
        raise ValueError("Coordinates must be finite Cartesian triples in angstrom.")
    if not math.isfinite(covalent_scale) or covalent_scale <= 0:
        raise ValueError("covalent_scale must be finite and positive.")
    atoms = tuple(ViewerAtom(i, label, element, point, ELEMENT_COLORS.get(element, "#7b8794"),
                             DISPLAY_RADII.get(element, 0.2))
                  for i, (label, element, point) in enumerate(zip(labels, elements, coordinates), 1))
    pairs = set()
    if bonds is not None:
        for left, right in bonds:
            if not (1 <= left <= len(atoms) and 1 <= right <= len(atoms)) or left == right:
                raise ValueError(f"Invalid bond endpoint: {left}-{right}.")
            pairs.add(tuple(sorted((left, right))))
    elif infer_bonds:
        candidates = []
        for i, left in enumerate(atoms):
            for right in atoms[i + 1:]:
                distance = math.dist(left.coordinates, right.coordinates)
                radius = COVALENT_RADII.get(left.element, 0.77) + COVALENT_RADII.get(right.element, 0.77)
                if 0.1 < distance <= covalent_scale * radius:
                    candidates.append((distance, left.index, right.index))
        # A close hydrogen bond must not become a second covalent bond to H.
        occupied_h = set()
        for _, left, right in sorted(candidates):
            hydrogens = {i for i in (left, right) if atoms[i - 1].element == "H"}
            if not hydrogens & occupied_h:
                pairs.add((left, right))
                occupied_h.update(hydrogens)
    kind = "explicit" if bonds is not None else "inferred"
    return ViewerMolecule(title, atoms, tuple(ViewerBond(a, b, kind) for a, b in sorted(pairs)), (), (), source)


def _auto_labels(elements):
    counts = Counter()
    result = []
    for element in elements:
        counts[element] += 1
        result.append(f"{element}{counts[element]}")
    return result


def _pdb(path, frame, infer_bonds, covalent_scale):
    lines = path.read_text().splitlines()
    model_lines = []
    current = []
    has_models = any(line.startswith("MODEL ") for line in lines)
    for line in lines:
        if line.startswith("MODEL "):
            current = []
        elif line.startswith("ENDMDL"):
            model_lines.append(current)
            current = []
        elif line.startswith(("ATOM  ", "HETATM")):
            current.append(line)
    if current:
        model_lines.append(current)
    if not has_models:
        model_lines = [[line for line in lines if line.startswith(("ATOM  ", "HETATM"))]]
    try:
        atom_lines = model_lines[frame]
    except IndexError as error:
        raise ValueError(f"PDB frame {frame} is unavailable.") from error
    if any(line[16:17].strip() not in ("", "A") for line in atom_lines):
        raise ValueError("PDB contains alternate locations; select one conformation before viewing.")
    labels = [line[12:16].strip() for line in atom_lines]
    if any(not label for label in labels):
        raise ValueError("PDB atom names must be nonempty.")
    if len(set(labels)) != len(labels):
        labels = [f"{line[12:16].strip()}@{line[21:22].strip() or '_'}:{line[22:27].strip()}"
                  for line in atom_lines]
        if len(set(labels)) != len(labels):
            # CP2 full-cell exports reuse both atom names and residue IDs.
            labels = [f'{label}:{i}' for i, label in enumerate(labels, 1)]
    from ase.data import atomic_numbers
    elements = []
    for line in atom_lines:
        element = line[76:78].strip().capitalize()
        if not element:
            import re
            letters = re.sub('[^A-Za-z]', '', line[12:16].strip())
            # PDB atom names beginning with a space encode one-letter elements.
            element = letters[0].upper() if line[12:13] == ' ' or line[12:13].isdigit() else letters[:2].capitalize()
            if element not in atomic_numbers:
                element = letters[0].upper()
        if element not in atomic_numbers:
            raise ValueError(f"Unknown PDB element {element!r}.")
        elements.append(element)
    xyz = [(float(line[30:38]), float(line[38:46]), float(line[46:54])) for line in atom_lines]
    serial_to_index = {int(line[6:11]): i for i, line in enumerate(atom_lines, 1)}
    bonds = set()
    for line in lines:
        if line.startswith('CONECT'):
            numbers = [int(line[i:i+5]) for i in range(6, len(line), 5) if line[i:i+5].strip()]
            for other in numbers[1:]:
                if numbers[0] in serial_to_index and other in serial_to_index:
                    a, b = serial_to_index[numbers[0]], serial_to_index[other]
                    if a != b:
                        bonds.add(tuple(sorted((a, b))))
    mol = molecule_from_coordinates(labels, elements, xyz, title=path.stem, source=str(path),
                                    bonds=bonds or None, infer_bonds=infer_bonds, covalent_scale=covalent_scale)
    crystal = None
    if any(line.startswith('CRYST1') for line in lines):
        from ..crystal_structure import CrystalStructure, AtomRecord
        crystal = CrystalStructure.from_file(path)
        crystal.atoms = [AtomRecord(a.label, a.element, a.coordinates) for a in mol.atoms]
    return StructureInput(mol, crystal=crystal)


def read_structure(path: str | Path, *, frame=0, infer_bonds=True, covalent_scale=1.25) -> StructureInput:
    path = Path(path).resolve()
    if frame < 0:
        raise ValueError("Frame indices must be nonnegative (zero-based).")
    suffix = path.suffix.lower()
    if suffix in ('.pdb', '.ent'):
        return _pdb(path, frame, infer_bonds, covalent_scale)
    if suffix in ('.cif', '.res'):
        if frame:
            raise ValueError("CIF/RES input currently supports its first structure only.")
        from ..crystal_structure import CrystalStructure
        structure = CrystalStructure.from_file(path)
        notices = []
        if suffix == '.cif':
            report = CrystalStructure.inspect_cif_unit_cell_expansion(path)
            if report.partial_occupancies:
                notices.append("Partial occupancies are present; sites are displayed without an occupancy model.")
        mol = molecule_from_coordinates([a.label for a in structure.atoms], [a.element for a in structure.atoms],
                                       [a.coordinates for a in structure.atoms], title=path.stem, source=str(path),
                                       infer_bonds=infer_bonds, covalent_scale=covalent_scale)
        return StructureInput(mol, crystal=structure, notices=tuple(notices))
    if suffix in ('.mol', '.sdf'):
        from rdkit import Chem
        if suffix == '.mol':
            if frame:
                raise ValueError("MOL contains one molecule; use frame 0.")
            molecule = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=False)
        else:
            supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
            if frame >= len(supplier):
                raise ValueError(f"SDF frame {frame} is unavailable.")
            molecule = supplier[frame]
        if molecule is None or not molecule.GetNumConformers():
            raise ValueError(f"Could not read molecular coordinates from {path}.")
        symbols = [a.GetSymbol() for a in molecule.GetAtoms()]
        mol = molecule_from_coordinates(_auto_labels(symbols), symbols, molecule.GetConformer().GetPositions(),
                                       title=path.stem, source=str(path),
                                       bonds=[(b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1) for b in molecule.GetBonds()])
        return StructureInput(mol, labels_explicit=False)
    # ASE supplies XYZ, extended XYZ, POSCAR/CONTCAR and other coordinate formats.
    from ase.io import read
    try:
        atoms = read(path, index=frame)
    except (IndexError, StopIteration) as error:
        raise ValueError(f"Frame {frame} is unavailable in {path}.") from error
    elements = atoms.get_chemical_symbols()
    labels = atoms.arrays.get('labels')
    if labels is None:
        labels = atoms.arrays.get('atom_labels')
    explicit = labels is not None
    labels = list(map(str, labels)) if explicit else _auto_labels(elements)
    mol = molecule_from_coordinates(labels, elements, atoms.positions, title=path.stem, source=str(path),
                                    infer_bonds=infer_bonds, covalent_scale=covalent_scale)
    crystal = None
    if atoms.cell.rank == 3 and any(atoms.pbc):
        from ..crystal_structure import CrystalStructure, AtomRecord
        crystal = CrystalStructure([AtomRecord(a.label, a.element, a.coordinates) for a in mol.atoms],
                                   tuple(map(float, atoms.cell.cellpar())),
                                   lattice_matrix=tuple(map(tuple, atoms.cell.array)), explict_unit_cell=True,
                                   name=path.stem)
    return StructureInput(mol, explicit, crystal)
