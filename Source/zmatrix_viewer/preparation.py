"""Prepare viewer scenes with explicit identity mapping and proper-rotation alignment."""
from __future__ import annotations

from dataclasses import replace
from hashlib import sha256
import json
import math
from pathlib import Path
import re

from .geometry import build_viewer_molecule, measure_angle_degrees, measure_dihedral_degrees
from .model import ViewerCoordinate, ViewerDocument, ViewerMolecule, ViewerScene, ZMatrixAtom, ZMatrixDocument
from .parser import load_zmatrix, parse_zmatrix_text
from .readers import StructureInput, molecule_from_coordinates, read_structure


def read_dofs(source) -> dict[str, dict]:
    """Read a JSON list/object or bnd/ang/dih definitions (e.g. CP2 input.in)."""
    if source is None:
        return {}
    if isinstance(source, dict):
        result = {str(k).lower(): (dict(v) if isinstance(v, dict) else {}) for k, v in source.items()}
    elif isinstance(source, (tuple, list, set, frozenset)):
        result = {str(name).lower(): {} for name in source}
    else:
        text = Path(source).read_text()
        if text.lstrip().startswith(('[', '{')):
            return read_dofs(json.loads(text))
        result = {}
        for line in text.splitlines():
            if line.lower().startswith('# independent_dofs:'):
                result.update(read_dofs(line.split(':', 1)[1].split()))
                continue
            fields = line.split()
            if fields and re.fullmatch(r'(bnd|ang|dih)\d+', fields[0], re.I):
                bounds = {}
                if len(fields) >= 3:
                    bounds = dict(zip(('lower', 'upper'), (float(v.replace('D', 'E').replace('d', 'e')) for v in fields[1:3])))
                result[fields[0].lower()] = bounds
    for name, bounds in result.items():
        if not re.fullmatch(r'(bnd|ang|dih)\d+', name):
            raise ValueError(f"Unknown internal-coordinate name {name!r}; use bndN, angN or dihN.")
        if set(bounds) - {'lower', 'upper'}:
            raise ValueError(f"Unknown bound fields for {name}; use lower and upper.")
        for key, bound in bounds.items():
            bounds[key] = float(bound)
            if not math.isfinite(bounds[key]):
                raise ValueError(f"Nonfinite bound for {name}.")
        if 'lower' in bounds and 'upper' in bounds and bounds['lower'] > bounds['upper']:
            raise ValueError(f"Reversed bounds for {name}; use an unwrapped periodic interval.")
    return result


def _load_topology(path):
    """Numeric # ZMAT v1, or CP2 label-only topology evaluated on supplied coordinates."""
    path = Path(path)
    text = path.read_text()
    lines = [line.split() for line in text.splitlines() if line.strip() and not line.lstrip().startswith('#')]
    labelled = bool(lines and (lines[0][0].lower().startswith('z-matrix') or
                    (len(lines[0]) == 1 and re.fullmatch(r'[A-Za-z]+\d+', lines[0][0]))))
    if not labelled:
        return load_zmatrix(path), False
    lines = [line for line in lines if not line[0].lower().startswith('z-matrix')]
    indices = {line[0]: i for i, line in enumerate(lines, 1)}
    if len(indices) != len(lines):
        raise ValueError("Label-only Z-matrix contains duplicate labels.")
    atoms = []
    for i, fields in enumerate(lines, 1):
        if len(fields) != min(i, 4):
            raise ValueError(f"Label-only Z-matrix row {i} needs {min(i,4)-1} references.")
        refs = []
        for label in fields[1:]:
            if label not in indices or indices[label] >= i:
                raise ValueError(f"Z-matrix row {i} reference {label} is not an earlier atom.")
            refs.append(indices[label])
        if len(refs) != len(set(refs)):
            raise ValueError(f"Z-matrix row {i} repeats a reference.")
        refs += [None] * (3-len(refs))
        element = re.match(r'[A-Za-z]+', fields[0])[0].capitalize()
        atoms.append(ZMatrixAtom(i, i, fields[0], element, refs[0], None, refs[1], None, refs[2], None))
    return ZMatrixDocument(path.stem, tuple(atoms), frozenset(), str(path)), True


def _load_input(path, **kwargs):
    path = Path(path)
    if path.suffix.lower() in ('.zmat', '.zmatrix'):
        if kwargs.get('frame', 0) != 0:
            raise ValueError("Numeric Z-matrix files contain one conformation; use frame 0.")
        return StructureInput(build_viewer_molecule(load_zmatrix(path),
                              infer_bonds=kwargs.get('infer_bonds', True), covalent_scale=kwargs.get('covalent_scale', 1.25)))
    return read_structure(path, **kwargs)


def _select_molecule(data, index, *, infer_bonds=True, covalent_scale=1.25):
    if data.crystal is None:
        if index is not None:
            raise ValueError("Molecule selection requires a periodic CIF/RES/PDB structure.")
        return data
    molecules = data.crystal.detect_molecules()
    if index is None:
        # Unwrap all supplied atoms, preserving their identity and original row order.
        by_label = {a.label: a for group in molecules for a in group}
        records = [by_label[a.label] for a in data.molecule.atoms]
    else:
        if not 1 <= index <= len(molecules):
            raise ValueError(f"Molecule {index} is unavailable; found {len(molecules)} components.")
        records = molecules[index-1]
    # Preserve supplied PDB connectivity; infer on the unwrapped coordinates otherwise.
    new_ids = {a.label: i for i, a in enumerate(records, 1)}
    explicit = [b for b in data.molecule.bonds if b.kind == 'explicit']
    pairs = [(new_ids[data.molecule.atoms[b.left-1].label], new_ids[data.molecule.atoms[b.right-1].label])
             for b in explicit if data.molecule.atoms[b.left-1].label in new_ids and data.molecule.atoms[b.right-1].label in new_ids]
    mol = molecule_from_coordinates([a.label for a in records], [a.element for a in records],
                                    [a.coordinates for a in records], title=data.molecule.title, source=data.molecule.source_name,
                                    bonds=pairs if explicit else None, infer_bonds=infer_bonds, covalent_scale=covalent_scale)
    return replace(data, molecule=mol)


def map_atoms(target, source: ViewerMolecule, atom_map=None, *, labels_explicit=True):
    """Return target-order source indices. Explicit maps run target label -> source label/index."""
    target = tuple(target)
    if atom_map is not None and not isinstance(atom_map, dict):
        atom_map = json.loads(Path(atom_map).read_text())
    by_label = {a.label: a.index for a in source.atoms}
    target_labels = {a.label for a in target}
    if atom_map is not None:
        if set(atom_map) != target_labels:
            raise ValueError("Atom map must contain exactly every reference/Z-matrix label.")
        indices = []
        for atom in target:
            value = atom_map[atom.label]
            if isinstance(value, bool):
                raise ValueError("Atom map indices are 1-based integers, not booleans.")
            index = value if isinstance(value, int) else by_label.get(value)
            if index is None or not 1 <= index <= len(source.atoms):
                raise ValueError(f"Unknown source atom {value!r} for {atom.label}.")
            indices.append(index)
        method = "explicit atom map"
    elif labels_explicit and target_labels <= set(by_label):
        indices = [by_label[a.label] for a in target]
        method = "exact atom labels"
    elif not labels_explicit and len(target) == len(source.atoms):
        indices = list(range(1, len(target)+1))
        method = "file atom order (unlabelled input)"
    else:
        raise ValueError("Atom identities differ. Supply an explicit atom map (target label -> source label or 1-based index).")
    if len(indices) != len(set(indices)):
        raise ValueError("Atom mapping must be one-to-one.")
    for atom, index in zip(target, indices):
        if atom.element != source.atoms[index-1].element:
            raise ValueError(f"Element mismatch for {atom.label}; provide the correct atom map.")
    return indices, method


def evaluate_topology(topology, source, indices, *, infer_bonds=True, covalent_scale=1.25):
    coordinates = tuple(source.atoms[i-1].coordinates for i in indices)
    atoms = []
    for atom in topology.atoms:
        i, j, k, l = atom.row_index, atom.bond_to, atom.angle_to, atom.dihedral_to
        distance = math.dist(coordinates[i-1], coordinates[j-1]) if j else None
        angle = measure_angle_degrees(coordinates, i, j, k) if k else None
        dihedral = measure_dihedral_degrees(coordinates, i, j, k, l) if l else None
        if any(value is not None and not math.isfinite(value) for value in (distance, angle, dihedral)):
            raise ValueError(f"Degenerate geometry at Z-matrix row {i}; its internal coordinates are undefined.")
        atoms.append(replace(atom, bond_length=distance, angle_degrees=angle, dihedral_degrees=dihedral))
    evaluated = replace(topology, atoms=tuple(atoms))
    # Parser checks the same numeric interchange contract used by GaussianInputBuilder.
    evaluated = parse_zmatrix_text(numeric_zmatrix(evaluated), source_name=source.source_name)
    built = build_viewer_molecule(evaluated, infer_bonds=infer_bonds, covalent_scale=covalent_scale)
    # Rendering uses the experimental coordinates, not a reconstructed or optimized conformation.
    built = replace(built, atoms=tuple(replace(atom, coordinates=point) for atom, point in zip(built.atoms, coordinates)))
    return evaluated, built


def numeric_zmatrix(document):
    lines = ['# ZMAT v1', '# title: ' + document.title, '# labels: ' + ' '.join(a.label for a in document.atoms)]
    if document.explicit_bonds:
        lines.append('# bonds: ' + ' '.join(f'{a}-{b}' for a, b in sorted(document.explicit_bonds)))
    for atom in document.atoms:
        fields = [atom.element]
        for ref, value in [(atom.bond_to, atom.bond_length), (atom.angle_to, atom.angle_degrees), (atom.dihedral_to, atom.dihedral_degrees)]:
            if ref is not None:
                if value is None:
                    raise ValueError("A label-only Z-matrix needs a coordinate file before it can be displayed.")
                fields.extend((str(ref), format(value, '.17g')))
        lines.append(' '.join(fields))
    return '\n'.join(lines) + '\n'


def align_comparison(reference, comparison):
    """Equal-weight all-atom Kabsch fit; determinant +1 prevents reflection."""
    import numpy as np
    if len(reference.atoms) != len(comparison.atoms):
        raise ValueError("Overlay requires equal atom counts after mapping.")
    ref = np.array([a.coordinates for a in reference.atoms])
    comp = np.array([a.coordinates for a in comparison.atoms])
    r0, c0 = ref.mean(0), comp.mean(0)
    u, _, vt = np.linalg.svd((comp-c0).T @ (ref-r0))
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    fitted = (comp-c0) @ (u @ correction @ vt) + r0
    rmsd = float(np.sqrt(np.mean(np.sum((fitted-ref)**2, axis=1))))
    return replace(comparison, atoms=tuple(replace(a, coordinates=tuple(map(float, p))) for a,p in zip(comparison.atoms, fitted))), rmsd


def _coordinates(document, dofs, reference=None):
    result = []
    for atom in document.atoms:
        ids = (atom.row_index, atom.bond_to, atom.angle_to, atom.dihedral_to)
        for kind, length, value, unit in [('bnd',2,atom.bond_length,'Å'),('ang',3,atom.angle_degrees,'°'),('dih',4,atom.dihedral_degrees,'°')]:
            if value is None:
                continue
            name = f'{kind}{atom.row_index}'
            ref = None
            if reference:
                a = reference.atoms[atom.row_index-1]
                ref = {'bnd':a.bond_length,'ang':a.angle_degrees,'dih':a.dihedral_degrees}[kind]
            bounds = dofs.get(name, {})
            result.append(ViewerCoordinate(name, kind, ids[:length], tuple(document.atoms[i-1].label for i in ids[:length]),
                                          value, unit, name in dofs, ref, bounds.get('lower'), bounds.get('upper')))
    unknown = set(dofs) - {c.name for c in result}
    if unknown:
        raise ValueError(f"Independent DoFs absent from this topology: {', '.join(sorted(unknown))}.")
    return tuple(result)


def build_structure_view(path, *, reference=None, zmatrix=None, independent_dofs=None, dofs_file=None,
                         atom_map=None, crystal=False, frame=0, reference_frame=0, molecule_index=None,
                         reference_molecule_index=None, infer_bonds=True, covalent_scale=1.25, title=None):
    """Prepare any of the five molecular/crystal viewer workflows without changing input files."""
    # Bind provenance to the files now, even if callers change directory before export.
    path = Path(path).resolve()
    reference = Path(reference).resolve() if reference is not None else None
    zmatrix = Path(zmatrix).resolve() if zmatrix is not None else None
    dofs_file = Path(dofs_file).resolve() if dofs_file is not None else None
    if atom_map is not None and not isinstance(atom_map, dict):
        atom_map = Path(atom_map).resolve()
    sources = [path]
    dofs = read_dofs(independent_dofs)
    if dofs_file:
        dofs.update(read_dofs(dofs_file)); sources.append(Path(dofs_file))
    if atom_map is not None and not isinstance(atom_map, dict):
        sources.append(Path(atom_map))
    topology = None
    label_only = False
    if zmatrix:
        sources.append(Path(zmatrix))
        topology, label_only = _load_topology(zmatrix)
    elif path.suffix.lower() in ('.zmat', '.zmatrix'):
        topology, label_only = _load_topology(path)
    if label_only and not zmatrix:
        raise ValueError("For a label-only Z-matrix, pass the coordinate file as input and the topology with --zmatrix.")
    data = _load_input(path, frame=frame, infer_bonds=infer_bonds, covalent_scale=covalent_scale)
    notices = list(data.notices)
    full_cell = None
    cell = ()
    if crystal:
        if reference or zmatrix or molecule_index is not None or dofs or atom_map is not None:
            raise ValueError("Crystal cell view is separate from a molecular overlay/mapping; select a molecule for those workflows.")
        if data.crystal is None:
            raise ValueError("Crystal view needs cell and symmetry information (CIF, RES or periodic PDB).")
        structure = data.crystal
        if structure.explict_unit_cell:
            full_structure = structure
            structure = structure.reduce_to_asymmetric_unit()
            notices.append("Asymmetric unit inferred from the explicit cell using spglib (symprec 0.05 Å).")
        elif path.suffix.lower() == '.cif':
            from ..crystal_structure import CrystalStructure
            full_structure = CrystalStructure.expand_cif_to_unit_cell(path)
        else:
            full_structure = structure.expand_to_explicit_unit_cell()
        def cell_molecule(s, name):
            groups = s.detect_molecules()
            atoms = [atom for group in groups for atom in group]
            # Expanded sites may retain identical labels; qualify by molecular copy.
            from collections import Counter
            counts = Counter(a.label for a in atoms)
            labels = [f'{a.label}@{i}' if counts[a.label] > 1 else a.label
                      for i, group in enumerate(groups,1) for a in group]
            return molecule_from_coordinates(labels, [a.element for a in atoms], [a.coordinates for a in atoms],
                                             title=name, source=str(path), infer_bonds=infer_bonds, covalent_scale=covalent_scale)
        molecule = cell_molecule(structure, 'Asymmetric unit')
        full_cell = cell_molecule(full_structure, 'Full cell')
        cell = tuple(tuple(map(float, vector)) for vector in structure.cell.array)
        notices.append("Whole molecules are unwrapped across cell boundaries; some atoms may lie outside the cell outline.")
        return ViewerScene(title or path.stem, molecule, full_cell=full_cell, cell=cell, notices=tuple(notices),
                           source_paths=tuple(map(str,sources)), source_hashes=_hashes(sources))
    data = _select_molecule(data, molecule_index, infer_bonds=infer_bonds, covalent_scale=covalent_scale)
    molecule = data.molecule
    ref_molecule = None
    ref_document = None
    rmsd = None
    mapping = ()
    method = ""
    mapped_text = None
    evaluated = topology
    if topology:
        if independent_dofs is None and dofs_file is None:
            for line in Path(zmatrix or path).read_text().splitlines():
                if line.lower().startswith('# independent_dofs:'):
                    dofs.update(read_dofs(line.split(':',1)[1].split()))
        indices, method = map_atoms(topology.atoms, molecule, atom_map, labels_explicit=data.labels_explicit)
        mapping = tuple((a.label, molecule.atoms[i-1].label) for a,i in zip(topology.atoms, indices))
        evaluated, molecule = evaluate_topology(topology, molecule, indices, infer_bonds=infer_bonds, covalent_scale=covalent_scale)
        mapped_text = numeric_zmatrix(evaluated)
        if zmatrix and not label_only:
            ref_document = topology
            ref_molecule = build_viewer_molecule(topology, infer_bonds=infer_bonds, covalent_scale=covalent_scale)
    if reference:
        sources.append(Path(reference))
        ref_data = _select_molecule(_load_input(reference, frame=reference_frame, infer_bonds=infer_bonds,
                                              covalent_scale=covalent_scale), reference_molecule_index,
                                    infer_bonds=infer_bonds, covalent_scale=covalent_scale)
        notices.extend(ref_data.notices)
        ref_molecule = ref_data.molecule
        if topology:
            ids, _ = map_atoms(topology.atoms, ref_molecule, labels_explicit=ref_data.labels_explicit)
            ref_document, ref_molecule = evaluate_topology(topology, ref_molecule, ids, infer_bonds=infer_bonds, covalent_scale=covalent_scale)
        else:
            indices, method = map_atoms(ref_molecule.atoms, molecule, atom_map, labels_explicit=(data.labels_explicit and ref_data.labels_explicit))
            mapping = tuple((a.label, molecule.atoms[i-1].label) for a,i in zip(ref_molecule.atoms,indices))
            reverse = {old: new for new,old in enumerate(indices,1)}
            molecule = replace(molecule, atoms=tuple(replace(molecule.atoms[i-1], index=j, label=a.label)
                                for j,(a,i) in enumerate(zip(ref_molecule.atoms,indices),1)),
                                bonds=tuple(replace(b,left=reverse[b.left],right=reverse[b.right]) for b in molecule.bonds
                                            if b.left in reverse and b.right in reverse))
    if ref_molecule:
        molecule, rmsd = align_comparison(ref_molecule, molecule)
        notices.append("Overlay aligned by equal-weight all-atom translation and proper rotation; reflection is disabled.")
    if 'file atom order' in method:
        notices.append("This input has no persistent atom labels. Correspondence assumes unchanged atom order; provide an atom map if reordered.")
    coordinates = _coordinates(evaluated, dofs, ref_document) if evaluated else ()
    if dofs and not evaluated:
        raise ValueError("Independent DoFs require a numeric or label-only Z-matrix definition.")
    return ViewerScene(title or path.stem, molecule, ref_molecule, coordinates=coordinates, mapping=mapping,
                       rmsd=rmsd, mapping_method=method, notices=tuple(notices), source_paths=tuple(map(str,sources)),
                       source_hashes=_hashes(sources), mapped_zmatrix=mapped_text)


def _hashes(paths):
    return tuple((str(p), sha256(Path(p).read_bytes()).hexdigest()) for p in dict.fromkeys(paths))


def build_comparison_view(reference, comparison, **kwargs):
    return build_structure_view(comparison, reference=reference, **kwargs)


def build_mapping_view(experimental, zmatrix, **kwargs):
    return build_structure_view(experimental, zmatrix=zmatrix, **kwargs)


def build_crystal_view(path, **kwargs):
    return build_structure_view(path, crystal=True, **kwargs)


def build_gallery(scenes, *, title="CSPToolbox molecular viewer"):
    scenes = tuple(scenes)
    if not scenes:
        raise ValueError("A viewer document needs at least one scene.")
    return ViewerDocument(scenes, title)
