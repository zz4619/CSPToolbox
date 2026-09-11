"""One self-contained renderer for standalone documents and embedded fragments."""
from __future__ import annotations
from dataclasses import asdict
from html import escape
from importlib.resources import files
import json
import math
from pathlib import Path
import re
from uuid import uuid4
from .model import ViewerCoordinate, ViewerDocument, ViewerMolecule, ViewerScene


def molecule_payload(molecule: ViewerMolecule) -> dict:
    """Public full-precision payload: indices are 1-based, positions in angstrom."""
    ids = {a.index for a in molecule.atoms}
    if not ids or len(ids) != len(molecule.atoms):
        raise ValueError("Viewer atoms need nonempty, unique indices.")
    for atom in molecule.atoms:
        if len(atom.coordinates) != 3 or not all(math.isfinite(v) for v in atom.coordinates):
            raise ValueError("Viewer coordinates must be finite Cartesian triples.")
    for bond in molecule.bonds:
        if bond.left not in ids or bond.right not in ids or bond.left == bond.right:
            raise ValueError(f"Invalid viewer bond {bond.left}-{bond.right}.")
    for coordinate in molecule.dihedrals:
        if any(i not in ids for i in coordinate.atom_indices):
            raise ValueError("A dihedral refers to a nonexistent atom.")
    return {
        'title': molecule.title, 'source_name': molecule.source_name,
        'atoms': [dict(index=a.index, label=a.label, element=a.element, x=a.coordinates[0], y=a.coordinates[1],
                       z=a.coordinates[2], color=a.color, display_radius=a.display_radius) for a in molecule.atoms],
        'bonds': [asdict(b) for b in molecule.bonds], 'dihedrals': [asdict(d) for d in molecule.dihedrals],
        'warnings': list(molecule.warnings),
    }


_molecule_payload = molecule_payload  # Compatibility for existing private adapters.


def _document(value):
    if isinstance(value, ViewerDocument):
        document = value
    elif isinstance(value, ViewerScene):
        document = ViewerDocument((value,), value.title)
    elif isinstance(value, ViewerMolecule):
        coordinates = tuple(ViewerCoordinate(f'dih{d.row_index}', 'dih', d.atom_indices, d.atom_labels,
                                            d.value_degrees, '°') for d in value.dihedrals)
        document = ViewerDocument((ViewerScene(value.title, value, coordinates=coordinates),), value.title)
    else:
        raise TypeError("Expected ViewerMolecule, ViewerScene or ViewerDocument.")
    if not document.scenes:
        raise ValueError("A viewer needs at least one scene.")
    return document


def viewer_payload(value) -> dict:
    """Versioned multi-scene payload shared by all export modes."""
    document = _document(value)
    scenes = []
    for scene in document.scenes:
        data = asdict(scene)
        for key in ('molecule', 'reference', 'full_cell'):
            molecule = getattr(scene, key)
            data[key] = molecule_payload(molecule) if molecule else None
        ids = {a.index for a in scene.molecule.atoms}
        for c in scene.coordinates:
            if any(i not in ids for i in c.atom_indices) or not math.isfinite(c.value):
                raise ValueError(f"Invalid internal coordinate {c.name}.")
        if scene.cell and (len(scene.cell) != 3 or any(len(v) != 3 or not all(math.isfinite(x) for x in v) for v in scene.cell)):
            raise ValueError("Cell must be a finite 3×3 lattice matrix.")
        scenes.append(data)
    payload = dict(schema_version=document.schema_version, title=document.title, scenes=scenes)
    json.dumps(payload, allow_nan=False)
    return payload


def render_viewer_fragment(value, *, root_id=None) -> str:
    """Self-initializing fragment, requiring no host scripts, CSS, CDN or server."""
    root_id = root_id or 'csp-viewer-' + uuid4().hex
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]*', root_id):
        raise ValueError("root_id must be a safe HTML identifier.")
    assets = files(__package__).joinpath('assets')
    payload = json.dumps(viewer_payload(value), ensure_ascii=True, allow_nan=False).replace('<', '\\u003c')
    replacements = {'__ROOT_ID__': root_id, '__VIEWER_CSS__': assets.joinpath('viewer.css').read_text(encoding='utf-8'),
                    '__VIEWER_JS__': assets.joinpath('viewer.js').read_text(encoding='utf-8'), '__PAYLOAD__': payload}
    return re.sub('|'.join(map(re.escape, replacements)), lambda m: replacements[m[0]],
                  assets.joinpath('viewer.html').read_text(encoding='utf-8'))


def render_viewer_html(value) -> str:
    document = _document(value)
    return ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width,initial-scale=1">\n'
            '<title>' + escape(document.title) + '</title>\n'
            '<style>html{color-scheme:light dark}body{margin:0;padding:16px;background:light-dark(#ffffff,#171a20)}</style>\n'
            '</head>\n<body>\n' + render_viewer_fragment(document) + '\n</body>\n</html>\n')


def _check_destination(value, destination):
    path = Path(destination)
    sources = set()
    for scene in _document(value).scenes:
        sources.update(scene.source_paths)
        for molecule in (scene.molecule, scene.reference, scene.full_cell):
            if molecule and molecule.source_name:
                sources.add(molecule.source_name)
    for source in sources:
        other = Path(source)
        if path.resolve() == other.resolve() or (path.exists() and other.exists() and path.samefile(other)):
            raise ValueError("Viewer output must not overwrite an input structure or metadata file.")
    return path


def write_viewer_html(value, destination) -> Path:
    path = _check_destination(value, destination)
    rendered = render_viewer_html(value)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rendered, encoding='utf-8')
    return path


def write_viewer_fragment(value, destination, *, root_id=None) -> Path:
    path = _check_destination(value, destination)
    rendered = render_viewer_fragment(value, root_id=root_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rendered, encoding='utf-8')
    return path
