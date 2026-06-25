"""Interactive static viewer utilities for CSPToolbox Z-matrix files."""

from __future__ import annotations

from pathlib import Path

from .geometry import build_viewer_molecule, reconstruct_coordinates
from .html_export import render_viewer_html, write_viewer_html
from .model import (
    ViewerAtom,
    ViewerBond,
    ViewerDihedral,
    ViewerMolecule,
    ZMatrixAtom,
    ZMatrixDocument,
)
from .parser import load_zmatrix, parse_zmatrix_text

__all__ = [
    "ViewerAtom",
    "ViewerBond",
    "ViewerDihedral",
    "ViewerMolecule",
    "ZMatrixAtom",
    "ZMatrixDocument",
    "build_viewer_document",
    "build_zmatrix_viewer_document",
    "build_viewer_molecule",
    "load_zmatrix",
    "parse_zmatrix_text",
    "reconstruct_coordinates",
    "render_viewer_html",
    "render_zmatrix_viewer_html",
    "write_viewer_html",
    "write_zmatrix_viewer_html",
]


def build_viewer_document(
    zmat_path: str | Path,
    *,
    infer_bonds: bool = True,
    covalent_scale: float = 1.25,
) -> ViewerMolecule:
    """Load a Z-matrix file and return the complete viewer payload."""

    return build_viewer_molecule(
        load_zmatrix(zmat_path),
        infer_bonds=infer_bonds,
        covalent_scale=covalent_scale,
    )


build_zmatrix_viewer_document = build_viewer_document
render_zmatrix_viewer_html = render_viewer_html
write_zmatrix_viewer_html = write_viewer_html
