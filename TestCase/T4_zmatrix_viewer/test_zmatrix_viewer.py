"""Tests for the static Z-matrix viewer."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.zmatrix_viewer import build_viewer_molecule, load_zmatrix, parse_zmatrix_text  # noqa: E402
from Source.zmatrix_viewer.cli import main as viewer_main  # noqa: E402
from Source.zmatrix_viewer.geometry import (  # noqa: E402
    measure_angle_degrees,
    measure_dihedral_degrees,
    reconstruct_coordinates,
)
from Source.zmatrix_viewer.html_export import render_viewer_html, write_viewer_html  # noqa: E402


BRANCH_IMPROPER_ROOT = (
    PROJECT_ROOT
    / "TestCase"
    / "T3_zmatrix_generation"
    / "BranchImproper_API_Zmatrices"
)


BUTANE_ZMAT = """# ZMAT v1
# title: butane_chain
C
C 1 1.5357408636
C 2 1.5400000000 1 111.3857279319
C 1 1.6527552753 2 109.8627296098 3 157.2945906352
"""


METHANOL_ZMAT = """# ZMAT v1
# title: methanol
C
O 1 1.4300000000
H 1 1.0547511555 2 121.4295656148
H 1 1.0547511555 3 117.1408687703 2 180.0000000000
H 1 1.0900000000 3 90.0000000000 4 -90.0000000000
H 2 0.8850988645 1 122.0740008753 5 -90.0000000000
"""


LABELLED_ZMAT = """# ZMAT v1
# title: labelled butane
C1 C
C2 C 1 1.5400000000
C3 C 2 1.5400000000 1 111.0000000000
C4 C 3 1.5400000000 2 111.0000000000 1 60.0000000000
"""


class ZMatrixViewerTests(unittest.TestCase):
    def test_current_numeric_zmat_autogenerates_labels(self) -> None:
        document = parse_zmatrix_text(METHANOL_ZMAT, source_name="methanol.zmat")

        self.assertEqual("methanol", document.title)
        self.assertEqual(
            ["C1", "O1", "H1", "H2", "H3", "H4"],
            [atom.label for atom in document.atoms],
        )
        self.assertEqual(3, sum(atom.has_dihedral for atom in document.atoms))

    def test_labelled_rows_preserve_atom_labels(self) -> None:
        document = parse_zmatrix_text(LABELLED_ZMAT)

        self.assertEqual(["C1", "C2", "C3", "C4"], [atom.label for atom in document.atoms])
        self.assertEqual("C4", document.atoms[-1].label)
        self.assertAlmostEqual(60.0, document.atoms[-1].dihedral_degrees)

    def test_reconstruction_preserves_internal_coordinates(self) -> None:
        document = parse_zmatrix_text(BUTANE_ZMAT)
        coords = reconstruct_coordinates(document.atoms)

        for atom in document.atoms:
            if atom.bond_to is not None:
                self.assertAlmostEqual(
                    atom.bond_length,
                    _distance(coords, atom.row_index, atom.bond_to),
                    places=8,
                )
            if atom.angle_to is not None:
                self.assertAlmostEqual(
                    atom.angle_degrees,
                    measure_angle_degrees(coords, atom.row_index, atom.bond_to, atom.angle_to),
                    places=8,
                )
            if atom.dihedral_to is not None:
                measured = measure_dihedral_degrees(
                    coords,
                    atom.row_index,
                    atom.bond_to,
                    atom.angle_to,
                    atom.dihedral_to,
                )
                self.assertAlmostEqual(
                    0.0,
                    _angle_delta(measured, atom.dihedral_degrees),
                    places=8,
                )

    def test_viewer_payload_and_html_include_dihedral_rows(self) -> None:
        molecule = build_viewer_molecule(parse_zmatrix_text(LABELLED_ZMAT))
        html = render_viewer_html(molecule)

        self.assertEqual(4, len(molecule.atoms))
        self.assertEqual(1, len(molecule.dihedrals))
        self.assertEqual(("C4", "C3", "C2", "C1"), molecule.dihedrals[0].atom_labels)
        self.assertIn('id="moleculeCanvas"', html)
        self.assertIn("window.__viewerReady = true", html)
        self.assertIn("labelled butane", html)
        self.assertIn("is-atom-hover", html)
        self.assertIn("function hitTestAtom", html)
        self.assertIn("function atomHoverDihedrals", html)

    def test_cli_writes_standalone_html(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            zmat_path = tmp_path / "methanol.zmat"
            html_path = tmp_path / "methanol_viewer.html"
            zmat_path.write_text(METHANOL_ZMAT, encoding="utf-8")

            exit_code = viewer_main([str(zmat_path), "--output", str(html_path)])

            self.assertEqual(0, exit_code)
            self.assertTrue(html_path.is_file())
            self.assertIn("Dihedral Angles", html_path.read_text(encoding="utf-8"))

    def test_branch_improper_api_zmatrices_are_viewer_compatible(self) -> None:
        paths = sorted(BRANCH_IMPROPER_ROOT.glob("*.zmat"))
        if not paths:
            self.skipTest(f"No branch-improper Z-matrix files found under {BRANCH_IMPROPER_ROOT}")

        for path in paths:
            with self.subTest(path=path.name):
                molecule = build_viewer_molecule(load_zmatrix(path))
                self.assertGreater(len(molecule.atoms), 3)
                self.assertEqual(len(molecule.atoms) - 3, len(molecule.dihedrals))
                self.assertEqual((), molecule.warnings)


def _distance(coords: tuple[tuple[float, float, float], ...], left: int, right: int) -> float:
    return math.dist(coords[left - 1], coords[right - 1])


def _angle_delta(left: float, right: float) -> float:
    return ((left - right + 180.0) % 360.0) - 180.0


if __name__ == "__main__":
    unittest.main()
