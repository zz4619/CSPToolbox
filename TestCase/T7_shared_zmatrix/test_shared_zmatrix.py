"""Regression tests for shared Z-matrix topologies and Gaussian handoff."""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import (  # noqa: E402
    AtomRecord,
    CrystalStructure,
    zmatrix_topology_signature,
)
from Source.gaussian_input import GaussianInputBuilder, GaussianSettings  # noqa: E402
from Source.zmatrix_viewer import load_zmatrix  # noqa: E402


def _structure(name: str, coordinates: list[tuple[float, float, float]]) -> CrystalStructure:
    atoms = [
        AtomRecord(f"C{index}", "C", coordinate)
        for index, coordinate in enumerate(coordinates, start=1)
    ]
    return CrystalStructure(
        atoms=atoms,
        cell_parameters=(30.0, 30.0, 30.0, 90.0, 90.0, 90.0),
        name=name,
        explict_unit_cell=True,
    )


class SharedZMatrixTests(unittest.TestCase):
    def setUp(self) -> None:
        self.reference = _structure(
            "reference",
            [
                (10.00, 10.00, 10.00),
                (11.54, 10.00, 10.00),
                (12.10, 11.43, 10.00),
                (13.42, 11.75, 10.70),
            ],
        )
        self.alternative = _structure(
            "alternative",
            [
                (10.00, 10.00, 10.00),
                (11.52, 10.08, 10.02),
                (12.18, 11.36, 10.28),
                (13.28, 11.88, 9.43),
            ],
        )

    def test_shared_topology_changes_only_coordinate_values(self) -> None:
        topology = self.reference.generate_zmatrices()[0]
        evaluated = self.alternative.apply_zmatrix_topology(topology)

        self.assertEqual(
            zmatrix_topology_signature(topology),
            zmatrix_topology_signature(evaluated),
        )
        self.assertEqual(evaluated.warnings, [])
        self.assertTrue(
            any(
                left.dihedral_degrees != right.dihedral_degrees
                for left, right in zip(topology.entries, evaluated.entries)
                if left.dihedral_degrees is not None
            )
        )

    def test_saved_zmat_is_the_gaussian_geometry_source(self) -> None:
        topology = self.reference.generate_zmatrices()[0]
        evaluated = self.alternative.apply_zmatrix_topology(topology)
        builder = GaussianInputBuilder()

        with tempfile.TemporaryDirectory() as temporary_directory:
            zmat_path = Path(temporary_directory) / "molecule.zmat"
            zmat_path.write_text(
                builder.render_zmat_text("shared topology", evaluated),
                encoding="utf-8",
            )
            document = load_zmatrix(zmat_path)
            gaussian_text = builder.render_com_from_zmat_file(
                zmat_path,
                GaussianSettings(pcm_eps=None),
            )

        self.assertEqual(
            [atom.label for atom in document.atoms],
            evaluated.ordered_atom_labels,
        )
        self.assertNotIn("Variables:", gaussian_text)
        self.assertNotIn("Constants:", gaussian_text)
        self.assertIn(
            f"bnd2={document.atoms[1].bond_length:.6f}",
            gaussian_text,
        )

    def test_gaussian_builder_can_freeze_selected_internal_coordinates(self) -> None:
        topology = self.reference.generate_zmatrices()[0]
        builder = GaussianInputBuilder()
        frozen_value = topology.entries[3].dihedral_degrees
        self.assertIsNotNone(frozen_value)
        gaussian_text = builder.render_com_text(
            topology,
            GaussianSettings(
                pcm_eps=None,
                fixed_internal_coordinates=frozenset({"dih4"}),
            ),
        )

        self.assertNotIn("Variables:", gaussian_text)
        self.assertNotIn("Constants:", gaussian_text)
        self.assertNotIn("\ndih4=", gaussian_text)
        frozen_entry = topology.entries[3]
        self.assertIn(
            f"C {frozen_entry.bond_to} bnd4 "
            f"{frozen_entry.angle_to} ang4 "
            f"{frozen_entry.dihedral_to} {frozen_value:.6f}",
            gaussian_text,
        )

        with self.assertRaisesRegex(ValueError, "Unknown fixed internal coordinate"):
            builder.render_com_text(
                topology,
                GaussianSettings(
                    pcm_eps=None,
                    fixed_internal_coordinates=frozenset({"dih99"}),
                ),
            )


if __name__ == "__main__":
    unittest.main()
