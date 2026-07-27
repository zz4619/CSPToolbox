"""Focused tests for RDKit conformer generation and identity checks."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from csptoolbox.conformer_generation import (  # noqa: E402
    ConformerGenerationSettings,
    generate_conformer_ensemble,
    validate_ensemble_identity,
    write_conformer_artifacts,
)


CHLOROTHIAZIDE_SMILES = "C1=C2C(=CC(=C1Cl)S(=O)(=O)N)S(=O)(=O)N=CN2"
CHLOROTHIAZIDE_INCHI_KEY = "JBMKAUGHUNFTOL-UHFFFAOYSA-N"


class ConformerGenerationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.ensemble = generate_conformer_ensemble(
            CHLOROTHIAZIDE_SMILES,
            settings=ConformerGenerationSettings(
                num_conformers=20,
                random_seed=20260727,
                embedding_prune_rms_ang=0.25,
                deduplicate_rms_ang=0.35,
                max_iterations=1000,
                num_threads=1,
            ),
        )

    def test_chlorothiazide_identity_is_preserved(self) -> None:
        validate_ensemble_identity(
            self.ensemble,
            expected_formula="C7H6ClN3O4S2",
            expected_charge=0,
            expected_inchi_key=CHLOROTHIAZIDE_INCHI_KEY,
        )
        self.assertEqual(23, self.ensemble.molecule.GetNumAtoms())
        self.assertEqual(17, self.ensemble.molecule.GetNumHeavyAtoms())

    def test_selected_conformers_are_converged_and_energy_ranked(self) -> None:
        selected = self.ensemble.selected_records
        self.assertGreaterEqual(len(selected), 1)
        self.assertTrue(all(record.converged for record in selected))
        self.assertEqual(0.0, self.ensemble.best_record.relative_energy_kcal_mol)
        self.assertEqual(
            sorted(record.energy_kcal_mol for record in selected),
            [record.energy_kcal_mol for record in selected],
        )

    def test_identity_mismatch_fails_explicitly(self) -> None:
        with self.assertRaisesRegex(ValueError, "formal charge"):
            validate_ensemble_identity(self.ensemble, expected_charge=1)

    def test_artifacts_include_explicit_hydrogens_and_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            artifacts = write_conformer_artifacts(
                self.ensemble,
                temporary_directory,
                molecule_name="chlorothiazide",
            )
            xyz_lines = artifacts.lowest_energy_xyz_path.read_text(
                encoding="utf-8"
            ).splitlines()
            metadata = json.loads(
                artifacts.metadata_json_path.read_text(encoding="utf-8")
            )

            self.assertEqual("23", xyz_lines[0])
            self.assertEqual(25, len(xyz_lines))
            self.assertTrue(artifacts.selected_sdf_path.is_file())
            self.assertTrue(artifacts.manifest_csv_path.is_file())
            self.assertEqual("C7H6ClN3O4S2", metadata["formula"])
            self.assertEqual(0, metadata["formal_charge"])
            self.assertEqual(
                len(self.ensemble.selected_records),
                metadata["selected_conformer_count"],
            )


if __name__ == "__main__":
    unittest.main()
