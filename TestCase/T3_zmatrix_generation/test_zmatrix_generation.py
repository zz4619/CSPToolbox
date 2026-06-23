"""Smoke tests for CSPToolbox Z-matrix generation."""

from __future__ import annotations

import unittest

from zmatrix_generation_helpers import (
    classify_dihedral_references,
    first_zmatrix,
    sample_cases,
    validate_zmatrix,
)


class ZMatrixGenerationTests(unittest.TestCase):
    def test_synthetic_cases_generate_valid_zmatrices(self) -> None:
        for case in sample_cases():
            with self.subTest(case=case.name):
                zmatrix = first_zmatrix(case.structure)
                errors = validate_zmatrix(case.structure, zmatrix)
                self.assertEqual([], errors)
                self.assertEqual(
                    len(case.structure.atoms),
                    len(zmatrix.entries),
                )
                self.assertEqual(
                    len(case.structure.atoms),
                    len(zmatrix.ordered_atom_labels),
                )

    def test_expected_improper_labels_are_improper(self) -> None:
        for case in sample_cases():
            if not case.expected_improper_labels:
                continue
            with self.subTest(case=case.name):
                zmatrix = first_zmatrix(case.structure)
                classifications = classify_dihedral_references(case.structure, zmatrix)
                by_label = {item.atom_label: item.kind for item in classifications}
                for label in case.expected_improper_labels:
                    self.assertEqual("improper", by_label.get(label), label)

    def test_expected_proper_label_sets_are_present(self) -> None:
        for case in sample_cases():
            if not case.expected_proper_label_sets:
                continue
            with self.subTest(case=case.name):
                zmatrix = first_zmatrix(case.structure)
                classifications = classify_dihedral_references(case.structure, zmatrix)
                proper_label_sets = {
                    item.label_set for item in classifications if item.kind == "proper"
                }
                for expected in case.expected_proper_label_sets:
                    self.assertIn(expected, proper_label_sets)


if __name__ == "__main__":
    unittest.main()

