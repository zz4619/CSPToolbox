"""Regression tests for the engine-neutral Z-matrix atom matcher."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys
import unittest

import networkx as nx
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.zmatrix_mapping import (  # noqa: E402
    MAPPING_METHOD_VERSION,
    InternalCoordinateKey,
    MappingCandidate,
    MappingScoreSettings,
    ZMatrixSite,
    match_zmatrix_atoms,
    proper_rotation_rmsd,
    select_mapping_candidate,
)


def _candidate(
    *, score: float, rmsd: float, mapping: tuple[int, ...]
) -> MappingCandidate:
    return MappingCandidate(
        canonical_to_experimental=tuple(
            (f"A{position}", index)
            for position, index in enumerate(mapping, start=1)
        ),
        gross_bond_angle_mismatches=0,
        internal_coordinate_score=score,
        fixed_pair_all_atom_rmsd_angstrom=rmsd,
        heavy_atom_rmsd_angstrom=rmsd,
        bond_rms_delta_angstrom=0.0,
        bond_max_abs_delta_angstrom=0.0,
        angle_rms_delta_degrees=0.0,
        angle_max_abs_delta_degrees=0.0,
        rigid_torsion_rms_delta_degrees=0.0,
        rigid_torsion_max_abs_delta_degrees=0.0,
        reference_values=(),
        experimental_values=(),
    )


class ZMatrixMappingTests(unittest.TestCase):
    def test_hydrogen_permutation_uses_fixed_pair_all_atom_rmsd(self) -> None:
        sites = (
            ZMatrixSite(1, "C1", "C"),
            ZMatrixSite(2, "H1", "H", bond_to=1),
            ZMatrixSite(3, "H2", "H", bond_to=1, angle_to=2),
        )
        reference_graph = nx.Graph()
        experimental_graph = nx.Graph()
        for graph in (reference_graph, experimental_graph):
            graph.add_node(0, element="C")
            graph.add_node(1, element="H")
            graph.add_node(2, element="H")
            graph.add_edges_from(((0, 1), (0, 2)))
        reference_coordinates = {
            0: (0.0, 0.0, 0.0),
            1: (1.00, 0.00, 0.00),
            2: (-0.25, 0.95, 0.00),
        }
        # Experimental atom indices 1 and 2 deliberately occupy the opposite
        # template positions.  Heavy-atom RMSD is identical for both mappings.
        experimental_coordinates = {
            0: (4.0, -2.0, 1.0),
            1: (3.75, -1.05, 1.0),
            2: (5.00, -2.00, 1.0),
        }
        flexible = frozenset(
            {
                InternalCoordinateKey("bond", 2),
                InternalCoordinateKey("bond", 3),
                InternalCoordinateKey("angle", 3),
            }
        )

        decision = match_zmatrix_atoms(
            sites=sites,
            reference_graph=reference_graph,
            reference_coordinates=reference_coordinates,
            canonical_to_reference={"C1": 0, "H1": 1, "H2": 2},
            experimental_graph=experimental_graph,
            experimental_coordinates=experimental_coordinates,
            flexible_coordinates=flexible,
        )

        self.assertEqual(MAPPING_METHOD_VERSION, decision.method_version)
        self.assertEqual(2, decision.primary_tie_count)
        self.assertEqual(1, decision.final_tie_count)
        self.assertEqual("fixed_pair_all_atom_rmsd", decision.selection_reason)
        self.assertEqual(
            {"C1": 0, "H1": 2, "H2": 1},
            decision.selected.as_index_mapping(),
        )
        self.assertAlmostEqual(
            0.0, decision.selected.fixed_pair_all_atom_rmsd_angstrom, places=12
        )

    def test_nonplanar_reflection_cannot_be_removed_by_proper_rotation(self) -> None:
        target = np.asarray(
            [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
            ]
        )
        reflected = target.copy()
        reflected[:, 0] *= -1.0

        self.assertGreater(proper_rotation_rmsd(reflected, target), 0.40)

    def test_primary_score_wins_before_fixed_pair_rmsd(self) -> None:
        lower_primary = _candidate(score=0.10, rmsd=1.0, mapping=(1, 0))
        lower_rmsd = _candidate(score=0.20, rmsd=0.0, mapping=(0, 1))

        decision = select_mapping_candidate((lower_rmsd, lower_primary))

        self.assertEqual(lower_primary, decision.selected)
        self.assertEqual("unique_primary_score", decision.selection_reason)

    def test_rmsd_is_used_only_inside_primary_tolerance(self) -> None:
        settings = MappingScoreSettings(primary_score_atol=1.0e-8)
        first = _candidate(score=0.1, rmsd=0.5, mapping=(1, 0))
        within_tolerance = _candidate(
            score=0.1 + 0.5e-8, rmsd=0.1, mapping=(0, 1)
        )

        decision = select_mapping_candidate(
            (first, within_tolerance), settings=settings
        )

        self.assertEqual(within_tolerance, decision.selected)
        self.assertEqual("fixed_pair_all_atom_rmsd", decision.selection_reason)

    def test_unresolved_tie_is_explicit_and_deterministic(self) -> None:
        first = _candidate(score=0.1, rmsd=0.2, mapping=(1, 0))
        second = replace(first, canonical_to_experimental=(("A1", 0), ("A2", 1)))

        decision = select_mapping_candidate((first, second))

        self.assertTrue(decision.mapping_ambiguous)
        self.assertEqual(2, decision.final_tie_count)
        self.assertEqual("deterministic_unresolved_tie", decision.selection_reason)
        self.assertEqual((0, 1), decision.selected.deterministic_key)


if __name__ == "__main__":
    unittest.main()
