"""Geometry invariance and exact-index regression tests; no CCDC dependency."""
import unittest
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.spatial.distance import cdist

from Source.pdd_descriptor import (
    PDDDescriptor, calculate_pdd_from_arrays, pdd_descriptor_distance, pdd_to_amd,
    _earth_movers_distance, _transport_constraints,
)
from Source.pdd_clustering import (
    PDDClusterRecord, PDDClusterSettings, cluster_pdd, threshold_agreement,
)


def descriptor(values, weights=None, elements=None):
    values = np.atleast_2d(values).astype(float)
    n, k = values.shape
    return PDDDescriptor('fixture', k, np.ones(n) / n if weights is None else np.asarray(weights),
                         values, tuple(elements or ['C'] * n), True, False, 0.0)


class PeriodicDescriptorTests(unittest.TestCase):
    def test_skew_equivalent_basis(self):
        reduced = np.diag([5., .5, 5.])
        skew = np.array([[5., 0., 0.], [500., .5, 0.], [0., 0., 5.]])
        a = calculate_pdd_from_arrays([[0, 0, 0]], reduced, ['C'])
        b = calculate_pdd_from_arrays([[0, 0, 0]], skew, ['C'])
        np.testing.assert_allclose(a.distances, b.distances, atol=1e-12)
        self.assertAlmostEqual(a.distances[0, 0], .5)

    def test_supercell_and_permutation_rotation_translation(self):
        cell = np.array([[5., 0., 0.], [1., 4., 0.], [.4, .7, 6.]])
        xyz = np.array([[.5, 1., .8], [2., 1., 1.5], [1., 2., 2.]])
        elements = ['C', 'H', 'O']
        a = calculate_pdd_from_arrays(xyz, cell, elements, k=30)
        supercell = cell.copy()
        supercell[0] *= 2
        b = calculate_pdd_from_arrays(np.vstack([xyz, xyz + cell[0]]), supercell, elements * 2, k=30)
        self.assertLess(pdd_descriptor_distance(a, b), 1e-10)
        rotation, _ = np.linalg.qr(np.random.default_rng(4).normal(size=(3, 3)))
        c = calculate_pdd_from_arrays((xyz @ rotation + [8, -3, 7])[::-1], cell @ rotation, elements[::-1], k=30)
        self.assertLess(pdd_descriptor_distance(a, c), 1e-10)

    def test_certified_neighbours_against_independent_cloud(self):
        rng = np.random.default_rng(91)
        cell = np.array([[4., 0, 0], [1.9, 3., 0], [-1., 1., 4.]])
        xyz = rng.random((6, 3)) @ cell
        pdd = calculate_pdd_from_arrays(xyz, cell, ['C'] * 6, k=25, collapse=False, lexsort=False)
        from itertools import product
        cloud = np.vstack([xyz + np.asarray(shift) @ cell for shift in product(range(-4, 5), repeat=3)])
        reference = np.sort(cdist(xyz, cloud), axis=1)[:, 1:26]
        np.testing.assert_allclose(pdd.distances, reference, atol=1e-12)

    def test_composition_and_invalid_inputs(self):
        a = descriptor([[1, 2], [1, 3]], elements=['C', 'H'])
        b = descriptor([[1, 2], [1, 3]], elements=['C', 'O'])
        with self.assertRaises(ValueError):
            pdd_descriptor_distance(a, b)
        with self.assertRaises(ValueError):
            calculate_pdd_from_arrays([[0, 0, 0]], np.zeros((3, 3)), ['C'])
        with self.assertRaises(ValueError):
            pdd_to_amd(descriptor([[np.nan, 2]]))

    def test_assignment_and_general_transport_match_lp(self):
        rng = np.random.default_rng(11)
        for a, b in [(np.ones(5) / 5, np.ones(5) / 5),
                     (np.array([.2, .3, .5]), np.array([.1, .2, .3, .4])),
                     (np.array([1.]), np.array([.2, .8]))]:
            cost = rng.random((len(a), len(b)))
            expected = linprog(cost.ravel(), A_eq=_transport_constraints(len(a), len(b)),
                               b_eq=np.r_[a, b], bounds=(0, None), method='highs').fun
            self.assertAlmostEqual(_earth_movers_distance(a, b, cost), expected, places=12)

    def test_amd_bound_and_loss_of_distribution_information(self):
        a, b = descriptor([[2.], [4.]]), descriptor([[3.], [3.]])
        self.assertEqual(float(abs(pdd_to_amd(a) - pdd_to_amd(b))[0]), 0.)
        self.assertAlmostEqual(pdd_descriptor_distance(a, b), 1.)

    def test_recorded_nicotinamide_reference_pairs(self):
        # Recompute geometry descriptors; CCDC labels are archived provenance,
        # not inferred from EMD and not re-evaluated without a CSD installation.
        import json
        from Source.crystal_structure import CrystalStructure
        from Source.pdd_descriptor import calculate_pdd
        root = Path(__file__).parent / 'fixtures'
        labels = json.loads((root / 'reference_pairs.json').read_text())
        for kind in ('positive', 'negative'):
            pair = labels[kind]
            descriptors = [calculate_pdd(CrystalStructure.from_file(root / f'nicotinamide_z4_{i}.res')
                                         .expand_to_explicit_unit_cell()) for i in (pair['a'], pair['b'])]
            self.assertAlmostEqual(pdd_descriptor_distance(*descriptors), pair['pdd'], places=8)
        controlled = labels['controlled_positive']
        structures = [root / f'nicotinamide_z4_{controlled["base_id"]}.res',
                      root / f'{controlled["control_id"]}.res']
        descriptors = [calculate_pdd(CrystalStructure.from_file(path).expand_to_explicit_unit_cell())
                       for path in structures]
        distance = pdd_descriptor_distance(*descriptors)
        self.assertAlmostEqual(distance, controlled['pdd'], places=8)
        self.assertLess(labels['negative']['pdd'], distance)


class ClusterTests(unittest.TestCase):
    def test_prepared_comparisons_match_general_pdd(self):
        from dataclasses import replace
        from Source.pdd_clustering import _prepare_comparison, _prepared_distance
        rng = np.random.default_rng(53)
        for wa, wb, ea, eb in [
            ([.25]*4, [.25]*4, ['C','C','H','H'], ['H','C','H','C']),
            ([.2,.3,.2,.3], [.1,.4,.15,.35], ['C','C','H','H'], ['C','C','H','H']),
            ([.25]*4, [1/6]*6, ['C','C','H','H'], ['C']*3+['H']*3),
        ]:
            a=descriptor(np.sort(rng.random((len(wa),20)),axis=1),weights=wa,elements=ea)
            b=descriptor(np.sort(rng.random((len(wb),20)),axis=1),weights=wb,elements=eb)
            for typed in [True,False]:
                left,right=replace(a,typed=typed),replace(b,typed=typed)
                self.assertAlmostEqual(_prepared_distance(_prepare_comparison(left),_prepare_comparison(right)),
                                       pdd_descriptor_distance(left,right),places=12)

    def test_index_equals_exhaustive_through_block_merges(self):
        rng = np.random.default_rng(7)
        records = []
        for i in range(260):
            # Some separated families, some duplicates, and nontrivial row masses.
            family = i % 113
            d = descriptor(np.sort(rng.random((3, 20)) * .001 + family * .2 + np.arange(20), axis=1),
                           weights=[.2, .3, .5])
            records.append(PDDClusterRecord(str(i), float(i), 1000., d))
        settings = PDDClusterSettings(.01, energy_tolerance=None, density_tolerance=None, block_size=8)
        indexed = cluster_pdd(records[::-1], settings)
        exhaustive = cluster_pdd(records, settings, indexed=False)
        self.assertEqual(indexed.assignments, exhaustive.assignments)
        self.assertEqual(len(indexed.representatives), 113)
        self.assertLess(indexed.statistics['amd_checks'], exhaustive.statistics['amd_checks'] / 5)
        self.assertGreater(indexed.statistics['index_builds'], 2)

    def test_threshold_boundary_and_no_transitive_chaining(self):
        records = [PDDClusterRecord(label, i, 1000, descriptor([[value]]))
                   for i, (label, value) in enumerate([('A', 1.), ('B', 1.125), ('C', 1.25)])]
        result = cluster_pdd(records, PDDClusterSettings(.125, energy_tolerance=None, block_size=1))
        self.assertEqual(result.assignments, {'A': 'A', 'B': 'A', 'C': 'C'})

    def test_emd_rejects_amd_collision(self):
        a, b = descriptor([[2.], [4.]]), descriptor([[3.], [3.]])
        result = cluster_pdd([PDDClusterRecord('A', 0, 1000, a), PDDClusterRecord('B', .1, 1000, b)],
                             PDDClusterSettings(.1))
        self.assertEqual(len(result.representatives), 2)
        self.assertEqual(result.statistics['emd_calls'], 1)

    def test_gates_and_confirmation_and_composition(self):
        d = descriptor([[1., 2.]])
        records = [PDDClusterRecord('A', 0, 1000, d), PDDClusterRecord('B', 1, 1000, d)]
        self.assertEqual(len(cluster_pdd(records, PDDClusterSettings(.01)).representatives), 2)
        settings = PDDClusterSettings(.01, energy_tolerance=None)
        self.assertEqual(len(cluster_pdd(records, settings, confirm=lambda a, b: False).representatives), 2)
        records[1] = PDDClusterRecord('B', .1, 1000, descriptor([[1., 2.]], elements=['O']))
        self.assertEqual(len(cluster_pdd(records, settings).representatives), 2)

    def test_zero_radius_and_different_element_fractions(self):
        a = descriptor([[1., 2.], [1., 2.]], weights=[.5, .5], elements=['C', 'H'])
        b = descriptor([[1., 2.], [1., 2.]], weights=[.25, .75], elements=['C', 'H'])
        records = [PDDClusterRecord('A', 0., 1000., a),
                   PDDClusterRecord('B', .1, 1000., a),
                   PDDClusterRecord('C', .2, 1000., b)]
        settings = PDDClusterSettings(0., numeric_slack=0., block_size=1)
        self.assertEqual(cluster_pdd(records, settings).assignments, {'A':'A', 'B':'A', 'C':'C'})
        self.assertEqual(cluster_pdd(records, settings).assignments,
                         cluster_pdd(records, settings, indexed=False).assignments)

    def test_no_calibration_positive_does_not_claim_recall(self):
        scores = threshold_agreement([.1, .3], [False, False], [.2])
        self.assertEqual(scores[0]['false_merge'], 1)
        self.assertIsNone(scores[0]['recall'])

    def test_empty_duplicate_id_and_bad_settings(self):
        self.assertEqual(cluster_pdd([], PDDClusterSettings(.1)).assignments, {})
        record = PDDClusterRecord('A', 0, 1000, descriptor([[1.]]))
        with self.assertRaises(ValueError):
            cluster_pdd([record, record], PDDClusterSettings(.1))
        with self.assertRaises(ValueError):
            PDDClusterSettings(-.1)

    def test_manifest_cache_and_public_namespace(self):
        from csptoolbox import load_pdd_cluster_manifest, cluster_pdd as public_cluster
        from Source.CLI_scripts.cluster_pdd import main
        import json
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            structure = root / 'cell.res'
            structure.write_text('TITL test\nCELL .71073 5 5 5 90 90 90\nLATT -1\nSFAC C H\nC1 1 .1 .2 .3\nH1 2 .3 .2 .3\nEND\n')
            manifest = root / 'manifest.csv'
            manifest.write_text('id,path,energy,density\nA,cell.res,-1,1000\nB,cell.res,-.9,1000\n')
            cache = root / 'cache.npz'
            records, identity = load_pdd_cluster_manifest(manifest, k=12, cache=cache)
            reloaded, cached_identity = load_pdd_cluster_manifest(manifest, k=12, cache=cache)
            self.assertEqual(identity, cached_identity)
            self.assertEqual(public_cluster(records, PDDClusterSettings(.001)).assignments,
                             public_cluster(reloaded, PDDClusterSettings(.001)).assignments)
            main([str(manifest), '--k', '12', '--cache', str(cache), '--threshold', '.001',
                  '--output', str(root / 'clusters.json')])
            self.assertEqual(json.loads((root / 'clusters.json').read_text())['assignments'], {'A': 'A', 'B': 'A'})
            structure.write_text(structure.read_text().replace('.3 .2 .3', '.4 .2 .3'))
            with self.assertRaises(ValueError):
                load_pdd_cluster_manifest(manifest, k=12, cache=cache)


if __name__ == '__main__':
    unittest.main()
