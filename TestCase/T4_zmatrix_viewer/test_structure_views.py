"""Scientific identity, geometry, file-safety and export regressions for all workflows."""
from dataclasses import replace
import json
import math
import os
from pathlib import Path
import re
import tempfile
import unittest

import numpy as np
from csptoolbox.zmatrix_viewer import (
    build_structure_view, build_comparison_view, build_mapping_view, build_crystal_view,
    build_gallery, build_viewer_document, load_zmatrix, parse_zmatrix_text,
    reconstruct_coordinates, render_viewer_fragment, render_viewer_html, write_viewer_html,
    write_viewer_fragment, viewer_payload, read_dofs,
)
from Source.zmatrix_viewer.cli import main
from Source.zmatrix_viewer.preparation import align_comparison, map_atoms
from Source.zmatrix_viewer.readers import molecule_from_coordinates, read_structure

FIXTURES = Path(__file__).parent / 'fixtures'


def xyz(path, elements, positions):
    path.write_text(f'{len(elements)}\nfixture\n' + ''.join(
        e + ' ' + ' '.join(format(float(x), '.17g') for x in p) + '\n' for e, p in zip(elements, positions)))
    return path


def pdb(path, labels, elements, positions, header=''):
    path.write_text(header + ''.join(
        f'HETATM{i:5d} {label:>4s} MOL A   1    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00          {e:>2s}\n'
        for i, (label, e, p) in enumerate(zip(labels, elements, positions), 1)))
    return path


def points(molecule):
    return np.array([a.coordinates for a in molecule.atoms])


def distances(molecule):
    p = points(molecule)
    return np.linalg.norm(p[:, None] - p[None, :], axis=-1)


class StructureViewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_plain_xyz_and_pdb_need_no_zmatrix_or_cell(self):
        coordinates = [(0, 0, 0), (1.43, 0, 0), (1.7, .9, .2)]
        x = xyz(self.root/'mol.xyz', ['C', 'O', 'H'], coordinates)
        p = pdb(self.root/'mol.pdb', ['C7', 'O9', 'H2'], ['C', 'O', 'H'], coordinates)
        p.write_text(p.read_text() + 'CONECT    1    2\nCONECT    2    3\n')
        for path in [x, p]:
            s = build_structure_view(path)
            np.testing.assert_allclose(points(s.molecule), coordinates, atol=0, rtol=0)
            self.assertFalse(s.coordinates)
            self.assertIsNone(s.full_cell)
        s = build_structure_view(p, infer_bonds=False)
        self.assertEqual(['C7', 'O9', 'H2'], [a.label for a in s.molecule.atoms])
        self.assertEqual([(1, 2), (2, 3)], [(b.left, b.right) for b in s.molecule.bonds])

    def test_numeric_zmatrix_and_explicit_dof_metadata(self):
        path = FIXTURES/'chain.zmat'
        s = build_structure_view(path)
        self.assertEqual(['dih4'], [c.name for c in s.coordinates if c.independent])
        s = build_structure_view(path, independent_dofs={'bnd2': {'lower': '1.4', 'upper': '1.6'}, 'ang3': {}})
        self.assertEqual(['bnd2', 'ang3'], [c.name for c in s.coordinates if c.independent])
        self.assertEqual(1.4, s.coordinates[0].lower)
        inp = self.root/'input.in'
        inp.write_text('header\ndih4 -60.0 300.0\n')
        self.assertEqual({'dih4': {'lower': -60., 'upper': 300.}}, read_dofs(inp))
        with self.assertRaisesRegex(ValueError, 'absent'):
            build_structure_view(path, independent_dofs=['dih90'])

    def test_pdb_repeated_sites_and_multiple_models(self):
        path = pdb(self.root/'repeated.pdb', ['C1','O1','C1','O1'], ['C','O','C','O'],
                   [(0,0,0),(1.4,0,0),(5,0,0),(6.4,0,0)])
        s = build_structure_view(path)
        self.assertEqual(4, len(set(a.label for a in s.molecule.atoms)))
        frame = pdb(self.root/'frame.pdb', ['C1','O1'], ['C','O'], [(0,0,0),(1.4,0,0)]).read_text()
        path.write_text('MODEL        1\n'+frame+'ENDMDL\nMODEL        2\n'+frame.replace('1.400','1.600')+'ENDMDL\n')
        self.assertAlmostEqual(1.6, distances(build_structure_view(path, frame=1).molecule)[0,1])
        with self.assertRaisesRegex(ValueError, 'unavailable'):
            build_structure_view(path, frame=2)
        path = pdb(self.root/'hydrogen.pdb', ['1HG1'], [''], [(0,0,0)])
        self.assertEqual('H', build_structure_view(path).molecule.atoms[0].element)

    def test_sdf_connectivity_and_xyz_frames(self):
        from rdkit import Chem
        molecule = Chem.MolFromSmiles('CO')
        conformer = Chem.Conformer(2)
        conformer.SetAtomPosition(0, (0, 0, 0)); conformer.SetAtomPosition(1, (1.4, 0, 0))
        molecule.AddConformer(conformer)
        path = self.root/'molecule.sdf'
        with Chem.SDWriter(str(path)) as writer:
            writer.write(molecule)
        s = build_structure_view(path, infer_bonds=False)
        self.assertEqual(1, len(s.molecule.bonds))
        self.assertAlmostEqual(1.4, distances(s.molecule)[0,1])
        path = self.root/'frames.xyz'
        path.write_text('1\nfirst\nC 0 0 0\n1\nsecond\nC 1 2 3\n')
        self.assertEqual((1.,2.,3.), build_structure_view(path, frame=1).molecule.atoms[0].coordinates)

    def test_overlay_proper_rotation_and_chirality(self):
        p = np.array([(0, 0, 0), (1.2, 0, 0), (0, 1.3, 0), (0, 0, 1.4)])
        ref = molecule_from_coordinates(['C1','N1','O1','F1'], ['C','N','O','F'], p)
        rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        comp = molecule_from_coordinates(['C1','N1','O1','F1'], ['C','N','O','F'], p @ rotation + [5, -3, 2])
        fitted, rmsd = align_comparison(ref, comp)
        self.assertLess(rmsd, 1e-14)
        np.testing.assert_allclose(distances(comp), distances(fitted), atol=1e-14, rtol=0)
        mirrored = replace(comp, atoms=tuple(replace(a, coordinates=(-a.coordinates[0], *a.coordinates[1:])) for a in comp.atoms))
        self.assertGreater(align_comparison(ref, mirrored)[1], .5)

    def test_xyz_overlay_notices_order_assumption(self):
        ref = xyz(self.root/'before.xyz', ['C','O','H'], [(0,0,0),(1.4,0,0),(1.4,.9,0)])
        comp = xyz(self.root/'after.xyz', ['C','O','H'], [(5,1,0),(6.4,1,0),(6.4,1.9,0)])
        s = build_comparison_view(ref, comp)
        self.assertLess(s.rmsd, 1e-14)
        self.assertIn('file atom order', s.mapping_method)
        self.assertTrue(any('unchanged atom order' in n for n in s.notices))

    def test_mapping_reordered_pdb_and_numeric_roundtrip(self):
        template = build_viewer_document(FIXTURES/'chain.zmat')
        order = [2, 0, 3, 1]
        path = pdb(self.root/'exp.pdb', [f'X{i+1}' for i in range(4)], ['C']*4,
                   points(template)[order] + [3, 5, -2])
        with self.assertRaisesRegex(ValueError, 'identities differ'):
            build_mapping_view(path, FIXTURES/'chain.zmat')
        mapping = {f'C{i+1}': f'X{order.index(i)+1}' for i in range(4)}
        s = build_mapping_view(path, FIXTURES/'chain.zmat', atom_map=mapping)
        original = points(read_structure(path).molecule)[np.argsort(order)]
        np.testing.assert_allclose(distances(s.molecule), np.linalg.norm(original[:,None]-original[None,:],axis=-1), atol=1e-14, rtol=0)
        evaluated = parse_zmatrix_text(s.mapped_zmatrix)
        reconstructed = np.array(reconstruct_coordinates(evaluated.atoms))
        np.testing.assert_allclose(np.linalg.norm(reconstructed[:,None]-reconstructed[None,:],axis=-1), distances(s.molecule), atol=1e-12, rtol=0)
        self.assertEqual(('C1','X2'), s.mapping[0])
        self.assertAlmostEqual(60, next(c.value for c in s.coordinates if c.name=='dih4'), delta=.05)

    def test_mapping_validation_and_label_only_topology(self):
        path = pdb(self.root/'mol.pdb', ['C1','O1','H1'], ['C','O','H'], [(0,0,0),(1.4,0,0),(1.4,.9,0)])
        topology = self.root/'Zmatrix'
        topology.write_text('Z-matrix for molecule 1\nC1\nO1 C1\nH1 O1 C1\n')
        scene = build_mapping_view(path, topology, independent_dofs=['ang3'])
        self.assertIsNone(scene.reference)
        self.assertAlmostEqual(90, scene.coordinates[-1].value)
        mol = scene.molecule
        with self.assertRaisesRegex(ValueError, 'one-to-one'):
            map_atoms(mol.atoms, mol, {'C1':1, 'O1':1, 'H1':3})
        with self.assertRaisesRegex(ValueError, 'Element mismatch'):
            map_atoms(mol.atoms, mol, {'C1':2, 'O1':1, 'H1':3})

    def test_crystal_asu_and_full_cell(self):
        s = build_crystal_view(FIXTURES/'pair.cif')
        self.assertEqual(2, len(s.molecule.atoms))
        self.assertEqual(4, len(s.full_cell.atoms))
        self.assertEqual(3, len(s.cell))
        self.assertAlmostEqual(1.4, distances(s.molecule)[0,1])
        self.assertEqual(4, len(set(a.label for a in s.full_cell.atoms)))

    def test_periodic_boundary_unwrap(self):
        path = self.root/'boundary.cif'
        path.write_text((FIXTURES/'pair.cif').read_text().replace('0.10 0.20', '0.94 0.20').replace('0.24 0.20', '0.08 0.20'))
        s = build_structure_view(path)
        self.assertAlmostEqual(1.4, distances(s.molecule)[0,1])
        self.assertFalse(build_structure_view(path, infer_bonds=False).molecule.bonds)

    def test_crystal_special_position_is_not_duplicated(self):
        path = self.root/'special.cif'
        path.write_text((FIXTURES/'pair.cif').read_text().replace('C1 C 0.10 0.20 0.30 1', 'C1 C 0 0 0 1').replace('O1 O 0.24 0.20 0.30 1', ''))
        s = build_crystal_view(path)
        self.assertEqual(1, len(s.full_cell.atoms))

    def test_metadata_and_reference_files_cannot_be_overwritten(self):
        dofs = self.root/'dofs.json'; dofs.write_text('["dih4"]')
        ref = self.root/'ref.zmat'; ref.write_text((FIXTURES/'chain.zmat').read_text())
        scene = build_structure_view(FIXTURES/'chain.zmat', reference=ref, dofs_file=dofs)
        for path in [dofs, ref]:
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                write_viewer_html(scene, path)

    def test_relative_input_protection_survives_changed_directory(self):
        source = self.root/'relative.zmat'
        source.write_text((FIXTURES/'chain.zmat').read_text())
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            scenes = [build_structure_view('relative.zmat'), build_viewer_document('relative.zmat')]
        finally:
            os.chdir(previous)
        for scene in scenes:
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                write_viewer_html(scene, source)

    def test_reject_invalid_bonds_and_source_overwrites(self):
        text = (FIXTURES/'chain.zmat').read_text()
        for bond in ['1-999', '0-2']:
            with self.subTest(bond=bond), self.assertRaises(ValueError):
                parse_zmatrix_text(text + '# bonds: ' + bond + '\n')
        source = self.root/'input.zmat'
        source.write_text(text)
        scene = build_structure_view(source)
        symlink = self.root/'symlink.html'; symlink.symlink_to(source)
        hardlink = self.root/'hardlink.html'; os.link(source, hardlink)
        for target in [source, symlink, hardlink]:
            for writer in [write_viewer_html, write_viewer_fragment]:
                with self.subTest(target=target.name), self.assertRaisesRegex(ValueError, 'overwrite'):
                    writer(scene, target)
        self.assertEqual(2, main([str(source), '-o', str(source)]))
        self.assertEqual(text, source.read_text())

    def test_safe_portable_multi_instance_payload(self):
        scene = build_structure_view(FIXTURES/'chain.zmat', title='</script><b>__ROOT_ID__ & molecule')
        gallery = build_gallery([scene, scene])
        self.assertEqual(2, len(viewer_payload(gallery)['scenes']))
        fragments = [render_viewer_fragment(gallery) for _ in range(2)]
        self.assertNotEqual(re.search(r'id="([^"]+)"',fragments[0])[1], re.search(r'id="([^"]+)"',fragments[1])[1])
        for fragment in fragments:
            payload = json.loads(re.search(r'<script type="application/json" data-role="data">(.*?)</script>', fragment, re.S)[1])
            self.assertEqual(scene.title, payload['scenes'][0]['title'])
            self.assertFalse(re.search(r'<script[^>]+src=', fragment))
            self.assertNotIn('window.__viewerReady', fragment)
        self.assertTrue(render_viewer_html(gallery).startswith('<!doctype html>'))


if __name__ == '__main__':
    unittest.main()
