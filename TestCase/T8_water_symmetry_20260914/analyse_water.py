#!/usr/bin/env python3
"""Reconstruct archived rigid P1 waters and run existing CSPToolbox symmetry APIs.

This study script does not minimize, round/symmetrize originals, or edit library
behaviour. Native CP2 final records have four decimal places; their coordinate
uncertainty must be considered when interpreting very tight symmetry tolerances.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.cell import Cell
from ase.geometry import find_mic
from scipy.optimize import linear_sum_assignment
import spglib

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Source.crystal_structure import AtomRecord, CrystalStructure

TOLERANCES = [0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rotation(phi, theta, psi):
    """CP2 rot_matrix_from_Euler_angles.f90, applied as A @ body_vector."""
    cp, ct, cs = np.cos([phi, theta, psi])
    sp, st, ss = np.sin([phi, theta, psi])
    return np.array([[cp*cs-sp*ct*ss, sp*cs+cp*ct*ss, st*ss],
                     [-cp*ss-sp*ct*cs, -sp*ss+cp*ct*cs, st*cs],
                     [sp*st, -cp*st, ct]])


def lattice(parameters):
    """CP2 column-vector lattice, returned as ASE/CSPToolbox row vectors."""
    a, b, c, alpha, beta, gamma = parameters
    ca, cb, cg = np.cos([alpha, beta, gamma])
    sg = np.sin(gamma)
    w2 = 1-ca*ca-cb*cb-cg*cg+2*ca*cb*cg
    assert w2 > 0 and sg > 0
    cell = np.array([[a, 0, 0], [b*cg, b*sg, 0],
                     [c*cb, c*(ca-cb*cg)/sg, c*np.sqrt(w2)/sg]])
    cellpar = list(parameters[:3]) + list(np.rad2deg(parameters[3:]))
    assert np.max(np.abs(cell-Cell.fromcellpar(cellpar).array)) < 1e-12
    return cell, cellpar


def body_from_lam(path):
    rows = path.read_text().splitlines()
    values = {r.split()[0]: float(r.split()[1]) for r in rows if
              r.startswith(('bnd2 ', 'bnd3 ', 'ang3 '))}
    assert set(values) == {'bnd2', 'bnd3', 'ang3'}
    assert any(r.startswith('H2  3') and '1 bnd3  2 ang3' in r for r in rows)
    angle = np.deg2rad(values['ang3'])
    body = np.array([[0., 0., 0.], [values['bnd2'], 0., 0.],
                     [values['bnd3']*np.cos(angle), values['bnd3']*np.sin(angle), 0.]])
    return body, values


def crystal(name, cell, cellpar, coords, symbols):
    return CrystalStructure(atoms=[AtomRecord(f'{s}{i+1}', s, tuple(x))
                            for i, (s, x) in enumerate(zip(symbols, coords))],
                            cell_parameters=tuple(cellpar), lattice_matrix=tuple(map(tuple, cell)),
                            space_group='P 1', name=name, explict_unit_cell=True)


def reconstruct(parameters, body, name):
    cell, cellpar = lattice(parameters[:6])
    coords = []
    for mol in np.asarray(parameters[6:]).reshape(16, 6):
        # zposnom denotes fractional first-atom origin, despite historical COM labels.
        r = rotation(*mol[3:])
        assert np.max(np.abs(r.T @ r-np.eye(3))) < 1e-14
        coords.extend(body @ r.T + mol[:3] @ cell)
    return crystal(name, cell, cellpar, coords, ['O', 'H', 'H']*16)


def match_max_distance(first, second, symbols, cell):
    distances = []
    for element in sorted(set(symbols)):
        inds = [i for i,s in enumerate(symbols) if s == element]
        a, b = np.array(first)[inds], np.array(second)[inds]
        delta = (a[:, None, :]-b[None, :, :]).reshape(-1, 3)
        _, norms = find_mic(delta, cell, pbc=True)
        costs = norms.reshape(len(a), len(b))
        ii, jj = linear_sum_assignment(costs)
        distances.extend(costs[ii, jj].tolist())
    return max(distances), float(np.sqrt(np.mean(np.square(distances))))


def symmetry_check(structure, tolerance, out):
    detection = structure.detect_space_group_symmetry(symprec=tolerance)
    reduced = structure.reduce_to_asymmetric_unit(symprec=tolerance)
    # Exercise the existing reduction AND expansion APIs, not only spglib discovery.
    expanded = reduced.expand_to_explicit_unit_cell()
    symbols = [a.element for a in structure.atoms]
    expanded_symbols = [a.element for a in expanded.atoms]
    roundtrip = {'same_element_counts': Counter(symbols) == Counter(expanded_symbols)}
    if roundtrip['same_element_counts']:
        xyz1 = np.array([a.coordinates for a in structure.atoms])
        xyz2 = np.array([a.coordinates for a in expanded.atoms])
        # Reorder species blocks without assigning molecular identities.
        order1 = sorted(range(len(symbols)), key=lambda i: symbols[i])
        order2 = sorted(range(len(symbols)), key=lambda i: expanded_symbols[i])
        maximum, rms = match_max_distance(xyz1[order1], xyz2[order2],
            [symbols[i] for i in order1], structure.cell.array)
        roundtrip.update(maximum_distance_a=maximum, rms_distance_a=rms)
    cell = structure._spglib_cell()
    dataset = spglib.get_symmetry_dataset(cell, symprec=tolerance)
    ops = []
    frac = cell[1]
    for rotation_matrix, translation in zip(dataset.rotations, dataset.translations):
        transformed = (frac @ rotation_matrix.T + translation) @ cell[0]
        maximum, rms = match_max_distance(frac @ cell[0], transformed, symbols, cell[0])
        ops.append(dict(rotation=rotation_matrix.tolist(), translation=translation.tolist(),
                        maximum_distance_a=maximum, rms_distance_a=rms))
    primitive = spglib.find_primitive(cell, symprec=tolerance)
    primitive_count = len(primitive[1]) if primitive is not None else None
    if out is not None:
        reduced.to_file(out)
    return dict(tolerance_a=tolerance, space_group=detection.symbol,
        space_group_number=detection.number, hall_number=detection.hall_number,
        atomic_orbits=list(detection.equivalent_atoms), asu_atoms=len(reduced.atoms),
        asu_element_counts=dict(Counter(a.element for a in reduced.atoms)),
        symmetry_operations=ops, primitive_atoms=primitive_count,
        primitive_volume_a3=float(abs(np.linalg.det(primitive[0]))) if primitive is not None else None,
        reduction_roundtrip=roundtrip)


def write_high_precision_cif(structure, path):
    """Coordinate-only P1 CIF; no inferred symmetry or experimental statistics."""
    xyz = np.array([a.coordinates for a in structure.atoms])
    frac = xyz @ np.linalg.inv(structure.cell.array)
    lines = ['data_'+structure.name, "_symmetry_space_group_name_H-M 'P 1'",
             '_space_group_IT_number 1']
    for key, value in zip(['length_a','length_b','length_c','angle_alpha','angle_beta','angle_gamma'],
                          structure.cell_parameters):
        lines.append('_cell_'+key+' '+format(value,'.12f'))
    lines += ['loop_', '_space_group_symop_operation_xyz', "'x,y,z'", 'loop_',
              '_atom_site_label', '_atom_site_type_symbol', '_atom_site_fract_x',
              '_atom_site_fract_y', '_atom_site_fract_z', '_atom_site_occupancy']
    for atom, point in zip(structure.atoms, frac):
        lines.append(f'{atom.label} {atom.element} '+' '.join(format(x,'.12f') for x in point)+' 1')
    path.write_text('\n'.join(lines)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    body, dimensions = body_from_lam(args.run/'rigid_lam_intra')
    finals = [list(map(float,line.split())) for line in (args.run/'crystals.out').read_text().splitlines()]
    starts = [list(map(float,line.split())) for line in (args.run/'starting_crystals.out').read_text().splitlines()]
    assert len(finals) == len(starts) == 4
    report = dict(schema=1, csptoolbox_revision='a40c80c759a8007aaa1afd13e880e2fee4c1dbf1',
        spglib_version=spglib.__version__, input_sha256={p.name:sha(p) for p in
        [args.run/'rigid_lam_intra',args.run/'crystals.out',args.run/'starting_crystals.out']},
        body_geometry=dimensions, structures=[],
        limitations=['Final native parameters are printed to four decimal places; starts to six.',
            'Geometry reconstruction is independent Python implementation of pinned CP2 rigid geometry.',
            'No energy evaluation, minimization or physical phase identification is performed.',
            'Identical elements are interchangeable for geometric symmetry; model charges on H2/H3 differ slightly.'])
    for row, (final, start) in enumerate(zip(finals, starts), 1):
        prefix = 6+int(final[5])+int(final[4])
        assert int(final[0]) == int(start[-2]) == 1
        assert int(final[4]) == 16 and int(final[1]) == 0
        assert len(final) == prefix+12+102 and len(start) == 105
        for state, parameters in [('start',start[1:-2]),('final',final[prefix+12:])]:
            name = f'water16_{row:02d}_{state}'
            structure = reconstruct(parameters, body, name)
            out = args.output/name
            out.mkdir()
            write_high_precision_cif(structure, out/'original.cif')
            xyz = np.array([a.coordinates for a in structure.atoms])
            np.savez(out/'geometry.npz', cell=structure.cell.array, positions=xyz,
                     symbols=np.array([a.element for a in structure.atoms]))
            entry = dict(name=name, row=row, state=state, native_final_energy_kj_mol=final[prefix],
                         cell_parameters=list(structure.cell_parameters), full_atoms=48, variants={})
            for variant, indices, symbols in [
                    ('full_OH', list(range(48)), ['O','H','H']*16),
                    ('oxygen_only', list(range(0,48,3)), ['O']*16),
                    ('hydrogen_proxy', list(range(48)), ['O','C','C']*16)]:
                model = crystal(name+'_'+variant,structure.cell.array,structure.cell_parameters,xyz[indices],symbols)
                write_high_precision_cif(model,out/(variant+'.cif'))
                entry['variants'][variant] = [symmetry_check(model,tol,out/(variant+'_asu_'+str(tol)+'.res'))
                                             for tol in TOLERANCES]
            assert [(r['space_group_number'],r['atomic_orbits']) for r in entry['variants']['full_OH']] == [
                   (r['space_group_number'],r['atomic_orbits']) for r in entry['variants']['hydrogen_proxy']]
            report['structures'].append(entry)
    (args.output/'symmetry_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in s.items() if k not in ['variants','cell_parameters']} |
        {'groups':{v:[(r['tolerance_a'],r['space_group'],r['asu_atoms'],r['primitive_atoms']) for r in rs]
                   for v,rs in s['variants'].items()}} for s in report['structures']],indent=2))


if __name__ == '__main__':
    main()
