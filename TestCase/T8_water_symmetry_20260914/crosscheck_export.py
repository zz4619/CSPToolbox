#!/usr/bin/env python3
"""Compare all eight archived water geometries with the native CP2 exporter.

Use direct Cartesian differences against the native ASU in original atom order.
Check the separate native full-cell export for integer lattice shifts only.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from analyse_water import body_from_lam, reconstruct, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, default=Path(__file__).with_name('inputs'))
    parser.add_argument('--reference', type=Path,
                        default=Path(__file__).with_name('results')/'native_export_all.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    reference = json.loads(args.reference.read_text())
    body, _ = body_from_lam(args.run/'rigid_lam_intra')
    finals = [list(map(float, line.split())) for line in
              (args.run/'crystals.out').read_text().splitlines()]
    starts = [list(map(float, line.split())) for line in
              (args.run/'starting_crystals.out').read_text().splitlines()]
    assert len(finals) == len(starts) == 4
    assert len(reference['records']) == 8
    assert {(r['row'], r['state']) for r in reference['records']} == {
        (i, state) for i in range(1, 5) for state in ('start', 'final')}
    comparisons = []
    for record in reference['records']:
        final, start = finals[record['row']-1], starts[record['row']-1]
        prefix = 6+int(final[5])+int(final[4])
        assert int(final[0]) == int(start[-2]) == 1
        assert int(final[1]) == record['native_status'] == 0
        assert int(final[4]) == 16
        assert len(final) == prefix+12+102 and len(start) == 105
        params = start[1:-2] if record['state'] == 'start' else final[prefix+12:]
        structure = reconstruct(params, body, 'native_export_crosscheck')
        xyz = np.array([atom.coordinates for atom in structure.atoms])
        native_xyz = np.asarray(record['xyz'])
        native_cell = np.asarray(record['cell'])
        assert xyz.shape == native_xyz.shape == (48, 3)
        assert native_cell.shape == (3, 3)
        assert record['asu_labels'] == record['full_labels'] == ['O1', 'H2', 'H2']*16
        atom_error = float(np.max(np.abs(xyz-native_xyz)))
        cell_error = float(np.max(np.abs(structure.cell.array-native_cell)))
        assert atom_error < 1e-11 and cell_error < 1e-11
        full_xyz = np.asarray(record['full_xyz'])
        assert full_xyz.shape == (48, 3)
        shift = (xyz-full_xyz) @ np.linalg.inv(native_cell)
        residual = (shift-np.rint(shift)) @ native_cell
        full_error = float(np.max(np.linalg.norm(residual, axis=1)))
        assert full_error < 1e-11
        comparisons.append(dict(row=record['row'], state=record['state'], atoms=48,
                                maximum_asu_cartesian_component_error_a=atom_error,
                                maximum_cell_component_error_a=cell_error,
                                maximum_full_cell_periodic_atom_error_a=full_error))
    report = dict(reference_sha256=sha(args.reference),
                  native_export_source_sha256=reference['source_sha256'],
                  all_checked=True, alignment=False, asu_periodic_wrapping=False,
                  full_cell_check='integer lattice translation only',
                  atom_permutation=False, comparisons=comparisons)
    text = json.dumps(report, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
