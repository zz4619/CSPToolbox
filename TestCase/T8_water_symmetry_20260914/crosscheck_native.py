#!/usr/bin/env python3
"""Check independent reconstruction against frozen native CP2 replay coordinates."""
import argparse
import json
from pathlib import Path

import numpy as np
from ase.geometry import find_mic

from analyse_water import body_from_lam, reconstruct, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, default=Path(__file__).with_name('inputs'))
    parser.add_argument('--reference', type=Path, default=Path(__file__).with_name('results')/'native_selected.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    data = json.loads(args.reference.read_text())
    body, _ = body_from_lam(args.run/'rigid_lam_intra')
    starts = [np.array(list(map(float,s.split()))[1:-2]) for s in
              (args.run/'starting_crystals.out').read_text().splitlines()]
    results = []
    for trial in data['trials']:
        if trial['reason'] != 'accepted':
            continue
        g = trial['geometry']
        params = list(g['lengths_angles'])
        for i in range(1,17):
            params += g['asm'][str(i)]['zposnom']+g['asm'][str(i)]['euler']
        structure = reconstruct(params,body,'native_reference')
        native = np.array([a['xyz'] for i in range(1,17) for a in g['cellmol'][str(i)]['atoms']])
        xyz = np.array([a.coordinates for a in structure.atoms])
        cell = np.array(g['cell_columns'])
        _, distance = find_mic(xyz-native,cell,pbc=True)
        diffs = [np.max(np.abs(np.asarray(params)-x)) for x in starts]
        row = int(np.argmin(diffs))
        printed = reconstruct(starts[row],body,'printed_reference')
        pxyz = np.array([a.coordinates for a in printed.atoms])
        _, pdistance = find_mic(pxyz-native,cell,pbc=True)
        result = dict(native_sample_id=trial['sample_id'],native_rank=trial['rank'],start_row=row+1,
            max_parameter_rounding=float(diffs[row]),maximum_cell_difference_a=float(np.max(abs(cell-structure.cell.array))),
            maximum_atom_error_a=float(max(distance)),printed_start_maximum_atom_error_a=float(max(pdistance)))
        assert result['maximum_atom_error_a'] < 1e-11, result
        assert result['maximum_cell_difference_a'] < 1e-11, result
        assert result['max_parameter_rounding'] < 5.01e-7, result
        assert result['printed_start_maximum_atom_error_a'] < 2e-5, result
        results.append(result)
    assert len(results) == 2
    report = dict(reference_sha256=sha(args.reference),all_checked=True,comparisons=results)
    text = json.dumps(report,indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
