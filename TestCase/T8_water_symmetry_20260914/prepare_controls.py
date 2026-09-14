#!/usr/bin/env python3
"""Known-symmetry controls: inversion-related waters and H-only broken symmetry."""
import argparse
import json
from pathlib import Path

import numpy as np
from ase.cell import Cell

from analyse_water import crystal, rotation, symmetry_check, write_high_precision_cif


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    cellpar = [12., 14., 16., 83., 91., 103.]
    cell = Cell.fromcellpar(cellpar).array
    angle = np.deg2rad(104.)
    body = np.array([[0.,0.,0.],[0.96,0.,0.],[0.96*np.cos(angle),0.96*np.sin(angle),0.]])
    rng = np.random.default_rng(20260914)
    first = []
    for _ in range(8):
        origin = rng.uniform(.1,.85,3) @ cell
        first.extend(body @ rotation(*rng.uniform(0,6,3)).T+origin)
    first = np.asarray(first)
    positions = np.concatenate([first,-first],axis=0)
    broken = positions.copy()
    broken[1:3] = (broken[1:3]-broken[0]) @ rotation(.73,.48,.25).T+broken[0]
    report = []
    for name,xyz in [('inversion_control',positions),('hydrogen_broken_control',broken)]:
        folder = args.output/name
        folder.mkdir()
        result = {'name':name,'variants':{}}
        for variant,inds,symbols in [('full_OH',list(range(48)),['O','H','H']*16),
                                   ('oxygen_only',list(range(0,48,3)),['O']*16),
                                   ('hydrogen_proxy',list(range(48)),['O','C','C']*16)]:
            structure = crystal(name+'_'+variant,cell,cellpar,xyz[inds],symbols)
            write_high_precision_cif(structure,folder/(variant+'.cif'))
            result['variants'][variant] = symmetry_check(structure,.01,None)
        expected = {'full_OH':2 if name=='inversion_control' else 1,
                    'oxygen_only':2,'hydrogen_proxy':2 if name=='inversion_control' else 1}
        assert all(result['variants'][k]['space_group_number']==v for k,v in expected.items())
        result['expected_space_group_numbers'] = expected
        # Apply the oxygen inversion directly to original species-labelled coordinates.
        full = crystal(name,cell,cellpar,xyz,['O','H','H']*16)
        from analyse_water import match_max_distance
        maximum, rms = match_max_distance(xyz,-xyz,['O','H','H']*16,cell)
        result['oxygen_inversion_applied_to_original_OH'] = {'maximum_distance_a':maximum,'rms_distance_a':rms}
        report.append(result)
    (args.output/'controls.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in r.items() if k!='variants'} for r in report],indent=2))


if __name__ == '__main__':
    main()
