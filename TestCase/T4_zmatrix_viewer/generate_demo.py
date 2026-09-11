"""Generate five workflow examples and a two-instance page for browser QA.

Run from the repository root. Optional --cp2-root uses the existing Salicylic_Acid
CP2 numerical-test inputs/results; it never modifies those scientific fixtures.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from csptoolbox.zmatrix_viewer import (
    build_structure_view, build_comparison_view, build_crystal_view, build_mapping_view,
    build_gallery, write_viewer_html, render_viewer_fragment,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--cp2-root', type=Path, help='CSPImperial/test/CP2 directory (optional)')
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    fixtures = Path(__file__).parent/'fixtures'
    if args.cp2_root:
        source = args.cp2_root.resolve()/'Salicylic_Acid'
        initial = source/'expcrys.pdb'
        final = source/'result'/'minimise_final_structure.pdb'
        topology = source/'Zmatrix'
        dofs = source/'input.in'
        with (source/'result/provenance/cp2_atom_mapping.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        # The archived map pins CP2 site order. Select the first whole molecule
        # from the full-cell output and explicitly map its renamed atom sites.
        atom_map = {row['cp2_label']: int(row['cp2_site_index_1based']) for row in rows}
        (out/'demo-cp2-atom-map.json').write_text(json.dumps(atom_map, indent=2)+'\n')
        numeric_text = build_mapping_view(initial, topology, dofs_file=dofs).mapped_zmatrix
        numeric = out/'demo-initial.zmat'
        numeric.write_text(numeric_text)
        scenes = [
            build_structure_view(initial, title='1 · Salicylic acid — PDB without a Z-matrix'),
            build_structure_view(numeric, dofs_file=dofs, title='2 · Z-matrix and independent DoFs'),
            build_comparison_view(initial, final, zmatrix=topology, dofs_file=dofs,
                                  molecule_index=1, atom_map=out/'demo-cp2-atom-map.json',
                                  title='3 · CP2 minimization — grey experimental reference'),
            build_crystal_view(initial, title='4 · Experimental crystal — asymmetric unit / full cell'),
            build_mapping_view(initial, topology, dofs_file=dofs,
                               title='5 · Experimental atoms mapped to the CP2 Z-matrix'),
        ]
    else:
        numeric = fixtures/'chain.zmat'
        mol = build_structure_view(numeric).molecule
        initial = out/'demo-chain.xyz'
        initial.write_text('4\nSynthetic geometry fixture, not a scientific result\n' + ''.join(
            a.element + ' ' + ' '.join(format(v, '.17g') for v in a.coordinates) + '\n' for a in mol.atoms))
        final = out/'demo-rotated.zmat'
        final.write_text(numeric.read_text().replace('1 60.0', '1 85.0'))
        scenes = [
            build_structure_view(initial, title='1 · XYZ without a Z-matrix'),
            build_structure_view(numeric, title='2 · Z-matrix and independent DoFs'),
            build_comparison_view(numeric, final, title='3 · Two conformations — grey reference'),
            build_crystal_view(fixtures/'pair.cif', title='4 · Synthetic P -1 cell'),
            build_mapping_view(initial, numeric, title='5 · Coordinate-to-Z-matrix mapping'),
        ]
    gallery = build_gallery(scenes, title='CSPToolbox — five viewer workflows')
    write_viewer_html(gallery, out/'viewer-workflows.html')
    # Ordinary page with no host runtime: two viewers must not share IDs/state.
    (out/'two-viewers.html').write_text('<!doctype html><html><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        '<style>body{margin:16px}#first-viewer{color-scheme:light;margin-bottom:32px}'
        '#second-viewer{color-scheme:dark}</style></head><body>' +
        render_viewer_fragment(scenes[1], root_id='first-viewer') +
        render_viewer_fragment(scenes[2], root_id='second-viewer') + '</body></html>')
    print(out/'viewer-workflows.html')
    print(out/'two-viewers.html')


if __name__ == '__main__':
    main()
