# Portable structure viewer

The five workflows share one Python scene model and one renderer. Exported HTML
contains its data, CSS and JavaScript: no CDN, Python server or Codex runtime is
needed to open the finished file. Use a modern browser with Canvas 2D support.

Install from the repository root with `python -m pip install -e .`. Commands below
use `csp-view`; `csp-zmat-viewer` is an alias. Without installing the console script,
use `python -m Source.CLI_scripts.zmatrix_viewer` from the repository root.

## 1. Coordinates without a Z-matrix

```bash
csp-view molecule.xyz -o molecule.html
csp-view molecule.pdb -o molecule.html
csp-view conformers.sdf --frame 2 -o third-conformer.html
```

XYZ/extXYZ, PDB/ENT (including PDB without `CRYST1`), MOL/SDF, CIF/RES and ASE-readable
coordinate formats are accepted. `--frame` is zero-based for multi-frame XYZ,
PDB models and SDF records. CIF/RES currently read the first structure only.
No Z-matrix is invented. Coordinates are displayed directly, except that whole
periodic molecules are unwrapped across cell boundaries.

PDB names are retained. Repeated names are qualified by chain/residue, then site
number if necessary (including CP2 cell exports). Alternate locations other than
blank/A are rejected: supply one resolved conformation. MOL/SDF and XYZ normally
have no persistent atom identities, so display labels are generated in file order.

Connectivity comes from PDB `CONECT`, MOL/SDF bonds, or a covalent-distance estimate.
Distance inference limits each H to one contact; it is a display aid, not bond-order
perception. `--no-infer-bonds` retains only supplied and Z-matrix construction links.
`--covalent-scale` controls display bond inference (default 1.25). Periodic component
detection uses the existing `CrystalStructure.detect_molecules` convention (default
1.20, with the existing hydrogen relaxation).

## 2. Z-matrix and independent DoFs

```bash
csp-view molecule.zmat --independent-dofs bnd2 ang3 dih8 dih10
csp-view molecule.zmat --dofs-file input.in
```

Numeric `.zmat`/`.zmatrix` files use CSPToolbox's `# ZMAT v1` convention. Atom labels,
references, signed torsions and explicit `# bonds:` comments are preserved. The
sidebar initially shows the supplied independent DoFs; switch to all dihedrals or
all coordinates. Clicking a row highlights its atom chain and shows its bounds.
These are viewing controls, not geometry-editing sliders.

DoF names use **1-based Z-matrix row numbers**: `bndN`, `angN`, `dihN`. Supply a list,
a JSON object, a CP2 `input.in` containing rows such as `dih8 -60 300`, or a numeric
Z-matrix comment `# independent_dofs: dih8 dih10`. Explicit arguments/files override
the comment. A JSON definition can include `{"dih8": {"lower": -60, "upper": 300}}`.
Bounds describe the source definition; coordinates are not clamped to them.
Use unwrapped intervals for periodic ranges. Unknown coordinates are rejected.

## 3. Before/after conformer comparison

```bash
csp-view after.pdb --reference before.pdb -o comparison.html
csp-view after.xyz --reference before.xyz --zmatrix known.zmat -o comparison.html
```

The reference is translucent grey; the comparison uses element colours. Alignment
uses equal-weight, all-atom least squares with translation and a proper rotation
(determinant +1). It never reflects or deforms either conformation. The reported
RMSD is in Å. With a shared topology, the sidebar also shows reference values and
differences; torsion differences wrap to [-180, 180).

Identity matching uses exact labels. For unlabelled inputs, equal atom counts and
matching element sequences permit a **file-order assumption**, displayed as a
notice. Reordered same-element atoms cannot be detected from XYZ alone: supply
`--atom-map` when order has changed. There is no automatic graph/symmetry atom
permutation or heavy-atom-only fit.

For a periodic output containing multiple molecules, select a whole component
with `--molecule-index 1`; the analogous reference option is
`--reference-molecule-index`. These component indices are 1-based. If atom names
have changed (as in CP2 output), supply a map derived from the calculation's
provenance; do not assume label spelling survived minimization.

## 4. Crystal asymmetric unit and full cell

```bash
csp-view experiment.cif --crystal -o crystal.html
csp-view experiment.res --crystal -o crystal.html
```

The Cell selector switches from the supplied asymmetric unit to its full unit
cell. CIF expansion uses the file's symmetry operations and special-position
deduplication. RES uses the existing LATT/SYMM expansion; periodic PDB uses its
declared space group. The lattice outline is shown in both states. Molecules are
unwrapped intact, so some atoms can extend outside that outline.

If the reader explicitly identifies an input as an already expanded cell, the
viewer uses the existing spglib reduction to infer the ASU (`symprec=0.05 Å`) and
reports this inference. An unmarked PDB with missing symmetry cannot establish
its crystallographic ASU. Prefer CIF/RES with reliable symmetry metadata.

Partial occupancies are reported but have no occupancy-dependent rendering model.
Disorder, extended bonded networks, very large cells and ambiguous symmetry need
separate scientific review. This is a molecular-crystal viewer, not a crystallographic
refinement tool. Crystal mode cannot be combined with a conformer overlay or
Z-matrix mapping; select a molecule for those workflows.

## 5. Experimental coordinates mapped to a known topology

```bash
csp-view experiment.pdb --zmatrix known.zmat --atom-map atom-map.json -o mapping.html
csp-view experiment.pdb --zmatrix Zmatrix --dofs-file input.in -o mapping.html
```

`--zmatrix` accepts a numeric file or CP2 label-only rows headed by
`Z-matrix for molecule 1`. Label-only topologies obtain all distances, angles and
torsions from the supplied experimental coordinates. Numeric topologies can also
provide a grey reference conformation. With label-only topology there is no
reference geometry to overlay.

The map runs **reference/Z-matrix label → input label or 1-based atom index**:

```json
{"C1": "C7", "O1": "O9", "H1": "H2"}
```

Include every topology label exactly once. Values must identify distinct atoms
with matching elements. With `--molecule-index`, numeric map indices refer to the
selected component in its preserved input order. Exact labels can select a subset
of a larger input; an explicit map can do the same. The Atom mapping table exposes
the chosen correspondence. When a separate `--reference` is supplied with a
topology, its labels must match that topology (or it must retain unlabelled order).

The viewer measures the supplied geometry; it does not optimize it. Save mapped
Z-matrix downloads the evaluated numeric definition with full-precision values.
The displayed geometry can be rigidly aligned to a reference, but the saved internal
coordinates retain the experimental geometry. Degenerate/undefined internal
coordinates are rejected.

## Python and embedding

```python
from csptoolbox.zmatrix_viewer import (
    build_structure_view, build_comparison_view, build_mapping_view,
    build_crystal_view, build_gallery, write_viewer_html, render_viewer_fragment,
)

scene = build_structure_view("molecule.zmat", independent_dofs=["dih8", "dih10"])
write_viewer_html(scene, "molecule.html")
gallery = build_gallery([scene, build_crystal_view("experiment.cif")])
write_viewer_html(gallery, "gallery.html")
fragment = render_viewer_fragment(scene)  # unique root ID; includes its runtime
```

Embedding must preserve and execute the inline scripts. A host that strips scripts
or inserts HTML via `innerHTML` without executing scripts will not initialize the
viewer. Use the standalone export in that case. Multiple fragments are independent;
if manually specifying `root_id`, keep IDs unique. A mounted root exposes
`root.cspViewer.getState()` and `destroy()` for integrations. Regenerate historical
HTML to receive fixes; exports do not update themselves.

Legacy `Source.zmatrix_viewer` imports and `build_viewer_document` still work; the
latter returns the historical `ViewerMolecule`. New builders return `ViewerScene`;
`build_gallery` returns `ViewerDocument`. `viewer_payload` exposes the shared versioned
data model. Assets are packaged under `Source/zmatrix_viewer/assets/`; do not fork
them into independent project renderers.

Python calculations and serialized positions use full floating-point precision.
Displayed distances use four decimals and angles two decimals. Export records
source paths and SHA-256 hashes. Writers reject output paths that alias any input
structure, topology, DoF file or map, including symlinks and hard links.

Controls: drag or arrow keys rotate; wheel or +/- zoom; Reset view (or Home while
the canvas is focused) resets the camera; Clear selection or Escape clears the
coordinate selection. Labels, H atoms and reference visibility can be toggled.
Atoms belonging to a selected coordinate remain visible when other H atoms are hidden.

## Verification and examples

```bash
python -m unittest discover -s TestCase/T4_zmatrix_viewer -p 'test_*.py' -v
python TestCase/T4_zmatrix_viewer/generate_demo.py --output-dir Codex_workspace/viewer_demo
python -m pip wheel --no-deps --no-build-isolation . -w /tmp/csptoolbox-wheel
```

The demo creates five scenes and `two-viewers.html` for integration checks. Synthetic
fixtures are software tests, not experimental or optimized results. In the multi-repo
CSP sandbox, add `--cp2-root ../CSPImperial/test/CP2` to generate a real salicylic-acid
example using archived experimental and minimized PDBs and their CP2 mapping.
This option requires those external fixtures; the unit tests do not.

For browser QA, open the output in Chrome. Exercise scene selection, a DoF row,
all-coordinate filtering, Clear selection, camera reset, reference/labels/H toggles,
ASU/full cell, mapping-table expansion and mapped-Z-matrix download. On the two-viewer
page, actions in one instance must leave the other unchanged. Check normal and narrow
viewports, a dark theme and the browser console. Numeric tests alone cannot verify
button wiring. An optional loopback preview server is
`python -m http.server 8766 --bind 127.0.0.1 --directory Codex_workspace/viewer_demo`.

No automated browser runner is configured in this repo; record which interactions
were actually exercised. The current geometry tests cover identity, proper-rotation
alignment, chirality preservation, periodic unwrapping, symmetry expansion, numeric
round trips, shared-parser compatibility and source-file protection.
