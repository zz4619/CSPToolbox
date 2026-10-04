# CSPToolbox

CSPToolbox is the reusable Python package for common crystal structure
prediction workflow utilities. The current cleanup separates importable library
code from one-off command-line scripts and keeps project-specific analysis in
`CSP-personal`.

## Source Layout

The importable modules live under `Source/`:

- `crystal_structure.py`: shared crystal/molecule data structures, CIF parsing,
  unit-cell expansion, reduction, symmetry checks, and z-matrix helpers.
- `gaussian_input.py`: Gaussian input builders and job artifact helpers.
- `csorm_input.py`: CSORM input builders and symmetry sanity checks.
- `csofm_input.py`: CSOFM input builders and Gaussian final-energy parsing.
- `mie_typing.py`: FIT/Mie atom typing, inter-spec parsing, and validation.
- `pdd_descriptor.py`: pointwise distance distribution descriptors and distance
  comparisons.
- `vasp_input.py`: VASP input builders for CIF and CONTCAR sources, including
  TPSS/PBE0 presets and default INCAR template paths.
- `vasp_results.py`: VASP result parsing for `vasp.out`, `OUTCAR`, `CONTCAR`,
  calculation health/status classification, and system summaries.
- `vasp_file_manifest.py`: reusable file manifest and tarball helpers for
  collecting selected VASP output files.
- `zmatrix_viewer/`: numeric `# ZMAT v1` parsing, Cartesian reconstruction,
  viewer molecule payload construction, and standalone HTML export.

`Source/__init__.py` lazily exports the main classes and functions so lightweight
tools can import CSPToolbox without immediately importing the full scientific
stack.

The `csptoolbox/` package is a compatibility namespace that re-exports selected
modules from `Source/`.

## PyZMAT Compatibility

`GaussianInputBuilder.render_zmat_text()` writes the numeric `# ZMAT v1`
interchange format. These files include stable atom labels and keep
Gaussian/CSPToolbox 1-based atom references on disk. PyZMAT converts those
references to its internal 0-based convention when loaded with:

```python
from pyzmat import ZMatrix

zmat = ZMatrix.load_from_csp_zmat("molecule.zmat")
```

Gaussian `.com` files use a separate syntax contract. Build them from a saved
and validated `# ZMAT v1` file with
`GaussianInputBuilder.write_com_from_zmat_file()`. The Gaussian symbolic
Z-matrix contains assignments such as `bnd2=1.234567` directly after the atom
rows and does not add a `Variables:` label. By default every internal
coordinate is variable. To hold selected coordinates fixed, pass their names
through `GaussianSettings.fixed_internal_coordinates`; the builder writes each
fixed value as a numeric literal in its Z-matrix atom row, omits its symbolic
assignment, and rejects unknown coordinate names. This avoids relying on a
Gaussian `Constants:` block for geometry symbols.

## CLI Scripts

Workflow scripts now live under `Source/CLI_scripts/`. These scripts are meant
to be runnable entry points around the library modules, not the primary location
for reusable logic.

Installed console commands defined in `pyproject.toml`:

- `csp-write-expanded-cifs`
- `csp-vasp-cif-inputs`
- `csp-vasp-pbe0-inputs`
- `csp-vasp-summary`
- `csp-vasp-manifest`
- `csp-zmat-viewer`
- `csp-view` (the same CLI, including coordinate files, overlays and crystals)

The CLI directory also contains plotting, rendering, CSORM/CSOFM generation,
Gaussian generation, and VASP summary scripts that were moved out of the top
level of `Source/` during cleanup.

## Molecular and Crystal Viewer

`csp-view` writes a standalone HTML viewer with its JavaScript, styling and
structure data included. It supports ordinary coordinate files without a
Z-matrix, Z-matrices with independent DoFs, grey-reference conformer overlays,
asymmetric-unit/full-cell views, and experimental-to-Z-matrix atom mappings.
`csp-zmat-viewer` remains available as an alias.

```bash
csp-view molecule.xyz --output molecule_viewer.html
csp-view molecule.zmat --independent-dofs dih8 dih10 dih12
csp-view after.pdb --reference before.pdb --output comparison.html
csp-view experiment.cif --crystal --output crystal.html
csp-view experiment.pdb --zmatrix Zmatrix --dofs-file input.in
```

See [the viewer guide](docs/VIEWER.md) for all five workflows, mapping direction,
crystal conventions, the Python API and browser verification. Independent DoFs
are supplied explicitly; the viewer does not infer which coordinates are free.

New integrations should use `csptoolbox.zmatrix_viewer` and
`render_viewer_fragment` rather than copying a renderer into a project-specific
script. The standalone page and fragment use the same packaged assets and
per-instance controls. Existing `Source.zmatrix_viewer` imports remain supported;
historical generated HTML and external project helpers are left in place.

## Cleanup Notes

Recent cleanup work moved reusable VASP functionality from personal scripts into
`Source/vasp_input.py`, `Source/vasp_results.py`, and
`Source/vasp_file_manifest.py`. Command-line wrappers were moved into
`Source/CLI_scripts/`.

The HA pair/ddU analysis was intentionally left out of CSPToolbox because it is
project-specific analysis. It now lives in:

`/Users/zianzhan/Desktop/CSP_sandbox/CSP-personal/Results/HA_ddU_analysis`

Older personal workflow scripts from `CSP-personal/2_VASP` were archived in:

`/Users/zianzhan/Desktop/CSP_sandbox/CSP-personal/legacy_scripts/2_VASP`

## Development

Install in editable mode from the repository root:

```bash
python -m pip install -e /Users/zianzhan/Desktop/CSP_sandbox/CSPToolbox
```

Run a lightweight syntax check:

```bash
python -m py_compile Source/*.py Source/CLI_scripts/*.py csptoolbox/*.py
```

Viewer and shared Z-matrix regression checks, from this repository root:

```bash
python -m unittest discover -s TestCase/T4_zmatrix_viewer -p 'test_*.py' -v
python -m unittest discover -s TestCase/T3_zmatrix_generation -p 'test_*.py' -v
python -m unittest discover -s TestCase/T7_shared_zmatrix -p 'test_*.py' -v
```

See [AGENTS.md](AGENTS.md) for scope boundaries, other discovered test commands,
packaging checks and completion criteria. No repository-wide lint, formatter,
type-checker or CI command is currently configured.
