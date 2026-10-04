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
- `pdd_clustering.py`: exact AMD indexing, PDD filtering and optional COMPACK
  confirmation; see [the clustering guide](docs/PDD_CLUSTERING.md) for calibration,
  cache provenance, scientific limits and the `csp-cluster-pdd` command.
- `pdd_clustering_parallel.py`: checkpointed batch clustering of uniform-row PDD
  array datasets, retaining the serial energy/identifier ordering.
- `vasp_input.py`: VASP input builders for CIF and CONTCAR sources, including
  TPSS/PBE0 presets and default INCAR template paths.
- `vasp_results.py`: VASP result parsing for `vasp.out`, `OUTCAR`, `CONTCAR`,
  calculation health/status classification, and system summaries.
- `vasp_file_manifest.py`: reusable file manifest and tarball helpers for
  collecting selected VASP output files.
- `zmatrix_viewer/`: numeric `# ZMAT v1` parsing, Cartesian reconstruction,
  viewer molecule payload construction, and standalone HTML export.
- `cp2_local_min.py`: prepares experimental-structure local-minimisation jobs
  for CrystalPredictor2 and collects their machine-readable result status.


For the distinction between landscape clustering and experimental matching,
read the [PDD → RMSD20 workflow](docs/PDD_RMSD20_WORKFLOW.md). It documents the
H-excluded screening settings, exact RMSD20 acceptance criteria, completed
calibration and remaining recall checks. The CLI's optional `--compack` uses
historical RMSD15 defaults and is not the RMSD20 study workflow.

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
- `csp-cluster-pdd` (AMD lookup, PDD filtering and optional COMPACK confirmation)
- `csp-zmat-viewer`
- `csp-view` (the same CLI, including coordinate files, overlays and crystals)
- `csp-cp2-local-min`
- `csp-cp2-local-min-batch`

The CLI directory also contains plotting, rendering, CSORM/CSOFM generation,
Gaussian generation, and VASP summary scripts that were moved out of the top
level of `Source/` during cleanup.

## CrystalPredictor2 Local Minimisation

`csp-cp2-local-min prepare` reorders an experimental structure into native CP2
order, writes fixed-column `expcrys.pdb`, and stages `input.in`, the unchanged
`potential.in`, and the native LAM database from the system's `5_Globalsearch`
directory (with `5_GlobSrch` accepted as a legacy name). The LAM supplies site
order/types, charges, and indexed references; unique labels and a second
connectivity check come from the system `Zmatrix`. A landscape structure whose
labels/order already match CP2 should be supplied as the mapping reference,
normally the closest RMSD_1/COMPACK match.

For a one-TYPE system, `--auto-single-type-occurrences` detects repeated
molecular components in the experimental asymmetric unit and changes only the
occurrence-count record in the case-local `input.in`. The manifest records both
input hashes and the old/new counts. This is required to prepare experimental
Z'>1 cases from a Z'=1 global-search input; `potential.in` and the LAM database
remain byte-identical to their sources.

```bash
csp-cp2-local-min prepare SYSTEM EXPERIMENTAL.res Local_Min_CP2/REFCODE \
  --reference 1=LANDSCAPE_REFERENCE.res --space-group P21/C \
  --auto-single-type-occurrences --allow-unvalidated \
  --cp2-executable /path/to/Minimise --stage-mode copy
# PBS_O_WORKDIR is the bundle contract: submit from inside the prepared folder.
cd Local_Min_CP2/REFCODE && qsub run_cp2_local_min.pbs
csp-cp2-local-min status Local_Min_CP2/REFCODE
```

Each prepared job contains a JSON provenance manifest, a TSV atom-mapping
table, and a versioned, engine-neutral `zmatrix_mapping.json`. The mapping
artifact records the CP2 global-search template, authoritative Z-matrix,
selected label-coordinate pairs, complete template/experimental internal
coordinates, fixed-pair RMSD, selection reason, and source hashes. A later
CSO-FM adapter can therefore consume exactly the same permutation after
verifying the Z-matrix hash; it must not rematch the atoms. When an executable
is supplied, the prepared job also contains a self-contained, single-core CX3
PBS bundle. The executable is copied byte-for-byte into the bundle; the runner
uses job-ID-keyed node-local scratch and copies results back. The supplied branch
binary needs only the CX3 production-tools and MKL runtime modules plus the
private departmental Linux NAG Kusari licence file at
`$HOME/.nag/ChemEngDept-nag_keys-2026-linux.txt`. Keep `$HOME/.nag` readable
only by its owner (mode `700`) and the licence file at mode `600`; create or
update that file directly on CX3, outside this repository, and do not share it
outside Chemical Engineering. The
generated runner exports `NAG_KUSARI_FILE` and stops before CP2 if the file is
missing or unreadable. For a permitted short licence preflight such as
`klcheck`, export
`NAG_KUSARI_FILE="$HOME/.nag/ChemEngDept-nag_keys-2026-linux.txt"` first. Run the
actual minimisation through the cluster scheduler. A different execution-host
path can be supplied with `--nag-license-file`, but licence contents are never
staged or written to provenance. The library does not submit jobs. The status
command reads the structured
`CP2_LOCAL_MIN_RESULT_V1` record, optimizer information, final energies, and
output structure.

The `--reference` inputs are CP2 global-search structures, not structures
polished with CSO-FM. Atom mapping is energy-model independent: graph-valid
candidates are ranked by a continuous, dimensionless RMS score against all
available coordinates in that CP2 template Z-matrix. The default scales are
0.05 A for bonds, 5 degrees for angles, 10 degrees for rigid torsions, and
15 degrees for independent torsions. Torsions use shortest-periodic
differences. The larger independent-torsion scale allows genuine conformational
motion while still distinguishing symmetry-related oxygen or hydrogen
permutations.
Candidates tied within the recorded primary-score tolerance are compared by
fixed-correspondence, all-atom RMSD after translation and a proper rotation
only. Atom pairs remain locked; COMPACK rematching, reflection, and torsional
optimization are not part of this final check. A remaining tie is recorded as
ambiguous before deterministic atom indices are used for reproducible output.

The 0.25 A bond and 15 degree angle "gross mismatch" counts are diagnostics,
not hard candidate-selection classes or automatic validation failures. They
remain in the mapping artifact for review, but a candidate is not discarded at
either boundary; the continuous score avoids discontinuous permutation changes.

Runnable-job validation is conservative: truncated, ambiguous, assumed-order,
reflected, and out-of-LAM-domain mappings do not pass automatic validation.
Rigid and independent torsions are compared to the CP2 global-search template
with shortest-periodic angular differences during atom assignment. After
selection, CP2 separately checks each independent torsion's shortest-periodic
distance to the full interval declared in `input.in`. An out-of-domain result
therefore remains a validation failure; the mapping score does not force a
candidate into the domain. The manifest explicitly records the score scales,
separate rigid- and independent-torsion diagnostics, the diagnostic-only gross
mismatch count, and that reflection was disabled.
The override remains available for other explicitly audited
mapping failures. The single-core runner also fixes OpenMP and MKL to one
thread.

Prepared runnable cases can be combined into one sequential, single-core PBS
job. Each case remains isolated and receives a unique node-local scratch
directory. The parent runner records failures and timeouts but continues to the
next experimental structure; rerunning the batch skips cases that already have
a successful structured CP2 result.

For inventory-driven preparation with `prepare_cp2_local_min_batch.py`, an
optional `exclude_reason` column keeps an in-scope experimental form in the
preparation-status table while preventing creation of its case directory and
excluding it from the PBS case manifest. Excluded rows are reported separately
and do not count as preparation failures. Use a specific scientific reason,
for example when a racemic crystal is incompatible with an enantiopure LAM
model; do not delete the inventory row or overwrite historical results.

```bash
csp-cp2-local-min assemble-batch Local_Min_CP2_All \
  --case-dir Local_Min_CP2_All/cases/SystemA/FORM01 \
  --case-dir Local_Min_CP2_All/cases/SystemA/FORM02 \
  --pbs-walltime 24:00:00 --per-case-timeout 15m
cd Local_Min_CP2_All && qsub run_all_cp2_local_min.pbs
```

The batch assembly step does not submit the job. It requires every listed case
to have the same executable checksum and PBS runtime settings, and writes a TSV
case manifest plus per-case and aggregate status files.

The supported pilot scope is single-component Z'=1 with an explicit CP2
global-search reference and nonambiguous, nontruncated mapping. Z'>1 and
multicomponent inputs are gated behind `--allow-unvalidated` because no
validation result set is available. Optimizer convergence is reported
separately and is not called a confirmed local minimum without Hessian or
perturbation evidence. TODOs are to validate Z'>1/multicomponent cases, make
absolute LAM references portable, and determine whether every experimental
structure reaches a genuine local minimum on the CP2 PES.

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
