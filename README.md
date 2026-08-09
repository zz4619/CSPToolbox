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
- `cp2_local_min.py`: prepares experimental-structure local-minimisation jobs
  for CrystalPredictor2 and collects their machine-readable result status.

`Source/__init__.py` lazily exports the main classes and functions so lightweight
tools can import CSPToolbox without immediately importing the full scientific
stack.

The `csptoolbox/` package is a compatibility namespace that re-exports selected
modules from `Source/`.

## PyZMAT Compatibility

`Source.gaussian_input.GaussianInputBuilder` writes Gaussian symbolic
Z-matrix inputs in a PyZMAT-readable form. The generated `.com` files now include
an explicit `Variables:` section, and internal-coordinate values are written as
Gaussian-style assignments such as `bnd2=1.234567`.

For standalone Z-matrix files, `GaussianInputBuilder.render_zmat_text()` writes
the numeric `# ZMAT v1` format. These files keep Gaussian/CSPToolbox 1-based
atom references on disk. PyZMAT converts those references to its internal
0-based convention when loaded with:

```python
from pyzmat import ZMatrix

zmat = ZMatrix.load_from_csp_zmat("molecule.zmat")
```

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
- `csp-cp2-local-min`

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
table, and—when an executable is supplied—a self-contained, single-core CX3 PBS
bundle. The executable is copied byte-for-byte into the bundle; the runner uses
job-ID-keyed node-local scratch and copies results back. The supplied branch
binary needs only the CX3 production-tools and MKL runtime modules plus the
private NAG Kusari licence file at `$HOME/.nag/license.dat`. Keep `$HOME/.nag`
readable only by its owner (mode `700`) and the licence file at mode `600`;
create or update that file directly on CX3, outside this repository. The
generated runner exports `NAG_KUSARI_FILE` and stops before CP2 if the file is
missing or unreadable. For a permitted short licence preflight such as
`klcheck`, export `NAG_KUSARI_FILE="$HOME/.nag/license.dat"` first. Run the
actual minimisation through the cluster scheduler. A different execution-host
path can be supplied with `--nag-license-file`, but licence contents are never
staged or written to provenance. The library does not submit jobs. The status
command reads the structured
`CP2_LOCAL_MIN_RESULT_V1` record, optimizer information, final energies, and
output structure.

Runnable-job validation is conservative: truncated, ambiguous, assumed-order,
grossly mismatched, and globally inverted mappings do not pass automatic
validation. Globally inverted mappings remain inspectable in a preparation-only
bundle, but cannot be made runnable even with `--allow-unvalidated-mapping`:
CSPToolbox does not yet prove that inversion is achiral or symmetry-equivalent.
When the heavy-atom torsion anchors are planar within 0.05 A, the matcher uses
one canonical same-handed convention because the torsion sign is not a useful
orientation discriminator. Otherwise, an inverted torsion sign is accepted as
`planar_inversion_equivalent` only when both complete mapped molecules are
planar within 0.05 A; inversion within that molecular plane is then equivalent
to a proper 180-degree rotation.
The override remains available for the other explicitly audited mapping
failures. The single-core runner also fixes OpenMP and MKL to one thread.

Prepared runnable cases can be combined into one sequential, single-core PBS
job. Each case remains isolated and receives a unique node-local scratch
directory. The parent runner records failures and timeouts but continues to the
next experimental structure; rerunning the batch skips cases that already have
a successful structured CP2 result.

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

The supported pilot scope is single-component Z'=1 with an explicit landscape
reference and nonambiguous, nontruncated mapping. Z'>1 and multicomponent inputs
are gated behind `--allow-unvalidated` because no validation result set is
available. Optimizer convergence is reported separately and is not called a
confirmed local minimum without Hessian or perturbation evidence. TODOs are to
validate Z'>1/multicomponent cases, add a flexible-torsion-bound precheck, make
absolute LAM references portable, and determine whether every experimental
structure reaches a genuine local minimum on the CP2 PES.

## Z-Matrix Viewer

`csp-zmat-viewer` writes a standalone interactive HTML file from a numeric
`# ZMAT v1` file. The output has no Python server or external JavaScript
runtime requirement, so it can be produced on an HPC filesystem and opened
locally.

```bash
csp-zmat-viewer molecule.zmat --output molecule_viewer.html
```

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
