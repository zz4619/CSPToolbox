# CSPToolbox contributor guidance

## Scope and orientation

This is a Python library for organic crystal-structure workflows. Its sibling
repositories (including CSPImperial and Delta-learning) are separate projects.
Check `git status --short` before changes and preserve unrelated edits. Do not
reset, stage or rewrite other work. Read the root README and closer guidance first.

Read [docs/MAINTENANCE.md](docs/MAINTENANCE.md) before changing code. It covers
where code goes, the coding rules, testing and git workflow.

Reusable implementations live in `Source/`. The crystal-structure core is the
`Source/crystal/` package, split by topic; `Source/crystal_structure.py` only
re-exports it, so add nothing there. `csptoolbox/` holds one forwarding module per
`Source` module. Keep both import paths working; T12 checks them, and
CSPImperial imports `Source.*`.

Put CLI argument handling in `Source/CLI_scripts/`, with reusable logic in the
library. Register installed commands in `pyproject.toml` rather than inventing
wrapper locations. Templates are under `Template/`; test cases and fixtures are
under `TestCase/` (index: `TestCase/README.md`). Do not treat the untracked
`MACE-Multipoles/` project, a separate repository, as part of this change.

Files exchanged with PyZMAT, ml-lams and CSPImperial are listed in
[docs/FORMATS.md](docs/FORMATS.md). Changing one changes another repository's input.

For viewer tasks read [docs/VIEWER.md](docs/VIEWER.md). `model.py`, `readers.py` and
`preparation.py` prepare scenes; `html_export.py` serializes them and uses the single
renderer under `Source/zmatrix_viewer/assets/`. The standalone and embedded exports
must share that renderer. Gaussian generation imports the viewer's numeric parser,
so parser changes require shared Z-matrix/Gaussian compatibility checks.

## Setup and discovered commands

Python >=3.10 and the dependencies in `pyproject.toml` are required. From this repo:

```bash
python -m pip install -e ".[dev]"
python -m pytest                    # every fast test, all cases (about 20 s)
python -m pytest -m slow            # full CE755 T1/T2 regressions (about 1 min)
python -m pytest TestCase/T11_input_snapshots   # one case
ruff check .                        # syntax errors and pyflakes findings
python -m pip wheel --no-deps --no-build-isolation . -w /tmp/csptoolbox-wheel
git diff --check
```

Run the cases relevant to the change. Before finishing any change to library
code, also run the full fast suite. Known scientific defects are strict expected
failures tagged `KI-n` (docs/KNOWN_ISSUES.md). If a change fixes one, remove its
marker and move the entry to *Resolved*. If a change alters generated input
files, regenerate the T11 snapshots with `CSPTOOLBOX_UPDATE_SNAPSHOTS=1` and
review the diff. CI (GitHub Actions) runs ruff and the whole suite on Python 3.10
and 3.12.

CCDC is an optional, separately installed/licensed integration. Do not invent a
pip package or cluster setup command for it. This workspace has used the `csp_310`
Conda environment, which lacks pytest and ruff; install the `dev` extra into an
environment you own, not into a shared one. `TestCase/conftest.py` keeps
matplotlib headless and redirects an unwritable `MPLCONFIGDIR`. No formatter or
type checker is configured; do not reformat files you are not otherwise changing.
The wheel build requires the setuptools/wheel build dependencies to be available;
`--no-build-isolation` avoids implicitly downloading another build environment.
T12 checks that the viewer HTML/CSS/JS assets are declared as package data.

For descriptor/clustering work read `docs/PDD_CLUSTERING.md`. Preserve full-cell
periodicity and typed element fractions. AMD shortlisting must not change results
relative to exhaustive PDD comparison. A PDD threshold is not interchangeable with
COMPACK RMSD; calibrate false merges and missed duplicates on independent pairs.

## Scientific and file conventions

- Preserve atom labels/order, element identity, units, chirality, Z-matrix references
  and numerical input precision. Document any deliberate change in these contracts.
- Crystal coordinates use Å; internal angles use degrees. Distinguish ASU from full
  cell and periodic unwrapping from a change in molecular conformation.
- DoF names refer to 1-based Z-matrix rows. Never infer the independent set from
  every available torsion. Use the calculation's definitions.
- Conformer fitting must not silently allow reflection or atom permutations. Require
  a trustworthy explicit map when labels/order do not establish correspondence.
- Treat existing scientific inputs, Gaussian/VASP outputs, numerical baselines and
  generated Z-matrices as evidence. Do not regenerate or overwrite them casually.
- Write review HTML, transient inputs and screenshots to ignored `Codex_workspace/`
  or outside the repo. Do not hand-edit historical generated HTML to fix a renderer.
- Preserve data templates and shell submission scripts unless the task requires them.
  Viewer work needs no HPC access, Gaussian run, Fortran edit or cluster submission.
- Do not place credentials, private keys or authentication commands in documentation.

## Done when

The requested workflows work through the public API and CLI, relevant regression
tests pass, and new assets survive packaging. For UI changes, exercise real controls
in a browser (including two simultaneous instances and a narrow viewport), inspect
the rendering and console, and report any unavailable browser check explicitly.
Maintain source-file protection and geometry/identity invariants. Update documented
commands and examples to match the implementation. Leave a scoped diff with unrelated
work preserved. The final response should identify changed areas/files, verification
results, material limitations/blockers and any recommended follow-up. Do not commit
or publish unless the user requests it or has already authorized that action.
