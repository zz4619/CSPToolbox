# CSPToolbox contributor guidance

## Scope and orientation

This is a Python library for organic crystal-structure workflows. Its sibling
repositories (including CSPImperial and Delta-learning) are separate projects.
Check `git status --short` before changes and preserve unrelated edits. Do not
reset, stage or rewrite other work. Read the root README and closer guidance first.

Reusable implementations live in `Source/`; the public `csptoolbox` namespace
forwards exports lazily through `Source/__init__.py`. Keep both import paths working.
Put CLI argument handling in `Source/CLI_scripts/`, with reusable logic in the
library. Register installed commands in `pyproject.toml` rather than inventing
wrapper locations. Templates are under `Template/`; fixture/study folders are
under `TestCase/`. Do not treat the untracked `MACE-Multipoles/` project as part of
this change without checking its own scope.

For viewer tasks read [docs/VIEWER.md](docs/VIEWER.md). `model.py`, `readers.py` and
`preparation.py` prepare scenes; `html_export.py` serializes them and uses the single
renderer under `Source/zmatrix_viewer/assets/`. The standalone and embedded exports
must share that renderer. Gaussian generation imports the viewer's numeric parser,
so parser changes require shared Z-matrix/Gaussian compatibility checks.

## Setup and discovered commands

Python >=3.10 and the dependencies in `pyproject.toml` are required. From this repo:

```bash
python -m pip install -e .
python -m Source.CLI_scripts.zmatrix_viewer --help
python -m unittest discover -s TestCase/T4_zmatrix_viewer -p 'test_*.py' -v
python -m unittest discover -s TestCase/T3_zmatrix_generation -p 'test_*.py' -v
python -m unittest discover -s TestCase/T7_shared_zmatrix -p 'test_*.py' -v
python -m unittest discover -s TestCase/T5_vasp_results -p 'test_*.py' -v
python -m unittest discover -s TestCase/T6_conformer_generation -p 'test_*.py' -v
python -m pip wheel --no-deps --no-build-isolation . -w /tmp/csptoolbox-wheel
git diff --check
```

Run the suites relevant to the change, not expensive unrelated scientific jobs.
T3 discovery adds its test directory for local helper imports. T7 may be locally
present but not yet committed; check its availability in a fresh clone. CCDC is an
optional, separately installed/licensed integration. Do not invent a pip package
or cluster setup command for it. This workspace has used the `csp_310` Conda
environment, but select an available environment by verifying its dependencies.
If matplotlib has no writable cache, set `MPLCONFIGDIR` to a task-owned temporary
directory rather than changing global user settings.

No repository-wide formatter, lint command, type-checker configuration or CI
workflow was found. Report these as unconfigured rather than claiming they pass.
The wheel build requires the setuptools/wheel build dependencies to be available;
`--no-build-isolation` avoids implicitly downloading another build environment.
Verify viewer HTML/CSS/JS assets are included in a wheel after packaging edits.

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
