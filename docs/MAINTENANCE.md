# Maintaining CSPToolbox

This guide is for anyone changing the code, human or agent. [AGENTS.md](../AGENTS.md)
summarises it for coding agents; scientific conventions are listed there too.

## Where things live

| Path | Contents | Rule |
|---|---|---|
| `Source/crystal/` | crystal-structure core, split by topic (see its `__init__.py`) | algorithms as functions; `structure.py` holds the API class |
| `Source/<topic>.py` | input builders and parsers for one program or method (`csorm_input`, `csofm_input`, `gaussian_input`, `vasp_*`, `mie_typing`, `pdd_*`, `cp2_local_min`, `zmatrix_mapping`, `conformer_generation`) | reuse `Source.crystal`; do not re-implement it |
| `Source/zmatrix_viewer/` | `# ZMAT v1` parser and HTML viewer | read `docs/VIEWER.md` first |
| `Source/CLI_scripts/` | command-line entry points | argument handling only; logic goes in the library |
| `Source/crystal_structure.py` | compatibility path for `Source.crystal` | keep; do not add code here |
| `csptoolbox/` | public import path: one forwarding file per `Source` module | T12 fails if one is missing |
| `Template/` | job templates, potentials, INCARs | change only when the job itself must change |
| `TestCase/` | test cases `T<n>_<topic>/` and shared `testsupport/` | see [TestCase/README.md](../TestCase/README.md) |
| `docs/` | user guides, [FORMATS.md](FORMATS.md), [KNOWN_ISSUES.md](KNOWN_ISSUES.md), this file | |

## Everyday commands

```bash
python -m pip install -e ".[dev]"     # editable install with pytest and ruff
python -m pytest                      # all fast tests (about 20 s)
python -m pytest -m slow              # full T1/T2 dataset regressions (about 1 min)
python -m pytest TestCase/T11_input_snapshots   # one case
ruff check .                          # syntax errors and pyflakes findings
CSPTOOLBOX_UPDATE_SNAPSHOTS=1 python -m pytest TestCase/T11_input_snapshots  # after an intended output change
```

CI (`.github/workflows/tests.yml`) runs `ruff check .` and every test, slow
included, on Python 3.10 and 3.12 for each push to `main` and each pull request.
Keep `main` green.

## Writing code

1. **No import-time side effects.**
   - Importing a module must not change global state: no `matplotlib.use`, no
     environment changes, no file I/O.
   - Import slow or optional packages (pyplot, `ase.io`, `ccdc`) inside the
     function that needs them.
   - T12 enforces this for `Source` and the crystal core.
2. **No personal or cluster paths in library code.** Take paths as arguments or
   settings. Cluster paths belong in `Template/` job scripts.
3. **One home per concept.**
   - Before writing a reader, a bond rule or a symmetry helper, check
     `Source.crystal`. Duplicates drift apart; KI-11 is an example.
   - When duplication is unavoidable (for example, mirroring a Fortran rule),
     name its source and test the two against each other.
4. **Module layout.**
   - Docstring first, then imports, constants and types, public functions, and
     finally private helpers.
   - A leading underscore means "used only in this module". A function used by
     another module is public.
   - Split a module that mixes topics or grows past roughly 800 lines.
5. **Readable signatures.**
   - Type hints on public functions, plus a docstring stating units (Å, degrees,
     Hartree, eV) and conventions (0- or 1-based, fractional or Cartesian).
   - Prefer keyword-only arguments for options.
6. **Fail loudly on scientific data.**
   - Raise `ValueError` naming the offending value; never silently drop, repair
     or reorder atoms.
   - Use `warnings.warn` only for conditions a caller may knowingly accept.
7. **Deterministic output.** Generated files must be byte-stable for the same
   input. Seed any randomness explicitly (see `conformer_generation`).
8. **Compatibility.**
   - `Source.*` and `csptoolbox.*` are both public, and CSPImperial imports
     `Source.*`.
   - When you move code, leave a re-export (as `crystal_structure.py` does).
   - Remove a public name only after a CHANGELOG entry announced the
     deprecation in an earlier version.

## Testing

| Kind | Example | Use it for |
|---|---|---|
| Unit | `T10_crystal_units` | one function's contract, edge cases, error messages |
| Dataset regression | `T1`, `T2` (marked `slow` for the full set) | behaviour over the 755 CE755 structures |
| Approval snapshot | `T11_input_snapshots` | every generated program input, compared as text |
| Contract | `T12_package_contract` | import paths, side effects, console scripts, packaged assets |
| Known issue | `T13_known_issues` | a confirmed defect, as a strict expected failure (`KI-n`) |
| Workflow | `T3`–`T9` | end-to-end behaviour of one workflow |

**Adding a test case**
- Create `TestCase/T<next free number>_<topic>/`. Never reuse a number;
  outside records cite these paths.
- Add a `README.md` if the case has data. It should state the data's source and
  how it was made.
- Register the case in `TestCase/README.md`.
- Use pytest style and the helpers in `TestCase/testsupport/`.
- Mark anything slower than a few seconds as `@pytest.mark.slow`, and give it a
  fast smoke subset.

**Changing behaviour on purpose**
1. Run the tests and read the failures.
2. Regenerate the snapshots that change and review `git diff` line by line.
3. Commit code, snapshots and docs together.

**Fixing a known issue**
1. Make the `KI-n` test pass.
2. Remove its `xfail` marker.
3. Move the entry to *Resolved* in `docs/KNOWN_ISSUES.md`.

**Pure refactors** must leave every snapshot unchanged, and `pytest -m slow`
must still pass. For large moves, also compare a before/after characterization
of the outputs over T1/T2, as done for the `Source.crystal` split (see CHANGELOG
0.2.0).

## Scientific data and evidence

- Fixtures under `TestCase/` are evidence. Do not regenerate or overwrite them
  to make a test pass. Add new data alongside, with provenance.
- Atom labels and order, element identity, units, chirality, Z-matrix
  references and numerical precision are contracts; see AGENTS.md. A deliberate
  change to one needs a CHANGELOG entry and, if it is an interchange format, an
  update to `docs/FORMATS.md`.
- Study outputs belong in `CSP_studies/`, not in this repository.

## Interchange formats

The files listed in [FORMATS.md](FORMATS.md) are read by other repositories,
and a change to one is a change to their input. Follow
[Changing a format](FORMATS.md#changing-a-format).

## Git workflow

- **Branches.**
  - Start one branch per task from `main`.
  - Merge with `--no-ff` once tests pass, then delete the branch.
  - Keep study-only branches unmerged as archives.
- **Commits.**
  - Prefix the subject with a conventional type and optional scope:
    `feat(cp2):`, `fix(symmetry):`, `refactor:`, `test:`, `docs:`, `build:`,
    `ci:`, `chore:`.
  - Write the subject as an imperative sentence of at most 72 characters; the
    body says why.
  - Keep one logical change per commit.
- **Releases.**
  - Update `CHANGELOG.md` under *Unreleased* with every user-visible change.
  - When CSPImperial needs a stable point, set the version in `pyproject.toml`,
    move *Unreleased* under that version, and tag `v<version>`.

## Engineering backlog

Ordered by value. Science fixes are tracked separately in
[KNOWN_ISSUES.md](KNOWN_ISSUES.md).

1. **`cp2_local_min.py` (3.7k lines)** has its own RES/PDB readers, SHELX parser,
   bond graph and geometry helpers. Move it onto `Source.crystal` behind the
   T9 tests and split it into preparation, mapping I/O, PBS batch and status
   modules.
2. **`vasp_results.py` (1.1k lines)** mixes parsers and summaries; split it the
   same way.
3. **Project-specific plotting scripts** in `Source/CLI_scripts/plot_*.py` and
   `generate_csorm_*_from_latest_contcar.py` default to paths on one
   workstation. Either make the paths required arguments or move the scripts
   to `CSP-personal`.
4. **No type checker.** Add mypy in non-strict mode on `Source/crystal` first.
5. **Naming kept for compatibility:**
   - the `Source` package name;
   - the `explict_unit_cell` field, which is written into files;
   - `DUMMEY_SYSTEM_NAME` in templates.

   Rename only with a deprecation period.
6. **`requirements.txt`** is a snapshot of one conda environment; `pyproject.toml`
   is authoritative.
