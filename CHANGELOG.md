# Changelog

User-visible changes, newest first. Add an *Unreleased* entry with every
change; see [docs/MAINTENANCE.md](docs/MAINTENANCE.md#git-workflow).

## Unreleased (0.2.0.dev0)

### Added
- CP2 experimental local-minimisation workflow: `cp2_local_min`,
  `zmatrix_mapping`, and the `csp-cp2-local-min` and `csp-cp2-local-min-batch`
  commands. An `exclude_reason` inventory column keeps a form out of a batch
  with a recorded reason.
- Indexed PDD clustering (serial and parallel), the `csp-cluster-pdd` command,
  and typed PDD comparison by element fractions.
- `zmatrix_topology_signature()` and `CrystalStructure.apply_zmatrix_topology()`,
  which measure one shared Z-matrix topology on other conformers.
- `parse_outcar_forces()` with per-block force statistics.
- `CrystalStructure.to_ase_atoms()`.
- `csptoolbox.crystal`, `.csofm_input`, `.csorm_input`, `.gaussian_input`,
  `.mie_typing` and `.zmatrix_mapping` import paths.
- Tests:
  - one `python -m pytest` command for every case;
  - T1/T2 dataset regressions over CE755;
  - input-file snapshots (T11);
  - package contract tests (T12);
  - known-issue specifications (T13);
  - crystal-core unit tests (T10).
- CI on Python 3.10 and 3.12, running ruff and the full suite.
- Documentation: `docs/FORMATS.md`, `docs/KNOWN_ISSUES.md`,
  `docs/MAINTENANCE.md`, `TestCase/README.md`.

### Changed
- **Gaussian `.com` format.** Files are written from the saved `# ZMAT v1`
  file. Assignments follow the atom rows without a `Variables:` heading, and
  fixed coordinates (`GaussianSettings.fixed_internal_coordinates`) are written
  as numbers in the atom row. `# ZMAT v1` gains a `# labels:` line. See
  `docs/FORMATS.md` and KI-7: PyZMAT does not yet read fixed coordinates in this
  form.
- **Reduction to the asymmetric unit.** One atom is still kept per symmetry
  orbit, but atoms are now chosen from whole molecules where possible.
- **`crystal_structure.py` split into the `Source.crystal` package**, by topic.
  - `Source.crystal_structure` still exports every public name.
  - Module-private helpers (formerly `_name` functions in that file) now live in
    `Source.crystal.<module>`, some under public names; none was imported
    outside the file.
  - All 19,634 characterization outputs over CE755 are bit-identical before and
    after.
- Importing CSPToolbox no longer calls `matplotlib.use("Agg")` and no longer
  imports pyplot or `ase.io`. Images use the matplotlib Figure API and are
  pixel-identical.
- Test case `T8_cp2_local_min_input` is renamed `T9_cp2_local_min_input`.
  `T8_pdd_clustering` keeps its name because study records cite it.
- Gas-phase VASP run scripts accept any MPI rank count. CSO-RM job templates
  look for the binary under `$HOME`.

### Fixed
- CIF symmetry operations are no longer evaluated with `eval()`, which a
  crafted CIF could exploit. Results are bit-identical.

### Removed
- `pandas` from the dependencies; nothing imports it.
- `csptoolbox/zmatrix_viewer.py`, which the `csptoolbox/zmatrix_viewer/`
  package always shadowed.

## 0.1.0

Initial package: crystal structures and file formats, CSO-RM/CSO-FM/Gaussian/
VASP input builders, Mie typing, PDD descriptors, the Z-matrix viewer and RDKit
conformer generation.
