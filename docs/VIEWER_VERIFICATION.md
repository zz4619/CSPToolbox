# Viewer verification — 2026-09-11

Verified locally in the CSP sandbox using
`/opt/anaconda3/envs/csp_310/bin/python -B`. Initial verification used the working
tree; commit preparation also verified a clean export of the staged files, as
described below. This report does not imply a release or HPC calculation.

## Automated checks

| Check | Result |
| --- | --- |
| `python -m unittest discover -s TestCase/T4_zmatrix_viewer -p 'test_*.py' -v` | 21 tests passed, including the existing 94 BranchImproper Z-matrix fixtures |
| `python -m unittest discover -s TestCase/T3_zmatrix_generation -p 'test_*.py' -v` | 4 tests passed |
| `python -m unittest discover -s TestCase/T7_shared_zmatrix -p 'test_*.py' -v` | 3 tests passed; shared topology and Gaussian input compatibility |
| `python -m Source.CLI_scripts.zmatrix_viewer --help` | New CLI options available |
| `python -m pip wheel --no-deps --no-build-isolation . -w /tmp/csptoolbox-wheel` | Wheel built; three renderer assets and both viewer entry points present; imported and rendered from an extracted wheel |
| `node --check` on the emitted inline runtime | Passed |
| `git diff --check` | Passed |

Numeric checks cover Cartesian/pair-distance preservation, signed-torsion round
trips, one-to-one element-preserving maps, translated/rotated conformer alignment,
rejection of reflected fits, PDB models/repeated labels, SDF connectivity, periodic
unwrapping, special-position deduplication and source protection through path,
symlink, hard-link and changed-working-directory aliases.

During commit preparation, the 21 viewer tests and 4 generation tests also passed
against a clean export of the exact staged files and their tracked dependencies.
This excludes the unrelated local Gaussian, VASP and crystal-structure edits.
The three T7 tests above were additional working-tree checks; that local test
directory is outside these viewer commits.

## Chrome interaction checks

Used the ordinary Chrome browser against a loopback preview of the generated
files, without a Codex page runtime. Verified:

- Five-scene selection and the separate three-molecule CP2 gallery.
- Independent-DoF button selection updates its pressed state, value, interval and
  highlighted chain. Salicylic acid has 3 independent DoFs, 13 dihedrals and 42
  internal coordinates; the filter returned those counts.
- Clear selection clears pressed rows. Keyboard rotation/zoom changes the camera
  description; Reset view restores its original orientation and 100% zoom.
- Labels/H/reference toggles work; controls in one embedded viewer leave the other
  viewer's selection and visibility settings unchanged.
- The CP2 experimental/full-cell toggle changes salicylic acid from 16 to 64 atoms.
- The experimental/comparison overlay renders a faint grey reference and reports
  aligned all-atom RMSD 0.0367 Å for the archived salicylic-acid PDB pair.
- The atom-mapping table expands. Save mapped Z-matrix produced a downloaded file
  identical to the scene's numeric Z-matrix payload (1006 bytes).
- Normal and narrow layouts are legible; the narrow document had no horizontal
  overflow. A second embedded viewer rendered correctly with a dark colour scheme.
- No warning/error entries were returned by the Chrome console check.

The browser connector's download-event wait timed out, but the actual downloaded
file was present and its complete contents matched the expected payload. Direct
`file://` opening and physical mouse drag/wheel gestures were not automated; the
camera was exercised through its keyboard controls. Other browsers were not tested.
No automated browser runner is configured; the reproducible example generator and
manual checklist are in [VIEWER.md](VIEWER.md).

## Real data and boundaries

Generated review examples under ignored `Codex_workspace/viewer_review/`:
`viewer-workflows.html`, `two-viewers.html` and `cp2-test-molecules.html`.
The CP2 gallery preserves 16/15/19 atoms and three independent DoFs for salicylic
acid, nicotinamide and vanillin. Their source data are in the sibling
`CSPImperial/test/CP2` tree. A separate real SALIAC20 RES expansion also returned
16 ASU atoms and 64 full-cell atoms.

The CP2 comparison uses the archived site-order mapping and selects one whole
molecule from the repeated full-cell output; no geometry optimization or Gaussian
calculation was performed. The standalone synthetic demo also generated successfully.

Six unrelated previously modified source/template files retained their original
SHA-256 hashes. Existing user changes in README and `Source/__init__.py` were retained
alongside the viewer documentation/exports. Historical exports and external project
helpers were not rewritten. New exports should use the shared public API.

Remaining limits: unlabelled atom order must be trustworthy or explicitly mapped;
independent DoFs must be supplied; disorder/extended networks and ambiguous crystal
symmetry need separate validation. There is no implementation blocker for the five
supported workflows. Recommended follow-up is to regenerate older project viewers
through the new API and add a browser runner if unattended CI coverage is wanted.
