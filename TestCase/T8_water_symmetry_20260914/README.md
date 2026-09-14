# Water Z′=16 symmetry validation, 14 September 2026

All four returned water minima retain **P1**, with **16 inequivalent waters**.
Their starting geometries are also P1. Neither CSPToolbox nor PLATON finds a
higher-symmetry oxygen framework or a smaller primitive cell within the tested
criteria. These are coordinate-based conclusions for the saved structures; they
do not establish experimental phase identity, thermodynamic stability or search
completeness.

| Native record | Final energy, kJ/mol water | Final full O/H group | Final oxygen group | Primitive-cell waters |
| --- | ---: | --- | --- | ---: |
| 1 | -34.1442 | P1 | P1 | 16 |
| 2 | -33.3346 | P1 | P1 | 16 |
| 3 | -33.7482 | P1 | P1 | 16 |
| 4 | -33.3319 | P1 | P1 | 16 |

The CP2 search ran in P1 and timed out with four converged returns; it is not a
completed 10,000-assignment search. The native row order pairs the two input
tables. The starting table records the geometry passed to minimization, after
screening/repair; its leading energy field is the returned final energy, not an
initial-state energy. Never infer an optimization energy change from that field.

## Sources and reconstruction

- Frozen CPU reference: CSPImperial branch `codex/cp2-high-z-screening-benchmarks`,
  reference commit `f5fd0a11e15afdf6946c6da5c50d25f1622f20cd`, scientific source
  `2fb1f4c0f2a585ea2eb055c2cd04cda2b2fe6aaf`, job `4050226.pbs-7`.
- Archive entry: `CSP_studies/cp2_high_z_screening_20260914_licence_fix/CPU_BASELINE_ARCHIVE_README.md`
  relative to the multi-repository workspace. Original archive files are unchanged.
- CSPToolbox library: clean commit `a40c80c759a8007aaa1afd13e880e2fee4c1dbf1`;
  this study changes no library or CP2 behaviour and excludes the unrelated dirty
  work in the original CSPToolbox checkout.
- Interpreter used: `/opt/anaconda3/envs/csp_310/bin/python`, Python 3.12,
  NumPy 2.4.1, ASE 3.27.0 and spglib 2.7.0. These dependencies were verified locally.

`analyse_water.py` reconstructs the three-atom rigid reference using CP2's
`build_molecule.f90`, `rot_matrix_from_Euler_angles.f90` and
`unit_cell_module.f90` conventions. Stored `zposnom` is the **fractional origin of
the first atom**, not the molecular centre of mass. Cell and Euler angles in the
native tables are radians; exported CIF angles are degrees. The cell matrix in
the NumPy data uses row vectors. Molecule/atom order and the two slightly different
O–H lengths are retained. Geometric symmetry allows same-element H permutations;
the source model has slightly unequal H charges, so coordinate symmetry alone
would not prove Hamiltonian symmetry under an H exchange.

Final native coordinates are printed to four decimal places and starting
parameters to six. Exporting more digits does not recover lost precision.
All **eight start/final structures** were independently exported through CP2's
native geometry routines. Direct comparisons against the native ASU coordinates
in the original atom order agree
within **4.45e-15 Å** for Cartesian components and **1.78e-15 Å** for cell
components, without alignment, wrapping or permutation. The compact native
coordinates, source hashes and results are retained in
`results/native_export_all.json`, `results/native_export_provenance.json` and
`results/native_export_crosscheck.json`; `crosscheck_export.py` repeats this check.
In P1 the ASU contains all 48 atoms. CP2's separate full-cell output moves some
molecules by whole lattice vectors; the same script checks that these positions
also agree modulo integer lattice translations, with no fitting or permutation.

Additionally, the unrounded rigid reconstruction matches two accepted native
replay coordinate sets with maximum periodic atom error **4.45e-15 Å**. Those
native replays correspond to starting rows 1 and 2. Reading the rounded starting
table instead raises the maximum difference to **8.69e-6 Å**, consistent with
printing precision. The selected native reference and its original source hash
are retained in `results/native_selected.json`; `crosscheck_native.py` repeats
this independent check. These are geometry checks, not minimization reruns.

## CSPToolbox reduction

The study calls the existing `detect_space_group_symmetry`,
`reduce_to_asymmetric_unit` and `expand_to_explicit_unit_cell` APIs. It also uses
spglib's primitive-cell search on the same geometry. Tests cover 0.0001, 0.001,
0.005, 0.01, 0.02, 0.05, 0.1 and 0.2 Å at spglib's automatic angular tolerance,
for each of eight geometries. Full O/H, oxygen-only and the auxiliary point-label
control all return P1 at every setting. The primitive cell retains 48 O/H atoms,
or 16 atoms for the oxygen-only framework.

Every returned operation is applied to the original coordinates and checked by
species-preserving periodic one-to-one assignment. Here the identity is the only
operation found. All 192 reduction/expansion checks preserve element counts.
The largest round-trip displacement is **1.02e-7 Å**: CSPToolbox's existing
expansion canonicalizes fractional coordinates to eight decimal places. This
measured rounding is far below the tested tolerances; no library fix was made.

`results/structures/` contains original full-cell CIFs, oxygen-only CIFs and
the existing API's ASU RES output at 0.05 Å. Since the group stays P1, the full
O/H ASU still contains all 16 waters. No symmetrized coordinates are substituted
for the originals.

## Independent PLATON check

The official [PLATON download](https://www.platonsoft.nl/xraysoft/Mac-OSX/platon/)
provides `platon_special.f`, a command-line build without X11. Its source banner
is PLATON-2026, internal version `200626`; the downloaded source is unchanged.
The complete build provenance and hashes are in `results/installation.json`.

PLATON is installed on CX3 at:

```text
$HOME/CSP_studies/water_z16_symmetry_20260914/software/platon/platon
```

The build command there was:

```bash
gfortran -O0 -w -o platon platon_special.f
```

GNU Fortran 8.5.0 compiled it successfully. A Mac attempt installed Homebrew GCC
16.1.0, but the Mac lacks the Apple SDK and could not link PLATON. The working
installation and all reported PLATON runs are on CX3; local native PLATON remains
unavailable until an Apple SDK is installed. No graphical interface is required.

For each structure, `run_platon.py` passes this instruction through standard input
to `platon -o input.cif`, followed by `END`:

```text
CALC ADDSYM EXACT angle distance distance distance
```

The [ADDSYM manual](https://platonsoft.nl/platon/pl000401.html) defines the angle
criterion and separate rotational, inversion and translational distance criteria.
The matrix uses angles 0.01, 0.1 and 1 degree and distances 0.001, 0.005, 0.01,
0.02, 0.05, 0.1 and 0.2 Å. These settings are documented separately from spglib;
equal numeric tolerances do not make the algorithms equivalent. `EXACT` means
zero permitted unmatched atoms within the specified tolerances, not zero
coordinate deviation.

**Ordinary ADDSYM excludes hydrogen atoms, including in EXACT mode.** The
downloaded source explicitly implements this in `PLA060`, and every ordinary
water run reports 16 included atoms. Therefore ordinary PLATON validates the
oxygen framework only. An auxiliary copy labels H positions as C and retains O
positions as O, allowing ADDSYM to process all 48 labelled points. Those files
are geometric tests, never chemical models or energy inputs. The proxy alone
does not establish full molecular symmetry.

All **504 water checks** ended normally and returned P1. The **126 control
checks** also passed:

| Constructed control | Full O/H via CSPToolbox | Ordinary PLATON (16 O atoms) | PLATON O/C point proxy (48 points) |
| --- | --- | --- | --- |
| 8 waters plus exact inversion partners | P-1 | P-1 | P-1 |
| Same O lattice, one water's H positions rotated | P1 | P-1 | P1 |

Applying the oxygen inversion directly to the original O/H positions gives
maximum residual 0 Å in the positive control and **0.952 Å** in the H-broken
control. This demonstrates why framework symmetry cannot be promoted to full
O/H symmetry. Actual water records have no additional proposed framework
operations to test. The controls are mathematical coordinate fixtures, not
energetically validated water crystals.

The first control run exposed the analysis parser's separate PLATON output format
for a detected symmetry increase. The parser was corrected and all controls were
rerun. The initial log is preserved in the external evidence archive; PLATON's
source and the native structures were not changed.

## Reproduce and verify

From the repository root, using an interpreter with the project dependencies:

```bash
python TestCase/T8_water_symmetry_20260914/crosscheck_export.py
python TestCase/T8_water_symmetry_20260914/crosscheck_native.py
MPLCONFIGDIR=/tmp/csp-water-symmetry-mpl python TestCase/T8_water_symmetry_20260914/analyse_water.py \
  --run TestCase/T8_water_symmetry_20260914/inputs --output /tmp/water-symmetry-fresh
MPLCONFIGDIR=/tmp/csp-water-symmetry-mpl python TestCase/T8_water_symmetry_20260914/prepare_controls.py \
  --output /tmp/water-symmetry-controls-fresh
python TestCase/T8_water_symmetry_20260914/run_platon.py \
  --platon /absolute/path/to/platon --inputs /tmp/water-symmetry-fresh \
  --output /tmp/water-platon-fresh
python TestCase/T8_water_symmetry_20260914/run_platon.py \
  --platon /absolute/path/to/platon --inputs /tmp/water-symmetry-controls-fresh \
  --output /tmp/water-platon-controls-fresh
(cd TestCase/T8_water_symmetry_20260914/inputs && shasum -a 256 -c SHA256SUMS)
(cd TestCase/T8_water_symmetry_20260914/results && shasum -a 256 -c SHA256SUMS)
git diff --check
```

All output paths must be fresh. The full PLATON logs, control outputs, inputs,
compiler output and executable are preserved in the external workspace study
`CSP_studies/water_z16_symmetry_20260914/platon_evidence.tar.gz`, also present
under the same study name on CX3. SHA-256:
`cf90ba4ecd42f2097fba947e76e59c6f574fa673c47cf3dbb5b43df7e2e4e820`.
The local copy's hash, 504 water results, 126 controls and included-atom counts
were verified. Compact JSON results and coordinate files are committed here;
the third-party PLATON source/executable stay in the external study.

Done when the input hashes, native geometry comparison, reduction round trips,
PLATON normal termination, tolerance matrix and positive/H-broken controls pass.
Do not interpret search SG labels, PLATON's non-H result, or a loose-tolerance
pseudosymmetry alone as a confirmed full-crystal symmetry increase.
