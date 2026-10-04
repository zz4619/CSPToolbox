# Interchange formats

This page lists the files that CSPToolbox writes or reads to exchange data with
PyZMAT, ml-lams, CSPImperial (CSO-RM, CSO-FM, CP2) and external codes. Each
entry says who writes the file and who reads it.

A change to any format below is a change to another repository's input. Follow
[Changing a format](#changing-a-format).

## Z-matrices

All the Z-matrix formats share one geometric convention:

- Row *i* places atom *i* by a bond to an earlier row, an angle to a second
  earlier row, and a dihedral to a third.
- The dihedral *i–b–a–d* uses the IUPAC sign convention.
- Units are Å and degrees.

CSPToolbox and PyZMAT reproduce identical dihedrals from the same file;
`Source/crystal/zmatrix_builder.py` and `Source/zmatrix_viewer/geometry.py`
implement the convention.

The formats differ in four ways: where the numbers sit, how references are
counted, whether atom labels are kept, and how fixed coordinates are written.

### A. `# ZMAT v1` (CSPToolbox numeric Z-matrix): the master copy

```text
# ZMAT v1
# title: SUCANH12_Mol1
# labels: C1_1 C2_1 C4_1 O1_1 O3_1 H3_1 H4_1 H1_1 H2_1 C3_1 O2_1
# bonds: 1-2 2-3 3-4 3-5 2-6 2-7 1-8 1-9 1-10 10-11
# torsion 0: 4-1-2-3 bond 4-3
C
C 1 1.5246077796
C 2 1.5022783103 1 104.3396878843
O 3 1.3865489907 2 110.2380948966 1 -1.6481554794
...
```

- Numbers sit in the atom rows, to 10 decimals.
- References count from 1.
- Labels, bonds and torsions are listed in `#` header lines. A
  `# torsion k: i-d-a-b bond i-b` line lists row *i*, then its dihedral, angle
  and bond references. This is the reverse of the row's *i–b–a–d* order, so
  read it carefully.
- The parser also accepts rows of the form `label element ...`.
- Written by `GaussianInputBuilder.render_zmat_text()`.
- Read by `Source.zmatrix_viewer.load_zmatrix()`; the `.com` writer, the viewer
  and CSPImperial Stage 01/01b all read through it.
- PyZMAT reads it with `ZMatrix.load_from_csp_zmat()`. That reader is not
  committed upstream yet, and it ignores the label, bond and torsion headers.

### B. Gaussian `.com` with named variables

Atom rows refer to names (`bnd3`, `ang3`, `dih4`), and the numbers follow the
rows in a separate block. Four versions are in use:

| Version | Written by | Numbers block | Fixed coordinates | Decimals |
|---|---|---|---|---|
| B1 | CSPToolbox before 991a225 | `Variables:` heading, then `bnd3=1.09`; `EPS=` listed under `Constants:` | none | 6 |
| B2 | CSPToolbox from 991a225 (`render_com_text`) | `bnd3=1.09` with no heading; `EPS=` in its own section after a blank line | written as a number directly in the atom row | 6 |
| B3 | CSPImperial Stage 01 (`_render_ml_lams_source`), as ml-lams input only | `Variables:` heading, then `bnd3 1.09` | none | 10 |
| B4 | PyZMAT `save_gaussian_com` (ml-lams Gaussian single point) | `Variables:` heading, then tab-separated values; dihedrals in (−180, 180] | listed under `Constants:` | 6 |

The only reader is PyZMAT's `ParseUtils.parse_gaussian_input`, which feeds
ml-lams:

- The committed PyZMAT requires the `Variables:` heading.
- An uncommitted local version also accepts B2, but it turns B2's inline fixed
  values into `None`, which is KI-7.

Treat `.com` as a Gaussian input, not as the interchange file. Write it from a
saved `# ZMAT v1` file with `write_com_from_zmat_file()`. Its 6 decimals lose
precision relative to `# ZMAT v1`.

### C. ORCA `* gzmt` block with `%geom` constraints (PyZMAT and ml-lams)

Numbers sit in the atom rows, and references in the `gzmt` block count from 1.
Constraint lines such as `{ B 3 0 C }` count atoms **from 0**. PyZMAT writes
these files (`save_orca_*`) and reads them back (`load_from_orca_input`,
`load_from_orca_output`).

### D. PyZMAT JSON

Used by ml-lams' Hessian stage. It stores the Z-matrix, 0-based connectivity,
energy, forces and Hessian, and is written and read only by PyZMAT.

### E. Label `Zmatrix` (CSPImperial Fortran)

Connectivity only, given as atom labels; there are no values.

```text
Z-matrix for molecule 1
C1
C2    C1
C3    C1    C2
O1    C3    C1    C2
```

- Written by `CSOFMInputBuilder.render_csofm_zmatrix_text()` and by CSPImperial
  runs; one is stored in each CP2 global-search directory.
- Read by CSO-FM (`main_SetupUnitCell_flexible`, `main_LAM_Integrity`),
  CrystalSetup (`extract_molecules`) and `cp2_local_min.parse_cp2_canonical_zmatrix`.

### Who uses which

| Format | CSPToolbox | PyZMAT | ml-lams | CSPImperial |
|---|---|---|---|---|
| A `# ZMAT v1` | writes, reads | reads (not committed upstream) | — | Stages 01/01b, through CSPToolbox |
| B `.com` | writes B2 | reads B1–B4, writes B4 | reads through PyZMAT | writes B3; strips B1's headings |
| C ORCA | — | writes, reads | uses | — |
| D JSON | — | writes, reads | uses | — |
| E label `Zmatrix` | writes (CSO-FM), reads (CP2) | — | — | Fortran reads |

## Crystal structures and program inputs

### SHELX `.res` (CSO-RM, CSO-FM, round trips)

Written by `CrystalStructure.to_file(..., fmt="res")` (`Source/crystal/res.py`):

| Record | Content |
|---|---|
| `TITL` | structure name |
| `REM SPACE_GROUP` | Hermann-Mauguin label |
| `REM EXPLICT_UNIT_CELL` | `TRUE` when the atoms fill the cell |
| `CELL` | 1.54184 (placeholder wavelength), then a, b, c, α, β, γ to 10 decimals |
| `LATT` | centring code 1–7 (P, I, R, F, A, B, C); positive means centrosymmetric with the inversion centre at the origin |
| `SYMM` | the remaining operations; `rounding=True` snaps translations within 0.01 of 0, ¼, ⅓, ½, ⅔ or ¾ |
| `SFAC` | elements in order of first appearance |
| atom rows | `label sfac x y z 11.00000 0.05000` (fractional to 10 decimals; occupancy and Uiso are placeholders) |

A reduced structure keeps the `LATT`/`SYMM` records it was reduced with, and
reading a `.res` file back keeps them as well. See KI-1, KI-3 and KI-4 for
current defects.

### CSO-RM job (`CSORMInputBuilder`)

| File | Content |
|---|---|
| `<name>.res` | molecules in descending molecular weight, atoms in detection order |
| `<name>.info` | `Num_Mols_Asm`, atoms per molecule, `exp_Ulatt`, `Pressure`, and one `scal_disp`/`scal_elec`/`Charges`/`Multiplicity` value per molecule |
| `CSO_RM.input` and `.inter` | copied from `Template/CSORM_input/` (`DUMMEY_SYSTEM_NAME`/`SYSTEM_NAME` and `RD filename` substituted) |

Atom labels must type as FIT/Mie sites present in the `.inter` file
(`mie_typing`).

### CSO-FM job (`CSOFMInputBuilder`)

| File | Content |
|---|---|
| `<name>.res` | atoms ordered by each molecule's Z-matrix |
| `mol.input` | `Nmol_asm`, `Nmol_unique`, then one block per species with charge, multiplicity, gas-phase energy (Hartree) and LAM database paths (see KI-6) |
| `Zmatrix` | format E |

### CP2 local minimisation (`cp2_local_min`)

`expcrys.pdb` uses fixed columns in CP2's atom order, alongside a copy of the
system's `input.in`. `zmatrix_mapping.json` records the atom mapping, scores and
provenance (`MAPPING_METHOD_VERSION`). The CP2 LAM database and the system-level
`Zmatrix` are authoritative and are read, never written.

### VASP (`VaspInputBuilder`)

`POSCAR` (Cartesian Å, atoms grouped by element in order of first appearance),
`KPOINTS`, `INCAR`
from `Template/VASP_input/`, `POTCAR` concatenated from the configured POTCAR
library, and the run script.

## Changing a format

1. Update this page and the writer in the same change.
2. Regenerate and review the T11 snapshots
   (`CSPTOOLBOX_UPDATE_SNAPSHOTS=1 python -m pytest TestCase/T11_input_snapshots`).
3. List the consumers from the tables above and check each one:
   - CSPImperial `workflow/stages/01_molecular_model` pins a CSPToolbox commit.
   - PyZMAT and ml-lams are maintained separately (Yixuan Huang); agree the
     change there before relying on it.
4. Record the change in `CHANGELOG.md` under the format's name.
