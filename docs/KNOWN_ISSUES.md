# Known scientific issues

This file lists scientific defects that have been confirmed but not yet fixed.
Each has an ID (`KI-n`), used by the tests that pin it.

A test that pins an issue states the correct behaviour and is marked as a
strict expected failure. When a fix lands, that test passes unexpectedly
("XPASS"), and pytest counts that as a failure. In the same change, remove the
test's marker and move the entry to *Resolved*.

Found in the review of 4 October 2026. "Reproduced" means a test or script
showed the defect; "from code" means it was found by reading the code and has
not been run.

## Open

### KI-1: SHELX `LATT` comes from the space-group symbol, not the cell

- **Symptom:** A full cell is reduced to the asymmetric unit and re-expanded.
  If the cell is in a non-standard centred setting, the re-expanded cell has
  the wrong number of atoms:
  - I2/a cells labelled C2/c give 1.5× the atoms.
  - C-1 gives half the atoms.
  - R-3 on rhombohedral axes gives 3× the atoms.
- **Cause:** `shelx_symmetry_records_from_symmetry` takes the centring letter
  from spglib's standard symbol. The operations, however, are in the input
  cell's setting.
- **Affected:** `reduce_to_asymmetric_unit()` and therefore CSO-RM `.res`
  files. In CE755 this hits 9 structures: AQUSEM, AWACAE, CLPSUL02, CRYSEN,
  CUBANE, DERBEH, LELVIK, MEYCEY, SALMID04.
- **Evidence:** reproduced; pinned by `TestCase/T2_unitcell_reduction` (strict
  expected failures).
- **Fix direction:** derive `LATT` from the pure translations among the
  operations.

### KI-2: Atoms on special positions are not merged on re-expansion

- **Symptom:** CEKGUU01 (P-3) and TETDAM03 (P6₃/m) re-expand with extra atoms.
- **Cause:** `expand_shelx_atoms` merges copies only when they agree to 1e-8 in
  fractional coordinates. Symmetry detection uses `symprec=0.05`, so atoms
  slightly off a special position are not merged.
- **Evidence:** reproduced; pinned by `TestCase/T2_unitcell_reduction`.

### KI-3: Positive `LATT` is written when the inversion centre is off the origin

- **Symptom:** In origin-choice-1 settings, the written operations generate the
  wrong positions: `I 41/a` gives 24 positions instead of 16, and `P 4/n`,
  `P n n n` and `F d d d` are wrong as well.
- **Cause:** SHELX `LATT>0` means "add -x,-y,-z". When the inversion is at, for
  example, (1/4,1/4,1/4), all operations must be written out with `LATT<0`
  instead.
- **Affected:** any structure whose symmetry is looked up by symbol. Lookup
  picks the first Hall setting, which is origin choice 1 for these groups. Also
  affects detection on input cells in origin choice 1.
- **Evidence:** reproduced; pinned by `TestCase/T13_known_issues`.

### KI-4: CSO-FM `.res` loses the cell setting

- **Symptom:** A P2₁/n crystal is written with P2₁/c operations, which describes
  a different crystal. For ABEROP the result is 3.45 Å away from the real
  structure.
- **Cause:** `CSOFMInputBuilder.write_job_from_crystal` rebuilds the structure
  without `shelx_latt_value`, `symmetry_operations` and `hall_number`. The
  writer then looks the operations up from the standard symbol.
- **Affected:** every CSO-FM input written for a cell that is not in the
  standard setting, including those from
  `Source/CLI_scripts/generate_csofm_from_contcar.py` (CONTCAR → asymmetric
  unit → CSO-FM). 152 of the 755 CE755 structures are deposited in P2₁/n.
- **Evidence:** reproduced; pinned by `TestCase/T13_known_issues`.

### KI-5: A multi-block CIF mixes blocks

- **Symptom:** The atoms come from the first data block, but the cell comes from
  the last one.
- **Cause:** `_parse_cif_blocks` merges every block's data items into one
  dictionary.
- **Evidence:** reproduced; pinned by `TestCase/T13_known_issues`.

### KI-6: CSO-FM `mol.input` counts are wrong for Z′ > 1

- **Symptom:** `Nmol_asm` is set to the number of unique species, and every
  species is written with `Number of molecules 1`. CSO-FM expects the total
  number of molecules in the asymmetric unit, and per-species counts that sum
  to it.
- **Affected:** any Z′ > 1 structure, including a dihydrate's two waters. The
  QICVAZ01 snapshot in T11 records the current output.
- **Evidence:** reproduced; pinned by `TestCase/T13_known_issues`.

### KI-7: Fixed Z-matrix coordinates in `.com` files are lost in PyZMAT

- **Symptom:** Since 991a225, fixed coordinates are written as numbers directly
  in the atom row.
  - PyZMAT's parser reads such a number as an unknown variable name and turns
    it into `None`. It also builds constraints only from a `Constants:` block,
    so the constraints are lost.
  - The committed PyZMAT also needs a `Variables:` heading.
- **Affected:** the handoff from CSPToolbox `.com` files to PyZMAT and ml-lams.
  This is a cross-repository contract; see [FORMATS.md](FORMATS.md).
- **Evidence:** reproduced with a script during the review. There is no test in
  this repository, because PyZMAT is not a dependency.

### KI-8: The gas-phase energy is accepted without checks

- **Symptom:** `read_gaussian_final_energy` returns the last `SCF Done` value,
  even from a job that ended in an error or did not converge.
- **Related issues:**
  - One Gaussian energy is used for every non-water molecule.
  - The water energy is hard-coded (`DEFAULT_WATER_GAS_PHASE_ENERGY`) with no
    record of its level of theory.
- **Evidence:** reproduced for the termination check; pinned by
  `TestCase/T13_known_issues`.

### KI-9: Partial occupancies are kept silently

- **Symptom:** `expand_cif_to_unit_cell` keeps every disorder component, so
  overlapping atoms reach molecule detection. Only
  `inspect_cif_unit_cell_expansion` reports the occupancies.
- **Evidence:** reproduced; pinned by `TestCase/T13_known_issues`.

### KI-10: CSO-FM `Zmatrix` may not match the `.res` atom order (from code)

`render_csofm_zmatrix_text` runs a second Z-matrix generation on the reordered
atoms. Z-matrix generation depends on geometry and atom index, so it is not
guaranteed to match the order used to write `.res`.

### KI-11: The Mie typing bond table differs from CSO-RM (from code)

`mie_typing.CSORM_BOND_DISTANCE` uses 1.5 Å for C–O, while CSO-RM
(`CSP_FUNCTIONS_mod.f90`) uses 1.62 Å. Molecule detection elsewhere uses a
third rule: ASE covalent radii.

### KI-12: One charge is written for every molecule in CSO-RM `.info` (from code)

`CSORMSettings.charge` is written for each molecule, which is wrong for salts.

### KI-13: The PDB writer truncates atom labels to 4 characters (from code)

Distinct labels can collide (for example `C10_12` and `C10_13` both become
`C10_`).

### KI-14: `CrystalStructure.from_file` ignores a CIF's symmetry operations (from code)

Expanding such a structure looks the setting up from the Hermann-Mauguin
symbol, which can choose the wrong setting (see KI-3). In contrast,
`expand_cif_to_unit_cell` uses the CIF's own operations.

### KI-15: Z-matrices of symmetry copies use the first graph isomorphism (from code)

`generate_zmatrices` maps the template onto each copy of a molecule with the
first isomorphism networkx returns. For atoms that are topologically equivalent
(for example methyl hydrogens), this can permute which atom gets which row.

## Resolved

| ID | Issue | Resolution |
|---|---|---|
| — | Typed PDD comparison required equal atom counts, so equivalent cells of different size could not be compared | Element fractions are compared since the PDD clustering merge (a151a9b) |
| — | CIF symmetry operations were evaluated with `eval()` | Replaced with an arithmetic-only evaluator that gives bit-identical results (641e90e) |
