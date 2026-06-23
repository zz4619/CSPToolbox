# T2 Unit-Cell Reduction

This test case checks the reverse direction of `T1_unitcell_expansion`.

T1 validates:

```text
CCDC asymmetric CIF -> CSPToolbox explicit unit cell
```

T2 validates:

```text
explicit P1 full-unit-cell CIF -> CSPToolbox asymmetric unit -> explicit unit cell
```

## Test Data

- CCDC asymmetric experimental CIFs:
  `Experimental/`

- CCDC Python API full-unit-cell CIFs:
  `Experimental_FullUnitCell/`

There are `755` CIF files in each folder. The filename stem is the CSD
refcode and is used to pair structures between the two folders.

The full-unit-cell CIFs were generated with the CCDC Python API using:

```python
packed = crystal.packing(inclusion="UniqueIncluded")
```

They are explicit P1 unit-cell CIFs and include the CSPToolbox metadata tag:

```text
_csptoolbox_explict_unit_cell true
```

## Important Oracle Detail

Asymmetric-unit reduction is not unique. A valid reduced asymmetric unit does
not have to match the original CCDC asymmetric CIF atom-for-atom.

For example, structures with atoms on special positions can reduce to fewer
representative atoms than the CCDC asymmetric CIF, while still expanding back
to the same full unit cell. Therefore T2 must validate reduction by
round-tripping through expansion, not by direct comparison to `Experimental/`.

## Workflow

For each paired refcode:

1. Read `Experimental_FullUnitCell/REFCODE.cif` as the full-cell input.
2. Reduce it with `CrystalStructure.reduce_to_asymmetric_unit(symprec=0.05)`.
3. Check that the reduced structure is not marked as an explicit unit cell.
4. Check that space-group and SHELX symmetry metadata were detected.
5. Re-expand the reduced structure with `expand_to_explicit_unit_cell()`.
6. Compare the re-expanded structure against the original full-cell CIF.
7. Independently expand `Experimental/REFCODE.cif` with
   `CrystalStructure.expand_cif_to_unit_cell()`.
8. Compare the re-expanded reduced structure against this T1-style expansion.

## Comparison Metrics

The helper module compares structures using:

- atom-count equality
- element-count equality
- lattice-metric equality
- typed PDD distance with `k=100`, Chebyshev metric, and default row collapse
- same-element minimum-image Cartesian site matching

The same-element site matching deliberately ignores atom labels and atom order,
because the reduced asymmetric unit and its re-expansion may choose different
representatives from the original CCDC asymmetric CIF.

## Expected Result

For the full validation run:

- reduction succeeds for `755/755` full-cell CIFs
- each reduced structure contains at least one atom and no more atoms than the
  full-cell input
- each reduced structure stores detected space-group and symmetry metadata
- re-expansion restores the original full-cell atom count and composition
- typed PDD distances are numerically zero to practical tolerance
- same-element site matching gives maximum displacement below the selected
  tolerance

## Suggested Test Strategy

Keep normal test discovery fast by running a representative smoke set by
default. Suggested smoke refcodes:

```text
UMIQEO
ABALAS
DMSULO04
CYCYPR
TRIZIN01
```

These include ordinary symmetry, centrosymmetry, and special-position cases.

Run all `755` structures only when explicitly requested, for example through an
environment variable such as:

```bash
CSPTOOLBOX_FULL_T2=1 python -m unittest discover -s TestCase -p 'test*.py'
```

