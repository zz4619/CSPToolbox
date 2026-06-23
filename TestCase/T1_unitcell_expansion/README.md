# T1 Unit-Cell Expansion

This test case checks whether CSPToolbox expands CCDC asymmetric experimental CIFs to the same explicit P1 full unit cells as the CCDC Python API.

The data were taken from the CE755 experimental structure set.

## Test Data

- CCDC asymmetric experimental CIFs:
  `Experimental/`

- CCDC Python API full-unit-cell CIFs:
  `Experimental_FullUnitCell/`

There are `755` CIF files in each folder. The filename stem is the CSD refcode and is used to pair structures between the two folders.

The full-unit-cell reference CIFs were generated with the CCDC Python API by reading each CSD entry as a `ccdc.crystal.Crystal` and calling:

```python
packed = crystal.packing(inclusion="UniqueIncluded")
```

The resulting packed structure was then written as an explicit P1 CIF. These files therefore represent CCDC's `Crystal.packing(inclusion="UniqueIncluded")` full-cell expansion convention.

Relevant CCDC API signature:

```python
packing(box_dimensions=((0, 0, 0), (1, 1, 1)), inclusion="CentroidIncluded")
```

This returns a molecule which fills some multiple of the unit cell of the crystal. The atoms to include are specified by the `inclusion` argument:

- `CentroidIncluded`: whole molecules are included if their centroid is within the box dimensions.
- `AllAtomsIncluded`: whole molecules are included only if all atoms of the molecule lie within the box.
- `AnyAtomIncluded`: whole molecules are included if any atom of the molecule lies within the box.
- `OnlyAtomsIncluded`: all and only the atoms lying within the box are included.
- `UniqueIncluded`: whole molecules are included if their centroid is within the box dimensions and they contribute unique box cell positions.

This test uses `inclusion="UniqueIncluded"` with the default unit-cell box dimensions.

## Workflow

1. Read each CIF in `Experimental/` with CSPToolbox.
2. Normalize isotope labels `D` and `T` to `H` for PDD compatibility.
3. Expand atom sites to an explicit P1 unit cell using CSPToolbox.
4. Deduplicate near-coincident same-element symmetry images using the new CIF site-merge tolerance in CSPToolbox.
5. Compare each CSPToolbox-expanded CIF against the matching CIF in `Experimental_FullUnitCell/`.
6. Calculate typed PDD distances with `k=100`, Chebyshev metric, and default PDD row collapse.
7. Run an additional same-element site-matching check using minimum-image Cartesian distances.

No extra missing hydrogens should be added for this test. Adding hydrogens only on the CSPToolbox side would no longer test whether CSPToolbox expansion matches the existing CCDC full-cell export.

## Expected Result

- CSPToolbox expansion succeeds for `755/755` structures.
- Typed PDD succeeds for `755/755` structures.
- No composition mismatches are found.
- No lattice metric mismatches are found.
- No hydrogen-only differences are found.
- Typed PDD distances are numerically zero to practical tolerance.
- Same-element site matching gives maximum displacement below `1e-6 A`.

In the validation run used to create this test case:

- Mean typed PDD was `2.93e-10`.
- Maximum typed PDD was `8.54e-08` for `TRIZIN01`.
- Maximum matched atom displacement was `6.61e-08 A` for `CYCYPR`.
- No site-matching displacement exceeded `1e-6 A`.

## Notes

The important edge cases are structures where atoms lie on special positions but coordinates are stored as rounded decimals. Without a Cartesian site-merge tolerance, symmetry expansion can create near-duplicate atoms separated by only about `1e-4 A` or less. The CSPToolbox CIF expansion path should therefore use same-element, minimum-image site merging with the default `0.01 A` tolerance.

Conclusion from the validation run: with the near-duplicate site-merge tolerance, CSPToolbox unit-cell expansion matches the CCDC Python API full-unit-cell export for these CE755 experimental structures to numerical precision.
