# T3 Z-Matrix Generation

This test case exercises CSPToolbox Z-matrix generation on small synthetic
molecules with known bonding patterns.

The goal is to test both the current generator and the classification logic
needed for the planned user-guided method where selected dihedrals remain
proper torsions and all other torsions should be improper torsions.

## Test Data

The structures are generated programmatically in
`zmatrix_generation_helpers.py` so the test case is small and deterministic.

Current synthetic cases:

- `butane_chain`: four carbon atoms in a non-linear bonded chain. This should
  produce at least one proper heavy-atom torsion over `C1-C2-C3-C4`.
- `methanol`: one methyl group attached to oxygen. This checks the current
  special improper handling for methyl hydrogens.

## What Is Checked

The test suite validates that:

- each synthetic molecule produces exactly one Z-matrix,
- every atom appears exactly once,
- every Z-matrix reference points to an earlier row,
- bond lengths, angles, and dihedrals are finite where required,
- generated references can be classified as proper, improper, or fallback
  using the molecule bonding graph,
- methanol methyl hydrogens after the first hydrogen use improper references,
- the butane heavy-atom chain produces a proper torsion.

## Run Tests

From the CSPToolbox repository root:

```bash
/opt/anaconda3/envs/csp_310/bin/python -B -m unittest discover \
  -s TestCase/T3_zmatrix_generation -p 'test*.py'
```

## Generate Report Files

The report runner writes example `.zmat` files and summary reports under this
same folder:

```bash
/opt/anaconda3/envs/csp_310/bin/python -B \
  TestCase/T3_zmatrix_generation/run_zmatrix_generation_report.py
```

Generated output:

- `Generated/*.zmat`
- `t3_zmatrix_generation_report.csv`
- `t3_zmatrix_generation_report.md`

## Export Labelled CONTCAR Diagrams

Use HA-pair TPSS converged VASP structures as visual starting points for
selecting user-defined proper dihedral angles:

```bash
/opt/anaconda3/envs/csp_310/bin/python -B \
  TestCase/T3_zmatrix_generation/export_contcar_labelled_diagrams.py --draw-box
```

The script reads `SYSTEM/CONTCAR` files from:

```text
/Users/zianzhan/Desktop/CSP_sandbox/CSP-personal/Results/VASP_HA_pair_TPSS_converged_structure
```

Generated output is written to:

```text
TestCase/T3_zmatrix_generation/HA_pair_TPSS_CONTCAR_labelled_diagrams
```

Output includes:

- `unit_cells/*_unit_cell.cif`
- `unit_cell_diagrams/*_unit_cell_labelled.png`
- `molecule_diagrams/*_molecule_NN_labelled.png`
- `labelled_diagram_manifest.csv`
