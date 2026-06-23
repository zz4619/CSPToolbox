# T2 Unit-Cell Reduction Report

Generated from `755` paired T2 CIF files.

## Configuration

- symprec: `0.05`
- PDD: typed, `k=100`, Chebyshev metric, default row collapse
- PDD tolerance: `1e-06`
- same-element site-matching tolerance: `1e-06 A`
- elapsed time: `107.2 s`

## Summary

- paired refcodes: `755`
- passed: `744`
- failed or errored: `11`
- mean PDD to full-cell reference: `2.974e-10`
- max PDD to full-cell reference: `8.542e-08` (TRIZIN01)
- mean PDD to T1 expansion: `2.026e-15`
- max PDD to T1 expansion: `4.552e-15` (DMEGLY01)
- max site displacement to full-cell reference: `6.612e-08 A` (CYCYPR)
- max site displacement to T1 expansion: `6.890e-15 A` (NAPOAC01)
- mean reduced/full atom-count ratio: `0.279`
- smallest reduced/full atom-count ratio: `0.050` (TCYETY01)

## Most Common Detected Space Groups

- `P 2_1/c`: `391`
- `P -1`: `90`
- `P 2_12_12_1`: `77`
- `P bca`: `36`
- `P 2_1`: `28`
- `C 2/c`: `27`
- `P na2_1`: `20`
- `P ca2_1`: `16`
- `P nma`: `13`
- `P c`: `8`
- `C c`: `6`
- `P bcn`: `5`
- `P 2_1/m`: `3`
- `P 2/c`: `2`
- `P 2_12_12`: `2`

## Failures

- `AQUSEM`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.
- `AWACAE`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.
- `CEKGUU01`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.
- `CLPSUL02`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.
- `CRYSEN`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.
- `CUBANE`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.
- `DERBEH`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.
- `LELVIK`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.
- `MEYCEY`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.
- `SALMID04`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.
- `TETDAM03`: error - ValueError: PDD comparison requires matching element counts in the explicit unit cell.

## Failure Diagnostics

The 11 failures are composition round-trip failures, not small numerical mismatches. In each case, the reduced structure re-expands to a different atom count/composition than both the CCDC full-cell CIF and the independent T1-style expansion of the original asymmetric CIF.

| Refcode | Asym | Full | Reduced | Re-expanded | T1 Expanded | Reduced SG | Diagnosis |
|---|---:|---:|---:|---:|---:|---|---|
| AQUSEM | 21 | 168 | 21 | 252 | 168 | C 2/c | reexpanded_composition_mismatch |
| AWACAE | 16 | 64 | 16 | 32 | 64 | P -1 | reexpanded_composition_mismatch |
| CEKGUU01 | 13 | 26 | 5 | 30 | 26 | P -3 | reexpanded_composition_mismatch |
| CLPSUL02 | 25 | 100 | 13 | 152 | 100 | C 2/c | reexpanded_composition_mismatch |
| CRYSEN | 30 | 120 | 15 | 180 | 120 | C 2/c | reexpanded_composition_mismatch |
| CUBANE | 16 | 16 | 4 | 48 | 16 | R -3 | reexpanded_composition_mismatch |
| DERBEH | 16 | 64 | 8 | 96 | 64 | C 2/c | reexpanded_composition_mismatch |
| LELVIK | 18 | 144 | 18 | 216 | 144 | C 2/c | reexpanded_composition_mismatch |
| MEYCEY | 21 | 42 | 11 | 64 | 42 | C 2 | reexpanded_composition_mismatch |
| SALMID04 | 17 | 136 | 17 | 204 | 136 | C 2/c | reexpanded_composition_mismatch |
| TETDAM03 | 20 | 40 | 4 | 48 | 40 | P 6_3/m | reexpanded_composition_mismatch |

Detailed count/composition diagnostics are in `t2_unitcell_reduction_failure_diagnostics.csv`.

## Output Files

- CSV details: `t2_unitcell_reduction_report.csv`
- Failure diagnostics: `t2_unitcell_reduction_failure_diagnostics.csv`
- Markdown summary: `t2_unitcell_reduction_report.md`
