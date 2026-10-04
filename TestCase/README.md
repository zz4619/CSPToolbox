# Test cases

Each `T<n>_<topic>/` folder holds one case: its tests, plus any fixtures, with a
README when the case has data. `testsupport/` holds helpers shared by the cases,
which tests import as `testsupport.<module>`.

Run everything from the repository root:

```bash
python -m pytest            # every fast test
python -m pytest -m slow    # full-dataset regressions (T1/T2)
```

Unittest-style cases (T3–T9) also still run with
`python -m unittest discover -s TestCase/<case> -p 'test_*.py'`.

| Case | What it checks | Data | Tier |
|---|---|---|---|
| T1_unitcell_expansion | CIF asymmetric unit → P1 cell matches the CCDC-packed cell site by site | CE755: 755 CCDC CIF pairs | smoke + slow |
| T2_unitcell_reduction | Full cell → asymmetric unit → full cell round trip; KI-1/KI-2 are expected failures | CE755 pairs | smoke + slow |
| T3_zmatrix_generation | Z-matrix generation, branch-improper references, labelled diagrams | HA systems, API Z-matrices | fast |
| T4_zmatrix_viewer | `# ZMAT v1` parsing, structure views, HTML export | in-test fixtures | fast |
| T5_vasp_results | OUTCAR force parsing | synthetic OUTCAR | fast |
| T6_conformer_generation | Seeded RDKit conformer ensembles | SMILES | fast |
| T7_shared_zmatrix | One Z-matrix topology across conformers; `.com` handoff | in-test fixtures | fast |
| T8_pdd_clustering | PDD descriptors and indexed clustering (serial and parallel) | nicotinamide RES fixtures, benchmark summaries | fast |
| T9_cp2_local_min_input | CP2 experimental local-minimisation inputs, atom mapping, batches | synthetic CP2 trees | fast |
| T10_crystal_units | Unit tests for the crystal core (symmetry operations, elements, file formats) | none | fast |
| T11_input_snapshots | Approval snapshots of CSO-RM, CSO-FM, Gaussian and VASP inputs for 7 structures | `expected/` | fast |
| T12_package_contract | `Source`/`csptoolbox` import paths, no import side effects, console scripts, assets | none | fast |
| T13_known_issues | Known scientific issues as strict expected failures (docs/KNOWN_ISSUES.md) | small inline fixtures | fast |

Rules for adding a case are in [docs/MAINTENANCE.md](../docs/MAINTENANCE.md#testing).
Never reuse a number: study records outside this repository cite these paths.
The next free number is **T14**.
