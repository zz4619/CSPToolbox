# AMD/PDD clustering regression and calibration

Run from the repository root using the existing scientific Python environment:

```bash
python -m unittest discover -s TestCase/T8_pdd_clustering -p 'test_*.py' -v
```

The 16 tests cover periodic geometry, equivalent supercells, skew lattice bases,
the AMD lower bound, exact transport shortcuts, fixed representative membership,
metadata boundaries, index rebuilds, composition, cache identity and public API/CLI
operation, including equivalence of prepared and general weighted PDD comparisons.
No CSD licence is needed. Six small RES fixtures are generated study
structures, not CSD database exports. Their recorded CCDC labels in
`fixtures/reference_pairs.json` are provenance; local tests recompute PDD values
and do not claim to rerun COMPACK.

## Nicotinamide calibration, 22 September 2026

The sample contains 283 minimized Nicotinamide Z′=4 structures: the lowest 256
energies plus the lowest two in each of 26 SG from the original 200,000-result
prefix. Native Analyse exported all 283 without a coarse merge. Energies span
−117.4607 to −97.8083 kJ/mol. Coordinate export is still a native Analyse step.

The reference is AXOSOW's unchanged `Clustering_serial.py`, SHA256
`a0537108be1231e431a4b474c23350007dec3445021fa525a54b15d616f14e3a`, using the
installed CX3 `CSD/2025.3.1` module (CCDC API 3.6.1). Eligibility requires strict
energy/density differences below 1 kJ/mol and 20 kg/m³. Matching requires 15
molecules and RMSD <0.20 Å, with hydrogen and bond counts respected.

The new pipeline uses typed PDD, k=100, row-collapse tolerance 1e−4 Å, eight AMD
index columns, 64-representative index blocks and **PDD EMD ≤0.20 Å followed by
the same COMPACK confirmation**. The 0.20 Å EMD threshold was selected from a
calibration grid; its numerical equality to the RMSD threshold is coincidental.

| Natural sample result | Reference | AMD/PDD + COMPACK |
| --- | ---: | ---: |
| Structures | 283 | 283 |
| Clusters | 282 | 282 |
| Matching assignments | — | 283/283 |
| COMPACK comparisons | 5,213 | 6 |
| Full AMD checks | — | 5,214 |
| PDD transport comparisons | — | 5,213 |

The index avoids most of the 39,646 potential representative checks; metadata
gates provide much of this reduction. AMD alone is a weak filter on this sample.
PDD filtering provides the major reduction in expensive COMPACK calls. The real
duplicate, structure 26 → representative 18, has PDD distance 0.0535201 Å and
COMPACK RMSD15 0.0233 Å.

A local replay with recorded confirmation labels took 3.82 s indexed versus
3.81 s for the vectorized exhaustive candidate scan, with identical assignments.
At this sample size, index lookup itself provides no measured time advantage;
both paths perform the same 5,213 PDD comparisons. Larger-scale index performance
still needs measurement.

## Threshold evidence

Sixteen source crystals provide 64 controlled pairs. Each control preserves the
molecular conformation and cell but applies independent rigid rotations and
translations to the four ASU molecules. Translation amplitudes are 0.02, 0.06,
0.15 and 0.35 Å; rotation angles are amplitude/5 radians. Eight source families
are used for calibration and eight held-out families for validation. The split
was made before generating the controls. These are **not minimized structures**;
their inherited energy metadata was not recalculated.

At the selected 0.20 Å PDD threshold:

| Pair set | True positive | False positive | False negative | True negative |
| --- | ---: | ---: | ---: | ---: |
| Controlled calibration | 24 | 0 | 0 | 8 |
| Controlled held-out validation | 24 | 0 | 0 | 8 |
| Natural reference comparisons | 1 | 5 | 0 | 5,207 |

A natural nonmatch (16 vs 9) has PDD 0.182061 Å, below the 0.185374 Å of one
held-out accepted control. Therefore **no single PDD threshold reproduces all
these COMPACK decisions**. Retaining COMPACK rejects the five natural false
positives. The sample contains only one natural duplicate, so this does not
establish universal recall. Keep the threshold explicit and validate more
independent minimized duplicates before production deployment.

## Timing and memory

All CX3 comparison runs use one numerical thread on a shared login node.
`calibration_summary.json` records the final run's precise values.

| Measured phase | Wall time |
| --- | ---: |
| Native export of the 283-structure sample | 59.77 s |
| Reference COMPACK active process time, two segments | 713.86 s |
| New clustering core, including six COMPACK calls | 5.07 s |
| New application, including cache loading and CCDC setup | 14.15 s |
| New complete Python process, using cached descriptors | 22.90 s |
| Initial descriptor preparation and imports, local workstation | 38.40 s |

The reference was bounded at 600 s after 234 complete records, then continued
for 113.84 s. A copied restart log was adapted to the original script's legacy
parser; the source and original logs were unchanged. The 713.86 s total includes
the interrupted partial work and two process starts, but excludes the gap and a
brief failed restart-parser attempt. This is an observed comparison, not a
controlled end-to-end speedup ratio.

Peak process RSS was 200,708 KiB for the initial reference segment and
238,936 KiB for the final cached implementation. Descriptor preparation ran
locally because the existing CSD environment lacks the geometry dependencies;
no packages were installed. Cached CX3 timings **exclude descriptor generation**.

The original 200,000-record baseline received SIGTERM after 139,370 native
records and 7,648.12 s; its Python stage never ran. This pilot does not complete
or replace that full pipeline. A geometry exporter bypassing native coarse
clustering and streaming/sharded descriptors remain necessary for a large run.
See [the implementation guide](../../docs/PDD_CLUSTERING.md) for API contracts
and memory scaling.

## Full evidence

Local workspace and CX3 home both contain the study folder
`CSP_studies/nicotinamide_pdd_calibration_20260922`. Its README indexes the input
selection, source hashes, reference logs, restart receipt, pair labels, descriptor
caches and timing files. Raw study outputs are kept there; only the compact
summary and regression fixtures belong in this test directory.

## Subsequent Z′=1 benchmark

`ben_z1_summary.json` records the completed test of every first-representative
export from Ben's Nicotinamide Z′=1 search: **4,953 inputs and 2,777 saved COMPACK
clusters**, with 2,176 duplicate assignments. It preserves the four historical
analysis batches and runs PDD without COMPACK confirmation. PDD thresholds
0.05/0.10/0.15/0.20/0.30 Å produce 3,446/2,911/2,686/2,484/805 clusters.
The summary includes label-invariant agreement metrics, not just cluster counts.
No production threshold or definition of acceptable similarity was selected.

The full inputs, reproducible runner, checks and interpretation are in workspace
`CSP_studies/nicotinamide_z1_ben_pdd_20260922/README.md`. This test replaces the
second-stage comparison only: native CP2 had already reduced 548,435 raw minima
to these 4,953 representatives. It does not demonstrate performance or agreement
when bypassing that native reduction.

## Replacement of both stages on all raw Z′=1 minima

The subsequent full experiment bypasses both native fingerprint reduction and
COMPACK. All **548,435 raw minima** were reconstructed using the original 426-LAM
model. Reconstruction was checked against every one of the 4,953 archived first
representatives and ASU/full-cell PDD equivalence against all 26 sampled SGs.

| Threshold | Full clusters | Anchor pair precision | Anchor pair recall | Clustering wall time | Workers |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.10 Å | 2,921 | 93.46% | 83.38% | 747.03 s | 1 |
| 0.15 Å | 2,696 | 85.45% | 89.64% | 721.43 s | 2 |

The historical final count is 2,777. Descriptor generation took 302.58 s with six
workers. The differing worker counts prevent a direct wall-time comparison between
thresholds. Reference agreement concerns the 4,953 labelled anchors only: the old
logs do not preserve every raw membership. Similar counts alone are insufficient.

`ben_z1_raw_summary.json`, `ben_z1_raw_controls.json` and
`ben_z1_raw_provenance.json` preserve compact results and environment/source
identities. [Reproduction scripts](ben_z1_raw/README.md) contain the narrowly
validated raw-record/LAM adapter, serial and parallel runners, and audits.
The approximately 6.58-GB descriptor cache, original data, all raw assignments and
selected original CP2 records stay in the external study folder.

The geometry-only control reproduces all earlier exported-coordinate assignments
exactly when native-average metadata is retained. Actual raw energy/density values
and the larger representative pool change greedy clustering decisions, explaining
why calibration of the second stage alone does not establish full-pipeline
agreement. No production similarity threshold has been accepted.
