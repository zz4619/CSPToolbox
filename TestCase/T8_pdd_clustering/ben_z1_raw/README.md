# Reproduce the raw Ben Nicotinamide Z′=1 experiment

These are the study-runner snapshots, not a general CP2 input reader. Copy them
into a fresh **external** study directory containing the selectively downloaded
original files listed below; the scripts resolve inputs beside themselves.
Keep multi-GB descriptor arrays, raw inputs and generated output outside Git.
Use the checked-out CSPToolbox on `PYTHONPATH`, with NumPy/SciPy/ASE and the regular
project dependencies already available. No CSD installation is used.

Source on CX3:
`/rds/general/user/zz4619/home/03_Ben_CSP/Nicotinamide/5_GlobSrch`.
Required files, with their original directory layout:

- `input.in`, `crystals_stable.out`, `flexible_lam_intra`, `potential.in`;
- `Analyse_1` through `Analyse_4`: `crystals_stable.out`, `Analyse_log.out`,
  `Clustering_log.out`, and the unpacked `unique_pool.tar.gz` coordinate archive.

`cp2_groups.txt` is a static symmetry-table dump from the legacy CP2 source,
printed by the included Fortran helper. The original `global.f90` and
`space_group_module.f90` were compiled unchanged with the existing CX3 Intel
environment. Rebuilding the dump is not needed to rerun the Python benchmark.

After copying these scripts into the fresh data directory:

```bash
export PYTHONPATH=/path/to/CSPToolbox
export MPLCONFIGDIR=/tmp/csp-pdd-mpl
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python benchmark.py prepare
python -u benchmark.py run --thresholds .10 .15
python raw_benchmark.py validate
python -u raw_benchmark.py generate --workers 6
python -u raw_benchmark.py run --thresholds .10
python -u raw_parallel.py --threshold .15 --workers 2
python raw_controls.py
python summarize_raw.py
```

`benchmark.py` builds the archived-representative cache and reference mapping;
its exported-only results are separate from the raw replacement benchmark.
`raw_benchmark.py` reconstructs every one of the 548,435 successful raw minima,
without native contact-fingerprint clustering or COMPACK, and preserves the
historical four-batch comparison boundary. All outputs refuse overwrites.

Raw labels in `raw_labels_<threshold>_A<batch>.npy` are zero-based original
batch-line indexes of fixed representatives. Raw search generation IDs are in
column 1 of the original records and are not those line indexes.

The saved native logs do not give membership labels for every raw minimum.
Report full-run counts separately from pair-agreement metrics on the 4,953 saved
first-representative anchors. Neither agreement with those anchors nor a count
near 2,777 validates a universal PDD threshold. The 0.10/0.15 Å values are trials.
The parallel runner only overlaps independent historical batches; order within
each batch is unchanged. It checks for existing outputs before starting workers
and saves a per-batch checkpoint. Do not compare its two-worker wall time directly
with the single-worker trial as a threshold or implementation speedup.

This adapter specifically requires the historical 15-atom molecule, three torsions,
426 LAMs, SMOOTH=0.4, and no hydrogen shortening or fixed dependent coordinates.
It reproduces the historical broad smoothing cutoff. Unsupported systems must
receive a separately validated parser/reconstructor; do not relax assertions to
force another system through this example.
The reconstruction adapter checks the audited input, LAM and symmetry-table
hashes and rejects a different model before applying its specific conventions.

The workspace study README contains the complete audit, timing interpretation,
and reconstruction checks. Its compact completed summary is stored one directory
up as `ben_z1_raw_summary.json`; the large `.npy` arrays stay in the study folder.
