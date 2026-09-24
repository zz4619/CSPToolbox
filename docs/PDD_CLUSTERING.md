# Indexed AMD/PDD clustering

For experimental matching, see [PDD_RMSD20_WORKFLOW.md](PDD_RMSD20_WORKFLOW.md).
Its screening cutoff and packing settings differ from landscape clustering.

This implementation replaces a scan over every representative with exact AMD
radius lookup, then compares the surviving PDDs using Earth Mover's Distance
(EMD). It never merges structures merely because AMD matches. It is a Python
analysis workflow; no CP2 production Fortran is changed.

## Scientific contract

- Input geometry is a full periodic unit cell, in Å. The file-manifest reader
  expands asymmetric-unit CIF/RES/PDB input through `CrystalStructure` first.
- PDD uses the nearest `k` atom distances, including intramolecular distances;
  the default is 100. Typing distinguishes the **central element**, not neighbour
  element, atom label, bond topology or molecular identity. H atoms are included.
- EMD compares cached weighted rows. Typed comparisons require equal element
  fractions, so equivalent supercells can match. Different stoichiometries remain
  separate. Distances alone cannot distinguish mirror images.
- AMD is the weighted column mean. Its Chebyshev distance is a lower bound on
  the full Chebyshev-ground-cost EMD. Indexing selected AMD columns weakens that
  bound but cannot discard additional PDD matches. Queries are exact radius
  queries, not approximate or limited to a fixed number of neighbours.
- The index searches immutable KDTree blocks and a pending buffer. Equal-sized
  blocks merge on insertion; every representative is always covered. No fixed
  dimensional index is guaranteed to outperform a scan on all datasets.
  Enabled energy/density gates are also indexed as normalized range coordinates;
  the original unscaled values are checked again before EMD.
- Records are sorted by `(energy, identifier)`. Every member must directly match
  a fixed, real representative; representatives never drift and clusters do not
  merge by transitive chains. This intentionally differs from native Analyse's
  moving centroids. Space-group labels do not partition the search.
- Optional metadata gates use strict differences `<1 kJ/mol` and `<20 kg/m3`,
  matching AXOSOW's final clustering eligibility criteria. Metadata must have a
  consistent per-molecule/formula-unit normalization. These gates are policies,
  not mathematical bounds on geometric similarity.
- PDD acceptance uses `EMD <= threshold`. A small numerical margin is used only
  for candidate lookup; it does not enlarge the final acceptance threshold.
- Descriptors are validated once during AMD construction. Element rows and
  uniform-weight assignment information are prepared once for the current
  record and retained representatives. Repeated comparisons reuse that data;
  they use the same EMD, including the general transport solver for unequal
  masses. This avoids repeatedly validating and regrouping unchanged arrays.

The periodic neighbour search now reduces the lattice basis and explicitly
bounds all images within a trial radius using reciprocal-vector lengths. It
returns only once every atom's kth neighbour lies within that covered radius.
The previous integer-shell plateau rule missed short lattice vectors in skew
bases. Regressions cover a 100-fold shear, equivalent supercells, rotation,
translation, atom reordering and an independent finite-cloud reference.

## API and command

```python
from csptoolbox import (
    PDDClusterSettings, load_pdd_cluster_manifest, cluster_pdd,
)

records, provenance = load_pdd_cluster_manifest(
    "manifest.csv", k=100, cache="descriptors.npz",
)
settings = PDDClusterSettings(threshold=0.20)
result = cluster_pdd(records, settings)
reference = cluster_pdd(records, settings, indexed=False)
assert result.assignments == reference.assignments
```

The example threshold is a **trial**, not a universal duplicate criterion.
The core API also accepts `PDDClusterRecord` objects with precomputed descriptors
and an optional `confirm(record, representative)` callback. `make_compack_confirmation`
constructs the licensed CCDC callback with the original AXOSOW settings: 15
matched molecules, RMSD strictly below 0.20 Å, hydrogen/bond counts respected.
Use one callback per process; do not share it between concurrent threads.

Manifest columns:

```csv
id,path,energy,density
minimum_1,structures/1.res,-110.1234,1350.0
minimum_2,structures/2.res,-109.9876,1351.0
```

Paths are relative to the manifest. Identifiers must be unique. A cache records
input-file contents, manifest contents and descriptor settings/version; changed
inputs are rejected rather than silently reusing stale descriptors.

```bash
python -m Source.CLI_scripts.cluster_pdd manifest.csv \
  --threshold 0.20 --cache descriptors.npz --output clusters.json
```

The installed command is `csp-cluster-pdd`. The command above performs PDD-only
clustering, as required for replacing both historical stages. For a separate
compatibility experiment, `--compack` adds confirmation using an already installed,
licensed CSD Python API; the library never installs it. Use `--exhaustive` to validate index
equivalence. A zero energy/density tolerance on the CLI disables that metadata
gate; the API uses `None`. Descriptor construction and numerical thread counts
default to one worker. Control BLAS threads through the normal environment.

Output includes assignments, representative IDs, populations, EMD/confirmation
counts and separate descriptor versus clustering times. Input files and earlier
output reports are not overwritten. Cached comparison requires NumPy/SciPy;
geometry construction also requires the regular CSPToolbox dependencies.

## Calibration and limitations

Do not copy CP2's 0.20 Å contact tolerance or COMPACK's 0.20 Å RMSD directly into
EMD as an equivalent threshold. Compare labelled pairs, report false merges and
missed duplicates, and also compare full cluster assignments. Split controlled
perturbations by original structure before generating pairs to avoid leakage.
`threshold_agreement` supplies pair metrics and returns undefined recall when
there are no positive examples.

The initial Nicotinamide Z′=4/26-SG study uses 283 real structures (lowest 256
energies plus the lowest two per SG from a selected 200,000-result prefix).
It also uses 64 rigid-molecule packing perturbations from 16 original structures,
split into eight source families for calibration and eight for validation.
Perturbations are controls, not minimized search results; their inherited energy
metadata must not be interpreted as recalculated energies.

Those earlier controls show that a single PDD threshold cannot reproduce all
COMPACK decisions: some accepted pairs exceed 0.18 Å EMD while some rejected real
pairs are below 0.20 Å. That compatibility experiment used AMD/PDD as
screening and retains COMPACK confirmation. No universal threshold or guaranteed
COMPACK recall is claimed. See `TestCase/T8_pdd_clustering/README.md` for the
completed study evidence and measured limits.

The batch API accepts descriptors backed by memory-mapped arrays, but still
keeps Python records and AMD vectors for every input. The NPZ CLI cache loads
all distances into memory. A million Z′=4 structures with 60 rows and k=100 need
about 48 GB for the distance arrays alone. A fully streaming CLI is a future
extension; the search trees contain only representative coordinates. Do not
construct an all-pairs distance matrix for production datasets.

The generic CLI consumes coordinate files; it does not yet provide a general CP2
LAM/record reader. The Z′=4 pilot used native Analyse's exported RES files.
For the clarified target of replacing **both** native reduction and COMPACK, the
Ben Nicotinamide Z′=1 study now has a separate raw-record reconstruction adapter,
validated against all 4,953 archived first representatives. It consumes all
548,435 raw minima and memory-maps per-batch distance arrays. This adapter is
specific to that historical 15-atom, three-torsion model; it must not silently
be used for arbitrary molecules, Z′ or LAM formats. The reproducible study lives
at `CSP_studies/nicotinamide_z1_ben_pdd_20260922` in the parent workspace.
Cached clustering timings still must not be called full-pipeline speedups.
The completed raw benchmark yields 2,921 clusters at 0.10 Å and 2,696 at 0.15 Å,
against 2,777 historically; agreement is assessed on the 4,953 labelled anchors.
See [the full benchmark results](../TestCase/T8_pdd_clustering/README.md#replacement-of-both-stages-on-all-raw-z1-minima).

## Verification

### Parallel batches for uniform PDD arrays

`prepare_pdd_batch_store` and `cluster_pdd_batched` provide a separate array API
for large datasets with identical central-element order and equal row weights.
The NPY input has shape `(structures, rows, k)` and dtype float64. This path does
not accept mixed row layouts, collapsed/nonuniform weights, or COMPACK callbacks;
use the existing `cluster_pdd` API for those cases. The serial default is unchanged.

```python
from csptoolbox import prepare_pdd_batch_store, cluster_pdd_batched, PDDClusterSettings

# Energies, densities and unique string IDs follow the PDD file's raw row order.
# elements gives the central element of each uniformly weighted PDD row.
prepare_pdd_batch_store("pdd.npy", energies, densities, identifiers, elements,
                        "array_store", workers=16)
result = cluster_pdd_batched("array_store", PDDClusterSettings(0.10),
                            "clusters_010", workers=15, batch_size=1024)
```

Run multiprocessing entry points under `if __name__ == "__main__":`. Fifteen
comparison workers plus the coordinator fit a 16-core allocation. Set numerical
library threads to one before starting Python. The AMD preparation stage can use
16 workers; it is separate from the clustering pool.

Each energy/ID-ordered batch searches a frozen representative snapshot in
parallel. The coordinator then checks only unmatched rows against representatives
created earlier in the same batch. An old representative always precedes a new
one, so a snapshot match cannot be superseded by a new representative. Both
paths select the earliest passing representative and use the existing scalar
PDD comparison. Independent shard clustering followed only by representative
merging does not provide this guarantee.

Workers map the same read-only PDD/AMD arrays. Each has a bounded LRU comparison
cache (256 records by default) and a local index containing published representative
IDs. Large prepared descriptor dictionaries are not copied into every worker.
All input files must remain immutable for the lifetime of a run. Store metadata
and ordering arrays have content hashes; the PDD cache's path, size and mtime are
checked, and callers can also record their previously audited PDD checksum.

The implementation reports snapshot and within-batch candidate/EMD counts and
timings separately. Worker durations are summed work times, not elapsed wall time;
the two `current_run_*_wall_s` fields measure the coordinator's phase durations.
Measure PBS/cgroup job memory; adding per-process RSS can double-count shared pages.

`labels.npy` maps each original input row to a zero-based original representative
row. `representatives.npy` lists representatives in energy/ID order. A completed
result checks representative self-assignment, populations and energy ordering.
`checkpoint.npz` is atomically replaced after completed batches (every 120 s by
default). Calling the same function on that output resumes from its last saved
batch; worker and batch counts can change, but data/settings cannot. Completed
outputs are verified and reused. A concurrent-writer lock is removed on clean
exit; after a hard kill, remove a stale lock only after verifying that its owner
and scheduler job have stopped. Partially built AMD stores are not resumable.

The September 2026 Nicotinamide study uses existing reconstructed Z′=2/Z′=4
caches. Descriptor/validation jobs use 16 cores and 96 GB; the updated production
jobs request 64 cores and 160 GB, with serial label checks before full production.
Its scripts and live scheduler records are in the parent workspace at
`CSP_studies/cp2_nicotinamide_pdd_parallel_20260922`. Speedups must come from its
matched benchmarks, not from multiplying the worker count. Small datasets can
be slower because of process startup and synchronization.

The parallel regression suite checks chains crossing batches, equal-energy ID
ordering, typed descriptors, threshold/metadata boundaries, index block merges,
bounded-cache eviction and restarting with different worker/batch counts.
It compares all labels, representative order and PDD comparison counts with
the serial implementation.

```bash
python -m unittest discover -s TestCase/T8_pdd_clustering -p 'test_*.py' -v
python -m Source.CLI_scripts.cluster_pdd --help
python -m pip wheel --no-deps --no-build-isolation . -w /tmp/csptoolbox-pdd-wheel
git diff --check
```

The tests need no CCDC installation. Full confirmation benchmarks use the site's
existing CSD module. No repository-wide formatter/linter/type-checker is configured.
