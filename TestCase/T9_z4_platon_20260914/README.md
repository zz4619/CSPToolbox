# PLATON checks of all Z′=4 returns, 14 September 2026

This extends the [Z′=16 water validation](../T8_water_symmetry_20260914/README.md)
to **all 10,000 terminal records each** from the frozen Ice_Z4 and
Nicotinamide_Z4 CPU benchmarks. It includes unsuccessful minimizations and
retains their native row IDs, statuses and energies. It does not cluster,
minimize, alter the archived coordinates, or identify experimental phases.

The complete batch runs on CX3 as job `4056251.pbs-7`, using 32 independent
processes and a one-hour allocation. Scientific results will be added after
checking completeness and any exclusions.
At the session handoff on 14 September 2026, the job was still queued for CPUs.
The thread follow-up `finish-z4-platon-validation` is **paused** because the user
is closing the session; its saved interval is four hours. The submitted CX3 job
has not been cancelled. Resume from the [pending result checklist](TODO.md),
and read the [session report](SESSION_REPORT.md) for completed findings and limits.

## Geometry and settings

- Coordinates come from a separate diagnostic using the pinned native CP2
  geometry routines and each run's own molecular/LAM inputs.
- The diagnostic selects the unique matching setting from the run-specific
  `SpaceList`. Numeric international group IDs alone do not distinguish settings.
- CIF asymmetric units use the first symmetry copy of each molecule in the
  final native full-cell output. `build_unit_cell` restores fixed SG cell
  lengths/angles after the preliminary ASU calculation.
- Native symmetry operations must form a closed group, preserve the cell metric,
  and reconstruct every full-cell atom modulo lattice translations. These tests
  are applied before a structure is passed to PLATON.
- PLATON receives the original native symmetry operations and the complete ASU.
  Atom labels are unique; coordinates retain 12 decimal places in the generated
  CIF, without claiming to recover precision lost in CP2's printed poses.
- `CALC ADDSYM EXACT 0.1 d d d` is run for `d = 0.01, 0.05, 0.1 Å`.
  The metric angle criterion is 0.1 degree and unmatched atoms are not allowed.
- Ordinary ADDSYM omits H. A separate geometric copy changes H labels/types to F,
  which is absent from both systems, while retaining C/N/O identity. This makes
  the H positions available to the typed-coordinate symmetry search. This copy
  is never used as a chemical model or energy input.
- Included ASU atom counts and reduced **primitive** cell counts must agree with
  the expected values. PLATON's printed `Unitcell` count excludes centering
  translations; it is not the original conventional full-cell count.
- Native export failures, PLATON timeouts and atom-count changes are retained
  explicitly. A normally exiting PLATON process alone is not a passed check.

## Exporter correction found during the pilot

The original diagnostic searched the whole `SpaceSupported` table by group
number and retained its last match. That substituted P2/n for requested P2/c,
P21/n for requested P21/c, and an alternate group-33 setting. A numeric group ID
is insufficient to select the run's actual setting. The original CP2 search
selects its groups by name from `input.in`.

`fix_exporter.py` changes only the separate diagnostic to consult the run's
`SpaceList` and stop if the numeric ID is ambiguous or absent. It preserves the
original diagnostic source and executable and records both old/new hashes.
CP2's scientific source and the original searches are unchanged.

An apparent 0.29 Å O–O separation in Ice row 33 was an artifact of that incorrect
export. The corrected native geometry puts the pair about 3.88 Å apart. PLATON
retains all atoms and reports international group 13, in its conventional P2/n
setting. That output setting name is not evidence of a symmetry increase.

The corrected pilot covered one representative of all 26 input groups in each
system plus all saved nonzero-status categories: **59 structures, 354 PLATON
checks, zero export or PLATON failures**. Earlier `pilot/`, `pilot_asu/` and
`spf_probe/` artifacts are retained for the diagnostic audit and are superseded.
Only `pilot_corrected/` and the corrected production `runs/` are usable results.

Independent spglib checks on the actual O/H coordinates confirm the pilot's
three full-atom higher-symmetry suggestions at 0.05 Å. Only row 4 is converged:
it gives group 33 and geometric Z′=1. Rows 49 and 287 have nonzero CP2 statuses
and remain failed minimizations despite their geometric symmetry. This small
pilot is not a frequency estimate for the 20,000-return population.

## Files and reproduction

The external study is `CSP_studies/z4_platon_20260914` under this workspace and
the CX3 home directory. It holds frozen inputs/executables, pilot logs, per-chunk
archives, full row/check tables and job logs. PLATON source and installation
provenance remain in `CSP_studies/water_z16_symmetry_20260914`.

`native_geometry.py` is an unchanged, hashed snapshot of the parser and native
expansion validator from the parallel scientific-analysis task; its origin is in
`native_geometry_provenance.json`. It is not a replacement library implementation.
`results/initial_copy_manifest.json` records the initial freeze, before the
diagnostic fix. The original driver/binary are retained as `*_original.f90` and
`export_legacy.x`; corrected hashes are in `results/exporter_fix.json` and the
production run's protocol manifest.

On CX3, `prepare_remote.py` creates a fresh study from the existing frozen source
study. After uploading these scripts into `scripts/`, load the already configured
CSP environment and run `fix_exporter.py` once. Then use the Python interpreter
from the installed CSD module, as demonstrated in `run_z4_platon.pbs`:

```bash
python scripts/run_batch.py --study "$PWD" --pilot --workers 2 --chunk-size 4 \
  --timeout 10 --output-name pilot_corrected
qsub run_z4_platon.pbs
python scripts/collect_results.py --study "$PWD"
```

After downloading the compact results, use an interpreter with spglib to run
`interpret_results.py --results /absolute/path/to/results`. It compares
international group numbers and symmetry operations per unit volume, accounting
for the reported cell transformation. A proposed Z′ assumes general-position
multiplicity and is a geometric estimate. Report ambiguous/excluded checks
explicitly and verify full molecular symmetry before making physical claims.

These are study-specific commands, not a generic CP2 installation recipe.
`prepare_remote.py` requires the named earlier study and refuses to overwrite an
existing target. The original inputs contain the named SG settings and must be
preserved alongside the printed pose records.

Production uses chunks of 100 returns. Each complete chunk has compressed logs
and checksummed results. Resumption requires identical input, executable, script
and tolerance hashes; incomplete chunks are not counted as complete. Temporary
working directories are removed only after their evidence is archived.

Verification commands for these scripts, from the CSPToolbox repository root:

```bash
python -m py_compile TestCase/T9_z4_platon_20260914/*.py
bash -n TestCase/T9_z4_platon_20260914/run_z4_platon.pbs
git diff --check
```

Done when every native row is accounted for exactly once, all expected checks
are present for exported rows, exclusions are explained, and evidence hashes and
the source-setting/coordinate contracts are verified. Distinguish a change of
space-group setting from a change of international group or primitive-cell size.
