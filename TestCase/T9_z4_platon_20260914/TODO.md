# TODO: collect and validate the complete Z′=4 PLATON results

Handoff: 14 September 2026. Read [SESSION_REPORT.md](SESSION_REPORT.md) and
[README.md](README.md) first. Job `4056251.pbs-7` was queued at handoff.
The scheduled follow-up is paused; collection requires a future session.

## Completed checks

- [x] Z′=16 ice: four converged structures and their starts retain P1; native
  geometry validation, 192 reduction round trips, 504 PLATON structure checks
  and 126 control checks completed. See [T8](../T8_water_symmetry_20260914/README.md).
- [x] Correct the diagnostic's SG-setting selection using the native run's
  `SpaceList`; preserve original sources, inputs and superseded exports.
- [x] Corrected Z′=4 pilot: 59 structures, 354 passed checks, zero export/check
  failures, covering all 26 groups and every native status category.
- [x] Independently check the pilot's three full-O/H proposals at 0.05 Å. Only
  Ice row 4 is converged (group 33, geometric Z′=1); rows 49/287 remain failures.
- [x] Submit the full 20,000-return batch, freeze its protocol and provide a
  collector that checks complete row/check coverage and archive hashes.

## Pending result checks

- [ ] Refresh `/opt/pbs/bin/qstat -fx 4056251.pbs-7` on CX3. Check scheduler exit
  status, `job_exit_status.txt`, `cp2_z4_platon.o4056251`, `progress.jsonl` and
  `batch.stderr` in `$HOME/CSP_studies/z4_platon_20260914`. If PBS history has
  expired, inspect the saved files; disappearance from qstat is not success.
- [ ] Confirm frozen `runs/protocol.json` and `runs/scripts/` use the corrected
  exporter and match the expected input/executable/script hashes. Never combine
  obsolete pilots with production chunks.
- [ ] If incomplete or timed out, diagnose first. Resume the existing chunked
  run only after the original job is no longer active and only with the same
  protocol. Do not overwrite completed chunks or launch a duplicate active job.
- [ ] Run the existing `scripts/collect_results.py --study "$PWD"` from the
  CX3 study root. It must account for rows 1–10,000 exactly once in each system,
  six unique checks per exported row, one protocol per system and correct
  hashes for all compressed evidence archives.
- [ ] Investigate every export failure, timeout, changed atom count, abnormal
  termination and unresolved group name. A failed scientific minimization and
  a failed symmetry check are different categories. Preserve both statuses.
- [ ] Check the two systems share the intended protocol. Record native group
  closure, metric and full-cell expansion error bounds and every exclusion.
- [ ] Copy compact `results/`, the frozen protocol/scripts, scheduler/job logs
  and full chunk evidence to the local study. Record sizes and verify hashes
  after transfer. Keep large third-party binaries/raw archives outside Git.
- [ ] Run `interpret_results.py --results /absolute/path/to/results` with the
  verified spglib environment. Investigate all unresolved/lower-symmetry cases.
- [ ] Report counts at each tolerance, split by system and convergence status:
  unchanged symmetry, setting-only changes, smaller primitive cells, higher
  heavy-atom symmetry and higher full-atom geometric symmetry. Compare group
  numbers and operations per volume, not only printed symbols or cell sizes.
- [ ] Independently validate representative proposals on the original actual
  element-labelled full coordinates, including converged low-energy examples
  and heavy-atom/full-atom disagreements. Verify periodic species-preserving
  assignments, molecular integrity and primitive-cell multiplicity. The proxy
  alone is insufficient for a full molecular symmetry claim.
- [ ] Save row-level result tables and a final interpretation with tolerance
  sensitivity, unresolved cases and scientific limits. Keep failed-minimization
  counts separate; do not label approximate symmetry as an experimental phase.
- [ ] Update this checklist/report, verify evidence hashes and links, and
  commit/push scoped compact artifacts on `codex/water-z16-symmetry-validation`.
  Do not merge dev/main or alter CP2 scientific sources. Report completion to
  the user; leave the follow-up paused when finished.

## Commands and completion criteria

Use the existing [PBS script](run_z4_platon.pbs) for the CX3 environment and
single-threaded numerical-library settings if a resumption is required. It uses
the installed CSD Python for the producer; do not assume default CX3 Python has
NumPy. The collector itself uses the standard library. The interpreter at
`/opt/anaconda3/envs/csp_310/bin/python` has spglib for local interpretation.

From this repository root, the existing lightweight script checks are:

```bash
python -m py_compile TestCase/T9_z4_platon_20260914/*.py
bash -n TestCase/T9_z4_platon_20260914/run_z4_platon.pbs
git diff --check
```

Done when all 20,000 rows are accounted for, six checks exist for every valid
export, exclusions/unresolved results are explained, full evidence has verified
local and CX3 copies, representative physical-coordinate checks are documented,
and final scientific counts are reported. Running jobs, partial chunks or
successful process exits alone do not satisfy completion.
