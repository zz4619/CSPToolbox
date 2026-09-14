# CP2 symmetry validation: session report and handoff

Recorded 14 September 2026, approximately 22:01 UTC (23:01 BST).

The Z′=16 ice study is complete. The corrected Z′=4 pilot is complete, but the
full 20,000-return PLATON batch has not yet run. Its results must not be inferred
from the pilot. The remaining work is in [TODO.md](TODO.md).

## Completed work

### Z′=16 ice: all four converged returns remain P1

Both starting and final geometries of all four returned structures retain P1,
with 16 inequivalent waters and no smaller primitive cell in the tested range.
The search timed out after these four converged returns; this is not a complete
10,000-candidate search or evidence about all possible ice structures.

- Reconstructed all eight start/final geometries and checked them against the
  pinned native CP2 geometry routines. Maximum Cartesian component disagreement:
  4.45e-15 Å, before effects from the original printed precision.
- Applied existing CSPToolbox symmetry detection, ASU reduction and expansion.
  All 192 round-trip checks preserved element counts; the maximum displacement
  was 1.02e-7 Å from the existing coordinate canonicalization.
- Installed the official command-line PLATON build on CX3. All 504 water checks
  and 126 positive/H-broken control checks passed.
- Tested oxygen-framework symmetry separately from full O/H geometry. Ordinary
  PLATON ADDSYM excludes H; controlled auxiliary point labels made H positions
  visible to its geometric test. These copies were never used as energy inputs.

The [Z′=16 report](../T8_water_symmetry_20260914/README.md) contains tolerance
matrices, installation provenance, input hashes, controls and reproduction steps.
Local PLATON linking remains unavailable without an Apple SDK; the working CX3
installation is sufficient for this study.

### Z′=4: corrected native export and validated pilot

The diagnostic exporter originally selected the last `SpaceSupported` entry
matching an international group number. That substituted settings such as P2/n
for P2/c and P21/n for P21/c. We corrected only this separate diagnostic to use
the unique matching entry in the run-specific `SpaceList`, rejecting ambiguity.
The original CP2 searches used the named settings from their inputs; neither
their source nor archived coordinates was changed.

ASU coordinates now come from each molecule's first symmetry copy in the final
native full-cell output. Before PLATON, checks establish operation closure,
cell-metric preservation and exact periodic reconstruction of the full cell.
The apparent 0.29 Å O–O contact in the incorrect Ice row-33 export disappears:
the corrected pair separation is approximately 3.88 Å.

| Corrected pilot | Structures | PLATON checks | Export/check failures | Maximum native expansion error |
| --- | ---: | ---: | ---: | ---: |
| Ice Z′=4 | 28 | 168 | 0 | 8.03e-15 Å |
| Nicotinamide Z′=4 | 31 | 186 | 0 | 7.85e-14 Å |

Each system's pilot covers all 26 input groups and every saved status category.
Independent spglib checks using the actual O/H coordinates agree with three Ice
proposals at 0.05 Å:

| Ice row | CP2 status | Proposed international group | Geometric Z′ |
| --- | ---: | --- | ---: |
| 4 | 0, converged | 33, Pna2₁ | 1 |
| 49 | 6, unsuccessful | 9, Cc | 2 |
| 287 | -1, unsuccessful | 4, P2₁ | 1 |

Only row 4 is a converged minimum. These are pilot findings, not population
frequencies or experimental phase assignments. The H-position proxy is a search
aid; physical interpretation requires the actual element-labelled geometry.

## Full batch: submitted, result checks pending

- Job: `4056251.pbs-7`, `cp2_z4_platon`, queue `v1_medium24`.
- Resources: 32 independent processes, 64 GB, one-hour walltime.
- Latest scheduler check at this handoff: **Q**, insufficient available CPUs.
- Scope: all 10,000 terminal returns from each system, including unsuccessful
  minimizations, with original row IDs and statuses retained.
- Protocol: ordinary and H-position-proxy ADDSYM at 0.01, 0.05 and 0.1 Å;
  `CALC ADDSYM EXACT`, angular criterion 0.1 degree.
- Expected coverage: 120,000 checks if all 20,000 rows export successfully;
  otherwise six checks per valid export and explicit accounting for every
  excluded row. A normal PLATON exit alone does not establish a valid result.

The automation `finish-z4-platon-validation` is **paused** for session closure.
Its saved interval is four hours. The CX3 job remains submitted and can run
independently of this session. A future session must inspect it before submitting
anything; there is no need to rerun either CP2 global search.

## Locations, ownership and preserved evidence

This study belongs to CSPToolbox branch
`codex/water-z16-symmetry-validation`, in the isolated local worktree
`CSPToolbox_water_symmetry`. Implementation/results were committed and pushed as
`b53d757` (Z′=16) and `19a0d7a` (corrected Z′=4 pilot and full-batch scripts).
No production branch was merged and no library or CP2 application behaviour was
changed. This handoff adds documentation only.

External evidence roots, relative to the workspace locally and home on CX3:

- `CSP_studies/water_z16_symmetry_20260914`: installation and Z′=16 evidence.
- `CSP_studies/z4_platon_20260914`: frozen inputs, diagnostic, pilot, chunk logs
  and pending full-batch results. Local compact evidence is already present.

The complete Z′=16 evidence archive has SHA-256
`cf90ba4ecd42f2097fba947e76e59c6f574fa673c47cf3dbb5b43df7e2e4e820`.
The corrected Z′=4 pilot archive has SHA-256
`fa01a384a4272b8e43d21d1b5fe32f08eacaa86abec80861ff85c1afd3bca741`.
Earlier `pilot/`, `pilot_asu/` and `spf_probe/` exports are superseded. Use only
`pilot_corrected/` and the corrected, protocol-hashed `runs/` for scientific work.

The landscape/RMSD20 task was handed to dedicated thread
`01a0a1a6-b775-7e13-8eb3-54231f9c4247` and completed separately at CSPImperial
commit `04f71e8`. Its report is
`CSPImperial_high_z_analysis/doc/CP2_Z4_LANDSCAPES_20260914.md` relative to the
workspace. Do not duplicate that task: its oxygen-shell matches do not replace
the symmetry validation pending here.

## Verification and remaining limits

The repository contains checksums and machine-readable validation results for
both completed studies. At this handoff, all files in the T8/T9 result checksum
manifests passed `shasum -a 256 -c SHA256SUMS`; both local evidence archives
matched the hashes above. `git diff --check` passed, and all 17 relative file
links across the five touched documentation files resolved. No application
tests or scientific reruns were needed for documentation-only changes.

The only current execution blocker is the CX3 CPU queue. Population-wide
symmetry counts, full-batch exclusions, evidence completeness and independent
validation of representative full-batch proposals remain **unverified** until
the checklist is completed. Coordinate symmetry alone does not establish phase
identity, thermodynamic stability, search completeness or Hamiltonian symmetry
under exchange of differently parametrized hydrogen sites.
