# PDD screening and RMSD20 experimental matching

Reviewed 24 September 2026. These are two related workflows with different
thresholds: landscape clustering groups predicted minima; experimental matching
shortlists prediction/reference pairs and confirms them by packing superposition.
Do not use the landscape clustering threshold as an experimental-match cutoff.

## Current landscape clustering

`csp-cluster-pdd` reads `id,path,energy,density` from a CSV manifest, expands
crystals to full periodic cells, builds element-typed k=100 PDDs, uses AMD as
a rigorous lower-bound index, and evaluates shortlisted PDDs by normalized EMD.
The current manifest command includes hydrogen. All space groups are compared
together for a composition; space-group labels are not a rejection filter.
See [PDD_CLUSTERING.md](PDD_CLUSTERING.md) for the full contract and commands.

Nicotinamide experiments used PDD thresholds 0.10 and 0.15 angstrom, with pair
energy difference <1 kJ/mol and density difference <20 kg/m3. These metadata
gates are scientific policies, not geometric lower bounds. Representatives are
fixed structures selected in energy/identifier order, not moving averages.
`cluster_pdd_batched` preserves this ordering with parallel comparisons against
a frozen representative set and ordered reconciliation within each batch.

For the forthcoming online study, retain only successful finite minima with
`E_final <= E_reference_Z1 + 20 kJ/mol` for descriptor/clustering work. The
reference is preset from the known Z'=1 search minimum and must be recorded with
its model identity. Do include structures below that reference. Save raw starts,
results and failed minimizations separately. No descriptor is needed for an
out-of-window result. This is an analysis filter, not an initial-energy rejection
rule. Online arrival-order clusters will be provisional and are distinct from
the energy-ordered offline contract above.

## Experimental matching: PDD then RMSD20

1. Select eligible predicted structures/representatives and record their source
   IDs, energies, densities, space groups and Z'. Record experimental Z' too.
2. Remove H from both descriptor centres and periodic neighbour sets. Construct
   element-typed PDDs with k=100 (k=60 is a tested comparison). Geometry remains
   periodic; centre weights represent complete symmetry orbits.
3. Use PDD to shortlist candidate/reference pairs. Heavy-atom PDD <=0.30 angstrom
   is a **provisional experimental-screening setting**, not a universal match
   criterion. It recovered nine previously confirmed Z'=1 pairs in the existing
   nicotinamide calibration; higher-Z' recall is not established.
4. Confirm shortlisted pairs with the existing licensed CCDC PackingSimilarity
   API using the settings below. Do not install CCDC implicitly.

```python
from ccdc.crystal import PackingSimilarity

similarity = PackingSimilarity()
settings = dict(
    ignore_hydrogen_counts=False,
    ignore_hydrogen_positions=True,
    ignore_bond_counts=False,
    ignore_bond_types=True,
    allow_molecular_differences=False,
    distance_tolerance=0.40,
    angle_tolerance=40,
    packing_shell_size=20,
)
for name, value in settings.items():
    setattr(similarity.settings, name, value)
result = similarity.compare(reference_crystal, candidate_crystal)
accepted = (result is not None and result.nmatched_molecules == 20
            and result.rmsd < 0.80)
```

An RMSD for fewer than 20 matched molecules is a partial-shell result, not
RMSD20. Record errors/timeouts separately from non-matches. Persist the CCDC
version, actual settings, file hashes, matched count, RMSD, elapsed time and
decision for every pair. Use one API instance per process.

**The `csp-cluster-pdd --compack` option is different:** it retains historical
AXOSOW defaults of 15 molecules and RMSD <0.20 angstrom. Merely requesting
`shell_size=20` from its callback does not reproduce the full settings above.
The published CLI currently has no end-to-end heavy-PDD/RMSD20 matching command;
the completed matching runs use study-specific drivers. This README documents
the actual workflow rather than claiming that interface already exists.

## Evidence and remaining validation

The September 23 nicotinamide study compared all-atom k=100 with H-excluded
k=60/100. At PDD <=0.30 angstrom, all recovered the nine known Z'=1 pairs;
additional pairs rejected by the saved RMSD20 criterion numbered 5,101, 675 and
289 respectively. Six selected higher-Z' candidates were subsequently checked:
three Z'=2 candidates matched all 20 molecules; three Z'=4 candidates had partial
shells. These selected checks do not establish complete landscape coverage.

Study evidence remains outside Git under the workspace's
`CSP_studies/nicotinamide_heavy_pdd_20260923/` (README, RESULTS, verification,
report audit and `rmsd20_selected_20260923/`). Large descriptor caches and licensed
experimental inputs are not bundled into this documentation commit.

No PDD-to-RMSD20 rejection bound is proven. The remaining validation must measure
missed confirmed matches, including independent checks above the proposed PDD
cutoff. An energy-window or representative-only result is not a guarantee for
all raw structures. Cluster representatives should not silently substitute for
all members when claiming experimental-match recall.

## Verification

From the repository root with the documented Python dependencies:

```sh
python -m unittest discover -s TestCase/T8_pdd_clustering -p 'test_*.py' -v
python -m Source.CLI_scripts.cluster_pdd --help
git diff --check
```

These cover indexing/exhaustive equivalence, periodic geometry, uniform versus
weighted transport and batch/restart semantics. They do not exercise licensed
RMSD20; that requires the CSD environment and recorded reference-pair controls.
The broader experimental-matching recall audit remains deferred.
