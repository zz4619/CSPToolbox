"""Command-line interface for deterministic RDKit conformer generation."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from Source.conformer_generation import (
    ConformerGenerationSettings,
    generate_conformer_ensemble,
    validate_ensemble_identity,
    write_conformer_artifacts,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate an ETKDGv3 conformer ensemble, optimize it with an RDKit "
            "force field, and write a lowest-energy XYZ initial guess."
        )
    )
    parser.add_argument("smiles", help="SMILES fixing connectivity and protonation.")
    parser.add_argument("--name", default="molecule", help="Molecule/output name.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for XYZ, SDF, CSV, and JSON artifacts.",
    )
    parser.add_argument("--num-conformers", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260727)
    parser.add_argument("--embed-prune-rms", type=float, default=0.25)
    parser.add_argument("--deduplicate-rms", type=float, default=0.35)
    parser.add_argument(
        "--force-field",
        choices=("MMFF94", "MMFF94s", "UFF"),
        default="MMFF94s",
    )
    parser.add_argument("--max-iterations", type=int, default=1000)
    parser.add_argument("--num-threads", type=int, default=1)
    parser.add_argument("--expected-formula")
    parser.add_argument("--expected-charge", type=int)
    parser.add_argument("--expected-inchi-key")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    settings = ConformerGenerationSettings(
        num_conformers=args.num_conformers,
        random_seed=args.seed,
        embedding_prune_rms_ang=args.embed_prune_rms,
        deduplicate_rms_ang=args.deduplicate_rms,
        force_field=args.force_field,
        max_iterations=args.max_iterations,
        num_threads=args.num_threads,
    )
    ensemble = generate_conformer_ensemble(args.smiles, settings=settings)
    validate_ensemble_identity(
        ensemble,
        expected_formula=args.expected_formula,
        expected_charge=args.expected_charge,
        expected_inchi_key=args.expected_inchi_key,
    )
    artifacts = write_conformer_artifacts(
        ensemble,
        args.output_dir,
        molecule_name=args.name,
    )

    best = ensemble.best_record
    print(f"Name: {args.name}")
    print(f"Formula: {ensemble.formula}")
    print(f"Formal charge: {ensemble.formal_charge}")
    print(f"InChIKey: {ensemble.inchi_key}")
    print(f"Canonical SMILES: {ensemble.canonical_smiles}")
    print(f"Embedded conformers: {len(ensemble.records)}")
    print(f"Selected conformers: {len(ensemble.selected_records)}")
    print(
        f"Lowest {settings.force_field} energy: "
        f"{best.energy_kcal_mol:.8f} kcal/mol"
    )
    print(f"Lowest-energy XYZ: {artifacts.lowest_energy_xyz_path}")
    print(f"Selected-conformer SDF: {artifacts.selected_sdf_path}")
    print(f"Manifest CSV: {artifacts.manifest_csv_path}")
    print(f"Metadata JSON: {artifacts.metadata_json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
