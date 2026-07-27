"""Deterministic RDKit conformer generation for molecular CSP inputs.

The generated geometries are force-field-relaxed initial guesses. They are not
intended to replace the subsequent gas-phase quantum-chemical optimization.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path
import re

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors


@dataclass(frozen=True)
class ConformerGenerationSettings:
    """Settings for embedding, optimizing, and deduplicating conformers."""

    num_conformers: int = 100
    random_seed: int = 20260727
    embedding_prune_rms_ang: float = 0.25
    deduplicate_rms_ang: float = 0.35
    force_field: str = "MMFF94s"
    max_iterations: int = 1000
    num_threads: int = 1

    def __post_init__(self) -> None:
        if self.num_conformers < 1:
            raise ValueError("num_conformers must be at least 1.")
        if self.random_seed < 0:
            raise ValueError("random_seed must be non-negative.")
        if self.embedding_prune_rms_ang < 0:
            raise ValueError("embedding_prune_rms_ang must be non-negative.")
        if self.deduplicate_rms_ang < 0:
            raise ValueError("deduplicate_rms_ang must be non-negative.")
        if self.force_field not in {"MMFF94", "MMFF94s", "UFF"}:
            raise ValueError("force_field must be MMFF94, MMFF94s, or UFF.")
        if self.max_iterations < 1:
            raise ValueError("max_iterations must be at least 1.")
        if self.num_threads < 0:
            raise ValueError("num_threads must be non-negative.")


@dataclass(frozen=True)
class ConformerRecord:
    """Optimization and selection result for one embedded conformer."""

    conformer_id: int
    optimization_status: int
    converged: bool
    energy_kcal_mol: float
    relative_energy_kcal_mol: float
    selected: bool = False
    duplicate_of: int | None = None
    heavy_atom_rms_ang: float | None = None


@dataclass(frozen=True)
class ConformerEnsemble:
    """An RDKit molecule and its ranked conformer-generation provenance."""

    molecule: Chem.Mol
    input_smiles: str
    canonical_smiles: str
    formula: str
    formal_charge: int
    inchi_key: str
    settings: ConformerGenerationSettings
    records: tuple[ConformerRecord, ...]

    @property
    def selected_records(self) -> tuple[ConformerRecord, ...]:
        return tuple(record for record in self.records if record.selected)

    @property
    def best_record(self) -> ConformerRecord:
        selected = self.selected_records
        if not selected:
            raise RuntimeError("The conformer ensemble has no selected conformers.")
        return min(selected, key=lambda record: record.energy_kcal_mol)


@dataclass(frozen=True)
class ConformerArtifacts:
    """Files written for a conformer ensemble."""

    output_directory: Path
    lowest_energy_xyz_path: Path
    selected_sdf_path: Path
    manifest_csv_path: Path
    metadata_json_path: Path


def generate_conformer_ensemble(
    smiles: str,
    *,
    settings: ConformerGenerationSettings | None = None,
) -> ConformerEnsemble:
    """Generate a ranked, heavy-atom-RMSD-deduplicated conformer ensemble.

    The input SMILES fixes connectivity, formal charge, stereochemistry, and
    proton placement. This function deliberately does not enumerate tautomers
    or protonation states.
    """

    config = settings or ConformerGenerationSettings()
    parent = Chem.MolFromSmiles(smiles)
    if parent is None:
        raise ValueError(f"RDKit could not parse the supplied SMILES: {smiles!r}")

    Chem.AssignStereochemistry(parent, cleanIt=True, force=True)
    molecule = Chem.AddHs(parent)
    canonical_smiles = Chem.MolToSmiles(
        parent,
        canonical=True,
        isomericSmiles=True,
    )
    formula = rdMolDescriptors.CalcMolFormula(molecule)
    formal_charge = int(Chem.GetFormalCharge(molecule))
    inchi_key = Chem.MolToInchiKey(parent)

    embed_parameters = AllChem.ETKDGv3()
    embed_parameters.randomSeed = config.random_seed
    embed_parameters.pruneRmsThresh = config.embedding_prune_rms_ang
    embed_parameters.enforceChirality = True
    embed_parameters.useSmallRingTorsions = True
    embed_parameters.useMacrocycleTorsions = True
    embed_parameters.numThreads = config.num_threads

    conformer_ids = list(
        AllChem.EmbedMultipleConfs(
            molecule,
            numConfs=config.num_conformers,
            params=embed_parameters,
        )
    )
    if not conformer_ids:
        raise RuntimeError("RDKit ETKDGv3 did not generate any conformers.")

    optimization_results = _optimize_conformers(molecule, config)
    optimized_ids = [conformer.GetId() for conformer in molecule.GetConformers()]
    if len(optimization_results) != len(optimized_ids):
        raise RuntimeError(
            "RDKit returned an unexpected number of force-field optimization results."
        )

    preliminary_records = [
        ConformerRecord(
            conformer_id=conformer_id,
            optimization_status=int(status),
            converged=int(status) == 0,
            energy_kcal_mol=float(energy),
            relative_energy_kcal_mol=0.0,
        )
        for conformer_id, (status, energy) in zip(
            optimized_ids,
            optimization_results,
            strict=True,
        )
    ]
    converged_energies = [
        record.energy_kcal_mol for record in preliminary_records if record.converged
    ]
    if not converged_energies:
        raise RuntimeError(
            f"No {config.force_field} conformer optimization reached convergence."
        )

    minimum_energy = min(converged_energies)
    relative_records = [
        replace(
            record,
            relative_energy_kcal_mol=record.energy_kcal_mol - minimum_energy,
        )
        for record in preliminary_records
    ]
    records = _select_unique_conformers(
        molecule,
        relative_records,
        rms_threshold_ang=config.deduplicate_rms_ang,
    )

    return ConformerEnsemble(
        molecule=molecule,
        input_smiles=smiles,
        canonical_smiles=canonical_smiles,
        formula=formula,
        formal_charge=formal_charge,
        inchi_key=inchi_key,
        settings=config,
        records=tuple(records),
    )


def validate_ensemble_identity(
    ensemble: ConformerEnsemble,
    *,
    expected_formula: str | None = None,
    expected_charge: int | None = None,
    expected_inchi_key: str | None = None,
) -> None:
    """Fail explicitly when generated molecular identity is not the expected one."""

    errors: list[str] = []
    if expected_formula is not None and ensemble.formula != expected_formula:
        errors.append(
            f"formula is {ensemble.formula}, expected {expected_formula}"
        )
    if expected_charge is not None and ensemble.formal_charge != expected_charge:
        errors.append(
            f"formal charge is {ensemble.formal_charge}, expected {expected_charge}"
        )
    if (
        expected_inchi_key is not None
        and ensemble.inchi_key.upper() != expected_inchi_key.upper()
    ):
        errors.append(
            f"InChIKey is {ensemble.inchi_key}, expected {expected_inchi_key}"
        )
    if errors:
        raise ValueError("Molecular identity validation failed: " + "; ".join(errors))


def write_conformer_artifacts(
    ensemble: ConformerEnsemble,
    output_directory: str | Path,
    *,
    molecule_name: str = "molecule",
) -> ConformerArtifacts:
    """Write the selected ensemble and its lowest-energy Cartesian geometry."""

    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)
    safe_name = _safe_name(molecule_name)

    xyz_path = output_path / f"{safe_name}_lowest_{ensemble.settings.force_field.lower()}.xyz"
    sdf_path = output_path / f"{safe_name}_selected_conformers.sdf"
    csv_path = output_path / f"{safe_name}_conformer_manifest.csv"
    json_path = output_path / f"{safe_name}_conformer_metadata.json"

    xyz_path.write_text(
        _render_xyz(ensemble, ensemble.best_record, molecule_name),
        encoding="utf-8",
    )
    _write_selected_sdf(ensemble, sdf_path, molecule_name)
    _write_manifest_csv(ensemble, csv_path)
    json_path.write_text(
        json.dumps(_metadata_payload(ensemble, molecule_name, xyz_path, sdf_path), indent=2)
        + "\n",
        encoding="utf-8",
    )

    return ConformerArtifacts(
        output_directory=output_path,
        lowest_energy_xyz_path=xyz_path,
        selected_sdf_path=sdf_path,
        manifest_csv_path=csv_path,
        metadata_json_path=json_path,
    )


def _optimize_conformers(
    molecule: Chem.Mol,
    settings: ConformerGenerationSettings,
) -> list[tuple[int, float]]:
    if settings.force_field.startswith("MMFF"):
        if not AllChem.MMFFHasAllMoleculeParams(molecule):
            raise ValueError(
                f"RDKit does not have all {settings.force_field} parameters "
                "for the supplied molecule."
            )
        return list(
            AllChem.MMFFOptimizeMoleculeConfs(
                molecule,
                numThreads=settings.num_threads,
                maxIters=settings.max_iterations,
                mmffVariant=settings.force_field,
            )
        )

    if not AllChem.UFFHasAllMoleculeParams(molecule):
        raise ValueError(
            "RDKit does not have all UFF parameters for the supplied molecule."
        )
    return list(
        AllChem.UFFOptimizeMoleculeConfs(
            molecule,
            numThreads=settings.num_threads,
            maxIters=settings.max_iterations,
        )
    )


def _select_unique_conformers(
    molecule: Chem.Mol,
    records: list[ConformerRecord],
    *,
    rms_threshold_ang: float,
) -> list[ConformerRecord]:
    heavy_molecule = Chem.RemoveHs(molecule)
    selected_ids: list[int] = []
    selected_records: dict[int, ConformerRecord] = {}

    for record in sorted(
        records,
        key=lambda item: (
            not item.converged,
            item.energy_kcal_mol,
            item.conformer_id,
        ),
    ):
        if not record.converged:
            selected_records[record.conformer_id] = record
            continue

        duplicate_id: int | None = None
        duplicate_rms: float | None = None
        for selected_id in selected_ids:
            probe = Chem.Mol(heavy_molecule)
            rms = float(
                rdMolAlign.GetBestRMS(
                    probe,
                    heavy_molecule,
                    prbId=record.conformer_id,
                    refId=selected_id,
                )
            )
            if duplicate_rms is None or rms < duplicate_rms:
                duplicate_rms = rms
                duplicate_id = selected_id
            if rms < rms_threshold_ang:
                break

        if duplicate_rms is not None and duplicate_rms < rms_threshold_ang:
            selected_records[record.conformer_id] = replace(
                record,
                duplicate_of=duplicate_id,
                heavy_atom_rms_ang=duplicate_rms,
            )
            continue

        selected_ids.append(record.conformer_id)
        selected_records[record.conformer_id] = replace(
            record,
            selected=True,
            heavy_atom_rms_ang=duplicate_rms,
        )

    return [
        selected_records[record.conformer_id]
        for record in sorted(
            records,
            key=lambda item: (
                not item.converged,
                item.energy_kcal_mol,
                item.conformer_id,
            ),
        )
    ]


def _render_xyz(
    ensemble: ConformerEnsemble,
    record: ConformerRecord,
    molecule_name: str,
) -> str:
    conformer = ensemble.molecule.GetConformer(record.conformer_id)
    lines = [
        str(ensemble.molecule.GetNumAtoms()),
        (
            f"{molecule_name}; conformer_id={record.conformer_id}; "
            f"{ensemble.settings.force_field}_energy_kcal_mol="
            f"{record.energy_kcal_mol:.8f}"
        ),
    ]
    for atom in ensemble.molecule.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{atom.GetSymbol():<2} "
            f"{position.x: .10f} {position.y: .10f} {position.z: .10f}"
        )
    return "\n".join(lines) + "\n"


def _write_selected_sdf(
    ensemble: ConformerEnsemble,
    path: Path,
    molecule_name: str,
) -> None:
    writer = Chem.SDWriter(str(path))
    try:
        for rank, record in enumerate(ensemble.selected_records, start=1):
            record_molecule = Chem.Mol(ensemble.molecule)
            conformer = Chem.Conformer(
                ensemble.molecule.GetConformer(record.conformer_id)
            )
            record_molecule.RemoveAllConformers()
            record_molecule.AddConformer(conformer, assignId=True)
            record_molecule.SetProp("_Name", f"{molecule_name}_conformer_{rank}")
            record_molecule.SetIntProp("ConformerID", record.conformer_id)
            record_molecule.SetDoubleProp(
                "Energy_kcal_mol",
                record.energy_kcal_mol,
            )
            record_molecule.SetDoubleProp(
                "RelativeEnergy_kcal_mol",
                record.relative_energy_kcal_mol,
            )
            record_molecule.SetProp("ForceField", ensemble.settings.force_field)
            record_molecule.SetProp("CanonicalSMILES", ensemble.canonical_smiles)
            record_molecule.SetProp("Formula", ensemble.formula)
            record_molecule.SetProp("InChIKey", ensemble.inchi_key)
            writer.write(record_molecule)
    finally:
        writer.close()


def _write_manifest_csv(ensemble: ConformerEnsemble, path: Path) -> None:
    fieldnames = [
        "conformer_id",
        "optimization_status",
        "converged",
        "energy_kcal_mol",
        "relative_energy_kcal_mol",
        "selected",
        "duplicate_of",
        "heavy_atom_rms_ang",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in ensemble.records:
            writer.writerow(asdict(record))


def _metadata_payload(
    ensemble: ConformerEnsemble,
    molecule_name: str,
    xyz_path: Path,
    sdf_path: Path,
) -> dict[str, object]:
    return {
        "molecule_name": molecule_name,
        "input_smiles": ensemble.input_smiles,
        "canonical_smiles": ensemble.canonical_smiles,
        "formula": ensemble.formula,
        "formal_charge": ensemble.formal_charge,
        "inchi_key": ensemble.inchi_key,
        "atom_count": ensemble.molecule.GetNumAtoms(),
        "heavy_atom_count": ensemble.molecule.GetNumHeavyAtoms(),
        "settings": asdict(ensemble.settings),
        "embedded_conformer_count": len(ensemble.records),
        "selected_conformer_count": len(ensemble.selected_records),
        "best_conformer": asdict(ensemble.best_record),
        "lowest_energy_xyz": xyz_path.name,
        "selected_conformer_sdf": sdf_path.name,
        "energy_note": (
            "Force-field energies rank initial guesses only; they are not "
            "gas-phase quantum-chemical energies."
        ),
        "tautomer_note": (
            "No tautomer or protonation-state enumeration was performed. "
            "The supplied SMILES fixes molecular identity."
        ),
    }


def _safe_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return cleaned or "molecule"
