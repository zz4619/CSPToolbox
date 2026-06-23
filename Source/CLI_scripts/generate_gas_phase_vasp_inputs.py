"""Generate gas-phase VASP inputs for non-water molecules from CSOFM structures."""

from __future__ import annotations

import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Source.crystal_structure import AtomRecord, CrystalStructure
from Source.vasp_input import VaspInputBuilder


SOURCE_ROOT = ROOT / "TestCase" / "CSOFM_from_CONTCAR"
OUTPUT_ROOT = ROOT.parent / "Calculations" / "01_HA_vasp_gas_phase_energy"
SUMMARY_CSV = OUTPUT_ROOT / "summary.csv"
BOX_LENGTH = 20.0


def _is_water_molecule(molecule: list[AtomRecord]) -> bool:
    counts: dict[str, int] = {}
    for atom in molecule:
        counts[atom.element] = counts.get(atom.element, 0) + 1
    return counts == {"H": 2, "O": 1}


def _source_structure_path(system_dir: Path) -> Path:
    structure_path = system_dir / f"{system_dir.name}.res"
    if not structure_path.is_file():
        raise FileNotFoundError(f"Missing system .res file: {structure_path}")
    return structure_path


def _job_name(system_name: str, molecule_index: int, total_jobs: int) -> str:
    if total_jobs == 1:
        return system_name
    return f"{system_name}_Mol{molecule_index}"


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    builder = VaspInputBuilder()
    rows: list[dict[str, object]] = []

    systems = sorted(path for path in SOURCE_ROOT.iterdir() if path.is_dir())
    for system_dir in systems:
        system_name = system_dir.name
        try:
            structure = CrystalStructure.from_file(_source_structure_path(system_dir))
            groups = structure.deduplicate_molecules()
            non_water_groups = [group for group in groups if not _is_water_molecule(group.representative_molecule)]

            if not non_water_groups:
                rows.append(
                    {
                        "system_name": system_name,
                        "job_name": "",
                        "status": "skipped_water_only",
                        "source_structure": str(system_dir / f"{system_name}.res"),
                        "output_dir": "",
                        "atom_count": "",
                        "error": "",
                    }
                )
                continue

            for molecule_index, group in enumerate(non_water_groups, start=1):
                job_name = _job_name(system_name, molecule_index, len(non_water_groups))
                gas_phase_structure = structure.generate_gas_phase_vasp_structure(
                    group.representative_molecule,
                    box_length=BOX_LENGTH,
                    name=job_name,
                )
                output_dir = OUTPUT_ROOT / job_name
                artifacts = builder.write_gas_phase_job(
                    gas_phase_structure,
                    output_dir,
                    system_name=job_name,
                )
                rows.append(
                    {
                        "system_name": system_name,
                        "job_name": job_name,
                        "status": "ok",
                        "source_structure": str(system_dir / f"{system_name}.res"),
                        "output_dir": str(artifacts.output_dir),
                        "atom_count": len(gas_phase_structure.atoms),
                        "error": "",
                    }
                )
        except Exception as exc:
            rows.append(
                {
                    "system_name": system_name,
                    "job_name": "",
                    "status": "error",
                    "source_structure": str(system_dir / f"{system_name}.res"),
                    "output_dir": "",
                    "atom_count": "",
                    "error": str(exc),
                }
            )

    if rows:
        with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    print(f"systems={len(systems)}")
    print(f"jobs={sum(row['status'] == 'ok' for row in rows)}")
    print(f"skipped={sum(row['status'] == 'skipped_water_only' for row in rows)}")
    print(f"errors={sum(row['status'] == 'error' for row in rows)}")
    print(f"summary_csv={SUMMARY_CSV}")


if __name__ == "__main__":
    main()
