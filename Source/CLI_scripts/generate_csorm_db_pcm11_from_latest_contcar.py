from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source import CSORMInputBuilder, CSORMSettings, CrystalStructure, read_contcar_as_crystal


DEFAULT_RESULTS_CSV = PROJECT_ROOT.parent / "Calculations" / "01_results.csv"
DEFAULT_CONTCAR_ROOT = PROJECT_ROOT / "TestCase" / "VASP" / "01_HA_latest_contcar"
DEFAULT_OPTIMISED_CIF_ROOT = Path("/Users/zianzhan/Desktop/CSP_sandbox/CSP-personal/2_VASP/01_HA/Optimised_CIF")
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT.parent / "Calculations" / "01_CSORM_db_pcm11"
DEFAULT_SUMMARY_CSV = DEFAULT_OUTPUT_ROOT / "csorm_db_pcm11_generation_summary.csv"
DEFAULT_TEMPLATE_DIR = PROJECT_ROOT / "Template" / "CSORM_input"
DEFAULT_POTENTIAL_FILE = DEFAULT_TEMPLATE_DIR / "DB14_pcm11.inter"
DEFAULT_CSO_RM_TEMPLATE = DEFAULT_TEMPLATE_DIR / "CSO_RM.input"

SPECIAL_WYCKOFF_SYSTEMS = {
    "ACXMAC",
    "ASAXUO",
    "AVEMUK",
    "HXBIUR10",
    "MEDZIE",
    "THIOUR08",
}


def _find_latest_contcar_path(contcar_root: Path, system_name: str) -> Path | None:
    matches = sorted((contcar_root / system_name).rglob("CONTCAR"))
    if not matches:
        return None
    return matches[-1]


def _starting_poscar_path(optimised_cif_root: Path, system_name: str) -> Path:
    return optimised_cif_root / f"{system_name}_{system_name}_POSCAR.cif"


def _reduced_molecule_signature(structure: CrystalStructure) -> tuple[int, ...]:
    reduced = structure.reduce_to_asymmetric_unit(symprec=0.05)
    molecules = reduced.detect_molecules()
    return tuple(sorted((len(molecule) for molecule in molecules), reverse=True))


def _valid_numeric(text: str) -> bool:
    try:
        float(text)
    except Exception:
        return False
    return True


def main() -> int:
    results_csv = DEFAULT_RESULTS_CSV
    contcar_root = DEFAULT_CONTCAR_ROOT
    optimised_cif_root = DEFAULT_OPTIMISED_CIF_ROOT
    output_root = DEFAULT_OUTPUT_ROOT
    summary_csv = DEFAULT_SUMMARY_CSV

    output_root.mkdir(parents=True, exist_ok=True)

    for entry in output_root.iterdir():
        if entry.is_dir():
            shutil.rmtree(entry)

    rows: list[dict[str, object]] = []
    with results_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            system_name = row["identifier"].strip()
            notes = (row.get("notes") or "").strip()
            e_latt = (row.get("E_latt") or "").strip()

            summary_row: dict[str, object] = {
                "system_name": system_name,
                "status": "",
                "reason": "",
                "e_latt": e_latt,
                "latest_contcar_path": "",
                "starting_poscar_path": "",
                "starting_poscar_signature": "",
                "reduced_contcar_signature": "",
                "output_dir": "",
                "potential_file": str(DEFAULT_POTENTIAL_FILE),
                "cso_rm_template": str(DEFAULT_CSO_RM_TEMPLATE),
            }

            if not _valid_numeric(e_latt):
                summary_row["status"] = "skipped"
                summary_row["reason"] = "invalid_E_latt"
                rows.append(summary_row)
                continue
            if notes:
                summary_row["status"] = "skipped"
                summary_row["reason"] = "metadata_notes_present"
                rows.append(summary_row)
                continue
            if system_name in SPECIAL_WYCKOFF_SYSTEMS:
                summary_row["status"] = "unsafe"
                summary_row["reason"] = "special_wyckoff_site"
                rows.append(summary_row)
                continue

            contcar_path = _find_latest_contcar_path(contcar_root, system_name)
            if contcar_path is None:
                summary_row["status"] = "unsafe"
                summary_row["reason"] = "missing_latest_contcar"
                rows.append(summary_row)
                continue

            poscar_path = _starting_poscar_path(optimised_cif_root, system_name)
            if not poscar_path.is_file():
                summary_row["status"] = "unsafe"
                summary_row["reason"] = "missing_starting_poscar_cif"
                summary_row["latest_contcar_path"] = str(contcar_path)
                rows.append(summary_row)
                continue

            latest_structure = read_contcar_as_crystal(contcar_path, name=f"{system_name}_CONTCAR")
            reduced_structure = latest_structure.reduce_to_asymmetric_unit(symprec=0.05)
            starting_structure = CrystalStructure.from_file(poscar_path, fmt="cif")

            starting_signature = _reduced_molecule_signature(starting_structure)
            contcar_signature = tuple(
                sorted((len(molecule) for molecule in reduced_structure.detect_molecules()), reverse=True)
            )

            summary_row["latest_contcar_path"] = str(contcar_path)
            summary_row["starting_poscar_path"] = str(poscar_path)
            summary_row["starting_poscar_signature"] = ";".join(str(value) for value in starting_signature)
            summary_row["reduced_contcar_signature"] = ";".join(str(value) for value in contcar_signature)

            if contcar_signature != starting_signature:
                summary_row["status"] = "unsafe"
                summary_row["reason"] = "molecule_signature_mismatch"
                rows.append(summary_row)
                continue

            settings = CSORMSettings(
                cso_rm_template_path=DEFAULT_CSO_RM_TEMPLATE,
                potential_file_path=DEFAULT_POTENTIAL_FILE,
            )
            builder = CSORMInputBuilder(settings=settings)
            output_dir = output_root / system_name
            builder.write_job_from_crystal(reduced_structure, output_dir, system_name=system_name)
            summary_row["status"] = "generated"
            summary_row["reason"] = "ok"
            summary_row["output_dir"] = str(output_dir)
            rows.append(summary_row)

    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    generated = sum(1 for row in rows if row["status"] == "generated")
    unsafe = sum(1 for row in rows if row["status"] == "unsafe")
    skipped = sum(1 for row in rows if row["status"] == "skipped")
    print(f"generated={generated}")
    print(f"unsafe={unsafe}")
    print(f"skipped={skipped}")
    print(f"summary_csv={summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
