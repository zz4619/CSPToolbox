"""Batch-generate reduced-cell CSOFM inputs from testcase CONTCAR structures."""

from __future__ import annotations

import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Source.csofm_input import CSOFMInputBuilder
from Source.gaussian_input import GaussianInputBuilder
from Source.vasp_results import read_contcar_as_crystal


TESTCASE_ROOT = ROOT / "TestCase"
CONTCAR_ROOT = TESTCASE_ROOT / "VASP" / "01_HA_latest_contcar"
GAUSSIAN_LOG_ROOT = TESTCASE_ROOT / "Gaussian" / "hydrate_local_minimisation_gaussian_logs"
OUTPUT_ROOT = TESTCASE_ROOT / "CSOFM_from_CONTCAR"
SUMMARY_CSV = OUTPUT_ROOT / "summary.csv"


def _latest_contcar_path(system_root: Path) -> Path:
    contcars = sorted(
        system_root.rglob("CONTCAR"),
        key=lambda path: (len(path.relative_to(system_root).parts), path.as_posix()),
    )
    if not contcars:
        raise FileNotFoundError(f"No CONTCAR found under {system_root}")
    return contcars[-1]


def _gaussian_log_path(system_name: str) -> Path:
    log_path = GAUSSIAN_LOG_ROOT / system_name / f"{system_name}_1.log"
    if not log_path.is_file():
        raise FileNotFoundError(f"No Gaussian log found for {system_name}: {log_path}")
    return log_path


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    builder = CSOFMInputBuilder()
    builder.write_batch_submission_scripts(OUTPUT_ROOT)
    gaussian_builder = GaussianInputBuilder()
    rows: list[dict[str, object]] = []

    systems = sorted(path for path in CONTCAR_ROOT.iterdir() if path.is_dir())
    for index, system_root in enumerate(systems, start=1):
        system_name = system_root.name
        try:
            contcar_path = _latest_contcar_path(system_root)
            gaussian_log_path = _gaussian_log_path(system_name)

            full_structure = read_contcar_as_crystal(contcar_path, name=system_name)
            reduced_structure = full_structure.reduce_to_asymmetric_unit()
            molecules = reduced_structure.detect_molecules()
            zmatrices = reduced_structure.generate_zmatrices()

            output_dir = OUTPUT_ROOT / system_name
            output_dir.mkdir(parents=True, exist_ok=True)

            reduced_cif_path = output_dir / f"{system_name}_asymmetric_unit.cif"
            reduced_res_path = output_dir / f"{system_name}_asymmetric_unit.res"
            image_path = output_dir / f"{system_name}_unit_cell_molecules.png"
            legacy_zmat_path = output_dir / f"{system_name}_molecules.zmat"
            gaussian_zmat_paths = []

            if legacy_zmat_path.exists():
                legacy_zmat_path.unlink()

            reduced_structure.to_file(reduced_cif_path, fmt="cif")
            reduced_structure.to_file(reduced_res_path, fmt="res")
            reduced_structure.write_unit_cell_molecule_image(image_path, draw_box=True)
            for molecule_index, zmatrix in enumerate(zmatrices, start=1):
                zmat_path = output_dir / f"{system_name}_{molecule_index}.zmat"
                zmat_path.write_text(
                    gaussian_builder.render_zmat_text(f"{system_name} molecule {molecule_index}", zmatrix),
                    encoding="utf-8",
                )
                gaussian_zmat_paths.append(str(zmat_path))

            artifacts = builder.write_job_from_crystal(
                reduced_structure,
                gaussian_log_path,
                output_dir,
                system_name=system_name,
            )

            rows.append(
                {
                    "system_name": system_name,
                    "status": "ok",
                    "contcar_path": str(contcar_path),
                    "gaussian_log_path": str(gaussian_log_path),
                    "space_group_symbol": reduced_structure.space_group,
                    "molecule_count": len(molecules),
                    "zmatrix_count": len(zmatrices),
                    "reduced_atom_count": len(reduced_structure.atoms),
                    "reduced_cif_path": str(reduced_cif_path),
                    "reduced_res_path": str(reduced_res_path),
                    "image_path": str(image_path),
                    "gaussian_zmat_paths": ";".join(gaussian_zmat_paths),
                    "csofm_zmatrix_path": "" if artifacts.zmatrix_path is None else str(artifacts.zmatrix_path),
                    "csofm_dir": str(artifacts.output_dir),
                    "error": "",
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "system_name": system_name,
                    "status": "error",
                    "contcar_path": "",
                    "gaussian_log_path": "",
                    "space_group_symbol": "",
                    "molecule_count": "",
                    "zmatrix_count": "",
                    "reduced_atom_count": "",
                    "reduced_cif_path": "",
                    "reduced_res_path": "",
                    "image_path": "",
                    "gaussian_zmat_paths": "",
                    "csofm_zmatrix_path": "",
                    "csofm_dir": "",
                    "error": str(exc),
                }
            )
        if index % 10 == 0 or index == len(systems):
            print(f"processed={index}/{len(systems)}", flush=True)

    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"systems={len(rows)}")
    print(f"successes={sum(row['status'] == 'ok' for row in rows)}")
    print(f"errors={sum(row['status'] != 'ok' for row in rows)}")
    print(f"summary_csv={SUMMARY_CSV}")

if __name__ == "__main__":
    main()
