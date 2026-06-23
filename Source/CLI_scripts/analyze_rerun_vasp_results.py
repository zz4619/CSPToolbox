from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.vasp_results import VaspSystemParser


DEFAULT_INPUT_ROOT = PROJECT_ROOT.parent / "Calculations" / "03_Rerun"
DEFAULT_OUTPUT_CSV = DEFAULT_INPUT_ROOT / "vasp_summary.csv"
DEFAULT_IMAGE_DIR = DEFAULT_INPUT_ROOT / "molecule_images"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze rerun VASP results and render one molecule image per CONTCAR."
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help=f"Directory containing one subdirectory per rerun system. Default: {DEFAULT_INPUT_ROOT}",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=DEFAULT_OUTPUT_CSV,
        help=f"Destination CSV path. Default: {DEFAULT_OUTPUT_CSV}",
    )
    parser.add_argument(
        "--image-dir",
        type=Path,
        default=DEFAULT_IMAGE_DIR,
        help=f"Directory for one PNG per system. Default: {DEFAULT_IMAGE_DIR}",
    )
    return parser.parse_args()


def analyze_systems(input_root: Path, image_dir: Path) -> list[dict[str, object]]:
    parser = VaspSystemParser()
    image_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for system_dir in sorted(path for path in input_root.iterdir() if path.is_dir()):
        vasprun_path = system_dir / "vasprun.xml"
        contcar_path = system_dir / "CONTCAR"
        if not vasprun_path.is_file() or not contcar_path.is_file():
            continue
        try:
            system = parser.parse_system(system_dir)
            latest = system.latest_calculation
            if latest is None:
                continue

            crystal = system.read_latest_contcar_as_crystal()
            molecules = crystal.detect_molecules()
            image_path = image_dir / f"{system.system_name}.png"
            crystal.write_unit_cell_molecule_image(image_path, draw_box=True)

            rows.append(
                {
                    "system_name": system.system_name,
                    "calculation_count": system.calculation_count,
                    "status": latest.status,
                    "ionic_iterations": latest.ionic_iterations,
                    "total_electronic_iterations": latest.total_electronic_iterations,
                    "last_electronic_iterations": latest.last_electronic_iterations,
                    "final_energy_ev": latest.final_energy_ev,
                    "converged": latest.converged,
                    "converged_ionic": latest.converged_ionic,
                    "converged_electronic": latest.converged_electronic,
                    "cpu_time_seconds": latest.cpu_time_seconds,
                    "molecule_count": len(molecules),
                    "atom_count": len(crystal.atoms),
                    "image_path": str(image_path),
                    "contcar_path": str(contcar_path),
                    "vasprun_path": str(vasprun_path),
                    "parse_error": latest.parse_error or "",
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "system_name": system_dir.name,
                    "calculation_count": "",
                    "status": "error",
                    "ionic_iterations": "",
                    "total_electronic_iterations": "",
                    "last_electronic_iterations": "",
                    "final_energy_ev": "",
                    "converged": "",
                    "converged_ionic": "",
                    "converged_electronic": "",
                    "cpu_time_seconds": "",
                    "molecule_count": "",
                    "atom_count": "",
                    "image_path": "",
                    "contcar_path": str(contcar_path),
                    "vasprun_path": str(vasprun_path),
                    "parse_error": str(exc),
                }
            )

    return rows


def write_summary(rows: list[dict[str, object]], output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("No rerun VASP rows were generated.")

    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_root = args.input_root.resolve()
    output_csv = args.output_csv.resolve()
    image_dir = args.image_dir.resolve()

    if not input_root.is_dir():
        raise FileNotFoundError(f"Input root does not exist: {input_root}")

    rows = analyze_systems(input_root, image_dir)
    write_summary(rows, output_csv)

    print(f"systems={len(rows)}")
    print(f"summary_csv={output_csv}")
    print(f"image_dir={image_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
