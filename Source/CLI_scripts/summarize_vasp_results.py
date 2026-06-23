from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.vasp_results import VaspSystemParser, summarize_system_rows


DEFAULT_INPUT_ROOT = PROJECT_ROOT / "TestCase" / "VASP" / "01_HA_vasprun_xml"
DEFAULT_OUTPUT_CSV = PROJECT_ROOT / "TestCase" / "VASP" / "01_HA_vasp_summary.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize nested VASP vasprun.xml results for all systems into one CSV."
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help=f"Directory containing one subdirectory per system. Default: {DEFAULT_INPUT_ROOT}",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=DEFAULT_OUTPUT_CSV,
        help=f"Destination CSV path. Default: {DEFAULT_OUTPUT_CSV}",
    )
    return parser.parse_args()


def build_summary(input_root: Path) -> list[dict[str, object]]:
    parser = VaspSystemParser()
    rows: list[dict[str, object]] = []

    for system_dir in sorted(path for path in input_root.iterdir() if path.is_dir()):
        if not any(system_dir.rglob("vasprun.xml")):
            continue
        system = parser.parse_system(system_dir)
        system_rows = summarize_system_rows(system)
        for row_index, row in enumerate(system_rows):
            row["calculation_count"] = system.calculation_count
            row["is_latest_calculation"] = row_index == len(system_rows) - 1
            row["system_dir"] = str(system.system_dir)
            row["vasprun_path"] = str(system.calculations[row_index].vasprun_path)
            row["calculation_dir"] = str(system.calculations[row_index].calculation_dir)
        rows.extend(system_rows)

    return rows


def write_summary(rows: list[dict[str, object]], output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("No VASP result rows were generated.")

    fieldnames = [
        "system_name",
        "calculation_count",
        "calculation_index",
        "calculation_chain",
        "is_latest_calculation",
        "status",
        "ionic_iterations",
        "total_electronic_iterations",
        "last_electronic_iterations",
        "final_energy_ev",
        "converged",
        "converged_ionic",
        "converged_electronic",
        "cpu_time_seconds",
        "parse_error",
        "system_dir",
        "calculation_dir",
        "vasprun_path",
    ]

    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_root = args.input_root.resolve()
    output_csv = args.output_csv.resolve()

    if not input_root.is_dir():
        raise FileNotFoundError(f"Input root does not exist: {input_root}")

    rows = build_summary(input_root)
    write_summary(rows, output_csv)

    print(f"Wrote {len(rows)} rows to {output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
