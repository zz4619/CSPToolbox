from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import CrystalStructure
from Source.vasp_input import VaspInputBuilder


DEFAULT_PRIMARY_CIF_ROOT = PROJECT_ROOT / "TestCase" / "Structures" / "add_H_cif_file"
DEFAULT_FALLBACK_CIF_ROOT = PROJECT_ROOT / "TestCase" / "Structures" / "cif_file"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT.parent / "Calculations" / "01_Rerun"
DEFAULT_SUMMARY_NAME = "vasp_input_generation_summary.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate VASP input structures from CIF files, preferring one source "
            "directory and falling back to another when needed."
        )
    )
    parser.add_argument(
        "--primary-cif-root",
        type=Path,
        default=DEFAULT_PRIMARY_CIF_ROOT,
        help=f"Preferred CIF directory. Default: {DEFAULT_PRIMARY_CIF_ROOT}",
    )
    parser.add_argument(
        "--fallback-cif-root",
        type=Path,
        default=DEFAULT_FALLBACK_CIF_ROOT,
        help=f"Fallback CIF directory. Default: {DEFAULT_FALLBACK_CIF_ROOT}",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help=f"Directory for generated VASP inputs. Default: {DEFAULT_OUTPUT_ROOT}",
    )
    parser.add_argument(
        "--summary-name",
        default=DEFAULT_SUMMARY_NAME,
        help=f"Summary CSV filename written under the output root. Default: {DEFAULT_SUMMARY_NAME}",
    )
    parser.add_argument(
        "systems",
        nargs="+",
        help="System identifiers to generate.",
    )
    return parser.parse_args()


def _select_cif_path(system_name: str, primary_root: Path, fallback_root: Path) -> tuple[Path, str]:
    primary = primary_root / f"{system_name}.cif"
    if primary.is_file():
        return primary, "primary"

    fallback = fallback_root / f"{system_name}.cif"
    if fallback.is_file():
        return fallback, "fallback"

    raise FileNotFoundError(
        f"No CIF found for {system_name} in {primary_root} or {fallback_root}"
    )


def main() -> int:
    args = parse_args()
    builder = VaspInputBuilder()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / args.summary_name

    rows: list[dict[str, str]] = []
    for system_name in args.systems:
        row = {
            "system_name": system_name,
            "status": "",
            "source_priority": "",
            "cif_path": "",
            "output_dir": "",
            "atom_count": "",
            "space_group": "",
            "error": "",
        }
        try:
            cif_path, source_priority = _select_cif_path(
                system_name,
                args.primary_cif_root.resolve(),
                args.fallback_cif_root.resolve(),
            )
            structure = CrystalStructure.expand_cif_to_unit_cell(cif_path)
            artifacts = builder.write_job_from_structure(
                structure,
                output_root / system_name,
                system_name=system_name,
            )
            row.update(
                {
                    "status": "ok",
                    "source_priority": source_priority,
                    "cif_path": str(cif_path),
                    "output_dir": str(artifacts.output_dir),
                    "atom_count": str(len(structure.atoms)),
                    "space_group": structure.space_group,
                }
            )
            print(f"OK   {system_name}: {source_priority} -> {cif_path}")
        except Exception as exc:
            row.update(
                {
                    "status": "error",
                    "error": str(exc),
                }
            )
            print(f"FAIL {system_name}: {exc}")
        rows.append(row)

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    ok_count = sum(row["status"] == "ok" for row in rows)
    error_count = len(rows) - ok_count
    print(f"generated={ok_count}")
    print(f"errors={error_count}")
    print(f"summary_csv={summary_path}")
    return 0 if error_count == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
