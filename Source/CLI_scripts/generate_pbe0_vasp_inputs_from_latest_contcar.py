"""Generate PBE0 VASP inputs from final CONTCAR files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Source.vasp_input import VaspInputBuilder, VaspSettings


DEFAULT_SOURCE_ROOT = ROOT / "TestCase" / "VASP" / "01_HA_latest_contcar"
DEFAULT_OUTPUT_ROOT = ROOT / "TestCase" / "VASP" / "02_PBE0"
DEFAULT_INCAR_TEMPLATE = ROOT / "Template" / "VASP_input" / "PBE0_INCAR_1000eV_point_calc"
DEFAULT_KPOINT_DENSITY = 0.05
SUMMARY_NAME = "pbe0_input_generation_summary.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate PBE0 VASP input folders from latest final CONTCAR files."
    )
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--incar-template", type=Path, default=DEFAULT_INCAR_TEMPLATE)
    parser.add_argument(
        "--kpoint-density",
        type=float,
        default=DEFAULT_KPOINT_DENSITY,
        help="K-point spacing multiplier in units of 2*pi Angstrom^-1.",
    )
    return parser.parse_args()


def _latest_contcar(system_dir: Path) -> Path:
    matches = sorted(system_dir.rglob("CONTCAR"), key=lambda path: (len(path.parts), str(path)))
    if not matches:
        raise FileNotFoundError(f"No CONTCAR found below {system_dir}")
    return matches[-1]


def main() -> int:
    args = parse_args()
    source_root = args.source_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    builder = VaspInputBuilder(
        VaspSettings.for_preset(
            "pbe0_point",
            incar_template_path=args.incar_template.resolve(),
            kpoint_spacing_multiplier=args.kpoint_density,
        )
    )

    rows: list[dict[str, str]] = []
    for system_dir in sorted(path for path in source_root.iterdir() if path.is_dir()):
        if system_dir.resolve() == output_root:
            continue
        system_name = system_dir.name
        row = {
            "system_name": system_name,
            "status": "",
            "source_contcar": "",
            "output_dir": "",
            "atom_count": "",
            "kpoints_grid": "",
            "error": "",
            }
        try:
            contcar_path = _latest_contcar(system_dir)
            artifacts = builder.write_job_from_contcar(
                contcar_path,
                output_root / system_name,
                system_name=system_name,
            )
            kpoints_grid = artifacts.kpoints_path.read_text(encoding="utf-8").splitlines()[3].strip()
            row.update(
                {
                    "status": "ok",
                    "source_contcar": str(contcar_path),
                    "output_dir": str(artifacts.output_dir),
                    "atom_count": str(len(artifacts.structure.atoms)),
                    "kpoints_grid": kpoints_grid,
                }
            )
            print(f"OK   {system_name}: {kpoints_grid}")
        except Exception as exc:
            row.update({"status": "error", "error": str(exc)})
            print(f"FAIL {system_name}: {exc}")
        rows.append(row)

    summary_path = output_root / SUMMARY_NAME
    if rows:
        with summary_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    ok_count = sum(row["status"] == "ok" for row in rows)
    error_count = len(rows) - ok_count
    print(f"systems={len(rows)}")
    print(f"generated={ok_count}")
    print(f"errors={error_count}")
    print(f"summary_csv={summary_path}")
    return 0 if error_count == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
