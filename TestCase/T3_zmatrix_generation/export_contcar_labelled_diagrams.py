"""Export labelled molecular diagrams from HA-pair TPSS CONTCAR files."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import sys


CASE_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = CASE_ROOT.parents[1]
DEFAULT_INPUT_ROOT = (
    PROJECT_ROOT.parent
    / "CSP-personal"
    / "Results"
    / "VASP_HA_pair_TPSS_converged_structure"
)
DEFAULT_OUTPUT_ROOT = CASE_ROOT / "HA_pair_TPSS_CONTCAR_labelled_diagrams"

CACHE_ROOT = Path("/private/tmp/csptoolbox_t3_contcar_diagram_cache")
os.environ.setdefault("MPLCONFIGDIR", str(CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_ROOT / "xdg"))
for cache_dir in (Path(os.environ["MPLCONFIGDIR"]), Path(os.environ["XDG_CACHE_HOME"])):
    cache_dir.mkdir(parents=True, exist_ok=True)

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.vasp_results import read_contcar_as_crystal  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read HA-pair TPSS CONTCAR files and export labelled molecule diagrams "
            "for selecting proper dihedral definitions."
        )
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help=f"Root containing SYSTEM/CONTCAR files. Default: {DEFAULT_INPUT_ROOT}",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_ROOT}",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional maximum number of systems to process.",
    )
    parser.add_argument(
        "--draw-box",
        action="store_true",
        help="Draw the unit-cell box on full-cell diagrams.",
    )
    return parser.parse_args()


def contcar_paths(input_root: Path) -> list[Path]:
    return sorted(input_root.glob("*/CONTCAR"))


def export_diagrams(
    *,
    input_root: Path,
    output_root: Path,
    limit: int | None = None,
    draw_box: bool = False,
) -> list[dict[str, object]]:
    output_root.mkdir(parents=True, exist_ok=True)
    unit_cell_dir = output_root / "unit_cells"
    full_diagram_dir = output_root / "unit_cell_diagrams"
    molecule_diagram_dir = output_root / "molecule_diagrams"
    unit_cell_dir.mkdir(parents=True, exist_ok=True)
    full_diagram_dir.mkdir(parents=True, exist_ok=True)
    molecule_diagram_dir.mkdir(parents=True, exist_ok=True)

    paths = contcar_paths(input_root)
    if limit is not None:
        paths = paths[:limit]

    rows: list[dict[str, object]] = []
    for index, contcar_path in enumerate(paths, start=1):
        system_name = contcar_path.parent.name
        row: dict[str, object] = {
            "system_name": system_name,
            "status": "ok",
            "contcar_path": str(contcar_path),
            "atom_count": "",
            "molecule_count": "",
            "molecule_sizes": "",
            "unit_cell_cif": "",
            "unit_cell_diagram": "",
            "molecule_diagrams": "",
            "error": "",
        }

        try:
            structure = read_contcar_as_crystal(contcar_path, name=system_name)
            molecules = structure.detect_molecules()

            unit_cell_cif = unit_cell_dir / f"{system_name}_unit_cell.cif"
            unit_cell_diagram = full_diagram_dir / f"{system_name}_unit_cell_labelled.png"
            structure.to_file(unit_cell_cif, fmt="cif")
            structure.write_unit_cell_molecule_image(
                unit_cell_diagram,
                title=f"{system_name} unit cell",
                draw_box=draw_box,
            )

            molecule_paths: list[str] = []
            for molecule_index, molecule in enumerate(molecules, start=1):
                gas_structure = structure.generate_gas_phase_vasp_structure(
                    molecule,
                    name=f"{system_name}_molecule_{molecule_index:02d}",
                )
                molecule_path = (
                    molecule_diagram_dir
                    / f"{system_name}_molecule_{molecule_index:02d}_labelled.png"
                )
                gas_structure.write_unit_cell_molecule_image(
                    molecule_path,
                    title=f"{system_name} molecule {molecule_index}",
                    draw_box=False,
                )
                molecule_paths.append(str(molecule_path))

            row.update(
                {
                    "atom_count": len(structure.atoms),
                    "molecule_count": len(molecules),
                    "molecule_sizes": ";".join(str(len(molecule)) for molecule in molecules),
                    "unit_cell_cif": str(unit_cell_cif),
                    "unit_cell_diagram": str(unit_cell_diagram),
                    "molecule_diagrams": ";".join(molecule_paths),
                }
            )
        except Exception as exc:
            row["status"] = "error"
            row["error"] = str(exc)

        rows.append(row)
        if index % 10 == 0 or index == len(paths):
            print(f"processed={index}/{len(paths)}", flush=True)

    return rows


def write_manifest(output_root: Path, rows: list[dict[str, object]]) -> Path:
    manifest_path = output_root / "labelled_diagram_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return manifest_path


def main() -> int:
    args = parse_args()
    rows = export_diagrams(
        input_root=args.input_root.resolve(),
        output_root=args.output_root.resolve(),
        limit=args.limit,
        draw_box=args.draw_box,
    )
    if not rows:
        raise FileNotFoundError(f"No CONTCAR files found under {args.input_root.resolve()}")

    manifest_path = write_manifest(args.output_root.resolve(), rows)
    success_count = sum(row["status"] == "ok" for row in rows)
    error_count = len(rows) - success_count

    print(f"systems={len(rows)}")
    print(f"successes={success_count}")
    print(f"errors={error_count}")
    print(f"output_root={args.output_root.resolve()}")
    print(f"manifest={manifest_path}")

    return 1 if error_count else 0


if __name__ == "__main__":
    raise SystemExit(main())

