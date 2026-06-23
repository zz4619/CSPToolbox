from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import CrystalStructure


DEFAULT_INPUT_ROOT = Path(
    "/Users/zianzhan/Desktop/CSP_sandbox/CSP-personal/2_VASP/01_HA/Optimised_CIF"
)
DEFAULT_OUTPUT_DIR = DEFAULT_INPUT_ROOT.parent / "Optimised_CIF_molecule_images"
DEFAULT_PATTERN = "*_CONTCAR.cif"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read Optimised_CIF CONTCAR CIF files with CSPToolbox and render one "
            "PNG per detected molecule."
        )
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help=f"Directory containing Optimised_CIF files. Default: {DEFAULT_INPUT_ROOT}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for molecule-image outputs. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--pattern",
        default=DEFAULT_PATTERN,
        help=f"Glob used to select input CIF files. Default: {DEFAULT_PATTERN}",
    )
    parser.add_argument(
        "--reduce-to-asymmetric-unit",
        action="store_true",
        help="Reduce each CIF to the asymmetric unit before molecule detection.",
    )
    parser.add_argument(
        "--symprec",
        type=float,
        default=0.05,
        help="Symmetry tolerance used with --reduce-to-asymmetric-unit. Default: 0.05",
    )
    return parser.parse_args()


def _collect_input_files(input_root: Path, pattern: str) -> list[Path]:
    files = sorted(path for path in input_root.glob(pattern) if path.is_file())
    if not files:
        raise FileNotFoundError(f"No CIF files matching {pattern!r} found under {input_root}")
    return files


def _prepare_structure(
    cif_path: Path,
    *,
    reduce_to_asymmetric_unit: bool,
    symprec: float,
) -> CrystalStructure:
    structure = CrystalStructure.from_file(cif_path, fmt="cif")
    if reduce_to_asymmetric_unit:
        return structure.reduce_to_asymmetric_unit(symprec=symprec)
    return structure


def render_file(
    cif_path: Path,
    *,
    output_dir: Path,
    reduce_to_asymmetric_unit: bool,
    symprec: float,
) -> dict[str, object]:
    structure = _prepare_structure(
        cif_path,
        reduce_to_asymmetric_unit=reduce_to_asymmetric_unit,
        symprec=symprec,
    )
    molecules = structure.detect_molecules()
    if not molecules:
        raise ValueError("No molecules detected.")

    structure_output_dir = output_dir / cif_path.stem
    structure_output_dir.mkdir(parents=True, exist_ok=True)

    image_paths: list[str] = []
    for molecule_index, molecule in enumerate(molecules, start=1):
        molecule_name = f"{cif_path.stem}_mol{molecule_index:02d}"
        image_path = structure_output_dir / f"{molecule_name}.png"
        gas_phase = structure.generate_gas_phase_vasp_structure(
            molecule,
            name=molecule_name,
        )
        gas_phase.write_unit_cell_molecule_image(
            image_path,
            title=molecule_name,
            draw_box=False,
        )
        image_paths.append(str(image_path))

    return {
        "input_cif": str(cif_path),
        "structure_name": cif_path.stem,
        "reduced_to_asymmetric_unit": reduce_to_asymmetric_unit,
        "atom_count_used_for_detection": len(structure.atoms),
        "molecule_count": len(molecules),
        "molecule_sizes": ";".join(str(len(molecule)) for molecule in molecules),
        "output_dir": str(structure_output_dir),
        "image_paths": ";".join(image_paths),
        "error": "",
    }


def render_all(
    *,
    input_root: Path,
    output_dir: Path,
    pattern: str,
    reduce_to_asymmetric_unit: bool,
    symprec: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    input_files = _collect_input_files(input_root, pattern)
    output_dir.mkdir(parents=True, exist_ok=True)

    for index, cif_path in enumerate(input_files, start=1):
        try:
            row = render_file(
                cif_path,
                output_dir=output_dir,
                reduce_to_asymmetric_unit=reduce_to_asymmetric_unit,
                symprec=symprec,
            )
            rows.append(row)
            print(
                f"OK   {index:04d}/{len(input_files):04d} {cif_path.name} "
                f"molecules={row['molecule_count']}"
            )
        except Exception as error:
            rows.append(
                {
                    "input_cif": str(cif_path),
                    "structure_name": cif_path.stem,
                    "reduced_to_asymmetric_unit": reduce_to_asymmetric_unit,
                    "atom_count_used_for_detection": "",
                    "molecule_count": "",
                    "molecule_sizes": "",
                    "output_dir": "",
                    "image_paths": "",
                    "error": str(error),
                }
            )
            print(f"FAIL {index:04d}/{len(input_files):04d} {cif_path.name}: {error}")

    return rows


def write_summary(rows: list[dict[str, object]], output_dir: Path) -> Path:
    if not rows:
        raise ValueError("No rows to write.")

    summary_path = output_dir / "summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return summary_path


def main() -> int:
    args = parse_args()
    input_root = args.input_root.resolve()
    output_dir = args.output_dir.resolve()

    if not input_root.is_dir():
        raise FileNotFoundError(f"Input root does not exist: {input_root}")

    rows = render_all(
        input_root=input_root,
        output_dir=output_dir,
        pattern=args.pattern,
        reduce_to_asymmetric_unit=args.reduce_to_asymmetric_unit,
        symprec=args.symprec,
    )
    summary_path = write_summary(rows, output_dir)

    successes = sum(not row["error"] for row in rows)
    failures = len(rows) - successes
    print(f"files={len(rows)}")
    print(f"successes={successes}")
    print(f"failures={failures}")
    print(f"output_dir={output_dir}")
    print(f"summary_csv={summary_path}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
