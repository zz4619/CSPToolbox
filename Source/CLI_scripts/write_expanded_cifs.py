from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import CrystalStructure


DEFAULT_INPUT_ROOT = PROJECT_ROOT / "TestCase" / "Structures"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "TestCase" / "ExpandedCIF"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Expand testcase CIF files to full unit-cell structures and write them "
            "back out as P1 CIF files."
        )
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help=f"Directory to scan for input CIF files. Default: {DEFAULT_INPUT_ROOT}",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help=f"Directory for expanded CIF outputs. Default: {DEFAULT_OUTPUT_ROOT}",
    )
    return parser.parse_args()


def write_all(input_root: Path, output_root: Path) -> tuple[list[Path], list[tuple[Path, str]]]:
    cif_files = sorted(input_root.rglob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No CIF files found under {input_root}")

    successes: list[Path] = []
    failures: list[tuple[Path, str]] = []

    for cif_path in cif_files:
        relative_path = cif_path.relative_to(input_root)
        output_path = output_root / relative_path.parent / f"{cif_path.stem}_expanded.cif"
        try:
            expanded = CrystalStructure.from_cif_unit_cell(cif_path)
            output_structure = CrystalStructure(
                atoms=expanded.atoms,
                cell_parameters=expanded.cell_parameters,
                lattice_matrix=expanded.lattice_matrix,
                space_group="P 1",
                name=f"{expanded.name}_expanded",
                explict_unit_cell=True,
            )
            output_structure.to_file(output_path, fmt="cif")
            successes.append(relative_path)
            print(f"OK   {relative_path} -> {output_path}")
        except Exception as error:
            failures.append((relative_path, str(error)))
            print(f"FAIL {relative_path}: {error}")

    return successes, failures


def main() -> int:
    args = parse_args()
    successes, failures = write_all(
        input_root=args.input_root.resolve(),
        output_root=args.output_root.resolve(),
    )

    print("\nSummary")
    print(f"  success: {len(successes)}")
    print(f"  failed:  {len(failures)}")
    print(f"  output:  {args.output_root.resolve()}")

    if failures:
        print("Failures:")
        for path, message in failures[:20]:
            print(f"  {path}: {message}")
        if len(failures) > 20:
            print(f"  ... and {len(failures) - 20} more")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
