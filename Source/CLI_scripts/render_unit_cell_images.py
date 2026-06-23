from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import CrystalStructure


DEFAULT_INPUT_DIR = Path(
    "/Users/zianzhan/Desktop/CSP_sandbox/CSPToolbox/TestCase/Structures/pdb_file"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/zianzhan/Desktop/CSP_sandbox/CSPToolbox/TestCase/Structures/pdb_molecule_images"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse PDB files, detect molecules in the unit cell, and write one "
            "image per structure containing all detected molecules."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Directory containing input PDB files. Default: {DEFAULT_INPUT_DIR}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for PNG outputs. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--draw-box",
        action="store_true",
        help="Draw the projected unit-cell box behind the molecules.",
    )
    return parser.parse_args()


def render_all(input_dir: Path, output_dir: Path, draw_box: bool) -> tuple[list[str], list[tuple[str, str]]]:
    pdb_files = sorted(input_dir.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found under {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    successes: list[str] = []
    failures: list[tuple[str, str]] = []

    for pdb_path in pdb_files:
        output_path = output_dir / f"{pdb_path.stem}.png"
        try:
            structure = CrystalStructure.from_file(pdb_path, fmt="pdb")
            structure.detect_molecules()
            structure.write_unit_cell_molecule_image(
                output_path,
                draw_box=draw_box,
                title=pdb_path.stem,
            )
            successes.append(pdb_path.name)
            print(f"OK   {pdb_path.name} -> {output_path.name}")
        except Exception as error:
            failures.append((pdb_path.name, str(error)))
            print(f"FAIL {pdb_path.name}: {error}")

    return successes, failures


def main() -> int:
    args = parse_args()
    successes, failures = render_all(
        input_dir=args.input_dir.resolve(),
        output_dir=args.output_dir.resolve(),
        draw_box=args.draw_box,
    )

    print("\nSummary")
    print(f"  success: {len(successes)}")
    print(f"  failed:  {len(failures)}")

    if failures:
        print("Failures:")
        for name, message in failures[:20]:
            print(f"  {name}: {message}")
        if len(failures) > 20:
            print(f"  ... and {len(failures) - 20} more")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
