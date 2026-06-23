from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import CrystalStructure
from Source.gaussian_input import GaussianInputBuilder, GaussianSettings


DEFAULT_INPUT_ROOT = PROJECT_ROOT / "TestCase"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "TestCase" / "Gaussian"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse all CIF files below an input root, deduplicate unit-cell molecules, "
            "and write Gaussian .com files plus run scripts for each unique molecule."
        )
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help=f"Root directory to scan for CIF files. Default: {DEFAULT_INPUT_ROOT}",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help=f"Root directory for generated Gaussian job folders. Default: {DEFAULT_OUTPUT_ROOT}",
    )
    return parser.parse_args()


def render_all(
    input_root: Path,
    output_root: Path,
    settings: GaussianSettings | None = None,
) -> tuple[list[tuple[Path, int]], list[tuple[Path, str]]]:
    cif_files = sorted(input_root.rglob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No CIF files found under {input_root}")

    builder = GaussianInputBuilder()
    successes: list[tuple[Path, int]] = []
    failures: list[tuple[Path, str]] = []

    for cif_path in cif_files:
        relative_path = cif_path.relative_to(input_root)
        structure_output_dir = output_root / relative_path.parent / cif_path.stem
        try:
            structure = CrystalStructure.from_file(cif_path, fmt="cif")
            artifacts = builder.write_unique_jobs(
                structure=structure,
                output_dir=structure_output_dir,
                settings=settings,
            )
            successes.append((relative_path, len(artifacts)))
            print(
                f"OK   {relative_path} -> {len(artifacts)} unique molecule job(s) in "
                f"{structure_output_dir}"
            )
        except Exception as error:
            failures.append((relative_path, str(error)))
            print(f"FAIL {relative_path}: {error}")

    return successes, failures


def main() -> int:
    args = parse_args()
    successes, failures = render_all(
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
