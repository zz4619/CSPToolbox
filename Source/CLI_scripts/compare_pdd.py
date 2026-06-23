from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import CrystalStructure
from Source.pdd_descriptor import calculate_pdd, pdd_distance_breakdown
from Source.vasp_results import read_contcar_as_crystal


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare two explicit unit-cell crystal structures using the PDD descriptor."
    )
    parser.add_argument("structure_a", type=Path, help="Path to the first structure.")
    parser.add_argument("structure_b", type=Path, help="Path to the second structure.")
    parser.add_argument("--format-a", default=None, help="Optional format for structure A.")
    parser.add_argument("--format-b", default=None, help="Optional format for structure B.")
    parser.add_argument("--k", type=int, default=100, help="Number of nearest neighbours in the PDD.")
    parser.add_argument(
        "--metric",
        default="chebyshev",
        help="Distance metric passed to scipy.spatial.distance.cdist. Default: chebyshev",
    )
    parser.add_argument(
        "--show-pdd",
        action="store_true",
        help="Print the descriptor rows for both structures.",
    )
    parser.add_argument(
        "--geometry-only",
        action="store_true",
        help="Ignore element types and compare all PDD rows globally.",
    )
    return parser.parse_args()


def _read_structure(path: Path, fmt: str | None) -> CrystalStructure:
    normalized = fmt.lower() if fmt is not None else None
    if normalized in {"contcar", "poscar"} or (fmt is None and path.name in {"CONTCAR", "POSCAR"}):
        return read_contcar_as_crystal(path, name=path.stem or path.parent.name)
    return CrystalStructure.from_file(path, fmt=fmt)


def _print_descriptor(title: str, descriptor) -> None:
    print(title)
    for element, row in zip(descriptor.center_elements, descriptor.matrix, strict=True):
        values = " ".join(f"{value:.8f}" for value in row)
        print(f"  {element}  {values}")


def main() -> int:
    args = parse_args()
    structure_a = _read_structure(args.structure_a.resolve(), args.format_a)
    structure_b = _read_structure(args.structure_b.resolve(), args.format_b)
    total, by_element = pdd_distance_breakdown(
        structure_a,
        structure_b,
        k=args.k,
        typed=not args.geometry_only,
        metric=args.metric,
    )

    print(f"{structure_a.name}: explict_unit_cell={structure_a.explict_unit_cell}, atoms={len(structure_a.atoms)}")
    print(f"{structure_b.name}: explict_unit_cell={structure_b.explict_unit_cell}, atoms={len(structure_b.atoms)}")
    mode = "geometry-only" if args.geometry_only else "typed"
    print(f"PDD distance (k={args.k}, metric={args.metric}, mode={mode}): {total:.10f}")
    print("Breakdown:")
    for element, value in by_element.items():
        print(f"  {element}: {value:.10f}")

    if args.show_pdd:
        descriptor_a = calculate_pdd(structure_a, k=args.k, typed=not args.geometry_only)
        descriptor_b = calculate_pdd(structure_b, k=args.k, typed=not args.geometry_only)
        _print_descriptor(f"\n{structure_a.name} PDD", descriptor_a)
        _print_descriptor(f"\n{structure_b.name} PDD", descriptor_b)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
