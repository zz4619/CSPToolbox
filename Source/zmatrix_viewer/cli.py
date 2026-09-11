"""Command-line interface for portable structure viewers."""

from __future__ import annotations

import argparse
from pathlib import Path

from .preparation import build_structure_view
from .html_export import write_viewer_html, write_viewer_fragment


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a portable interactive molecular or crystal HTML viewer."
    )
    parser.add_argument("input", type=Path, help="XYZ, PDB, CIF, RES, MOL/SDF or numeric # ZMAT v1 file.")
    parser.add_argument("--reference", type=Path, help="Reference conformation, shown in translucent grey.")
    parser.add_argument("--zmatrix", type=Path, help="Known numeric or CP2 label-only topology to map onto the input.")
    parser.add_argument("--independent-dofs", nargs="+", help="Selected bndN/angN/dihN names, separated by spaces or commas.")
    parser.add_argument("--dofs-file", type=Path, help="JSON definitions or CP2 input.in with DoF bounds.")
    parser.add_argument("--atom-map", type=Path, help="JSON: reference/topology label -> input label or 1-based index.")
    parser.add_argument("--crystal", action="store_true", help="Show asymmetric unit with a full-cell toggle.")
    parser.add_argument("--frame", type=int, default=0, help="Zero-based input frame/model index.")
    parser.add_argument("--reference-frame", type=int, default=0, help="Zero-based reference frame/model index.")
    parser.add_argument("--molecule-index", type=int, help="1-based connected molecule from a periodic input.")
    parser.add_argument("--reference-molecule-index", type=int, help="1-based connected molecule from a periodic reference.")
    parser.add_argument("--fragment", action="store_true", help="Write a self-contained embeddable fragment instead of a full page.")
    parser.add_argument("--title", help="Viewer title.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output HTML path. Defaults to '<input>_viewer.html'.",
    )
    parser.add_argument(
        "--no-infer-bonds",
        action="store_true",
        help="Only draw explicit connectivity and Z-matrix construction bonds.",
    )
    parser.add_argument(
        "--covalent-scale",
        type=float,
        default=1.25,
        help="Covalent-radius multiplier used when inferring display bonds.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    output = args.output or args.input.with_name(f"{args.input.stem}_viewer.html")
    dofs = [name for group in args.independent_dofs for name in group.split(',')] if args.independent_dofs else None
    try:
        scene = build_structure_view(
            args.input, reference=args.reference, zmatrix=args.zmatrix, independent_dofs=dofs,
            dofs_file=args.dofs_file, atom_map=args.atom_map, crystal=args.crystal,
            frame=args.frame, reference_frame=args.reference_frame, molecule_index=args.molecule_index,
            reference_molecule_index=args.reference_molecule_index, title=args.title,
            infer_bonds=not args.no_infer_bonds, covalent_scale=args.covalent_scale,
        )
        writer = write_viewer_fragment if args.fragment else write_viewer_html
        html_path = writer(scene, output)
    except (OSError, ValueError) as exc:
        import sys
        print(f"csp-view: {exc}", file=sys.stderr)
        return 2
    print(f"viewer_html={html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
