"""Command-line interface for the static Z-matrix viewer."""

from __future__ import annotations

import argparse
from pathlib import Path

from . import build_viewer_document
from .html_export import write_viewer_html


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a standalone interactive HTML viewer for a CSPToolbox Z-matrix."
    )
    parser.add_argument("zmatrix", type=Path, help="Input # ZMAT v1 file.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output HTML path. Defaults to '<input>_viewer.html'.",
    )
    parser.add_argument(
        "--no-infer-bonds",
        action="store_true",
        help="Only draw Z-matrix construction bonds and explicit # bonds comments.",
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
    output = args.output or args.zmatrix.with_name(f"{args.zmatrix.stem}_viewer.html")
    molecule = build_viewer_document(
        args.zmatrix,
        infer_bonds=not args.no_infer_bonds,
        covalent_scale=args.covalent_scale,
    )
    html_path = write_viewer_html(molecule, output)
    print(f"viewer_html={html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

