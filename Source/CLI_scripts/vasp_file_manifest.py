"""Command-line wrapper for VASP file manifest helpers."""

from __future__ import annotations

import argparse
from pathlib import Path

from Source.vasp_file_manifest import (
    DEFAULT_VASP_FILENAMES,
    collect_vasp_manifest,
    create_tar_from_manifest,
    render_summary,
    summarize_manifest,
    write_manifest_csv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a dry-run-safe manifest of selected VASP files."
    )
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument(
        "--filenames",
        nargs="+",
        default=list(DEFAULT_VASP_FILENAMES),
        help="File basenames to collect.",
    )
    parser.add_argument("--output-csv", type=Path)
    parser.add_argument("--archive", type=Path, help="Optional tar.gz archive path.")
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually write the optional archive. CSV output is always written when requested.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    entries = collect_vasp_manifest(args.root, filenames=args.filenames)
    summary = summarize_manifest(entries, expected_filenames=args.filenames)
    print(render_summary(summary))

    if args.output_csv is not None:
        manifest_path = write_manifest_csv(entries, args.output_csv)
        print(f"manifest_csv={manifest_path}")

    if args.archive is not None:
        archive_path = create_tar_from_manifest(entries, args.archive, apply=args.apply)
        action = "archive" if args.apply else "dry_run_archive"
        print(f"{action}={archive_path}")
        if not args.apply:
            print("archive_written=false")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

