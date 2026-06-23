"""Safe file manifest helpers for VASP result collection."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import tarfile
from typing import Iterable


DEFAULT_VASP_FILENAMES = ("INCAR", "POSCAR", "CONTCAR", "vasprun.xml", "vasp.out")


@dataclass(frozen=True)
class VaspManifestEntry:
    root: Path
    path: Path
    relative_path: Path
    system_name: str
    calculation_chain: tuple[str, ...]
    filename: str
    size_bytes: int


@dataclass(frozen=True)
class VaspManifestSummary:
    root: Path
    file_count: int
    system_count: int
    calculation_dir_count: int
    counts_by_filename: dict[str, int]
    systems_missing_expected_files: dict[str, tuple[str, ...]]


def collect_vasp_manifest(
    root: str | Path,
    *,
    filenames: Iterable[str] = DEFAULT_VASP_FILENAMES,
) -> list[VaspManifestEntry]:
    """Collect selected VASP files under a calculation root.

    The root is expected to contain one directory per system, optionally with a
    nested non-branching calculation chain such as ``SYSTEM/1/2``. The function
    does not write, copy, or delete anything.
    """

    root_path = Path(root).resolve()
    if not root_path.is_dir():
        raise FileNotFoundError(f"Manifest root is not a directory: {root_path}")

    wanted = set(filenames)
    entries: list[VaspManifestEntry] = []
    for path in sorted(p for p in root_path.rglob("*") if p.is_file() and p.name in wanted):
        relative = path.relative_to(root_path)
        parts = relative.parts
        if not parts:
            continue
        system_name = parts[0]
        chain = tuple(parts[1:-1])
        entries.append(
            VaspManifestEntry(
                root=root_path,
                path=path,
                relative_path=relative,
                system_name=system_name,
                calculation_chain=chain,
                filename=path.name,
                size_bytes=path.stat().st_size,
            )
        )
    return entries


def summarize_manifest(
    entries: list[VaspManifestEntry],
    *,
    expected_filenames: Iterable[str] = DEFAULT_VASP_FILENAMES,
) -> VaspManifestSummary:
    if not entries:
        return VaspManifestSummary(
            root=Path("."),
            file_count=0,
            system_count=0,
            calculation_dir_count=0,
            counts_by_filename={name: 0 for name in expected_filenames},
            systems_missing_expected_files={},
        )

    root = entries[0].root
    expected = tuple(expected_filenames)
    counts = {name: 0 for name in expected}
    systems: set[str] = set()
    calculation_dirs: set[tuple[str, tuple[str, ...]]] = set()
    present_by_system: dict[str, set[str]] = {}

    for entry in entries:
        systems.add(entry.system_name)
        calculation_dirs.add((entry.system_name, entry.calculation_chain))
        present_by_system.setdefault(entry.system_name, set()).add(entry.filename)
        counts[entry.filename] = counts.get(entry.filename, 0) + 1

    missing = {
        system: tuple(name for name in expected if name not in present)
        for system, present in sorted(present_by_system.items())
        if any(name not in present for name in expected)
    }

    return VaspManifestSummary(
        root=root,
        file_count=len(entries),
        system_count=len(systems),
        calculation_dir_count=len(calculation_dirs),
        counts_by_filename=counts,
        systems_missing_expected_files=missing,
    )


def write_manifest_csv(entries: list[VaspManifestEntry], output_csv: str | Path) -> Path:
    """Write a manifest CSV. This is the only required write operation."""

    path = Path(output_csv)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "system_name",
        "calculation_chain",
        "filename",
        "relative_path",
        "size_bytes",
        "absolute_path",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for entry in entries:
            writer.writerow(
                {
                    "system_name": entry.system_name,
                    "calculation_chain": "/".join(entry.calculation_chain) or ".",
                    "filename": entry.filename,
                    "relative_path": str(entry.relative_path),
                    "size_bytes": entry.size_bytes,
                    "absolute_path": str(entry.path),
                }
            )
    return path


def create_tar_from_manifest(
    entries: list[VaspManifestEntry],
    archive_path: str | Path,
    *,
    apply: bool = False,
) -> Path:
    """Create a tar.gz archive from manifest entries.

    By default this is a dry run and only returns the intended archive path.
    Pass ``apply=True`` to write the archive.
    """

    archive = Path(archive_path)
    if not apply:
        return archive
    archive.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "w:gz") as handle:
        for entry in entries:
            handle.add(entry.path, arcname=str(entry.relative_path))
    return archive


def render_summary(summary: VaspManifestSummary) -> str:
    lines = [
        f"root={summary.root}",
        f"files={summary.file_count}",
        f"systems={summary.system_count}",
        f"calculation_dirs={summary.calculation_dir_count}",
    ]
    for filename, count in sorted(summary.counts_by_filename.items()):
        lines.append(f"{filename}={count}")
    if summary.systems_missing_expected_files:
        lines.append(f"systems_with_missing_expected_files={len(summary.systems_missing_expected_files)}")
        for system, missing in list(summary.systems_missing_expected_files.items())[:20]:
            lines.append(f"missing {system}: {','.join(missing)}")
        extra = len(summary.systems_missing_expected_files) - 20
        if extra > 0:
            lines.append(f"... {extra} more systems with missing files")
    else:
        lines.append("all_systems_have_expected_files=true")
    return "\n".join(lines)

