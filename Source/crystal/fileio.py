"""Helpers shared by the CIF, PDB and RES readers and writers."""

from __future__ import annotations

from pathlib import Path


def normalize_format(path: Path, fmt: str | None) -> str:
    """Return 'cif', 'pdb' or 'res' from ``fmt`` or the file suffix."""

    if fmt is not None:
        return fmt.lower()
    suffix = path.suffix.lower()
    if suffix == ".cif":
        return "cif"
    if suffix == ".pdb":
        return "pdb"
    if suffix == ".res":
        return "res"
    raise ValueError(f"Could not infer format from path: {path}")


def strip_quotes(value: str) -> str:
    """Strip whitespace and one layer of single or double quotes."""

    return value.strip().strip("'").strip('"')


def parse_bool_tag(value: str | None, *, default: bool) -> bool:
    """Parse a true/false tag written by CSPToolbox; ``default`` when absent or unknown."""

    if value is None:
        return default

    normalized = strip_quotes(value).strip().lower()
    if normalized in {"true", "t", "yes", "y", "1"}:
        return True
    if normalized in {"false", "f", "no", "n", "0"}:
        return False
    return default
