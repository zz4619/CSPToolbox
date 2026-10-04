"""Element-symbol normalisation and ordering."""

from __future__ import annotations

import re
from typing import Iterable

from ase.data import atomic_numbers


def guess_element_from_label(label: str) -> str:
    """Guess the element from the leading letters of an atom label."""

    match = re.match(r"([A-Za-z]+)", label)
    if not match:
        raise ValueError(f"Could not infer element from label: {label}")
    return normalize_element_symbol(match.group(1))


def normalize_element_symbol(symbol: str) -> str:
    """Return the canonical element symbol (D and T map to H)."""

    raw = re.sub(r"[^A-Za-z]", "", str(symbol).strip())
    if not raw:
        raise ValueError(f"Could not normalize element symbol: {symbol!r}")
    if raw.upper() in {"D", "T"}:
        return "H"

    candidates = [raw]
    title_case = raw[0].upper() + raw[1:].lower()
    upper_case = raw.upper()
    lower_case = raw.lower()
    first_letter = raw[0].upper()
    for candidate in (title_case, upper_case, lower_case, first_letter):
        if candidate not in candidates:
            candidates.append(candidate)

    for candidate in candidates:
        if candidate in atomic_numbers:
            return candidate

    raise ValueError(f"Unsupported element symbol: {symbol!r}")


def element_order(elements: Iterable[str]) -> list[str]:
    """Return the elements in order of first appearance."""

    ordered: list[str] = []
    for element in elements:
        if element not in ordered:
            ordered.append(element)
    return ordered
