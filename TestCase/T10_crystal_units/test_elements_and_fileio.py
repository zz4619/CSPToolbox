"""Unit tests for element and file-format helpers in Source.crystal."""

from __future__ import annotations

from pathlib import Path

import pytest

from Source.crystal.elements import element_order, guess_element_from_label, normalize_element_symbol
from Source.crystal.fileio import normalize_format, parse_bool_tag


@pytest.mark.parametrize(
    ("raw", "expected"),
    [("C", "C"), ("cl", "Cl"), ("CL", "Cl"), ("Br1", "Br"), ("D", "H"), ("T", "H"), ("n", "N")],
)
def test_normalize_element_symbol(raw, expected):
    assert normalize_element_symbol(raw) == expected


def test_label_guess_uses_leading_letters():
    assert guess_element_from_label("C12") == "C"
    assert guess_element_from_label("Cl3") == "Cl"


def test_element_order_keeps_first_appearance():
    assert element_order(["C", "H", "C", "N", "H"]) == ["C", "H", "N"]


@pytest.mark.parametrize(("name", "expected"), [("a.CIF", "cif"), ("b.pdb", "pdb"), ("c.res", "res")])
def test_format_from_suffix(name, expected):
    assert normalize_format(Path(name), None) == expected


def test_unknown_suffix_is_rejected():
    with pytest.raises(ValueError):
        normalize_format(Path("structure.xyz"), None)


@pytest.mark.parametrize(("value", "expected"), [("TRUE", True), ("'false'", False), ("1", True), ("maybe", None), (None, None)])
def test_parse_bool_tag(value, expected):
    assert parse_bool_tag(value, default=None) is expected
