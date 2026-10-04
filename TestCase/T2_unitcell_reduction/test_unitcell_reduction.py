"""T2: full cell -> asymmetric unit -> full cell must reproduce the crystal.

Structures that currently fail are listed in KNOWN_FAILURES with the issue that
explains them (docs/KNOWN_ISSUES.md). They are strict expected failures: when a
fix makes one pass, pytest reports XPASS as a failure, and the entry must then
be removed. The full set is marked ``slow``.
"""

from __future__ import annotations

import warnings

import pytest

from Source.crystal_structure import CrystalStructure
from testsupport.crystal_compare import max_same_element_site_displacement
from testsupport.paths import T2_ASYMMETRIC_DIR, T2_FULL_CELL_DIR, paired_refcodes

SITE_TOLERANCE_ANGSTROM = 1e-6
SMOKE_REFCODES = ("ABALAS", "AQUSEM", "CYCYPR", "DMSULO04", "TETDAM03", "TRIZIN01", "UMIQEO")
ALL_REFCODES = paired_refcodes(T2_ASYMMETRIC_DIR, T2_FULL_CELL_DIR)

_KI1 = "KI-1: SHELX LATT taken from the standard symbol, not the cell's centring"
_KI2 = "KI-2: special-position copies are not merged on re-expansion"
KNOWN_FAILURES = {
    "AQUSEM": _KI1,
    "AWACAE": _KI1,
    "CLPSUL02": _KI1,
    "CRYSEN": _KI1,
    "CUBANE": _KI1,
    "DERBEH": _KI1,
    "LELVIK": _KI1,
    "MEYCEY": _KI1,
    "SALMID04": _KI1,
    "CEKGUU01": _KI2,
    "TETDAM03": _KI2,
}


def _params(refcodes):
    return [
        pytest.param(
            refcode,
            marks=pytest.mark.xfail(reason=KNOWN_FAILURES[refcode], strict=True)
            if refcode in KNOWN_FAILURES
            else (),
        )
        for refcode in refcodes
    ]


def _check_round_trip(refcode: str) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = CrystalStructure.from_file(T2_FULL_CELL_DIR / f"{refcode}.cif")
        reduced = full.reduce_to_asymmetric_unit()
        reexpanded = reduced.expand_to_explicit_unit_cell()
    assert 1 <= len(reduced.atoms) <= len(full.atoms)
    assert reduced.hall_number is not None
    displacement = max_same_element_site_displacement(reexpanded, full)
    assert displacement <= SITE_TOLERANCE_ANGSTROM, f"{refcode}: {displacement:.3e} Å"


@pytest.mark.parametrize("refcode", _params(SMOKE_REFCODES))
def test_reduction_round_trip_smoke(refcode: str) -> None:
    _check_round_trip(refcode)


@pytest.mark.slow
@pytest.mark.parametrize("refcode", _params(ALL_REFCODES))
def test_reduction_round_trip_all(refcode: str) -> None:
    _check_round_trip(refcode)
