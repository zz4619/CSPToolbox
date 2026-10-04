"""T1: CIF asymmetric unit -> explicit P1 cell must match CCDC's packed cell.

Every CE755 structure is compared site by site with the CCDC Python API
expansion (see README.md). The default run checks a smoke subset; the full set
is marked ``slow``: ``python -m pytest -m slow TestCase/T1_unitcell_expansion``.
"""

from __future__ import annotations

import warnings

import pytest

from Source.crystal_structure import CrystalStructure
from testsupport.crystal_compare import max_same_element_site_displacement
from testsupport.paths import T1_ASYMMETRIC_DIR, T1_FULL_CELL_DIR, paired_refcodes

SITE_TOLERANCE_ANGSTROM = 1e-6
SMOKE_REFCODES = ("ABALAS", "AQUSEM", "CUBANE", "CYCYPR", "TRIZIN01", "UMIQEO")
ALL_REFCODES = paired_refcodes(T1_ASYMMETRIC_DIR, T1_FULL_CELL_DIR)


def _check_expansion(refcode: str) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expanded = CrystalStructure.expand_cif_to_unit_cell(T1_ASYMMETRIC_DIR / f"{refcode}.cif")
        reference = CrystalStructure.from_file(T1_FULL_CELL_DIR / f"{refcode}.cif")
    assert expanded.explict_unit_cell
    displacement = max_same_element_site_displacement(expanded, reference)
    assert displacement <= SITE_TOLERANCE_ANGSTROM, f"{refcode}: {displacement:.3e} Å"


@pytest.mark.parametrize("refcode", SMOKE_REFCODES)
def test_expansion_matches_ccdc_smoke(refcode: str) -> None:
    _check_expansion(refcode)


@pytest.mark.slow
@pytest.mark.parametrize("refcode", ALL_REFCODES)
def test_expansion_matches_ccdc_all(refcode: str) -> None:
    _check_expansion(refcode)
