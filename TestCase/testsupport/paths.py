"""Locations of the repository and of the shared structure datasets."""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
TESTCASE_ROOT = REPO_ROOT / "TestCase"

# CE755 experimental set: CCDC asymmetric-unit CIFs and CCDC-packed P1 full cells.
T1_ASYMMETRIC_DIR = TESTCASE_ROOT / "T1_unitcell_expansion" / "Experimental"
T1_FULL_CELL_DIR = TESTCASE_ROOT / "T1_unitcell_expansion" / "Experimental_FullUnitCell"
T2_ASYMMETRIC_DIR = TESTCASE_ROOT / "T2_unitcell_reduction" / "Experimental"
T2_FULL_CELL_DIR = TESTCASE_ROOT / "T2_unitcell_reduction" / "Experimental_FullUnitCell"


def paired_refcodes(asymmetric_dir: Path, full_cell_dir: Path) -> list[str]:
    """Return the sorted refcodes present in both folders, failing if unpaired."""

    asymmetric = {path.stem for path in asymmetric_dir.glob("*.cif")}
    full_cell = {path.stem for path in full_cell_dir.glob("*.cif")}
    if asymmetric != full_cell:
        raise AssertionError(
            f"Unpaired CIFs: only asymmetric={sorted(asymmetric - full_cell)[:10]}, "
            f"only full cell={sorted(full_cell - asymmetric)[:10]}"
        )
    return sorted(asymmetric)
