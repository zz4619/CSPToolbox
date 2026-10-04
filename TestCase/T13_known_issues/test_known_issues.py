"""T13: executable specifications of known scientific issues (docs/KNOWN_ISSUES.md).

Each test states the correct behaviour and is a strict expected failure until
the issue is fixed; ``raises=`` pins the expected failure, so an error in
the test itself is not mistaken for the known issue. When a fix lands, the test XPASSes, which pytest reports as
a failure: remove its xfail marker and update docs/KNOWN_ISSUES.md in the same
change. Passing control tests show the checks themselves are sound.
"""

from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
import pytest
import spglib

from Source.crystal.records import AtomRecord
from Source.crystal.symmetry import _space_group_hall_number, expand_shelx_atoms, shelx_symmetry_records
from Source.crystal_structure import CrystalStructure
from Source.csofm_input import CSOFMInputBuilder, read_gaussian_final_energy
from testsupport.crystal_compare import max_same_element_site_displacement
from testsupport.paths import T2_FULL_CELL_DIR

GENERAL_POSITION = np.array([0.1234, 0.2345, 0.3456])
FAKE_GAUSSIAN_LOG = " SCF Done:  E(RPBE1PBE) =  -512.123456789     A.U. after   12 cycles\n"


def _shelx_orbit(symbol: str) -> set[tuple[float, ...]]:
    lattice = np.diag([11.0, 11.0, 13.0])
    latt, operations = shelx_symmetry_records(symbol)
    atoms = expand_shelx_atoms(
        atoms=[AtomRecord("X1", "C", tuple(GENERAL_POSITION @ lattice))],
        cell=lattice,
        latt_value=latt,
        symmetry_operations=operations,
    )
    inverse = np.linalg.inv(lattice)
    return {tuple(np.round((np.asarray(a.coordinates) @ inverse) % 1.0, 6) % 1.0) for a in atoms}


def _database_orbit(symbol: str) -> set[tuple[float, ...]]:
    symmetry = spglib.get_symmetry_from_database(_space_group_hall_number(symbol))
    return {
        tuple(np.round((rotation @ GENERAL_POSITION + translation) % 1.0, 6) % 1.0)
        for rotation, translation in zip(symmetry["rotations"], symmetry["translations"])
    }


@pytest.mark.parametrize("symbol", ["P 21/c", "C 2/c", "P b c a", "P -1"])
def test_control_shelx_records_reproduce_the_orbit(symbol):
    assert _shelx_orbit(symbol) == _database_orbit(symbol)


@pytest.mark.xfail(strict=True, reason="KI-3: LATT>0 written when the inversion centre is off the origin", raises=AssertionError)
@pytest.mark.parametrize("symbol", ["I 41/a", "P 4/n", "P n n n", "F d d d"])
def test_ki3_shelx_records_reproduce_origin_choice_1_orbits(symbol):
    assert _shelx_orbit(symbol) == _database_orbit(symbol)


def _csofm_round_trip(refcode: str, tmp_path: Path) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = CrystalStructure.from_file(T2_FULL_CELL_DIR / f"{refcode}.cif")
        asu = full.reduce_to_asymmetric_unit()
    log = tmp_path / "gas.log"
    log.write_text(FAKE_GAUSSIAN_LOG, encoding="utf-8")
    CSOFMInputBuilder().write_job_from_crystal(asu, log, tmp_path / "csofm", system_name=refcode, hpc_workdir="/hpc")
    written = CrystalStructure.from_file(tmp_path / "csofm" / f"{refcode}.res").expand_to_explicit_unit_cell()
    return max_same_element_site_displacement(written, full)


def test_control_csofm_res_reproduces_a_p21c_crystal(tmp_path):
    assert _csofm_round_trip("ABALAS", tmp_path) < 1e-4


@pytest.mark.xfail(strict=True, reason="KI-4: CSO-FM .res drops LATT/SYMM; P21/n is written with P21/c operations", raises=AssertionError)
def test_ki4_csofm_res_reproduces_a_p21n_crystal(tmp_path):
    assert _csofm_round_trip("ABEROP", tmp_path) < 1e-4


TWO_BLOCK_CIF = """data_first
_cell_length_a 5.0
_cell_length_b 5.0
_cell_length_c 5.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.5 0.5 0.5
data_second
_cell_length_a 20.0
_cell_length_b 20.0
_cell_length_c 20.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
N9 N 0.1 0.1 0.1
"""


@pytest.mark.xfail(strict=True, reason="KI-5: a multi-block CIF mixes the first block's atoms with the last block's cell", raises=AssertionError)
def test_ki5_multi_block_cif_is_not_mixed(tmp_path):
    path = tmp_path / "two_blocks.cif"
    path.write_text(TWO_BLOCK_CIF, encoding="utf-8")
    structure = CrystalStructure.from_file(path)
    assert structure.cell_parameters[0] == 5.0 and structure.atoms[0].label == "C1"


@pytest.mark.xfail(strict=True, reason="KI-6: mol.input writes Nmol_asm = unique species and 'Number of molecules 1'", raises=AssertionError)
def test_ki6_mol_input_counts_every_molecule_for_z_prime_2(tmp_path):
    refcode = "QICVAZ01"  # P21/c, two molecules of one species in the asymmetric unit
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        asu = CrystalStructure.from_file(T2_FULL_CELL_DIR / f"{refcode}.cif").reduce_to_asymmetric_unit()
    assert len(asu.detect_molecules()) == 2
    log = tmp_path / "gas.log"
    log.write_text(FAKE_GAUSSIAN_LOG, encoding="utf-8")
    CSOFMInputBuilder().write_job_from_crystal(asu, log, tmp_path / "csofm", system_name=refcode, hpc_workdir="/hpc")
    fields = dict(
        line.split()[:2] for line in (tmp_path / "csofm" / "mol.input").read_text().splitlines()
        if line.startswith(("Nmol_asm", "Nmol_unique"))
    )
    assert fields == {"Nmol_asm": "2", "Nmol_unique": "1"}


@pytest.mark.xfail(strict=True, reason="KI-8: the gas-phase energy is read without checking Gaussian terminated normally", raises=pytest.fail.Exception)
def test_ki8_failed_gaussian_job_is_rejected(tmp_path):
    log = tmp_path / "failed.log"
    log.write_text(FAKE_GAUSSIAN_LOG + " Error termination via Lnk1e in l9999.exe\n", encoding="utf-8")
    with pytest.raises(ValueError):
        read_gaussian_final_energy(log)


DISORDERED_CIF = """data_disorder
_cell_length_a 10.0
_cell_length_b 10.0
_cell_length_c 10.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
C1 C 0.10 0.10 0.10 1.0
O1A O 0.20 0.10 0.10 0.6
O1B O 0.10 0.20 0.10 0.4
"""


@pytest.mark.xfail(strict=True, reason="KI-9: partial occupancies are kept silently by expand_cif_to_unit_cell", raises=pytest.fail.Exception)
def test_ki9_partial_occupancy_is_reported(tmp_path):
    path = tmp_path / "disorder.cif"
    path.write_text(DISORDERED_CIF, encoding="utf-8")
    with pytest.warns(UserWarning, match="(?i)occupanc"):
        CrystalStructure.expand_cif_to_unit_cell(path)
