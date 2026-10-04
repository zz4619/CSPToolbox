"""T11: approval tests for every generated program input.

For seven CE755 structures the test reduces the experimental full cell and
writes CSO-RM, CSO-FM, Gaussian and VASP inputs plus the structure files and
Mie atom types. Each generated text file must equal its reviewed copy under
``expected/<REFCODE>/``. Large files copied verbatim from ``Template/``
(potentials, POTCARs, CSO_FM.input) are not snapshotted.

The snapshots record current behaviour, including known issues such as KI-6
(Z'>1 ``mol.input`` counts, see QICVAZ01). After an intended change, regenerate
and review the diff (see ``testsupport/snapshots.py``).
"""

from __future__ import annotations

from pathlib import Path
import warnings

import pytest

from Source.crystal_structure import CrystalStructure
from Source.csofm_input import CSOFMInputBuilder
from Source.csorm_input import CSORMInputBuilder
from Source.gaussian_input import GaussianInputBuilder
from Source.mie_typing import detect_fit_atom_types
from Source.vasp_input import VaspInputBuilder
from testsupport.paths import T2_FULL_CELL_DIR
from testsupport.snapshots import assert_text_snapshot

EXPECTED_ROOT = Path(__file__).resolve().parent / "expected"
REFCODES = (
    "BOTYIT",    # P2_1/c, chlorine
    "FIRLAU",    # Pbca, fluorine
    "MALEHY11",  # P-1, N-H and O-H hydrogens
    "PURINE",    # Pna2_1
    "QICVAZ01",  # P2_1/c with Z' = 2
    "SUCANH12",  # P2_12_12_1
    "XAGWIM",    # C2/c, centred cell
)
SKIPPED_COPIES = {"POTCAR", "CSO_FM.input"}
FAKE_GAUSSIAN_LOG = " SCF Done:  E(RPBE1PBE) =  -512.123456789     A.U. after   12 cycles\n"


def _collect(directory: Path, prefix: str) -> dict[str, str]:
    files = {}
    for path in sorted(p for p in directory.rglob("*") if p.is_file()):
        if path.name in SKIPPED_COPIES or path.suffix == ".inter":
            continue
        files[f"{prefix}/{path.relative_to(directory).as_posix()}"] = path.read_text(encoding="utf-8")
    return files


def _attempt(prefix: str, write, directory: Path) -> dict[str, str]:
    """Snapshot the files a writer produces, or the error it raises."""

    try:
        write(directory)
    except Exception as error:  # noqa: BLE001 - a raised error is recorded behaviour
        return {f"{prefix}/ERROR.txt": f"{type(error).__name__}: {error}\n"}
    return _collect(directory, prefix)


def generate_outputs(refcode: str, work: Path) -> dict[str, str]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = CrystalStructure.from_file(T2_FULL_CELL_DIR / f"{refcode}.cif")
        full.explict_unit_cell = True
        asu = full.reduce_to_asymmetric_unit()

    outputs: dict[str, str] = {}
    structure_dir = work / "structure"
    structure_dir.mkdir()
    asu.to_file(structure_dir / "asu.res", rounding=True)
    asu.to_file(structure_dir / "asu.cif")
    asu.to_file(structure_dir / "asu.pdb")
    outputs.update(_collect(structure_dir, "structure"))
    outputs["structure/fit_types.txt"] = "".join(
        f"{t.atom_label} {t.element} {t.fit_type} {t.molecule_index} {t.bonded_parent_label}\n"
        for t in detect_fit_atom_types(asu)
    )

    log = work / "gas.log"
    log.write_text(FAKE_GAUSSIAN_LOG, encoding="utf-8")
    outputs.update(_attempt(
        "csorm", lambda d: CSORMInputBuilder().write_job_from_crystal(asu, d, system_name=refcode), work / "csorm"))
    outputs.update(_attempt(
        "csofm",
        lambda d: CSOFMInputBuilder().write_job_from_crystal(asu, log, d, system_name=refcode, hpc_workdir="/hpc/work"),
        work / "csofm",
    ))
    outputs.update(_attempt("gaussian", lambda d: GaussianInputBuilder().write_unique_jobs(asu, d), work / "gaussian"))
    outputs.update(_attempt(
        "vasp", lambda d: VaspInputBuilder().write_job_from_structure(full, d, system_name=refcode), work / "vasp"))
    return outputs


@pytest.mark.parametrize("refcode", REFCODES)
def test_generated_inputs_match_snapshots(refcode: str, tmp_path: Path) -> None:
    outputs = generate_outputs(refcode, tmp_path)
    expected_dir = EXPECTED_ROOT / refcode
    stale = {
        p.relative_to(expected_dir).as_posix() for p in expected_dir.rglob("*") if p.is_file()
    } - set(outputs) if expected_dir.is_dir() else set()
    assert not stale, f"Snapshots no longer produced (delete them if intended): {sorted(stale)}"
    for name, text in outputs.items():
        assert_text_snapshot(text, expected_dir / name)
