"""Generate CSORM input bundles from crystal structures."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import shutil

from ase.data import atomic_masses, atomic_numbers

from .crystal_structure import AtomRecord, CrystalStructure
from .mie_typing import validate_inter_atom_types


TEMPLATE_DIR = Path(__file__).resolve().parents[1] / "Template" / "CSORM_input"
DEFAULT_INFO_TEMPLATE = TEMPLATE_DIR / "SYSTEM_NAME.info"
DEFAULT_RES_TEMPLATE = TEMPLATE_DIR / "SYSTEM_NAME.res"
DEFAULT_CSO_RM_TEMPLATE = TEMPLATE_DIR / "CSO_RM.input"
DEFAULT_POTENTIAL_FILE = TEMPLATE_DIR / "bit_mp_mie_pcm11.inter"
DEFAULT_RUN_TEMPLATE = TEMPLATE_DIR / "runSingleCSORM.csh"


@dataclass(frozen=True)
class CSORMSettings:
    template_dir: Path = TEMPLATE_DIR
    info_template_path: Path = DEFAULT_INFO_TEMPLATE
    res_template_path: Path = DEFAULT_RES_TEMPLATE
    cso_rm_template_path: Path = DEFAULT_CSO_RM_TEMPLATE
    potential_file_path: Path = DEFAULT_POTENTIAL_FILE
    run_template_path: Path = DEFAULT_RUN_TEMPLATE
    pressure: float = 0.0
    exp_ulatt: float = 0.0
    scal_disp: float = 1.2
    scal_elec: float = 1.3
    charge: int = 0
    multiplicity: int = 1
    validate_atom_types: bool = True


@dataclass(frozen=True)
class CSORMJobArtifacts:
    system_name: str
    output_dir: Path
    structure: CrystalStructure
    molecule_atom_counts: list[int]
    res_path: Path
    info_path: Path
    cso_rm_input_path: Path
    potential_path: Path
    run_script_path: Path


class CSORMInputBuilder:
    """Build the static input files required for a CSORM lattice minimisation job."""

    def __init__(self, settings: CSORMSettings | None = None) -> None:
        self.settings = settings or CSORMSettings()

    def write_job_from_crystal(
        self,
        structure: CrystalStructure,
        output_dir: str | Path,
        *,
        system_name: str | None = None,
    ) -> CSORMJobArtifacts:
        _require_non_explict_unit_cell(
            structure,
            context="CSORM input generation requires a non-explict asymmetric-unit structure.",
        )
        job_name = system_name or structure.name
        destination = Path(output_dir)
        destination.mkdir(parents=True, exist_ok=True)

        ordered_structure, molecule_atom_counts = self.reorder_structure_for_csorm(structure, name=job_name)
        if self.settings.validate_atom_types:
            validate_inter_atom_types(ordered_structure, self.settings.potential_file_path)

        res_path = destination / f"{job_name}.res"
        info_path = destination / f"{job_name}.info"
        cso_rm_input_path = destination / self.settings.cso_rm_template_path.name
        potential_path = destination / self.settings.potential_file_path.name
        run_script_path = destination / self.settings.run_template_path.name

        ordered_structure.to_file(res_path, fmt="res", rounding=True)
        info_path.write_text(
            self.render_info_text(molecule_atom_counts),
            encoding="utf-8",
        )
        cso_rm_input_path.write_text(
            self.render_cso_rm_input_text(job_name),
            encoding="utf-8",
        )
        shutil.copyfile(self.settings.potential_file_path, potential_path)
        run_script_path.write_text(
            self.render_run_script_text(job_name),
            encoding="utf-8",
        )
        run_script_path.chmod(0o755)

        return CSORMJobArtifacts(
            system_name=job_name,
            output_dir=destination,
            structure=ordered_structure,
            molecule_atom_counts=molecule_atom_counts,
            res_path=res_path,
            info_path=info_path,
            cso_rm_input_path=cso_rm_input_path,
            potential_path=potential_path,
            run_script_path=run_script_path,
        )

    def reorder_structure_for_csorm(
        self,
        structure: CrystalStructure,
        *,
        name: str | None = None,
    ) -> tuple[CrystalStructure, list[int]]:
        molecules = structure.detect_molecules()
        ordered_molecules = [
            molecule
            for _, molecule in sorted(
                enumerate(molecules),
                key=lambda item: (-_molecular_weight(item[1]), item[0]),
            )
        ]
        molecule_atom_counts = [len(molecule) for molecule in ordered_molecules]

        reordered_atoms: list[AtomRecord] = []
        for molecule in ordered_molecules:
            reordered_atoms.extend(molecule)

        return (
            CrystalStructure(
                atoms=reordered_atoms,
                cell_parameters=structure.cell_parameters,
                lattice_matrix=structure.lattice_matrix,
                space_group=structure.space_group,
                hall_number=structure.hall_number,
                name=name or structure.name,
                explict_unit_cell=structure.explict_unit_cell,
                shelx_latt_value=structure.shelx_latt_value,
                symmetry_operations=structure.symmetry_operations,
            ),
            molecule_atom_counts,
        )

    def render_info_text(self, molecule_atom_counts: list[int]) -> str:
        if not molecule_atom_counts:
            raise ValueError("At least one detected molecule is required for CSORM input.")

        n_mols = len(molecule_atom_counts)
        return (
            f"Num_Mols_Asm {n_mols}\n"
            f"Num_Atoms {' '.join(str(value) for value in molecule_atom_counts)}\n"
            f"exp_Ulatt {self.settings.exp_ulatt}\n"
            f"Pressure {self.settings.pressure}\n"
            f"scal_disp {' '.join(str(self.settings.scal_disp) for _ in molecule_atom_counts)}\n"
            f"scal_elec {' '.join(str(self.settings.scal_elec) for _ in molecule_atom_counts)}\n"
            f"Charges {' '.join(str(self.settings.charge) for _ in molecule_atom_counts)}\n"
            f"Multiplicity {' '.join(str(self.settings.multiplicity) for _ in molecule_atom_counts)}\n"
        )

    def render_cso_rm_input_text(self, system_name: str) -> str:
        text = self.settings.cso_rm_template_path.read_text(encoding="utf-8")
        text = text.replace("DUMMEY_SYSTEM_NAME", system_name)
        text = text.replace("SYSTEM_NAME", system_name)
        text = re.sub(
            r"(^\s*RD filename:\s+).*$",
            rf"\1{self.settings.potential_file_path.name}",
            text,
            flags=re.MULTILINE,
        )
        return text

    def render_run_script_text(self, system_name: str) -> str:
        text = self.settings.run_template_path.read_text(encoding="utf-8")
        text = text.replace("DUMMY_SYSTEM_NAME", system_name)
        text = text.replace("SYSTEM_NAME", system_name)
        return text


def _molecular_weight(molecule: list[AtomRecord]) -> float:
    return float(
        sum(
            atomic_masses[atomic_numbers[atom.element]]
            for atom in molecule
        )
    )


def _require_non_explict_unit_cell(structure: CrystalStructure, *, context: str) -> None:
    if structure.explict_unit_cell:
        raise ValueError(context)
