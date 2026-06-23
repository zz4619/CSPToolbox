"""Generate CSOFM input bundles from crystal structures and Gaussian logs."""

from __future__ import annotations

from collections import Counter
import csv
from dataclasses import dataclass
from pathlib import Path
import re
import shutil

from ase.data import atomic_masses, atomic_numbers

from .crystal_structure import AtomRecord, CrystalStructure, MoleculeGroup
from .vasp_results import read_contcar_as_crystal


TEMPLATE_DIR = Path(__file__).resolve().parents[1] / "Template" / "CSOFM_input"
DEFAULT_CSO_FM_TEMPLATE = TEMPLATE_DIR / "CSO_FM.input"
DEFAULT_POTENTIAL_FILE = TEMPLATE_DIR / "bit_mp_mie_pcm11.inter"
DEFAULT_ZMATRIX_TEMPLATE = TEMPLATE_DIR / "Zmatrix"
DEFAULT_RUN_CSOFM_TEMPLATE = TEMPLATE_DIR / "runSingleCSOFM.csh"
DEFAULT_RUN_INTEGRITY_TEMPLATE = TEMPLATE_DIR / "runSingleIntegrity.csh"
DEFAULT_ENV_CSV = TEMPLATE_DIR / "csofm_env.csv"
DEFAULT_BATCH_LAM_INTEGRITY_SCRIPT = "batchRunLAMIntegrity.sh"
DEFAULT_BATCH_CSOFM_SCRIPT = "batchRunCSOFM.sh"

DEFAULT_WATER_NAME = "WATER"
DEFAULT_WATER_GAS_PHASE_ENERGY = -76.3688745965


@dataclass(frozen=True)
class CSOFMSettings:
    template_dir: Path = TEMPLATE_DIR
    cso_fm_template_path: Path = DEFAULT_CSO_FM_TEMPLATE
    potential_file_path: Path = DEFAULT_POTENTIAL_FILE
    zmatrix_template_path: Path = DEFAULT_ZMATRIX_TEMPLATE
    run_csofm_template_path: Path = DEFAULT_RUN_CSOFM_TEMPLATE
    run_integrity_template_path: Path = DEFAULT_RUN_INTEGRITY_TEMPLATE
    env_csv_path: Path = DEFAULT_ENV_CSV
    water_name: str = DEFAULT_WATER_NAME
    water_gas_phase_energy_hartree: float = DEFAULT_WATER_GAS_PHASE_ENERGY
    api_charge: int = 0
    api_multiplicity: int = 1
    water_charge: int = 0
    water_multiplicity: int = 1
    include_zmatrix: bool = True
    include_run_scripts: bool = True


@dataclass(frozen=True)
class CSOFMMoleculeDefinition:
    index: int
    molecule_name: str
    atom_count: int
    charge: int
    multiplicity: int
    gas_phase_energy_hartree: float
    global_intra: str
    global_elec: str
    is_water: bool
    signature: str


@dataclass(frozen=True)
class CSOFMJobArtifacts:
    system_name: str
    output_dir: Path
    structure: CrystalStructure
    res_path: Path
    mol_input_path: Path
    cso_fm_input_path: Path
    potential_path: Path
    zmatrix_path: Path | None
    run_csofm_path: Path | None
    run_integrity_path: Path | None
    gaussian_log_path: Path
    gaussian_final_energy_hartree: float
    molecule_definitions: list[CSOFMMoleculeDefinition]


class CSOFMInputBuilder:
    """Build the static input files required for a CSOFM local minimisation job."""

    def __init__(self, settings: CSOFMSettings | None = None) -> None:
        self.settings = settings or CSOFMSettings()

    def read_gaussian_final_energy(self, log_path: str | Path) -> float:
        return read_gaussian_final_energy(log_path)

    def default_hpc_workdir(self) -> str:
        return _read_hpc_workdir(self.settings.env_csv_path)

    def default_csofm_walltime(self) -> str:
        return _read_csofm_walltime(self.settings.env_csv_path)

    def default_lam_integrity_walltime(self) -> str:
        return _read_lam_integrity_walltime(self.settings.env_csv_path)

    def read_contcar_as_crystal(
        self,
        contcar_path: str | Path,
        *,
        name: str | None = None,
        space_group: str = "P 1",
    ) -> CrystalStructure:
        return read_contcar_as_crystal(contcar_path, name=name, space_group=space_group)

    def write_job_from_crystal(
        self,
        structure: CrystalStructure,
        gaussian_log_path: str | Path,
        output_dir: str | Path,
        *,
        system_name: str | None = None,
        hpc_workdir: str | Path | None = None,
        csofm_walltime: str | None = None,
        lam_integrity_walltime: str | None = None,
        global_job_dir: str | Path | None = None,
    ) -> CSOFMJobArtifacts:
        _require_non_explict_unit_cell(
            structure,
            context="CSOFM input generation requires a non-explict asymmetric-unit structure.",
        )
        job_name = system_name or structure.name
        destination = Path(output_dir)
        destination.mkdir(parents=True, exist_ok=True)

        ordered_structure = self.reorder_structure_for_csofm(structure)
        gaussian_log = Path(gaussian_log_path)
        gaussian_energy = self.read_gaussian_final_energy(gaussian_log)
        resolved_hpc_workdir = str(Path(hpc_workdir)) if hpc_workdir is not None else self.default_hpc_workdir()
        resolved_csofm_walltime = csofm_walltime or self.default_csofm_walltime()
        resolved_lam_integrity_walltime = (
            lam_integrity_walltime or self.default_lam_integrity_walltime()
        )
        resolved_global_dir = (
            str(Path(global_job_dir))
            if global_job_dir is not None
            else str(Path(resolved_hpc_workdir) / job_name)
        )
        use_dof = (destination / "dof.txt").is_file()
        molecule_definitions = self.build_molecule_definitions(
            structure=ordered_structure,
            system_name=job_name,
            gaussian_final_energy_hartree=gaussian_energy,
            global_job_dir=resolved_global_dir,
        )

        serializable_structure = CrystalStructure(
            atoms=ordered_structure.atoms,
            cell_parameters=ordered_structure.cell_parameters,
            lattice_matrix=ordered_structure.lattice_matrix,
            space_group=ordered_structure.space_group,
            name=job_name,
            explict_unit_cell=ordered_structure.explict_unit_cell,
        )

        res_path = destination / f"{job_name}.res"
        mol_input_path = destination / "mol.input"
        cso_fm_input_path = destination / "CSO_FM.input"
        potential_path = destination / self.settings.potential_file_path.name
        zmatrix_path = destination / self.settings.zmatrix_template_path.name if self.settings.include_zmatrix else None
        run_csofm_path = destination / self.settings.run_csofm_template_path.name if self.settings.include_run_scripts else None
        run_integrity_path = destination / self.settings.run_integrity_template_path.name if self.settings.include_run_scripts else None

        serializable_structure.to_file(res_path, fmt="res", rounding=True)
        mol_input_path.write_text(
            self.render_mol_input_text(molecule_definitions),
            encoding="utf-8",
        )
        shutil.copyfile(self.settings.cso_fm_template_path, cso_fm_input_path)
        shutil.copyfile(self.settings.potential_file_path, potential_path)

        if zmatrix_path is not None:
            zmatrix_path.write_text(
                self.render_csofm_zmatrix_text(serializable_structure),
                encoding="utf-8",
            )

        if run_csofm_path is not None and run_integrity_path is not None:
            run_csofm_path.write_text(
                self.render_run_script_text(
                    self.settings.run_csofm_template_path,
                    job_name,
                    resolved_hpc_workdir,
                    resolved_csofm_walltime,
                    use_dof,
                ),
                encoding="utf-8",
            )
            run_integrity_path.write_text(
                self.render_run_script_text(
                    self.settings.run_integrity_template_path,
                    job_name,
                    resolved_hpc_workdir,
                    resolved_lam_integrity_walltime,
                    use_dof,
                ),
                encoding="utf-8",
            )
            run_csofm_path.chmod(0o755)
            run_integrity_path.chmod(0o755)

        return CSOFMJobArtifacts(
            system_name=job_name,
            output_dir=destination,
            structure=serializable_structure,
            res_path=res_path,
            mol_input_path=mol_input_path,
            cso_fm_input_path=cso_fm_input_path,
            potential_path=potential_path,
            zmatrix_path=zmatrix_path,
            run_csofm_path=run_csofm_path,
            run_integrity_path=run_integrity_path,
            gaussian_log_path=gaussian_log,
            gaussian_final_energy_hartree=gaussian_energy,
            molecule_definitions=molecule_definitions,
        )

    def build_molecule_definitions(
        self,
        *,
        structure: CrystalStructure,
        system_name: str,
        gaussian_final_energy_hartree: float,
        global_job_dir: str,
    ) -> list[CSOFMMoleculeDefinition]:
        groups = structure.deduplicate_molecules()
        if not groups:
            raise ValueError(f"No molecules detected in structure {system_name}")

        definitions: list[CSOFMMoleculeDefinition] = []
        for index, group in enumerate(groups, start=1):
            molecule = group.representative_molecule
            is_water = _is_water_molecule(molecule)
            if is_water:
                molecule_name = self.settings.water_name
                gas_phase_energy = self.settings.water_gas_phase_energy_hartree
                charge = self.settings.water_charge
                multiplicity = self.settings.water_multiplicity
            else:
                molecule_name = system_name
                gas_phase_energy = gaussian_final_energy_hartree
                charge = self.settings.api_charge
                multiplicity = self.settings.api_multiplicity

            definitions.append(
                CSOFMMoleculeDefinition(
                    index=index,
                    molecule_name=molecule_name,
                    atom_count=len(molecule),
                    charge=charge,
                    multiplicity=multiplicity,
                    gas_phase_energy_hartree=gas_phase_energy,
                    global_intra=str(Path(global_job_dir) / f"intra_lam_global_{index}"),
                    global_elec=str(Path(global_job_dir) / f"elec_lam_global_{index}"),
                    is_water=is_water,
                    signature=group.signature,
                )
            )

        return definitions

    def render_mol_input_text(
        self,
        molecule_definitions: list[CSOFMMoleculeDefinition],
    ) -> str:
        if not molecule_definitions:
            raise ValueError("At least one molecule definition is required.")

        n_mol_asm = len(molecule_definitions)
        n_mol_unique = len(molecule_definitions)

        lines = [
            "############################################################################################",
            "#                                Molecule Definition File                                  #",
            "############################################################################################",
            "  ",
            f"Nmol_asm  {n_mol_asm}       #number of distinct molecules present in the .res file    ",
            f"Nmol_unique  {n_mol_unique}    #number of unique molecules present in the .res file",
            "",
            "############################################################################################",
            "#      Molecules entered in the order they are processed in the specified .res file        #",
            "#               Dashed lines (------) are used as dividers between entries                 #",
            "############################################################################################",
            "",
        ]

        for definition in molecule_definitions:
            lines.extend(
                [
                    "--------------------------------------------------------------------------------------------",
                    "Number of molecules  1",
                    f"Molecule Name        {definition.molecule_name}",
                    f"nAtoms               {definition.atom_count}",
                    f"Charge               {definition.charge}",
                    f"Multiplicity         {definition.multiplicity}",
                    "scal_disp            1.0",
                    "scal_elec            1.0",
                    f"Gas phase energy     {definition.gas_phase_energy_hartree:.10f}  Hartree  ",
                    "",
                    "# Database file names and paths ",
                    f"Local intra     './intra_lam_tmp_{definition.index}'                         #Should be referenced locally",
                    f"Global intra    '{definition.global_intra}'",
                    f"Local elec      './elec_lam_tmp_{definition.index}'                         #Should be referenced locally",
                    f"Global elec     '{definition.global_elec}'    #Should be referenced globally",
                ]
            )

        lines.append("--------------------------------------------------------------------------------------------")
        return "\n".join(lines) + "\n"

    def render_csofm_zmatrix_text(self, structure: CrystalStructure) -> str:
        zmatrices = structure.generate_zmatrices()
        lines: list[str] = []
        for index, zmatrix in enumerate(zmatrices, start=1):
            lines.append(f"Z-matrix for molecule {index}")
            for entry_index, entry in enumerate(zmatrix.entries, start=1):
                row = [entry.label]
                if entry.bond_to is not None:
                    row.append(zmatrix.ordered_atom_labels[entry.bond_to - 1])
                if entry.angle_to is not None:
                    row.append(zmatrix.ordered_atom_labels[entry.angle_to - 1])
                if entry.dihedral_to is not None:
                    row.append(zmatrix.ordered_atom_labels[entry.dihedral_to - 1])
                lines.append("    ".join(row))
        return "\n".join(lines) + "\n"

    def reorder_structure_for_csofm(self, structure: CrystalStructure) -> CrystalStructure:
        molecules = structure.detect_molecules()
        zmatrices = structure.generate_zmatrices()
        if len(molecules) != len(zmatrices):
            raise ValueError("Detected molecule count does not match generated Z-matrix count.")

        ordered_pairs = sorted(
            enumerate(zip(molecules, zmatrices, strict=True)),
            key=lambda item: (-_molecular_weight(item[1][0]), item[0]),
        )

        reordered_atoms: list[AtomRecord] = []
        for _, (molecule, zmatrix) in ordered_pairs:
            atoms_by_label = {atom.label: atom for atom in molecule}
            reordered_atoms.extend(atoms_by_label[label] for label in zmatrix.ordered_atom_labels)

        return CrystalStructure(
            atoms=reordered_atoms,
            cell_parameters=structure.cell_parameters,
            lattice_matrix=structure.lattice_matrix,
            space_group=structure.space_group,
            name=structure.name,
            explict_unit_cell=structure.explict_unit_cell,
        )

    def render_run_script_text(
        self,
        template_path: str | Path,
        system_name: str,
        hpc_workdir: str | Path,
        walltime: str,
        use_dof: bool,
    ) -> str:
        template_text = Path(template_path).read_text(encoding="utf-8")
        script = (
            template_text
            .replace("SYSTEM_NAME", system_name)
            .replace("HPC_WORKDIR", str(Path(hpc_workdir)))
        )
        script = re.sub(
            r"^#PBS -l walltime=.*$",
            f"#PBS -l walltime={walltime}",
            script,
            flags=re.MULTILINE,
        )
        if not use_dof:
            script = script.replace(" DOF >", " >")
            script = script.replace(' DOF >', ' >')
            script = script.replace(' DOF 2>&1', ' 2>&1')
            script = script.replace('" DOF >', '" >')
        return script

    def render_batch_submission_script_text(self, run_script_name: str) -> str:
        return (
            "#!/bin/bash\n"
            "set -eu\n"
            "\n"
            'SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)\n'
            'ROOT_DIR="${1:-$SCRIPT_DIR}"\n'
            "\n"
            'find "$ROOT_DIR" -mindepth 1 -maxdepth 1 -type d -print | sort | while IFS= read -r system_dir; do\n'
            f'    if [ ! -f "$system_dir/{run_script_name}" ]; then\n'
            "        continue\n"
            "    fi\n"
            '    system_name=$(basename "$system_dir")\n'
            f'    printf \'%s\\n\' "Submitting $system_name with qsub {run_script_name}"\n'
            "    (\n"
            '        cd "$system_dir"\n'
            f"        qsub {run_script_name}\n"
            "    )\n"
            "done\n"
        )

    def write_batch_submission_scripts(self, root_dir: str | Path) -> tuple[Path, Path]:
        destination = Path(root_dir)
        destination.mkdir(parents=True, exist_ok=True)

        lam_integrity_path = destination / DEFAULT_BATCH_LAM_INTEGRITY_SCRIPT
        csofm_path = destination / DEFAULT_BATCH_CSOFM_SCRIPT

        lam_integrity_path.write_text(
            self.render_batch_submission_script_text("runSingleIntegrity.csh"),
            encoding="utf-8",
        )
        csofm_path.write_text(
            self.render_batch_submission_script_text("runSingleCSOFM.csh"),
            encoding="utf-8",
        )
        lam_integrity_path.chmod(0o755)
        csofm_path.chmod(0o755)
        return lam_integrity_path, csofm_path


def read_gaussian_final_energy(log_path: str | Path) -> float:
    path = Path(log_path)
    if not path.is_file():
        raise FileNotFoundError(f"Gaussian log not found: {path}")

    energy: float | None = None
    pattern = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+([\-0-9.]+)")
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = pattern.search(line)
        if match:
            energy = float(match.group(1))

    if energy is None:
        raise ValueError(f"No 'SCF Done' energy found in Gaussian log: {path}")
    return energy


def _is_water_molecule(molecule: list[AtomRecord]) -> bool:
    if len(molecule) != 3:
        return False
    composition = Counter(atom.element for atom in molecule)
    return composition == Counter({"H": 2, "O": 1})


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


def _read_hpc_workdir(env_csv_path: str | Path) -> str:
    return _read_env_value(env_csv_path, "hpc_workdir")


def _read_csofm_walltime(env_csv_path: str | Path) -> str:
    return _read_env_value(env_csv_path, "CSOFM_walltime", fallback_tags=("walltime",))


def _read_lam_integrity_walltime(env_csv_path: str | Path) -> str:
    return _read_env_value(env_csv_path, "LAM_integrity_walltime", fallback_tags=("walltime",))


def _read_env_value(
    env_csv_path: str | Path,
    tag: str,
    *,
    fallback_tags: tuple[str, ...] = (),
) -> str:
    path = Path(env_csv_path)
    if not path.is_file():
        raise FileNotFoundError(f"CSOFM env csv not found: {path}")
    requested_tags = (tag, *fallback_tags)
    found_values: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            row_tag = row.get("tag")
            if row_tag in requested_tags:
                value = row.get("path", "").strip()
                if not value:
                    break
                found_values[row_tag] = value
    for requested_tag in requested_tags:
        if requested_tag in found_values:
            return found_values[requested_tag]
    raise ValueError(f"No {tag} entry found in {path}")
