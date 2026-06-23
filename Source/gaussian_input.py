"""Gaussian gas-phase input generation for deduplicated crystal molecules."""

from __future__ import annotations

from dataclasses import dataclass
import re
from pathlib import Path

from pymatgen.core import Molecule as PymatgenMolecule

from .crystal_structure import AtomRecord, CrystalStructure, MoleculeGroup, ZMatrixRepresentation


TEMPLATE_DIR = Path(__file__).resolve().parents[1] / "Template" / "Gaussian_input"
DEFAULT_COM_TEMPLATE = TEMPLATE_DIR / "template.com"
DEFAULT_RUN_TEMPLATE = TEMPLATE_DIR / "runGAUSSIAN_cx3.csh"


@dataclass(frozen=True)
class GaussianSettings:
    chk: str = "hess.chk"
    mem: str = "32000MB"
    nprocshared: int = 8
    functional: str = "PBE1PBE"
    basis_set: str = "6-311G(d,p)"
    charge: int = 0
    multiplicity: int = 1
    opt: str | None = "Z-matrix, CalcAll,MaxCycle=150,MaxStep=20"
    integral: str | None = "ultrafine"
    scf: str | None = "QC"
    pop: str | None = "hlygat"
    others: str | None = "nosymm"
    pcm_eps: float | None = 11.0
    title: str = "Gas-phase optimization from crystal"
    walltime: str = "6:00:00"
    pbs_select: str = "1:ncpus=16:mpiprocs=16:mem=40gb"
    array_range: str = "1-2"
    gaussian_command: str = "g16"
    formchk_command: str = "formchk"
    run_script_name: str = "runGAUSSIAN_cx3.csh"


@dataclass(frozen=True)
class GaussianJobArtifacts:
    job_name: str
    signature: str
    duplicate_count: int
    molecule_formula: str
    ordered_atom_labels: list[str]
    molecule: PymatgenMolecule
    job_directory: Path
    zmat_path: Path
    com_path: Path
    run_script_path: Path


class GaussianInputBuilder:
    """Write Gaussian inputs and run scripts for unique crystal molecules."""

    def __init__(
        self,
        *,
        com_template_path: str | Path = DEFAULT_COM_TEMPLATE,
        run_template_path: str | Path = DEFAULT_RUN_TEMPLATE,
    ) -> None:
        self.com_template_path = Path(com_template_path)
        self.run_template_path = Path(run_template_path)

    def write_unique_jobs(
        self,
        structure: CrystalStructure,
        output_dir: str | Path,
        settings: GaussianSettings | None = None,
    ) -> list[GaussianJobArtifacts]:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        config = settings or GaussianSettings()

        deduplicated = structure.deduplicate_molecules()
        zmatrices = structure.generate_zmatrices()
        artifacts: list[GaussianJobArtifacts] = []

        for group_index, group in enumerate(deduplicated, start=1):
            zmatrix = self._match_representative_zmatrix(group, zmatrices)
            ordered_atoms = self._ordered_atoms(group.representative_molecule, zmatrix)
            pmg_molecule = self._to_pymatgen_molecule(ordered_atoms)
            job_name = f"{_safe_job_name(structure.name)}_Mol{group_index}"
            job_directory = output_path / job_name
            job_directory.mkdir(parents=True, exist_ok=True)

            zmat_path = job_directory / f"{job_name}.zmat"
            com_path = job_directory / f"{job_name}.com"
            run_script_path = job_directory / config.run_script_name
            zmat_path.write_text(
                self.render_zmat_text(job_name, zmatrix),
                encoding="utf-8",
            )
            com_path.write_text(self.render_com_text(zmatrix, config), encoding="utf-8")
            run_script_path.write_text(
                self.render_run_script_text(job_name, config),
                encoding="utf-8",
            )
            run_script_path.chmod(0o755)

            artifacts.append(
                GaussianJobArtifacts(
                    job_name=job_name,
                    signature=group.signature,
                    duplicate_count=len(group.duplicate_molecules),
                    molecule_formula=pmg_molecule.composition.alphabetical_formula,
                    ordered_atom_labels=zmatrix.ordered_atom_labels,
                    molecule=pmg_molecule,
                    job_directory=job_directory,
                    zmat_path=zmat_path,
                    com_path=com_path,
                    run_script_path=run_script_path,
                )
            )

        return artifacts

    def render_com_text(
        self,
        zmatrix: ZMatrixRepresentation,
        settings: GaussianSettings,
    ) -> str:
        config = _normalized_settings(settings)
        lines = [
            f"%chk={config.chk}",
            f"%mem={config.mem}",
            f"%nprocshared={config.nprocshared}",
            _gaussian_route(config),
            "",
            config.title,
            "",
            f"{config.charge} {config.multiplicity}",
        ]
        variable_lines: list[str] = []

        for atom_position, entry in enumerate(zmatrix.entries, start=1):
            variable_index = atom_position
            if atom_position == 1:
                lines.append(entry.element)
                continue
            if entry.bond_to is None or entry.bond_length is None:
                raise ValueError(f"Missing bond definition for Z-matrix atom {atom_position}.")
            if atom_position == 2:
                lines.append(f"{entry.element} {entry.bond_to} bnd{variable_index}")
                variable_lines.append(f"bnd{variable_index}={entry.bond_length:.6f}")
                continue
            if entry.angle_to is None or entry.angle_degrees is None:
                raise ValueError(f"Missing angle definition for Z-matrix atom {atom_position}.")
            if atom_position == 3:
                lines.append(
                    f"{entry.element} {entry.bond_to} bnd{variable_index} "
                    f"{entry.angle_to} ang{variable_index}"
                )
                variable_lines.append(f"bnd{variable_index}={entry.bond_length:.6f}")
                variable_lines.append(f"ang{variable_index}={entry.angle_degrees:.6f}")
                continue
            if entry.dihedral_to is None or entry.dihedral_degrees is None:
                raise ValueError(
                    f"Missing dihedral definition for Z-matrix atom {atom_position}."
                )
            lines.append(
                f"{entry.element} {entry.bond_to} bnd{variable_index} "
                f"{entry.angle_to} ang{variable_index} "
                f"{entry.dihedral_to} dih{variable_index}"
            )
            variable_lines.append(f"bnd{variable_index}={entry.bond_length:.6f}")
            variable_lines.append(f"ang{variable_index}={entry.angle_degrees:.6f}")
            variable_lines.append(f"dih{variable_index}={entry.dihedral_degrees:.6f}")

        lines.append("")
        lines.append("Variables:")
        lines.extend(variable_lines)
        lines.append("Constants:")
        if config.pcm_eps is not None:
            lines.append(f"EPS={config.pcm_eps:.1f}")
        lines.extend(["", ""])
        return "\n".join(lines)

    def render_zmat_text(
        self,
        title: str,
        zmatrix: ZMatrixRepresentation,
    ) -> str:
        lines = [
            "# ZMAT v1",
            f"# title: {title}",
            "# bonds: " + " ".join(
                f"{entry.bond_to}-{index}"
                for index, entry in enumerate(zmatrix.entries, start=1)
                if entry.bond_to is not None
            ),
        ]
        torsion_index = 0
        for index, entry in enumerate(zmatrix.entries, start=1):
            if entry.dihedral_to is None or entry.angle_to is None or entry.bond_to is None:
                continue
            lines.append(
                f"# torsion {torsion_index}: "
                f"{index}-{entry.dihedral_to}-{entry.angle_to}-{entry.bond_to} "
                f"bond {index}-{entry.bond_to}"
            )
            torsion_index += 1

        for index, entry in enumerate(zmatrix.entries, start=1):
            row = [entry.element]
            if entry.bond_to is not None and entry.bond_length is not None:
                row.extend([str(entry.bond_to), f"{entry.bond_length:.10f}"])
            if entry.angle_to is not None and entry.angle_degrees is not None:
                row.extend([str(entry.angle_to), f"{entry.angle_degrees:.10f}"])
            if entry.dihedral_to is not None and entry.dihedral_degrees is not None:
                row.extend([str(entry.dihedral_to), f"{entry.dihedral_degrees:.10f}"])
            lines.append(" ".join(row))

        return "\n".join(lines) + "\n"

    def render_run_script_text(
        self,
        job_name: str,
        settings: GaussianSettings,
    ) -> str:
        config = _normalized_settings(settings)
        if self.run_template_path.exists():
            template = self.run_template_path.read_text(encoding="utf-8")
            script = template.replace("SYSTEM_NAME", job_name)
            script = re.sub(
                r"^#PBS -lwalltime=.*$",
                f"#PBS -lwalltime={config.walltime}",
                script,
                flags=re.MULTILINE,
            )
            script = re.sub(
                r"^#PBS -l select=.*$",
                f"#PBS -l select={config.pbs_select}",
                script,
                flags=re.MULTILINE,
            )
            script = re.sub(
                r"^#PBS -J .*?$",
                f"#PBS -J {config.array_range}",
                script,
                flags=re.MULTILINE,
            )
            script = re.sub(
                r"^g16\s+\$\{name\}\.com\s+\$\{name\}\.log\s*$",
                f"{config.gaussian_command} ${{name}}.com ${{name}}.log",
                script,
                flags=re.MULTILINE,
            )
            script = re.sub(
                r"^formchk -3 .*?$",
                (
                    f"{config.formchk_command} -3 {config.chk} "
                    f"{_formatted_checkpoint_name(config.chk)}"
                ),
                script,
                flags=re.MULTILINE,
            )
            return script.rstrip() + "\n"

        return (
            "#!/bin/sh\n"
            f"#PBS -lwalltime={config.walltime}\n"
            f"#PBS -l select={config.pbs_select}\n"
            f"#PBS -N {job_name}\n"
            f"#PBS -J {config.array_range}\n\n"
            "# Just adding this to bypass slow queues\n"
            "if [[ $PBS_ARRAY_INDEX -eq 2 ]];\n"
            "then\n"
            '    echo "Dummy array job"\n'
            "    exit\n"
            "fi\n\n"
            "#DEFINE directories\n"
            f"name={job_name}\n\n"
            "EPHEMERAL=~/../ephemeral/${PBS_JOBNAME}_${name}/\n"
            "INPUT=${PBS_O_WORKDIR}/\n"
            "OUTPUT=${PBS_O_WORKDIR}/\n\n"
            "#Load modules\n"
            "module load tools/prod\n"
            "module load Gaussian/16.C.02-AVX2\n\n"
            "# CREATE directories\n"
            "rm -rf $EPHEMERAL\n"
            "mkdir -p $EPHEMERAL\n\n"
            "#Execute\n"
            "cd $EPHEMERAL\n"
            "pwd\n"
            "cp ${INPUT}/${name}.com $EPHEMERAL/\n"
            f"{config.gaussian_command} ${{name}}.com ${{name}}.log \n"
            f"{config.formchk_command} -3 {config.chk} {_formatted_checkpoint_name(config.chk)}\n\n"
            "#Copy back\n"
            "cp $EPHEMERAL/* ${OUTPUT}/\n"
        )

    def _match_representative_zmatrix(
        self,
        group: MoleculeGroup,
        zmatrices: list[ZMatrixRepresentation],
    ) -> ZMatrixRepresentation:
        representative_labels = frozenset(atom.label for atom in group.representative_molecule)
        for zmatrix in zmatrices:
            if frozenset(zmatrix.ordered_atom_labels) == representative_labels:
                return zmatrix
        raise ValueError(
            "Could not match a Z-matrix representation to the deduplicated representative molecule."
        )

    @staticmethod
    def _ordered_atoms(
        molecule: list[AtomRecord],
        zmatrix: ZMatrixRepresentation,
    ) -> list[AtomRecord]:
        atoms_by_label = {atom.label: atom for atom in molecule}
        return [atoms_by_label[label] for label in zmatrix.ordered_atom_labels]

    @staticmethod
    def _to_pymatgen_molecule(molecule: list[AtomRecord]) -> PymatgenMolecule:
        species = [atom.element for atom in molecule]
        coordinates = [atom.coordinates for atom in molecule]
        return PymatgenMolecule(species, coordinates)


def _formatted_checkpoint_name(chk_name: str) -> str:
    chk_path = Path(chk_name)
    if chk_path.suffix.lower() == ".chk":
        return chk_path.with_suffix(".fchk").name
    return f"{chk_path.name}.fchk"


def _gaussian_route(settings: GaussianSettings) -> str:
    pieces = [f"{settings.functional}/{settings.basis_set}"]
    if settings.integral:
        pieces.append(f"int=({settings.integral})")
    if settings.others:
        pieces.append(settings.others)
    if settings.opt:
        pieces.append(f"opt=({settings.opt})")
    if settings.pop:
        pieces.append(f"pop=({settings.pop})")
    if settings.scf:
        pieces.append(f"scf=({settings.scf})")
    if settings.pcm_eps is not None:
        pieces.append("scrf=(PCM,Read)")
    return "#P " + " ".join(pieces)


def _normalized_settings(settings: GaussianSettings) -> GaussianSettings:
    return GaussianSettings(
        chk=settings.chk,
        mem=settings.mem,
        nprocshared=settings.nprocshared,
        functional=settings.functional,
        basis_set=settings.basis_set,
        charge=settings.charge,
        multiplicity=settings.multiplicity,
        opt=_normalize_optional_keyword(settings.opt),
        integral=_normalize_optional_keyword(settings.integral),
        scf=_normalize_optional_keyword(settings.scf),
        pop=_normalize_optional_keyword(settings.pop),
        others=_normalize_optional_keyword(settings.others),
        pcm_eps=settings.pcm_eps,
        title=settings.title,
        walltime=settings.walltime,
        pbs_select=settings.pbs_select,
        array_range=settings.array_range,
        gaussian_command=settings.gaussian_command,
        formchk_command=settings.formchk_command,
        run_script_name=settings.run_script_name,
    )


def _normalize_optional_keyword(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    if not stripped or stripped.lower() == "none":
        return None
    return stripped


def _safe_job_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return cleaned or "CrystalStructure"
