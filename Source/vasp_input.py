"""VASP input generation for periodic and gas-phase crystal structures."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from pathlib import Path
from typing import Iterable

from .crystal_structure import AtomRecord, CrystalStructure


MIN_DIST = 0.025 * 2 * math.pi
TEMPLATE_DIR = Path(__file__).resolve().parents[1] / "Template" / "VASP_input"
DEFAULT_TPSS_INCAR_TEMPLATE = TEMPLATE_DIR / "TPSS_INCAR_1000eV"
DEFAULT_PBE0_INCAR_TEMPLATE = TEMPLATE_DIR / "PBE0_INCAR_1000eV_point_calc"
DEFAULT_INCAR_TEMPLATE = DEFAULT_TPSS_INCAR_TEMPLATE
DEFAULT_RUN_TEMPLATE = TEMPLATE_DIR / "run_vasp.sh"
DEFAULT_POTCAR_LIB = TEMPLATE_DIR / "potpawPBE54"
VASP_PRESETS = {
    "tpss_relax": DEFAULT_TPSS_INCAR_TEMPLATE,
    "pbe0_point": DEFAULT_PBE0_INCAR_TEMPLATE,
}


@dataclass(frozen=True)
class VaspSettings:
    incar_template_path: Path = DEFAULT_INCAR_TEMPLATE
    run_template_path: Path = DEFAULT_RUN_TEMPLATE
    potcar_lib_path: Path = DEFAULT_POTCAR_LIB
    run_script_name: str = "run_vasp_1000eV.sh"
    kpoint_spacing_multiplier: float | None = None

    @classmethod
    def for_preset(
        cls,
        preset: str,
        **overrides,
    ) -> "VaspSettings":
        """Create settings for a named VASP setup.

        ``preset`` is intentionally small and explicit. It selects a template
        set only; callers can still override paths or k-point spacing.
        """

        if preset not in VASP_PRESETS:
            raise ValueError(
                f"Unknown VASP preset {preset!r}; expected one of {sorted(VASP_PRESETS)}"
            )
        settings = cls(incar_template_path=VASP_PRESETS[preset])
        if overrides:
            settings = replace(settings, **overrides)
        return settings


@dataclass(frozen=True)
class VaspJobArtifacts:
    system_name: str
    structure: CrystalStructure
    output_dir: Path
    poscar_path: Path
    kpoints_path: Path
    incar_path: Path
    potcar_path: Path
    run_script_path: Path


class VaspInputBuilder:
    """Generate VASP inputs from explicit-unit-cell crystal structures."""

    def __init__(self, settings: VaspSettings | None = None) -> None:
        self.settings = settings or VaspSettings()

    def write_job_from_structure(
        self,
        structure: CrystalStructure,
        output_dir: str | Path,
        *,
        system_name: str | None = None,
        kpoints_text: str | None = None,
        incar_text: str | None = None,
        run_script_text: str | None = None,
    ) -> VaspJobArtifacts:
        _require_explict_unit_cell(
            structure,
            context="VASP input generation requires an explict unit cell structure.",
        )
        job_name = system_name or structure.name
        destination = Path(output_dir)
        destination.mkdir(parents=True, exist_ok=True)

        poscar_path = destination / "POSCAR"
        kpoints_path = destination / "KPOINTS"
        incar_path = destination / "INCAR"
        potcar_path = destination / "POTCAR"
        run_script_path = destination / self.settings.run_script_name

        poscar_path.write_text(self.render_poscar_text(job_name, structure), encoding="utf-8")
        kpoints_path.write_text(
            kpoints_text or self.render_kpoints_text(structure.cell.array.tolist()),
            encoding="utf-8",
        )
        incar_path.write_text(
            incar_text or self.render_incar_text(job_name),
            encoding="utf-8",
        )
        potcar_path.write_text(self.build_potcar_text(structure.atoms), encoding="utf-8")
        run_script_path.write_text(
            run_script_text or self.render_run_script_text(job_name),
            encoding="utf-8",
        )
        run_script_path.chmod(0o755)

        return VaspJobArtifacts(
            system_name=job_name,
            structure=structure,
            output_dir=destination,
            poscar_path=poscar_path,
            kpoints_path=kpoints_path,
            incar_path=incar_path,
            potcar_path=potcar_path,
            run_script_path=run_script_path,
        )

    def write_job_from_cif(
        self,
        cif_path: str | Path,
        output_dir: str | Path,
        *,
        system_name: str | None = None,
        expand_unit_cell: bool = True,
    ) -> VaspJobArtifacts:
        """Generate a VASP job from a CIF file."""

        path = Path(cif_path)
        if expand_unit_cell:
            structure = CrystalStructure.from_cif_unit_cell(path)
        else:
            structure = CrystalStructure.from_file(path, fmt="cif")
            if not structure.explict_unit_cell:
                structure = structure.expand_to_explicit_unit_cell()
        return self.write_job_from_structure(
            structure,
            output_dir,
            system_name=system_name or path.stem,
        )

    def write_job_from_contcar(
        self,
        contcar_path: str | Path,
        output_dir: str | Path,
        *,
        system_name: str | None = None,
    ) -> VaspJobArtifacts:
        """Generate a VASP job from a POSCAR/CONTCAR-like file."""

        from .vasp_results import read_contcar_as_crystal

        path = Path(contcar_path)
        structure = read_contcar_as_crystal(
            path,
            name=system_name or path.parent.name,
        )
        return self.write_job_from_structure(
            structure,
            output_dir,
            system_name=system_name or structure.name,
        )

    def write_gas_phase_job(
        self,
        structure: CrystalStructure,
        output_dir: str | Path,
        *,
        system_name: str | None = None,
    ) -> VaspJobArtifacts:
        job_name = system_name or structure.name
        return self.write_job_from_structure(
            structure,
            output_dir,
            system_name=job_name,
            kpoints_text=self.render_gas_phase_kpoints_text(),
            incar_text=self.render_gas_phase_incar_text(job_name),
            run_script_text=self.render_gas_phase_run_script_text(job_name),
        )

    def render_poscar_text(self, system_name: str, structure: CrystalStructure) -> str:
        ordered_elements = _element_order(atom.element for atom in structure.atoms)
        counts = [sum(1 for atom in structure.atoms if atom.element == element) for element in ordered_elements]

        lines = [system_name, "1.00"]
        for vector in structure.cell.array.tolist():
            lines.append(f"{vector[0]:15.12f} {vector[1]:15.12f} {vector[2]:15.12f}")
        lines.append("  ".join(ordered_elements))
        lines.append("  ".join(str(count) for count in counts))
        lines.append("Cartesian")

        for element in ordered_elements:
            for atom in structure.atoms:
                if atom.element != element:
                    continue
                x, y, z = atom.coordinates
                lines.append(f"{x:15.12f} {y:15.12f} {z:15.12f}")

        return "\n".join(lines) + "\n"

    def render_kpoints_text(self, lattice: list[tuple[float, float, float]]) -> str:
        reciprocal = _reciprocal_lattice(lattice)
        spacing = (
            self.settings.kpoint_spacing_multiplier * 2 * math.pi
            if self.settings.kpoint_spacing_multiplier is not None
            else MIN_DIST
        )
        if spacing <= 0.0:
            raise ValueError("k-point spacing must be positive.")
        grid = [max(1, math.ceil(_norm(vector) / spacing)) for vector in reciprocal]
        lines = [
            "AutoKPOINTS",
            "0",
            "Gamma",
            f" {grid[0]} {grid[1]} {grid[2]}",
            " 0 0 0",
        ]
        return "\n".join(lines) + "\n"

    def render_incar_text(self, system_name: str) -> str:
        template_text = self.settings.incar_template_path.read_text(encoding="utf-8")
        return template_text.replace("Dummy_System", system_name)

    def render_run_script_text(self, system_name: str) -> str:
        template_text = self.settings.run_template_path.read_text(encoding="utf-8")
        return template_text.replace("Dummy_System", system_name)

    def render_gas_phase_kpoints_text(self) -> str:
        lines = [
            "AutoKPOINTS isol",
            "0",
            "Gamma",
            " 1 1 1",
            " 0 0 0",
        ]
        return "\n".join(lines) + "\n"

    def render_gas_phase_incar_text(self, system_name: str) -> str:
        template_text = self.render_incar_text(system_name)
        updated_lines: list[str] = []
        for line in template_text.splitlines():
            stripped = line.strip()
            if stripped.startswith("#SIGMA="):
                updated_lines.append("SIGMA= 0.05      # Width of smearing (eV) (Ignore for tetra)")
            elif stripped.startswith("ISMEAR=") or stripped.startswith("ISMEAR ="):
                updated_lines.append(
                    "ISMEAR= 0        # change between 0, -5; Type (0 = Gauss, -1 = Fermi, -2 = fixed INCAR, -3 = looped INCAR, -4 = tetra, -5 = Blochl corr. tetra)"
                )
            elif stripped.startswith("NSW=") or stripped.startswith("NSW ="):
                updated_lines.append(
                    "NSW= 0          # Number of steps allowed for relaxation (0 = single-point/no nuclei movement)"
                )
            elif stripped.startswith("IBRION=") or stripped.startswith("IBRION ="):
                updated_lines.append(
                    "IBRION= -1         # change according to convergence; -1 = single-point calculation, 0 = MD calculation, 1-8 = relaxation"
                )
            else:
                updated_lines.append(line)
        return "\n".join(updated_lines) + "\n"

    def render_gas_phase_run_script_text(self, system_name: str) -> str:
        script_text = self.render_run_script_text(system_name)
        script_text = script_text.replace("#$ -l h_rt=24:00:00", "#$ -l h_rt=4:00:00")
        script_text = script_text.replace("#$ -pe mpi 160", "#$ -pe mpi 40")
        return script_text

    def build_potcar_text(self, atoms: list[AtomRecord]) -> str:
        ordered_elements = _element_order(atom.element for atom in atoms)
        fragments: list[str] = []
        for element in ordered_elements:
            folder = "Cl" if element.upper() == "CL" else element
            potcar_path = self.settings.potcar_lib_path / folder / "POTCAR"
            if not potcar_path.exists():
                raise FileNotFoundError(f"Missing POTCAR for element {element}: {potcar_path}")
            fragments.append(potcar_path.read_text(encoding="utf-8"))
        return "".join(fragments)


def _cross(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
) -> tuple[float, float, float]:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _dot(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _require_explict_unit_cell(structure: CrystalStructure, *, context: str) -> None:
    if not structure.explict_unit_cell:
        raise ValueError(context)


def _norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(_dot(vector, vector))


def _reciprocal_lattice(
    lattice: list[tuple[float, float, float]],
) -> list[tuple[float, float, float]]:
    a_vec, b_vec, c_vec = lattice
    volume = _dot(a_vec, _cross(b_vec, c_vec))
    factor = 2 * math.pi / volume
    return [
        tuple(factor * value for value in _cross(b_vec, c_vec)),
        tuple(factor * value for value in _cross(c_vec, a_vec)),
        tuple(factor * value for value in _cross(a_vec, b_vec)),
    ]


def _element_order(elements: Iterable[str]) -> list[str]:
    ordered: list[str] = []
    for element in elements:
        if element not in ordered:
            ordered.append(element)
    return ordered
