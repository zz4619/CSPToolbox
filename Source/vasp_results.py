"""Parse nested VASP relaxation results for one system."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import TYPE_CHECKING
from xml.etree import ElementTree as ET

import numpy as np

if TYPE_CHECKING:
    from .crystal_structure import CrystalStructure


@dataclass(frozen=True)
class VaspCalculationResult:
    system_name: str
    calculation_index: int
    calculation_chain: tuple[int, ...]
    calculation_dir: Path
    vasprun_path: Path
    ionic_iterations: int | None
    total_electronic_iterations: int | None
    last_electronic_iterations: int | None
    final_energy_ev: float | None
    converged: bool | None
    converged_ionic: bool | None
    converged_electronic: bool | None
    cpu_time_seconds: float | None
    parse_error: str | None = None
    status: str = "unknown"

    @property
    def contcar_path(self) -> Path:
        return self.calculation_dir / "CONTCAR"

    def read_contcar_as_crystal(
        self,
        *,
        contcar_root: str | Path | None = None,
    ) -> CrystalStructure:
        return read_contcar_as_crystal(
            _resolve_contcar_path(self, contcar_root=contcar_root),
            name=f"{self.system_name}_calc{self.calculation_index}_CONTCAR",
        )


@dataclass(frozen=True)
class VaspSystemResult:
    system_name: str
    system_dir: Path
    calculation_count: int
    calculations: list[VaspCalculationResult]

    @property
    def latest_calculation(self) -> VaspCalculationResult | None:
        if not self.calculations:
            return None
        return self.calculations[-1]

    @property
    def latest_energy_ev(self) -> float | None:
        latest = self.latest_calculation
        if latest is None:
            return None
        return latest.final_energy_ev

    @property
    def latest_converged(self) -> bool | None:
        latest = self.latest_calculation
        if latest is None:
            return None
        return latest.converged

    def read_latest_contcar_as_crystal(
        self,
        *,
        contcar_root: str | Path | None = None,
    ) -> CrystalStructure:
        latest = self.latest_calculation
        if latest is None:
            raise ValueError(f"No calculations parsed for system {self.system_name}")
        return latest.read_contcar_as_crystal(contcar_root=contcar_root)


@dataclass(frozen=True)
class VaspOutSummary:
    path: Path
    exists: bool
    converged: bool | None
    final_free_energy_ev: float | None
    final_energy_sigma0_ev: float | None
    final_iteration_label: str | None
    fatal_markers: tuple[str, ...]


@dataclass(frozen=True)
class OutcarSummary:
    path: Path
    exists: bool
    reached_required_accuracy: bool | None
    has_general_timing: bool | None
    cpu_time_seconds: float | None
    max_memory_kb: float | None
    fatal_markers: tuple[str, ...]


@dataclass(frozen=True)
class OutcarForceAtom:
    atom_index: int
    element: str | None
    position: tuple[float, float, float]
    force: tuple[float, float, float]

    @property
    def force_norm_ev_ang(self) -> float:
        return float(np.linalg.norm(np.array(self.force, dtype=float)))


@dataclass(frozen=True)
class OutcarForceBlock:
    ionic_step: int
    atoms: tuple[OutcarForceAtom, ...]
    total_drift: tuple[float, float, float] | None

    @property
    def atom_count(self) -> int:
        return len(self.atoms)

    @property
    def positions_ang(self) -> np.ndarray:
        return np.array([atom.position for atom in self.atoms], dtype=float)

    @property
    def forces_ev_ang(self) -> np.ndarray:
        return np.array([atom.force for atom in self.atoms], dtype=float)

    @property
    def force_norms_ev_ang(self) -> np.ndarray:
        forces = self.forces_ev_ang
        if forces.size == 0:
            return np.array([], dtype=float)
        return np.linalg.norm(forces, axis=1)

    @property
    def max_force_norm_ev_ang(self) -> float | None:
        norms = self.force_norms_ev_ang
        if norms.size == 0:
            return None
        return float(np.max(norms))

    @property
    def rms_force_norm_ev_ang(self) -> float | None:
        norms = self.force_norms_ev_ang
        if norms.size == 0:
            return None
        return float(np.sqrt(np.mean(norms**2)))

    @property
    def mean_force_norm_ev_ang(self) -> float | None:
        norms = self.force_norms_ev_ang
        if norms.size == 0:
            return None
        return float(np.mean(norms))

    @property
    def p95_force_norm_ev_ang(self) -> float | None:
        norms = self.force_norms_ev_ang
        if norms.size == 0:
            return None
        return float(np.percentile(norms, 95))


@dataclass(frozen=True)
class OutcarForces:
    path: Path
    exists: bool
    nions: int | None
    ion_types: tuple[str, ...]
    ion_counts: tuple[int, ...]
    total_force_block_count: int
    blocks: tuple[OutcarForceBlock, ...]

    @property
    def force_block_count(self) -> int:
        return self.total_force_block_count

    @property
    def retained_force_block_count(self) -> int:
        return len(self.blocks)

    @property
    def final_block(self) -> OutcarForceBlock | None:
        if not self.blocks:
            return None
        return self.blocks[-1]


@dataclass(frozen=True)
class VaspCalculationHealth:
    calculation_dir: Path
    status: str
    fatal_markers: tuple[str, ...]
    vasp_out: VaspOutSummary
    outcar: OutcarSummary


class VaspSystemParser:
    """Parse one system directory containing nested VASP minimization runs."""

    def __init__(self, *, parse_potcar_file: bool = False) -> None:
        self.parse_potcar_file = parse_potcar_file

    def parse_system(self, system_dir: str | Path) -> VaspSystemResult:
        root = Path(system_dir)
        if not root.is_dir():
            raise FileNotFoundError(f"System directory not found: {root}")

        vasprun_paths = sorted(
            root.rglob("vasprun.xml"),
            key=lambda path: _calculation_sort_key(root, path),
        )
        if not vasprun_paths:
            raise FileNotFoundError(f"No vasprun.xml files found under {root}")

        calculations = [self.parse_calculation(root, path) for path in vasprun_paths]
        return VaspSystemResult(
            system_name=root.name,
            system_dir=root,
            calculation_count=len(calculations),
            calculations=calculations,
        )

    def parse_calculation(
        self,
        system_dir: str | Path,
        vasprun_path: str | Path,
    ) -> VaspCalculationResult:
        root = Path(system_dir)
        xml_path = Path(vasprun_path)
        calc_dir = xml_path.parent
        chain = _calculation_chain(root, xml_path)
        health = classify_vasp_calculation_dir(calc_dir)

        try:
            from pymatgen.io.vasp.outputs import Vasprun

            vasprun = Vasprun(
                str(xml_path),
                parse_dos=False,
                parse_eigen=False,
                parse_projected_eigen=False,
                parse_potcar_file=self.parse_potcar_file,
            )
        except Exception as error:
            return VaspCalculationResult(
                system_name=root.name,
                calculation_index=len(chain),
                calculation_chain=chain,
                calculation_dir=calc_dir,
                vasprun_path=xml_path,
                ionic_iterations=None,
                total_electronic_iterations=None,
                last_electronic_iterations=None,
                final_energy_ev=None,
                converged=None,
                converged_ionic=None,
                converged_electronic=None,
                cpu_time_seconds=health.outcar.cpu_time_seconds or _parse_cpu_time_seconds(calc_dir),
                parse_error=str(error),
                status=health.status,
            )

        ionic_steps = list(vasprun.ionic_steps)
        electronic_counts = [len(step.get("electronic_steps", [])) for step in ionic_steps]

        return VaspCalculationResult(
            system_name=root.name,
            calculation_index=len(chain),
            calculation_chain=chain,
            calculation_dir=calc_dir,
            vasprun_path=xml_path,
            ionic_iterations=int(vasprun.nionic_steps),
            total_electronic_iterations=sum(electronic_counts),
            last_electronic_iterations=electronic_counts[-1] if electronic_counts else 0,
            final_energy_ev=float(vasprun.final_energy),
            converged=bool(vasprun.converged),
            converged_ionic=bool(vasprun.converged_ionic),
            converged_electronic=bool(vasprun.converged_electronic),
            cpu_time_seconds=health.outcar.cpu_time_seconds or _parse_cpu_time_seconds(calc_dir),
            parse_error=None,
            status="completed_converged" if bool(vasprun.converged) else health.status,
        )

    def read_contcar_as_crystal(
        self,
        contcar_path: str | Path,
        *,
        name: str | None = None,
        space_group: str = "P 1",
    ) -> CrystalStructure:
        return read_contcar_as_crystal(contcar_path, name=name, space_group=space_group)


def _calculation_chain(system_dir: Path, vasprun_path: Path) -> tuple[int, ...]:
    relative_dir = vasprun_path.parent.relative_to(system_dir)
    if relative_dir == Path("."):
        return ()

    chain: list[int] = []
    for part in relative_dir.parts:
        if not part.isdigit():
            raise ValueError(
                f"Unexpected non-numeric minimization layer '{part}' under {system_dir}"
            )
        chain.append(int(part))
    return tuple(chain)


def _calculation_sort_key(system_dir: Path, vasprun_path: Path) -> tuple[int, tuple[int, ...]]:
    chain = _calculation_chain(system_dir, vasprun_path)
    return len(chain), chain


def _parse_cpu_time_seconds(calculation_dir: Path) -> float | None:
    vaspout_path = calculation_dir / "vaspout.xml"
    if not vaspout_path.exists():
        return None

    text = vaspout_path.read_text(encoding="utf-8", errors="ignore")
    explicit_patterns = (
        r"Total CPU time used \(sec\)\s*:\s*([0-9.+\-Ee]+)",
        r"cpu time\s*=\s*([0-9.+\-Ee]+)",
    )
    for pattern in explicit_patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return float(match.group(1))

    try:
        root = ET.fromstring(text)
    except ET.ParseError:
        return None

    for element in reversed(list(root.iter())):
        if element.tag != "time":
            continue
        if element.attrib.get("name", "").lower() not in {"total", "totalsc"}:
            continue
        parsed = _parse_time_pair(element.text or "")
        if parsed is not None:
            return parsed[0]
    return None


def parse_vasp_out(path: str | Path) -> VaspOutSummary:
    """Parse convergence and final energy markers from ``vasp.out``."""

    vasp_out_path = Path(path)
    if not vasp_out_path.exists():
        return VaspOutSummary(
            path=vasp_out_path,
            exists=False,
            converged=None,
            final_free_energy_ev=None,
            final_energy_sigma0_ev=None,
            final_iteration_label=None,
            fatal_markers=(),
        )

    text = vasp_out_path.read_text(encoding="utf-8", errors="ignore")
    converged = "reached required accuracy - stopping structural energy minimisation" in text
    fatal_markers = _detect_fatal_markers(text)

    final_free_energy_ev: float | None = None
    final_energy_sigma0_ev: float | None = None
    final_iteration_label: str | None = None
    for line in reversed(text.splitlines()):
        if "F=" not in line:
            continue
        free_match = re.search(r"\bF=\s*([+-]?[0-9.]+(?:[Ee][+-]?[0-9]+)?)", line)
        sigma_match = re.search(r"\bE0=\s*([+-]?[0-9.]+(?:[Ee][+-]?[0-9]+)?)", line)
        if free_match:
            final_free_energy_ev = float(free_match.group(1))
            final_iteration_label = line.split("F=", 1)[0].strip() or None
        if sigma_match:
            final_energy_sigma0_ev = float(sigma_match.group(1))
        break

    return VaspOutSummary(
        path=vasp_out_path,
        exists=True,
        converged=converged,
        final_free_energy_ev=final_free_energy_ev,
        final_energy_sigma0_ev=final_energy_sigma0_ev,
        final_iteration_label=final_iteration_label,
        fatal_markers=fatal_markers,
    )


def parse_outcar_summary(path: str | Path) -> OutcarSummary:
    """Parse completion, timing, memory, and fatal markers from ``OUTCAR``."""

    outcar_path = Path(path)
    if not outcar_path.exists():
        return OutcarSummary(
            path=outcar_path,
            exists=False,
            reached_required_accuracy=None,
            has_general_timing=None,
            cpu_time_seconds=None,
            max_memory_kb=None,
            fatal_markers=(),
        )

    text = outcar_path.read_text(encoding="utf-8", errors="ignore")
    cpu_time = _first_float_match(text, r"Total CPU time used \(sec\)\s*:\s*([0-9.+\-Ee]+)")
    max_memory = _first_float_match(text, r"Maximum memory used \(kb\)\s*:\s*([0-9.+\-Ee]+)")
    return OutcarSummary(
        path=outcar_path,
        exists=True,
        reached_required_accuracy="reached required accuracy - stopping structural energy minimisation" in text,
        has_general_timing="General timing and accounting" in text,
        cpu_time_seconds=cpu_time,
        max_memory_kb=max_memory,
        fatal_markers=_detect_fatal_markers(text),
    )


def parse_outcar_forces(
    path: str | Path,
    *,
    all_ionic_steps: bool = False,
) -> OutcarForces:
    """Parse atom-resolved force blocks from a VASP ``OUTCAR``.

    By default only the final ``POSITION / TOTAL-FORCE`` block is retained,
    which is the usual geometry-quality diagnostic. Set ``all_ionic_steps`` to
    retain every ionic step from a relaxation trajectory.
    """

    outcar_path = Path(path)
    if not outcar_path.exists():
        return OutcarForces(
            path=outcar_path,
            exists=False,
            nions=None,
            ion_types=(),
            ion_counts=(),
            total_force_block_count=0,
            blocks=(),
        )

    lines = outcar_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    nions = _parse_outcar_nions(lines)
    ion_types = _parse_outcar_ion_types(lines)
    ion_counts = _parse_outcar_ion_counts(lines)
    expanded_elements = _expand_ion_elements(ion_types, ion_counts)
    if nions is not None and len(expanded_elements) != nions:
        expanded_elements = ()

    blocks: list[OutcarForceBlock] = []
    block_index = 0
    line_index = 0
    while line_index < len(lines):
        if not _is_force_block_header(lines[line_index]):
            line_index += 1
            continue

        block_index += 1
        block, line_index = _parse_force_block(
            lines,
            line_index + 1,
            ionic_step=block_index,
            expanded_elements=expanded_elements,
        )
        if all_ionic_steps:
            blocks.append(block)
        else:
            blocks = [block]

    return OutcarForces(
        path=outcar_path,
        exists=True,
        nions=nions,
        ion_types=ion_types,
        ion_counts=ion_counts,
        total_force_block_count=block_index,
        blocks=tuple(blocks),
    )


def classify_vasp_calculation_dir(calculation_dir: str | Path) -> VaspCalculationHealth:
    """Classify one VASP calculation folder from text output files.

    The classifier is intentionally conservative: fatal markers win over
    completion markers, and complete-but-not-converged runs remain distinct.
    """

    directory = Path(calculation_dir)
    vasp_out = parse_vasp_out(directory / "vasp.out")
    outcar = parse_outcar_summary(directory / "OUTCAR")
    fatal_markers = tuple(dict.fromkeys((*vasp_out.fatal_markers, *outcar.fatal_markers)))

    if fatal_markers:
        status = "failed"
    elif outcar.has_general_timing and (outcar.reached_required_accuracy or vasp_out.converged):
        status = "completed_converged"
    elif outcar.has_general_timing:
        status = "completed_unconverged"
    elif vasp_out.converged:
        status = "completed_converged"
    elif vasp_out.final_free_energy_ev is not None:
        status = "incomplete_or_running"
    elif not vasp_out.exists and not outcar.exists and not (directory / "vasprun.xml").exists():
        status = "missing_output"
    else:
        status = "unknown"

    return VaspCalculationHealth(
        calculation_dir=directory,
        status=status,
        fatal_markers=fatal_markers,
        vasp_out=vasp_out,
        outcar=outcar,
    )


_FLOAT_PATTERN = r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?"
_FORCE_ROW_RE = re.compile(
    rf"^\s*({_FLOAT_PATTERN})\s+({_FLOAT_PATTERN})\s+({_FLOAT_PATTERN})\s+"
    rf"({_FLOAT_PATTERN})\s+({_FLOAT_PATTERN})\s+({_FLOAT_PATTERN})(?:\s+.*)?$"
)


def _parse_outcar_nions(lines: list[str]) -> int | None:
    for line in lines:
        match = re.search(r"\bNIONS\s*=\s*(\d+)", line)
        if match:
            return int(match.group(1))
    return None


def _parse_outcar_ion_types(lines: list[str]) -> tuple[str, ...]:
    ion_types: list[str] = []
    for line in lines:
        match = re.search(r"\bVRHFIN\s*=\s*([A-Za-z][A-Za-z]?)\s*:", line)
        if match:
            ion_types.append(_normalize_poscar_symbol(match.group(1)))
    return tuple(ion_types)


def _parse_outcar_ion_counts(lines: list[str]) -> tuple[int, ...]:
    for line in lines:
        if "ions per type" not in line:
            continue
        _, _, tail = line.partition("=")
        try:
            return tuple(int(token) for token in tail.split())
        except ValueError:
            return ()
    return ()


def _expand_ion_elements(
    ion_types: tuple[str, ...],
    ion_counts: tuple[int, ...],
) -> tuple[str, ...]:
    if len(ion_types) != len(ion_counts):
        return ()
    elements: list[str] = []
    for ion_type, count in zip(ion_types, ion_counts, strict=True):
        elements.extend([ion_type] * count)
    return tuple(elements)


def _is_force_block_header(line: str) -> bool:
    return "POSITION" in line and "TOTAL-FORCE" in line


def _parse_force_block(
    lines: list[str],
    line_index: int,
    *,
    ionic_step: int,
    expanded_elements: tuple[str, ...],
) -> tuple[OutcarForceBlock, int]:
    atoms: list[OutcarForceAtom] = []
    total_drift: tuple[float, float, float] | None = None
    index = line_index

    while index < len(lines):
        line = lines[index]
        if atoms and _is_separator_line(line):
            index += 1
            break
        match = _FORCE_ROW_RE.match(line)
        if match:
            values = tuple(float(match.group(number)) for number in range(1, 7))
            atom_index = len(atoms) + 1
            element = expanded_elements[atom_index - 1] if atom_index <= len(expanded_elements) else None
            atoms.append(
                OutcarForceAtom(
                    atom_index=atom_index,
                    element=element,
                    position=values[:3],
                    force=values[3:],
                )
            )
        index += 1

    for drift_index in range(index, min(index + 5, len(lines))):
        if "total drift" not in lines[drift_index]:
            continue
        values = [float(value) for value in re.findall(_FLOAT_PATTERN, lines[drift_index])]
        if len(values) >= 3:
            total_drift = tuple(values[-3:])
        index = drift_index + 1
        break

    return (
        OutcarForceBlock(
            ionic_step=ionic_step,
            atoms=tuple(atoms),
            total_drift=total_drift,
        ),
        index,
    )


def _is_separator_line(line: str) -> bool:
    stripped = line.strip()
    return bool(stripped) and set(stripped) == {"-"}


def _detect_fatal_markers(text: str) -> tuple[str, ...]:
    marker_patterns = {
        "ibzkpt_error": r"\bIBZKPT\b",
        "insufficient_memory": r"insufficient memory|out of memory|oom-kill",
        "killed": r"\bkilled\b",
        "vasp_fatal": r"VERY BAD NEWS|ZBRENT: fatal error|internal error",
        "error": r"^\s*ERROR\b",
    }
    found: list[str] = []
    for marker, pattern in marker_patterns.items():
        if re.search(pattern, text, flags=re.IGNORECASE | re.MULTILINE):
            found.append(marker)
    return tuple(found)


def _first_float_match(text: str, pattern: str) -> float | None:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def _parse_time_pair(text: str) -> tuple[float, float] | None:
    values = [token for token in text.split() if token]
    if len(values) < 2:
        return None
    try:
        return float(values[0]), float(values[1])
    except ValueError:
        return None


def _resolve_contcar_path(
    calculation: VaspCalculationResult,
    *,
    contcar_root: str | Path | None = None,
) -> Path:
    direct = calculation.contcar_path
    if direct.is_file():
        return direct

    if contcar_root is not None:
        root = Path(contcar_root)
        candidate = root / calculation.system_name
        for value in calculation.calculation_chain:
            candidate /= str(value)
        candidate /= "CONTCAR"
        if candidate.is_file():
            return candidate

    raise FileNotFoundError(f"CONTCAR file not found: {direct}")


def read_contcar_as_crystal(
    contcar_path: str | Path,
    *,
    name: str | None = None,
    space_group: str = "P 1",
) -> CrystalStructure:
    path = Path(contcar_path)
    if not path.is_file():
        raise FileNotFoundError(f"CONTCAR file not found: {path}")

    try:
        from pymatgen.io.vasp.inputs import Poscar

        structure = Poscar.from_file(path).structure
    except Exception:
        return _read_contcar_as_crystal_fallback(
            path,
            name=name,
            space_group=space_group,
        )

    return _crystal_from_pymatgen_structure(
        structure=structure,
        name=name or path.parent.name,
        space_group=space_group,
    )


def summarize_system_rows(system: VaspSystemResult) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for calc in system.calculations:
        rows.append(
            {
                "system_name": calc.system_name,
                "calculation_index": calc.calculation_index,
                "calculation_chain": ".".join(str(value) for value in calc.calculation_chain) or "0",
                "status": calc.status,
                "ionic_iterations": calc.ionic_iterations,
                "total_electronic_iterations": calc.total_electronic_iterations,
                "last_electronic_iterations": calc.last_electronic_iterations,
                "final_energy_ev": calc.final_energy_ev,
                "converged": calc.converged,
                "converged_ionic": calc.converged_ionic,
                "converged_electronic": calc.converged_electronic,
                "cpu_time_seconds": calc.cpu_time_seconds,
                "parse_error": calc.parse_error,
            }
        )
    return rows


def _crystal_from_pymatgen_structure(
    *,
    structure,
    name: str,
    space_group: str,
) -> CrystalStructure:
    from .crystal_structure import AtomRecord, CrystalStructure

    abc = tuple(float(value) for value in structure.lattice.abc)
    angles = tuple(float(value) for value in structure.lattice.angles)
    cell_parameters = abc + angles
    lattice_matrix = tuple(
        tuple(float(value) for value in vector)
        for vector in structure.lattice.matrix
    )

    label_counts: dict[str, int] = {}
    atoms: list[AtomRecord] = []
    for site in structure:
        element = site.specie.symbol
        label_counts[element] = label_counts.get(element, 0) + 1
        atoms.append(
            AtomRecord(
                label=f"{element}{label_counts[element]}",
                element=element,
                coordinates=tuple(float(value) for value in site.coords),
            )
        )

    return CrystalStructure(
        atoms=atoms,
        cell_parameters=cell_parameters,
        lattice_matrix=lattice_matrix,
        space_group=space_group,
        name=name,
        explict_unit_cell=True,
    )


def _read_contcar_as_crystal_fallback(
    path: Path,
    *,
    name: str | None,
    space_group: str,
) -> CrystalStructure:
    from .crystal_structure import AtomRecord, CrystalStructure

    lines = [line.rstrip() for line in path.read_text(encoding="utf-8", errors="ignore").splitlines() if line.strip()]
    if len(lines) < 8:
        raise ValueError(f"Invalid CONTCAR/POSCAR format: {path}")

    scale = float(lines[1].split()[0])
    lattice = []
    for index in range(2, 5):
        vector = [float(value) for value in lines[index].split()[:3]]
        lattice.append([scale * value for value in vector])

    raw_elements = lines[5].split()
    elements = [_normalize_poscar_symbol(token) for token in raw_elements]
    counts = [int(value) for value in lines[6].split()]
    if len(elements) != len(counts):
        raise ValueError(f"Element/count mismatch in {path}")

    coord_line_index = 7
    selective_dynamics = lines[coord_line_index].strip().lower().startswith("s")
    if selective_dynamics:
        coord_line_index += 1

    mode = lines[coord_line_index].strip().lower()
    direct_coordinates = mode.startswith("d")
    if not direct_coordinates and not mode.startswith("c"):
        raise ValueError(f"Unknown coordinate mode in {path}: {lines[coord_line_index]}")

    coordinate_start = coord_line_index + 1
    lattice_matrix = np.array(lattice, dtype=float)
    label_counts: dict[str, int] = {}
    atoms: list[AtomRecord] = []
    line_index = coordinate_start

    for element, count in zip(elements, counts, strict=True):
        for _ in range(count):
            tokens = lines[line_index].split()
            coords = np.array([float(value) for value in tokens[:3]], dtype=float)
            if direct_coordinates:
                cart = coords @ lattice_matrix
            else:
                cart = coords * scale
            label_counts[element] = label_counts.get(element, 0) + 1
            atoms.append(
                AtomRecord(
                    label=f"{element}{label_counts[element]}",
                    element=element,
                    coordinates=tuple(float(value) for value in cart),
                )
            )
            line_index += 1

    a_vec, b_vec, c_vec = lattice_matrix
    a = float(np.linalg.norm(a_vec))
    b = float(np.linalg.norm(b_vec))
    c = float(np.linalg.norm(c_vec))
    alpha = _angle_degrees(b_vec, c_vec)
    beta = _angle_degrees(a_vec, c_vec)
    gamma = _angle_degrees(a_vec, b_vec)

    return CrystalStructure(
        atoms=atoms,
        cell_parameters=(a, b, c, alpha, beta, gamma),
        lattice_matrix=tuple(tuple(float(value) for value in row) for row in lattice_matrix),
        space_group=space_group,
        name=name or path.parent.name,
        explict_unit_cell=True,
    )


def _normalize_poscar_symbol(token: str) -> str:
    token = token.strip()
    if not token:
        raise ValueError("Empty POSCAR species token.")
    if len(token) == 1:
        return token.upper()
    return token[0].upper() + token[1:].lower()


def _angle_degrees(left: np.ndarray, right: np.ndarray) -> float:
    numerator = float(np.dot(left, right))
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    cosine = max(-1.0, min(1.0, numerator / denominator))
    return float(np.degrees(np.arccos(cosine)))
