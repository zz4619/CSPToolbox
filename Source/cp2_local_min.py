"""Prepare and audit CrystalPredictor2 experimental local-minimisation jobs.

The CP2 molecular-information (LAM) database is authoritative for site order,
site types, charges, and indexed Z-matrix references.  Its site names are not
unique, so canonical atom labels come from the system-level ``Zmatrix`` file.
This module maps an experimental SHELX structure into that order, writes the
fixed-column ``expcrys.pdb`` read by CP2, stages unchanged CP2 inputs, and
records complete provenance.

The adapter's supported pilot scope is single-component Z'=1.  That is not a
claim that a mapping or an optimizer-converged structure is scientifically
validated.  Z'>1/multicomponent preparation and positive local-minimum evidence
remain explicit validation tasks.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import shlex
import shutil
from typing import Iterable, Mapping, Sequence

import networkx as nx
import numpy as np
from ase.cell import Cell
from ase.data import atomic_numbers, covalent_radii
from networkx.algorithms import isomorphism as nx_isomorphism


GLOBAL_SEARCH_DIRECTORY_NAMES = ("5_Globalsearch", "5_GlobSrch")
DEFAULT_COVALENT_SCALE = 1.20
HYDROGEN_COVALENT_SCALE = 1.30
DEFAULT_MAX_HEAVY_MAPPINGS = 100_000
STATUS_FILENAME = "cp2_local_min_status.json"
MANIFEST_FILENAME = "cp2_local_min_manifest.json"
MAPPING_FILENAME = "cp2_atom_mapping.tsv"
PBS_SCRIPT_FILENAME = "run_cp2_local_min.pbs"
DEFAULT_CX3_MODULES = (
    "tools/prod",
    "imkl/2022.1.0",
)
CP2_SUPPORTED_SPACE_GROUPS = (
    "P1",
    "P-1",
    "P21",
    "C2",
    "Pc",
    "Cm",
    "Cc",
    "P21/m",
    "C2/m",
    "P2/c",
    "P2/n",
    "P21/c",
    "P21/a",
    "P21/n",
    "C2/c",
    "P2221",
    "P21212",
    "P212121",
    "C2221",
    "Pca21",
    "Pmn21",
    "Pba2",
    "Pna21",
    "Pn21a",
    "Cmc21",
    "Aba2",
    "Fdd2",
    "Iba2",
    "Ima2",
    "Pnna",
    "Pccn",
    "Pbcm",
    "Pnnm",
    "Pmmn",
    "Pbcn",
    "Pbca",
    "Pnma",
    "Cmcm",
    "Cmca",
    "Fddd",
    "Ibam",
    "P41",
    "P43",
    "I-4",
    "P4/n",
    "P42/n",
    "I4/m",
    "I41/a",
    "P41212",
    "P43212",
    "P-421c",
    "I-42d",
    "P31",
    "P32",
    "R3",
    "P-3",
    "R-3",
    "P3121",
    "P3221",
    "R3c",
    "R-3c",
    "P61",
    "P63",
    "P63/m",
    "P6122",
    "P213",
    "Pa-3",
)


@dataclass(frozen=True)
class CP2Torsion:
    label: str
    lower_degrees: float
    upper_degrees: float


@dataclass(frozen=True)
class CP2MolecularType:
    index: int
    name: str
    occurrences: int
    atom_count: int
    lam_reference: str
    torsions: tuple[CP2Torsion, ...]


@dataclass(frozen=True)
class CP2InputDefinition:
    path: Path
    molecular_types: tuple[CP2MolecularType, ...]
    space_groups: tuple[str, ...]
    non_uniform_enabled: bool

    @property
    def molecular_type_count(self) -> int:
        return len(self.molecular_types)

    @property
    def asymmetric_unit_molecule_count(self) -> int:
        return sum(item.occurrences for item in self.molecular_types)

    @property
    def supported_scope(self) -> str:
        if self.molecular_type_count == 1 and self.asymmetric_unit_molecule_count == 1:
            return "supported_single_component_zprime1"
        return "experimental_unvalidated_zprime_gt1_or_multicomponent"

    @property
    def structurally_supported(self) -> bool:
        return self.supported_scope == "supported_single_component_zprime1"


@dataclass(frozen=True)
class CP2CanonicalSite:
    index: int
    label: str
    element: str
    bond_to: int | None
    angle_to: int | None
    dihedral_to: int | None


@dataclass(frozen=True)
class CP2CanonicalZMatrix:
    path: Path
    sites: tuple[CP2CanonicalSite, ...]


@dataclass(frozen=True)
class CP2LamSite:
    index: int
    label: str
    element: str
    site_type: str
    symmetry_type: int
    charge: float
    bond_to: int | None
    angle_to: int | None
    dihedral_to: int | None


@dataclass(frozen=True)
class CP2LamTopology:
    path: Path
    sites: tuple[CP2LamSite, ...]
    canonical_zmatrix_path: Path | None
    labels_generated: bool

    @property
    def labels(self) -> tuple[str, ...]:
        return tuple(site.label for site in self.sites)

    @property
    def elements(self) -> tuple[str, ...]:
        return tuple(site.element for site in self.sites)


@dataclass(frozen=True)
class StructureAtom:
    source_index: int
    label: str
    element: str
    fractional: tuple[float, float, float]
    cartesian: tuple[float, float, float]


@dataclass(frozen=True)
class StructureData:
    path: Path
    atoms: tuple[StructureAtom, ...]
    cell_parameters: tuple[float, float, float, float, float, float]
    space_group: str

    @property
    def cell(self) -> Cell:
        return Cell.fromcellpar(self.cell_parameters)


@dataclass(frozen=True)
class MappingMetrics:
    method: str
    heavy_mapping_count: int
    candidates_truncated: bool
    gross_bond_angle_mismatches: int
    internal_coordinate_score: float
    heavy_kabsch_rmsd_angstrom: float
    second_best_score_gap: float | None
    torsion_orientation: str
    mapping_ambiguous: bool
    reference_order_assumed: bool

    @property
    def validation_safe(self) -> bool:
        """Whether automatic runnable-job validation may accept this mapping."""

        return (
            not self.candidates_truncated
            and not self.mapping_ambiguous
            and not self.reference_order_assumed
            and self.gross_bond_angle_mismatches == 0
            and self.torsion_orientation == "same"
        )


@dataclass(frozen=True)
class ComponentMapping:
    molecular_type_index: int
    occurrence_index: int
    experimental_component_index: int
    canonical_to_experimental: tuple[tuple[str, int], ...]
    metrics: MappingMetrics

    def as_index_mapping(self) -> dict[str, int]:
        return dict(self.canonical_to_experimental)


@dataclass(frozen=True)
class CP2LocalMinArtifacts:
    output_dir: Path
    expcrys_pdb_path: Path
    manifest_path: Path
    mapping_tsv_path: Path
    input_path: Path
    potential_path: Path
    staged_lam_paths: tuple[Path, ...]
    staged_zmatrix_paths: tuple[Path, ...]
    executable_path: Path | None
    pbs_script_path: Path | None


@dataclass(frozen=True)
class CP2PBSSettings:
    walltime: str = "24:00:00"
    memory_gb: int = 8
    modules: tuple[str, ...] = DEFAULT_CX3_MODULES
    # The branch binary links NAG statically and has no dynamic MPI dependency;
    # only the Kusari licence location is needed in addition to MKL at runtime.
    # Keep this command configurable for differently installed CP2 builds.
    nagvars_command: str = (
        "export NAG_KUSARI_FILE=/sw-eb/software/NAGlib/license/license.dat"
    )
    queue: str | None = None


@dataclass(frozen=True)
class CP2LocalMinStatus:
    job_dir: str
    status: str
    finished: bool
    optimizer_converged: bool
    local_minimum_confirmed: bool
    local_minimum_evidence: str
    log_path: str | None
    log_candidates: tuple[str, ...]
    final_structure_path: str | None
    final_structure_parseable: bool
    exit_code: int | None
    result_schema: str | None
    optimizer_status: str | None
    optimizer_info: int | None
    result_values: dict[str, float]
    final_ifail: int | None
    ifail_history: tuple[int, ...]
    minimisation_attempts: int
    cpu_time_seconds: float | None
    energies_kj_mol: dict[str, dict[str, float | None]]
    failure_markers: tuple[str, ...]


@dataclass(frozen=True)
class _Connectivity:
    graph: nx.Graph
    shifts: dict[tuple[int, int], tuple[int, int, int]]
    components: tuple[tuple[int, ...], ...]
    unwrapped: dict[int, np.ndarray]


@dataclass(frozen=True)
class _ReferenceContext:
    molecular_type: CP2MolecularType
    topology: CP2LamTopology
    structure: StructureData
    component: tuple[int, ...]
    canonical_to_reference: dict[str, int]
    graph: nx.Graph
    unwrapped: dict[int, np.ndarray]
    order_assumed: bool


def discover_global_search_dir(system_dir: str | Path) -> Path:
    """Locate the CP2 global-search directory without silently choosing aliases."""

    root = Path(system_dir)
    candidates = [
        root / name for name in GLOBAL_SEARCH_DIRECTORY_NAMES if (root / name).is_dir()
    ]
    if len(candidates) > 1:
        raise ValueError(
            "Ambiguous CP2 global-search directories: "
            + ", ".join(str(path) for path in candidates)
        )
    if candidates:
        return candidates[0]
    expected = ", ".join(GLOBAL_SEARCH_DIRECTORY_NAMES)
    raise FileNotFoundError(
        f"No CP2 global-search directory under {root}; expected one of {expected}"
    )


def parse_cp2_input(path: str | Path) -> CP2InputDefinition:
    """Parse CP2 molecular TYPE blocks from ``input.in``."""

    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    header_index = None
    expected_count = None
    for index, raw in enumerate(lines):
        match = re.match(r"^\s*MOLECULAR\s+TYPES\s+(\d+)\b", raw, flags=re.IGNORECASE)
        if match:
            header_index = index
            expected_count = int(match.group(1))
            break
    if header_index is None or expected_count is None:
        raise ValueError(f"No 'MOLECULAR TYPES' block in {source}")

    molecular_types: list[CP2MolecularType] = []
    cursor = header_index + 1
    while len(molecular_types) < expected_count:
        cursor = _next_matching_line(
            lines, cursor, lambda text: text.upper().startswith("TYPE ")
        )
        type_tokens = lines[cursor].strip().split(maxsplit=1)
        if len(type_tokens) != 2 or not type_tokens[1].strip():
            raise ValueError(f"Invalid TYPE line in {source}: {lines[cursor]!r}")
        name = type_tokens[1].strip()
        cursor += 1

        occurrences, cursor = _read_next_integer(
            lines, cursor, source, "molecular occurrence count"
        )
        atom_count, cursor = _read_next_integer(
            lines, cursor, source, "molecular atom count"
        )
        lam_line, cursor = _read_next_data_line(
            lines, cursor, source, "LAM database reference"
        )
        lam_tokens = shlex.split(lam_line, comments=True)
        if not lam_tokens:
            raise ValueError(f"Empty LAM database reference in {source}")
        lam_reference = lam_tokens[0]
        torsion_count, cursor = _read_next_integer(
            lines, cursor, source, "flexible torsion count"
        )

        torsions: list[CP2Torsion] = []
        for _ in range(torsion_count):
            torsion_line, cursor = _read_next_data_line(
                lines, cursor, source, "flexible torsion"
            )
            tokens = torsion_line.split()
            if len(tokens) < 3:
                raise ValueError(
                    f"Invalid flexible torsion row in {source}: {torsion_line!r}"
                )
            torsions.append(
                CP2Torsion(
                    label=tokens[0],
                    lower_degrees=float(tokens[1].replace("D", "E").replace("d", "e")),
                    upper_degrees=float(tokens[2].replace("D", "E").replace("d", "e")),
                )
            )

        if occurrences < 1 or atom_count < 1:
            raise ValueError(
                f"Non-positive occurrence or atom count for TYPE {name!r} in {source}"
            )
        molecular_types.append(
            CP2MolecularType(
                index=len(molecular_types) + 1,
                name=name,
                occurrences=occurrences,
                atom_count=atom_count,
                lam_reference=lam_reference,
                torsions=tuple(torsions),
            )
        )

    space_groups: list[str] = []
    in_space_groups = False
    for raw in lines:
        stripped = raw.strip()
        upper = stripped.upper()
        if upper == "SPACE_GROUPS":
            in_space_groups = True
            continue
        if upper == "END_SPACE_GROUPS":
            in_space_groups = False
            continue
        if in_space_groups and stripped and set(stripped) != {"-"}:
            space_groups.extend(stripped.split())
    if not space_groups:
        raise ValueError(f"No SPACE_GROUPS entries found in {source}")

    return CP2InputDefinition(
        path=source,
        molecular_types=tuple(molecular_types),
        space_groups=tuple(space_groups),
        non_uniform_enabled=any(
            raw.strip().upper().startswith("NON_UNIFORM") for raw in lines
        ),
    )


def parse_cp2_canonical_zmatrix(
    path: str | Path, atom_count: int,
) -> CP2CanonicalZMatrix:
    """Parse the unique atom labels and textual references in system ``Zmatrix``."""

    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    header_index = next(
        (
            index
            for index, raw in enumerate(lines)
            if raw.strip().lower().startswith("z-matrix for molecule")
        ),
        None,
    )
    if header_index is None:
        raise ValueError(f"No 'Z-matrix for molecule' header in {source}")
    rows: list[list[str]] = []
    cursor = header_index + 1
    while len(rows) < atom_count:
        if cursor >= len(lines):
            raise ValueError(
                f"Canonical Zmatrix {source} ended before {atom_count} rows"
            )
        tokens = lines[cursor].split()
        cursor += 1
        if not tokens:
            continue
        expected = 1 + min(len(rows), 3)
        if len(tokens) != expected:
            raise ValueError(
                f"Invalid canonical Zmatrix row {len(rows) + 1} in {source}: "
                f"expected {expected} fields, found {len(tokens)}"
            )
        rows.append(tokens)

    labels = [tokens[0] for tokens in rows]
    if len(labels) != len(set(labels)):
        raise ValueError(f"Canonical Zmatrix labels are not unique in {source}")
    label_to_index = {label: index for index, label in enumerate(labels, start=1)}
    sites: list[CP2CanonicalSite] = []
    for row_index, tokens in enumerate(rows, start=1):
        references: list[int] = []
        for label in tokens[1:]:
            reference = label_to_index.get(label)
            if reference is None or reference >= row_index:
                raise ValueError(
                    f"Canonical Zmatrix row {row_index} has invalid reference {label!r}"
                )
            references.append(reference)
        padded = references + [None] * (3 - len(references))
        sites.append(
            CP2CanonicalSite(
                index=row_index,
                label=tokens[0],
                element=_element_from_site_label(tokens[0]),
                bond_to=padded[0],
                angle_to=padded[1],
                dihedral_to=padded[2],
            )
        )
    return CP2CanonicalZMatrix(path=source, sites=tuple(sites))


def parse_cp2_lam_topology(
    path: str | Path,
    atom_count: int,
    *,
    canonical_zmatrix_path: str | Path | None = None,
    allow_generated_labels: bool = False,
) -> CP2LamTopology:
    """Read CP2 LAM site order and combine it with unique canonical labels."""

    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    from_index = next(
        (index for index, raw in enumerate(lines) if raw.lstrip().startswith("From")),
        None,
    )
    if from_index is None:
        raise ValueError(f"No 'From' Z-matrix section in CP2 LAM database {source}")

    raw_sites: list[dict[str, object]] = []
    cursor = from_index + 1
    while len(raw_sites) < atom_count:
        if cursor >= len(lines):
            raise ValueError(
                f"LAM database {source} ended before {atom_count} sites were read"
            )
        raw = lines[cursor].strip()
        cursor += 1
        if not raw:
            continue
        tokens = raw.split()
        row_index = len(raw_sites) + 1
        required = 3 + 2 * min(row_index - 1, 3)
        if len(tokens) < required:
            raise ValueError(
                f"Invalid CP2 LAM Z-matrix row {row_index} in {source}: {raw!r}"
            )
        site_type = tokens[0]
        element = _element_from_site_label(site_type)
        try:
            symmetry_type = int(tokens[1])
            charge = _fortran_float(tokens[2])
            references = [
                int(tokens[offset]) for offset in (3, 5, 7) if offset < required
            ]
        except ValueError as error:
            raise ValueError(
                f"Invalid numeric field in CP2 LAM Z-matrix row {row_index}: {raw!r}"
            ) from error
        for reference in references:
            if not 1 <= reference < row_index:
                raise ValueError(
                    f"CP2 LAM row {row_index} has forward/invalid reference {reference}"
                )
        padded = references + [None] * (3 - len(references))
        raw_sites.append(
            {
                "index": row_index,
                "element": element,
                "site_type": site_type,
                "symmetry_type": symmetry_type,
                "charge": charge,
                "bond_to": padded[0],
                "angle_to": padded[1],
                "dihedral_to": padded[2],
            }
        )

    canonical = None
    labels_generated = False
    if canonical_zmatrix_path is not None:
        canonical = parse_cp2_canonical_zmatrix(canonical_zmatrix_path, atom_count)
        labels = [site.label for site in canonical.sites]
        for raw_site, canonical_site in zip(raw_sites, canonical.sites):
            if raw_site["element"] != canonical_site.element:
                raise ValueError(
                    f"Element mismatch at CP2 site {canonical_site.index}: LAM "
                    f"{raw_site['site_type']} vs Zmatrix {canonical_site.label}"
                )
            lam_references = (
                raw_site["bond_to"],
                raw_site["angle_to"],
                raw_site["dihedral_to"],
            )
            canonical_references = (
                canonical_site.bond_to,
                canonical_site.angle_to,
                canonical_site.dihedral_to,
            )
            if lam_references != canonical_references:
                raise ValueError(
                    f"Connectivity mismatch at CP2 site {canonical_site.index} between "
                    f"LAM {source} and canonical Zmatrix {canonical.path}"
                )
    elif allow_generated_labels:
        counts: Counter[str] = Counter()
        labels = []
        for raw_site in raw_sites:
            element = str(raw_site["element"])
            counts[element] += 1
            labels.append(f"{element}{counts[element]}")
        labels_generated = True
    else:
        raise ValueError(
            f"CP2 LAM site types are not unique atom labels; provide the system Zmatrix "
            f"for {source}, or explicitly allow generated per-element labels"
        )

    sites = tuple(
        CP2LamSite(label=label, **raw_site)
        for label, raw_site in zip(labels, raw_sites)
    )
    return CP2LamTopology(
        path=source,
        sites=sites,
        canonical_zmatrix_path=canonical.path if canonical is not None else None,
        labels_generated=labels_generated,
    )


def read_structure(
    path: str | Path,
    *,
    space_group: str | None = None,
    require_space_group: bool = True,
) -> StructureData:
    source = Path(path)
    suffix = source.suffix.lower()
    if suffix == ".res":
        return _read_res_structure(
            source, space_group=space_group, require_space_group=require_space_group,
        )
    if suffix == ".pdb":
        return _read_pdb_structure(
            source, space_group=space_group, require_space_group=require_space_group,
        )
    raise ValueError(f"Unsupported structure format for CP2 mapping: {source}")


def match_experimental_to_cp2(
    experimental: StructureData,
    input_definition: CP2InputDefinition,
    topologies: Mapping[int, CP2LamTopology],
    reference_paths: Mapping[int, str | Path] | None = None,
    *,
    max_heavy_mappings: int = DEFAULT_MAX_HEAVY_MAPPINGS,
) -> tuple[ComponentMapping, ...]:
    """Map experimental components into CP2 type/occurrence/site order."""

    if max_heavy_mappings < 1:
        raise ValueError("max_heavy_mappings must be positive")
    valid_type_indices = {item.index for item in input_definition.molecular_types}
    unknown_reference_indices = set(reference_paths or {}) - valid_type_indices
    if unknown_reference_indices:
        raise ValueError(
            "Reference paths supplied for unknown molecular TYPE indices: "
            + ", ".join(str(value) for value in sorted(unknown_reference_indices))
        )

    experimental_connectivity = _build_connectivity(experimental)
    components = experimental_connectivity.components
    expected_components = input_definition.asymmetric_unit_molecule_count
    if len(components) != expected_components:
        raise ValueError(
            f"Detected {len(components)} experimental molecular components; "
            f"CP2 input expects {expected_components}"
        )

    resolved_reference_paths = {
        int(key): Path(value) for key, value in (reference_paths or {}).items()
    }
    contexts: dict[int, _ReferenceContext] = {}
    for molecular_type in input_definition.molecular_types:
        topology = topologies[molecular_type.index]
        reference_path = resolved_reference_paths.get(molecular_type.index)
        if reference_path is None:
            contexts[molecular_type.index] = _identity_reference_context(
                molecular_type, topology, experimental, experimental_connectivity,
            )
        else:
            contexts[molecular_type.index] = _build_reference_context(
                molecular_type,
                topology,
                read_structure(reference_path, require_space_group=False),
            )

    candidates: dict[tuple[int, int], tuple[dict[str, int], MappingMetrics]] = {}
    for molecular_type in input_definition.molecular_types:
        context = contexts[molecular_type.index]
        for component_index, component in enumerate(components, start=1):
            if _component_formula(experimental, component) != _topology_formula(
                context.topology
            ):
                continue
            try:
                candidates[(molecular_type.index, component_index)] = _match_component(
                    context,
                    experimental,
                    experimental_connectivity,
                    component,
                    max_heavy_mappings=max_heavy_mappings,
                )
            except ValueError:
                continue

    assignment = _assign_components(input_definition, components, candidates)
    mappings: list[ComponentMapping] = []
    for molecular_type in input_definition.molecular_types:
        selected_components = sorted(
            assignment[molecular_type.index],
            key=lambda component_index: min(components[component_index - 1]),
        )
        for occurrence_index, component_index in enumerate(
            selected_components, start=1
        ):
            mapping, metrics = candidates[(molecular_type.index, component_index)]
            mappings.append(
                ComponentMapping(
                    molecular_type_index=molecular_type.index,
                    occurrence_index=occurrence_index,
                    experimental_component_index=component_index,
                    canonical_to_experimental=tuple(
                        (site.label, mapping[site.label])
                        for site in topologies[molecular_type.index].sites
                    ),
                    metrics=metrics,
                )
            )
    return tuple(mappings)


def prepare_cp2_local_min_inputs(
    system_dir: str | Path,
    experimental_structure: str | Path,
    output_dir: str | Path,
    *,
    reference_paths: Mapping[int, str | Path] | None = None,
    zmatrix_paths: Mapping[int, str | Path] | None = None,
    global_search_dir: str | Path | None = None,
    stage_mode: str = "symlink",
    space_group: str | None = None,
    system_name: str | None = None,
    refcode: str | None = None,
    compack_metadata: Mapping[str, object] | None = None,
    max_heavy_mappings: int = DEFAULT_MAX_HEAVY_MAPPINGS,
    allow_generated_labels: bool = False,
    allow_unvalidated: bool = False,
    allow_unvalidated_mapping: bool = False,
    cp2_executable: str | Path | None = None,
    pbs_settings: CP2PBSSettings | None = None,
) -> CP2LocalMinArtifacts:
    """Prepare a complete CP2 local-minimisation input directory."""

    if stage_mode not in {"symlink", "copy"}:
        raise ValueError("stage_mode must be 'symlink' or 'copy'")
    system_root = Path(system_dir)
    global_root = (
        Path(global_search_dir)
        if global_search_dir is not None
        else discover_global_search_dir(system_root)
    )
    input_source = global_root / "input.in"
    potential_source = global_root / "potential.in"
    for required in (input_source, potential_source):
        if not required.is_file():
            raise FileNotFoundError(f"Missing required CP2 input: {required}")

    input_definition = parse_cp2_input(input_source)
    if not input_definition.structurally_supported and not allow_unvalidated:
        raise ValueError(
            f"CP2 input scope is {input_definition.supported_scope}; Z'>1 and "
            "multicomponent preparation requires explicit allow_unvalidated=True"
        )
    executable_source = Path(cp2_executable) if cp2_executable is not None else None
    if executable_source is not None and not executable_source.is_file():
        raise FileNotFoundError(
            f"CP2 Minimise executable not found: {executable_source}"
        )

    resolved_zmatrix_paths = {
        int(key): Path(value) for key, value in (zmatrix_paths or {}).items()
    }
    valid_type_indices = {item.index for item in input_definition.molecular_types}
    unknown_zmatrix_indices = set(resolved_zmatrix_paths) - valid_type_indices
    if unknown_zmatrix_indices:
        raise ValueError(
            "Zmatrix paths supplied for unknown molecular TYPE indices: "
            + ", ".join(str(value) for value in sorted(unknown_zmatrix_indices))
        )
    if input_definition.molecular_type_count == 1 and 1 not in resolved_zmatrix_paths:
        zmatrix_candidates = [
            path
            for path in (system_root / "Zmatrix", global_root / "Zmatrix")
            if path.is_file()
        ]
        if len(zmatrix_candidates) == 2 and _sha256(zmatrix_candidates[0]) != _sha256(
            zmatrix_candidates[1]
        ):
            raise ValueError(
                "Ambiguous non-identical canonical Zmatrix files: "
                + ", ".join(str(path) for path in zmatrix_candidates)
            )
        if zmatrix_candidates:
            resolved_zmatrix_paths[1] = zmatrix_candidates[0]

    topologies: dict[int, CP2LamTopology] = {}
    lam_sources: dict[int, Path] = {}
    for molecular_type in input_definition.molecular_types:
        lam_source = _resolve_lam_source(global_root, molecular_type.lam_reference)
        if not lam_source.is_file():
            raise FileNotFoundError(
                f"Missing CP2 LAM database for TYPE {molecular_type.name!r}: {lam_source}"
            )
        lam_sources[molecular_type.index] = lam_source
        topologies[molecular_type.index] = parse_cp2_lam_topology(
            lam_source,
            molecular_type.atom_count,
            canonical_zmatrix_path=resolved_zmatrix_paths.get(molecular_type.index),
            allow_generated_labels=allow_generated_labels,
        )

    if input_definition.non_uniform_enabled:
        non_uniform_source = global_root / "NON_UNIFORM_LAM_RELEVANCE"
        if not non_uniform_source.is_file():
            raise FileNotFoundError(
                f"input.in enables NON_UNIFORM but {non_uniform_source} is missing"
            )
    else:
        non_uniform_source = None

    experimental = read_structure(experimental_structure, space_group=space_group)
    mappings = match_experimental_to_cp2(
        experimental,
        input_definition,
        topologies,
        reference_paths,
        max_heavy_mappings=max_heavy_mappings,
    )
    supplied_reference_indices = {int(index) for index in (reference_paths or {})}
    all_types_referenced = valid_type_indices.issubset(supplied_reference_indices)
    mapping_validated = (
        input_definition.structurally_supported
        and all_types_referenced
        and all(not topology.labels_generated for topology in topologies.values())
        and all(item.metrics.validation_safe for item in mappings)
    )
    if (
        executable_source is not None
        and not mapping_validated
        and not allow_unvalidated_mapping
    ):
        raise ValueError(
            "Refusing to generate a runnable PBS job for an unvalidated atom mapping; "
            "inspect a preparation without --cp2-executable or explicitly use "
            "allow_unvalidated_mapping=True"
        )
    connectivity = _build_connectivity(experimental)

    destination = Path(output_dir)
    if destination.exists() and not destination.is_dir():
        raise FileExistsError(
            f"CP2 local-min output path is not a directory: {destination}"
        )
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError(
            f"Refusing to overwrite non-empty CP2 local-min directory {destination}"
        )
    destination.mkdir(parents=True, exist_ok=True)

    staged_input = destination / "input.in"
    staged_potential = destination / "potential.in"
    _stage_file(input_source, staged_input, stage_mode)
    _stage_file(potential_source, staged_potential, stage_mode)
    if _sha256(potential_source) != _sha256(staged_potential):
        raise AssertionError(
            "Staged potential.in differs byte-for-byte from the global-search source"
        )

    staged_lam_paths: list[Path] = []
    staged_lam_by_source: dict[Path, Path] = {}
    for molecular_type in input_definition.molecular_types:
        source = lam_sources[molecular_type.index]
        if source in staged_lam_by_source:
            continue
        reference = Path(molecular_type.lam_reference)
        if reference.is_absolute():
            staged = source
        else:
            if ".." in reference.parts:
                raise ValueError(
                    f"Unsafe relative LAM reference {molecular_type.lam_reference!r}; "
                    "use an input-local path"
                )
            staged = destination / reference
            _stage_file(source, staged, stage_mode)
        staged_lam_by_source[source] = staged
        staged_lam_paths.append(staged)

    staged_zmatrix_paths: list[Path] = []
    staged_zmatrix_by_type: dict[int, Path] = {}
    for type_index, source in sorted(resolved_zmatrix_paths.items()):
        target = destination / (
            "Zmatrix" if type_index == 1 else f"Zmatrix_type_{type_index}"
        )
        _stage_file(source, target, stage_mode)
        staged_zmatrix_paths.append(target)
        staged_zmatrix_by_type[type_index] = target

    optional_assets: list[dict[str, object]] = []
    staged_non_uniform = None
    if non_uniform_source is not None:
        staged_non_uniform = destination / non_uniform_source.name
        _stage_file(non_uniform_source, staged_non_uniform, stage_mode)
        optional_assets.append(_file_provenance(non_uniform_source, staged_non_uniform))

    expcrys_path = destination / "expcrys.pdb"
    cp2_space_group = resolve_cp2_space_group(experimental.space_group)
    mapping_rows = _write_expcrys_pdb(
        expcrys_path,
        experimental,
        input_definition,
        topologies,
        mappings,
        connectivity,
        cp2_space_group=cp2_space_group,
    )
    mapping_path = destination / MAPPING_FILENAME
    _write_tsv(mapping_path, mapping_rows)

    reference_manifest = {
        str(index): _path_provenance(Path(path))
        for index, path in sorted((reference_paths or {}).items())
    }
    runtime_paths = [staged_input, staged_potential, expcrys_path]
    runtime_paths.extend(
        path for path in staged_lam_paths if path.is_relative_to(destination)
    )
    if staged_non_uniform is not None:
        runtime_paths.append(staged_non_uniform)
    staged_executable = None
    pbs_script_path = None
    effective_pbs_settings = pbs_settings or CP2PBSSettings()
    if executable_source is not None:
        staged_executable = destination / "cp2_Minimise"
        _stage_file(executable_source, staged_executable, "copy")
        staged_executable.chmod(executable_source.stat().st_mode & 0o777)
        if _sha256(executable_source) != _sha256(staged_executable):
            raise AssertionError("Staged CP2 executable is not byte-identical")
        pbs_script_path = destination / PBS_SCRIPT_FILENAME
        build_cp2_pbs_script(
            pbs_script_path,
            destination,
            staged_executable,
            runtime_paths,
            settings=effective_pbs_settings,
            job_name=f"cp2lm_{refcode or Path(experimental_structure).stem}",
        )

    manifest = {
        "schema_version": 1,
        "system_name": system_name or system_root.name,
        "refcode": refcode or Path(experimental_structure).stem,
        "supported_scope": input_definition.supported_scope,
        "structurally_supported": input_definition.structurally_supported,
        "mapping_validated": mapping_validated,
        "scientifically_validated": False,
        "validation_note": "Optimizer convergence and positive local-minimum evidence require separate validation.",
        "system_dir": str(system_root.resolve()),
        "global_search_dir": str(global_root.resolve()),
        "experimental_structure": _path_provenance(Path(experimental_structure)),
        "reference_structures": reference_manifest,
        "experimental_space_group_input": experimental.space_group,
        "cp2_space_group": cp2_space_group,
        "input": _file_provenance(input_source, staged_input),
        "potential": {
            **_file_provenance(potential_source, staged_potential),
            "byte_identical": True,
        },
        "molecular_types": [
            {
                **asdict(item),
                "torsions": [asdict(torsion) for torsion in item.torsions],
                "lam": _file_provenance(
                    lam_sources[item.index],
                    staged_lam_by_source[lam_sources[item.index]],
                ),
                "canonical_zmatrix": (
                    _file_provenance(
                        resolved_zmatrix_paths[item.index],
                        staged_zmatrix_by_type[item.index],
                    )
                    if item.index in resolved_zmatrix_paths
                    else None
                ),
                "canonical_labels_generated": topologies[item.index].labels_generated,
                "canonical_atom_order": list(topologies[item.index].labels),
                "authoritative_lam_site_types": [
                    site.site_type for site in topologies[item.index].sites
                ],
            }
            for item in input_definition.molecular_types
        ],
        "mappings": [
            {
                "molecular_type_index": item.molecular_type_index,
                "occurrence_index": item.occurrence_index,
                "experimental_component_index": item.experimental_component_index,
                "metrics": asdict(item.metrics),
            }
            for item in mappings
        ],
        "atom_mapping_tsv": str(mapping_path.resolve()),
        "expcrys_pdb": _path_provenance(expcrys_path),
        "optional_assets": optional_assets,
        "stage_mode": stage_mode,
        "allow_unvalidated": allow_unvalidated,
        "allow_unvalidated_mapping": allow_unvalidated_mapping,
        "allow_generated_labels": allow_generated_labels,
        "compack_metadata": dict(compack_metadata or {}),
        "execution": (
            {
                "cp2_executable": _path_provenance(executable_source),
                "staged_cp2_executable": _file_provenance(
                    executable_source, staged_executable,
                ),
                "pbs_script": _path_provenance(pbs_script_path),
                "pbs_settings": asdict(effective_pbs_settings),
                "submission_performed": False,
            }
            if executable_source is not None and pbs_script_path is not None
            else None
        ),
        "checks": {
            "component_count_matches_cp2": True,
            "atom_mapping_bijective": True,
            "atom_elements_match": True,
            "potential_bytes_preserved": True,
            "pdb_atom_count": len(mapping_rows),
        },
    }
    manifest_path = destination / MANIFEST_FILENAME
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    return CP2LocalMinArtifacts(
        output_dir=destination,
        expcrys_pdb_path=expcrys_path,
        manifest_path=manifest_path,
        mapping_tsv_path=mapping_path,
        input_path=staged_input,
        potential_path=staged_potential,
        staged_lam_paths=tuple(staged_lam_paths),
        staged_zmatrix_paths=tuple(staged_zmatrix_paths),
        executable_path=staged_executable,
        pbs_script_path=pbs_script_path,
    )


def collect_cp2_local_min_status(
    job_dir: str | Path, *, log_path: str | Path | None = None, write_json: bool = True,
) -> CP2LocalMinStatus:
    """Parse CP2 minimisation outputs and optionally write a status JSON file."""

    root = Path(job_dir)
    log_candidates = _find_log_candidates(root)
    selected_log = (
        Path(log_path)
        if log_path is not None
        else (log_candidates[0] if log_candidates else None)
    )
    text = ""
    if selected_log is not None and selected_log.is_file():
        text = selected_log.read_text(encoding="utf-8", errors="replace")

    ifail_history: list[int] = []
    for line in text.splitlines():
        if "IFAIL" not in line.upper():
            continue
        integers = re.findall(r"[-+]?\d+", line)
        if integers:
            ifail_history.append(int(integers[-1]))
    final_ifail = ifail_history[-1] if ifail_history else None
    (
        result_schema,
        optimizer_status,
        optimizer_info,
        result_values,
    ) = _parse_result_marker(text)

    exit_code = None
    exit_status_path = root / "exit_status.txt"
    if exit_status_path.is_file():
        tokens = exit_status_path.read_text(encoding="utf-8", errors="replace").split()
        if tokens:
            try:
                exit_code = int(tokens[0])
            except ValueError:
                exit_code = None

    final_candidates = [
        root / "minimise_final_structure.pdb",
        root / "minimize_final_structure.pdb",
    ]
    final_structure = next(
        (
            path
            for path in final_candidates
            if path.is_file() and path.stat().st_size > 0
        ),
        None,
    )
    final_structure_parseable = (
        _is_parseable_cp2_pdb(final_structure) if final_structure is not None else False
    )
    failure_markers = tuple(
        marker
        for marker in (
            "segmentation fault",
            "segmentation violation",
            "terminated abnormally",
            "energy failed",
            "forrtl: severe",
        )
        if marker in text.lower()
    )
    completion_marker = "MINIMISATION CPU TIME" in text.upper()

    effective_info = optimizer_info if optimizer_info is not None else final_ifail
    if result_schema is not None:
        if (
            optimizer_status in {"TORSION_BOUND", "OPTIMIZER_FAILED"}
            or (effective_info is not None and effective_info != 0)
            or (exit_code is not None and exit_code != 0)
        ):
            status = "failed"
        elif (
            optimizer_status == "OPTIMIZER_CONVERGED"
            and effective_info == 0
            and final_structure_parseable
        ):
            status = "optimizer_converged"
        else:
            status = "completed_unverified"
    elif final_ifail is not None and final_ifail != 0:
        status = "failed"
    elif exit_code is not None and exit_code != 0:
        status = "failed"
    elif failure_markers:
        status = "failed"
    elif final_ifail == 0 and final_structure_parseable and exit_code in {None, 0}:
        status = "completed_unverified"
    elif final_structure_parseable and exit_code == 0 and completion_marker:
        status = "completed_unverified"
    elif not text and exit_code is None and final_structure is None:
        status = "not_started"
    else:
        status = "incomplete"

    attempts = [
        int(value)
        for value in re.findall(
            r"Performing\s+Minimisation\s+attempt\s+(\d+)", text, flags=re.IGNORECASE
        )
    ]
    cpu_matches = re.findall(
        r"MINIMISATION\s+CPU\s+TIME:\s*([-+0-9.EeDd]+)\s*sec",
        text,
        flags=re.IGNORECASE,
    )
    cpu_time = _fortran_float(cpu_matches[-1]) if cpu_matches else None
    energies = _parse_energy_report(text)
    for result_key, energy_key in (
        ("utot_kj_mol_entity", "Utot"),
        ("uvdw_kj_mol_entity", "Uvdw"),
        ("uelec_kj_mol_entity", "Uelec"),
        ("uintra_kj_mol_entity", "Uintra"),
    ):
        if result_key in result_values:
            energies.setdefault(energy_key, {"initial": None, "final": None})[
                "final"
            ] = result_values[result_key]

    result = CP2LocalMinStatus(
        job_dir=str(root.resolve()),
        status=status,
        finished=status in {"optimizer_converged", "failed", "completed_unverified"},
        optimizer_converged=status == "optimizer_converged",
        local_minimum_confirmed=False,
        local_minimum_evidence="not_evaluated_hessian_or_perturbation",
        log_path=str(selected_log.resolve())
        if selected_log is not None and selected_log.exists()
        else None,
        log_candidates=tuple(str(path.resolve()) for path in log_candidates),
        final_structure_path=str(final_structure.resolve())
        if final_structure is not None
        else None,
        final_structure_parseable=final_structure_parseable,
        exit_code=exit_code,
        result_schema=result_schema,
        optimizer_status=optimizer_status,
        optimizer_info=optimizer_info,
        result_values=result_values,
        final_ifail=final_ifail,
        ifail_history=tuple(ifail_history),
        minimisation_attempts=max(attempts, default=0),
        cpu_time_seconds=cpu_time,
        energies_kj_mol=energies,
        failure_markers=failure_markers,
    )
    if write_json:
        (root / STATUS_FILENAME).write_text(
            json.dumps(asdict(result), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return result


def _read_res_structure(
    path: Path, *, space_group: str | None, require_space_group: bool,
) -> StructureData:
    lines = (
        path.read_text(encoding="utf-8", errors="replace")
        .replace("\r", "")
        .splitlines()
    )
    sfac: list[str] = []
    cell_parameters = None
    parsed_space_group = space_group
    atoms_raw: list[tuple[str, str, tuple[float, float, float]]] = []
    cursor = 0
    while cursor < len(lines):
        raw = lines[cursor]
        tokens = raw.split()
        cursor += 1
        if not tokens:
            continue
        keyword = tokens[0].upper()
        if keyword == "END":
            break
        if keyword == "HKLF":
            break
        if keyword == "PART":
            try:
                part_number = int(tokens[1]) if len(tokens) >= 2 else 0
            except ValueError as error:
                raise ValueError(f"Invalid PART record in {path}: {raw!r}") from error
            if part_number != 0:
                raise ValueError(
                    f"Disordered SHELX PART {part_number} is unsupported for CP2 mapping in {path}"
                )
            continue
        if keyword == "CELL":
            if len(tokens) < 8:
                raise ValueError(f"Invalid CELL record in {path}: {raw!r}")
            cell_parameters = tuple(_fortran_float(value) for value in tokens[2:8])
            continue
        if keyword == "SFAC":
            sfac = [_normalize_element(value) for value in tokens[1:]]
            continue
        if keyword == "REM" and len(tokens) >= 3 and tokens[1].upper() == "SPACE_GROUP":
            parsed_space_group = " ".join(tokens[2:])
            continue
        if keyword == "TITL" and parsed_space_group is None:
            match = re.search(r"\bin\s+([^\s]+)\s*$", raw, flags=re.IGNORECASE)
            if match:
                parsed_space_group = match.group(1)
            continue
        atom = _parse_res_atom_tokens(tokens, sfac)
        if atom is None:
            continue
        label, element, fractional = atom
        continuation = raw.rstrip()
        while continuation.endswith("="):
            if cursor >= len(lines):
                raise ValueError(
                    f"Unterminated SHELX continuation for atom {label} in {path}"
                )
            continuation = lines[cursor].rstrip()
            cursor += 1
            continuation_tokens = continuation.split()
            numeric = (
                continuation_tokens[:-1]
                if continuation_tokens and continuation_tokens[-1] == "="
                else continuation_tokens
            )
            try:
                for value in numeric:
                    _fortran_float(value)
            except ValueError as error:
                raise ValueError(
                    f"Non-numeric SHELX continuation for atom {label} in {path}: {continuation!r}"
                ) from error
        atoms_raw.append((label, element, fractional))

    if cell_parameters is None:
        raise ValueError(f"No CELL record found in {path}")
    if not atoms_raw:
        raise ValueError(f"No SHELX atom records found in {path}")
    if parsed_space_group is None:
        has_symmetry = any(
            raw.split() and raw.split()[0].upper() == "SYMM" for raw in lines
        )
        latt_values = [
            raw.split()[1]
            for raw in lines
            if len(raw.split()) >= 2 and raw.split()[0].upper() == "LATT"
        ]
        latt = latt_values[0] if len(latt_values) == 1 else None
        if not has_symmetry and latt == "-1":
            parsed_space_group = "P 1"
        elif not has_symmetry and latt == "1":
            parsed_space_group = "P -1"
        elif require_space_group:
            raise ValueError(
                f"Could not determine a unique space-group symbol from SHELX LATT/SYMM records "
                f"in {path}; pass space_group explicitly"
            )
        else:
            parsed_space_group = "UNKNOWN"

    cell = Cell.fromcellpar(cell_parameters)
    atoms = tuple(
        StructureAtom(
            source_index=index,
            label=label,
            element=element,
            fractional=tuple(float(value) for value in fractional),
            cartesian=tuple(
                float(value) for value in np.dot(np.asarray(fractional), cell.array)
            ),
        )
        for index, (label, element, fractional) in enumerate(atoms_raw)
    )
    labels = [atom.label for atom in atoms]
    if len(labels) != len(set(labels)):
        raise ValueError(f"Duplicate atom labels in {path}")
    return StructureData(
        path=path,
        atoms=atoms,
        cell_parameters=tuple(float(value) for value in cell_parameters),
        space_group=parsed_space_group,
    )


def _read_pdb_structure(
    path: Path, *, space_group: str | None, require_space_group: bool,
) -> StructureData:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    cryst1 = next((line for line in lines if line.startswith("CRYST1")), None)
    if cryst1 is None:
        raise ValueError(f"No CRYST1 record found in {path}")
    cell_parameters = (
        float(cryst1[6:15]),
        float(cryst1[15:24]),
        float(cryst1[24:33]),
        float(cryst1[33:40]),
        float(cryst1[40:47]),
        float(cryst1[47:54]),
    )
    parsed_space_group = space_group or cryst1[55:66].strip()
    if not parsed_space_group:
        if require_space_group:
            raise ValueError(f"No space-group symbol in CRYST1 record in {path}")
        parsed_space_group = "UNKNOWN"
    cell = Cell.fromcellpar(cell_parameters)
    atoms: list[StructureAtom] = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        label = line[12:16].strip() or f"ATOM{len(atoms) + 1}"
        raw_element = line[76:78].strip() if len(line) >= 78 else ""
        element = _normalize_element(raw_element or _element_from_site_label(label))
        cartesian = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        fractional_array = cell.scaled_positions(np.asarray([cartesian], dtype=float))[
            0
        ]
        atoms.append(
            StructureAtom(
                source_index=len(atoms),
                label=label,
                element=element,
                fractional=tuple(float(value) for value in fractional_array),
                cartesian=cartesian,
            )
        )
    if not atoms:
        raise ValueError(f"No PDB atoms found in {path}")
    return StructureData(
        path=path,
        atoms=tuple(atoms),
        cell_parameters=cell_parameters,
        space_group=parsed_space_group,
    )


def _parse_res_atom_tokens(
    tokens: Sequence[str], sfac: Sequence[str],
) -> tuple[str, str, tuple[float, float, float]] | None:
    if len(tokens) < 5 or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_'-]*", tokens[0]):
        return None
    if re.fullmatch(r"Q\d+", tokens[0], flags=re.IGNORECASE):
        return None
    try:
        fractional = tuple(_fortran_float(value) for value in tokens[2:5])
    except ValueError:
        return None
    try:
        sfac_index = int(tokens[1])
    except ValueError:
        try:
            element = _normalize_element(tokens[1])
        except ValueError:
            return None
    else:
        if not 1 <= sfac_index <= len(sfac):
            return None
        element = sfac[sfac_index - 1]
    return tokens[0], element, fractional  # type: ignore[return-value]


def _build_connectivity(structure: StructureData) -> _Connectivity:
    graph = nx.Graph()
    for atom in structure.atoms:
        graph.add_node(atom.source_index, element=atom.element)
    shifts: dict[tuple[int, int], tuple[int, int, int]] = {}
    cell = np.asarray(structure.cell.array, dtype=float)
    fractional = np.asarray([atom.fractional for atom in structure.atoms], dtype=float)
    radii = [
        float(covalent_radii[atomic_numbers[atom.element]]) for atom in structure.atoms
    ]

    candidate_edges: dict[tuple[int, int], tuple[float, tuple[int, int, int]]] = {}
    for left in range(len(structure.atoms)):
        for right in range(left + 1, len(structure.atoms)):
            delta = fractional[right] - fractional[left]
            base_shift = -np.rint(delta).astype(int)
            best_distance = float("inf")
            best_shift = (0, 0, 0)
            for offset in itertools.product((-1, 0, 1), repeat=3):
                shift_array = base_shift + np.asarray(offset, dtype=int)
                displacement = np.dot(delta + shift_array, cell)
                distance = float(np.linalg.norm(displacement))
                if distance < best_distance:
                    best_distance = distance
                    best_shift = tuple(int(value) for value in shift_array)
            scale = (
                HYDROGEN_COVALENT_SCALE
                if "H"
                in {structure.atoms[left].element, structure.atoms[right].element}
                else DEFAULT_COVALENT_SCALE
            )
            if best_distance <= scale * (radii[left] + radii[right]):
                candidate_edges[(left, right)] = (best_distance, best_shift)

    allowed = set(candidate_edges)
    for atom in structure.atoms:
        if atom.element != "H":
            continue
        incident = [edge for edge in allowed if atom.source_index in edge]
        if len(incident) > 1:
            keep = min(incident, key=lambda edge: (candidate_edges[edge][0], edge))
            allowed.difference_update(edge for edge in incident if edge != keep)

    for left, right in sorted(allowed):
        distance, shift = candidate_edges[(left, right)]
        graph.add_edge(left, right, distance=distance)
        shifts[(left, right)] = shift
        shifts[(right, left)] = tuple(-value for value in shift)

    components = tuple(
        tuple(sorted(component)) for component in nx.connected_components(graph)
    )
    components = tuple(sorted(components, key=lambda component: min(component)))
    unwrapped: dict[int, np.ndarray] = {}
    cartesian = np.asarray([atom.cartesian for atom in structure.atoms], dtype=float)
    for component in components:
        root = min(component)
        unwrapped[root] = cartesian[root].copy()
        stack = [root]
        while stack:
            current = stack.pop()
            for neighbor in graph.neighbors(current):
                if neighbor not in component or neighbor in unwrapped:
                    continue
                displacement = (
                    cartesian[neighbor]
                    + np.dot(shifts[(current, neighbor)], cell)
                    - cartesian[current]
                )
                unwrapped[neighbor] = unwrapped[current] + displacement
                stack.append(neighbor)
    return _Connectivity(
        graph=graph, shifts=shifts, components=components, unwrapped=unwrapped
    )


def _build_reference_context(
    molecular_type: CP2MolecularType,
    topology: CP2LamTopology,
    structure: StructureData,
) -> _ReferenceContext:
    connectivity = _build_connectivity(structure)
    formula = _topology_formula(topology)
    matching_components = [
        component
        for component in connectivity.components
        if _component_formula(structure, component) == formula
    ]
    if not matching_components:
        raise ValueError(
            f"Reference {structure.path} has no component matching TYPE {molecular_type.name!r}"
        )
    by_label = {atom.label: atom.source_index for atom in structure.atoms}
    labelled = [
        component
        for component in matching_components
        if set(topology.labels).issubset(
            {structure.atoms[index].label for index in component}
        )
    ]
    order_assumed = False
    if len(labelled) == 1:
        component = labelled[0]
        canonical_to_reference = {label: by_label[label] for label in topology.labels}
    elif len(matching_components) == 1:
        component = matching_components[0]
        ordered = sorted(component)
        observed_elements = tuple(structure.atoms[index].element for index in ordered)
        if observed_elements != topology.elements:
            raise ValueError(
                f"Reference labels do not match CP2 sites and atom order is incompatible in {structure.path}"
            )
        canonical_to_reference = dict(zip(topology.labels, ordered))
        order_assumed = True
    else:
        raise ValueError(
            f"Reference {structure.path} has multiple indistinguishable components for TYPE {molecular_type.name!r}"
        )
    return _ReferenceContext(
        molecular_type=molecular_type,
        topology=topology,
        structure=structure,
        component=component,
        canonical_to_reference=canonical_to_reference,
        graph=connectivity.graph.subgraph(component).copy(),
        unwrapped=connectivity.unwrapped,
        order_assumed=order_assumed,
    )


def _identity_reference_context(
    molecular_type: CP2MolecularType,
    topology: CP2LamTopology,
    experimental: StructureData,
    connectivity: _Connectivity,
) -> _ReferenceContext:
    labels = {atom.label: atom.source_index for atom in experimental.atoms}
    if not set(topology.labels).issubset(labels):
        raise ValueError(
            f"TYPE {molecular_type.name!r} requires a reference structure because experimental labels "
            "do not already match the authoritative CP2 LAM labels"
        )
    component_candidates = [
        component
        for component in connectivity.components
        if set(topology.labels).issubset(
            {experimental.atoms[index].label for index in component}
        )
    ]
    if len(component_candidates) != 1:
        raise ValueError(
            f"Could not identify a unique label-matched component for TYPE {molecular_type.name!r}"
        )
    component = component_candidates[0]
    return _ReferenceContext(
        molecular_type=molecular_type,
        topology=topology,
        structure=experimental,
        component=component,
        canonical_to_reference={label: labels[label] for label in topology.labels},
        graph=connectivity.graph.subgraph(component).copy(),
        unwrapped=connectivity.unwrapped,
        order_assumed=False,
    )


def _match_component(
    reference: _ReferenceContext,
    experimental: StructureData,
    experimental_connectivity: _Connectivity,
    component: tuple[int, ...],
    *,
    max_heavy_mappings: int,
) -> tuple[dict[str, int], MappingMetrics]:
    experimental_graph = experimental_connectivity.graph.subgraph(component).copy()
    reference_heavy = reference.graph.subgraph(
        [
            node
            for node in reference.graph
            if reference.structure.atoms[node].element != "H"
        ]
    ).copy()
    experimental_heavy = experimental_graph.subgraph(
        [node for node in experimental_graph if experimental.atoms[node].element != "H"]
    ).copy()
    matcher = nx_isomorphism.GraphMatcher(
        reference_heavy,
        experimental_heavy,
        node_match=lambda left, right: left["element"] == right["element"],
    )
    if not matcher.is_isomorphic():
        raise ValueError("Heavy-atom graphs are not isomorphic")

    candidates: list[tuple[tuple[object, ...], dict[str, int], dict[str, object]]] = []
    mapping_count = 0
    truncated = False
    for heavy_mapping in matcher.isomorphisms_iter():
        mapping_count += 1
        if mapping_count > max_heavy_mappings:
            truncated = True
            break
        try:
            full_reference_mapping, heavy_rmsd = _assign_hydrogens(
                reference, experimental, experimental_connectivity, heavy_mapping,
            )
        except ValueError:
            continue
        canonical_mapping = {
            canonical_label: full_reference_mapping[reference_index]
            for canonical_label, reference_index in reference.canonical_to_reference.items()
        }
        score = _score_zmatrix_mapping(
            reference, experimental, experimental_connectivity, canonical_mapping,
        )
        deterministic = tuple(
            canonical_mapping[site.label] for site in reference.topology.sites
        )
        rank = (
            score["gross_bond_angle_mismatches"],
            score["internal_coordinate_score"],
            heavy_rmsd,
            deterministic,
        )
        candidates.append(
            (rank, canonical_mapping, {**score, "heavy_rmsd": heavy_rmsd})
        )
    if not candidates:
        raise ValueError("No complete graph-preserving atom mapping was generated")
    candidates.sort(key=lambda item: item[0])
    best_rank, best_mapping, best_score = candidates[0]
    second_gap = None
    mapping_ambiguous = False
    if len(candidates) > 1:
        second_gap = float(candidates[1][2]["internal_coordinate_score"]) - float(
            best_score["internal_coordinate_score"]
        )
        second_score = candidates[1][2]
        mapping_ambiguous = (
            int(second_score["gross_bond_angle_mismatches"])
            == int(best_score["gross_bond_angle_mismatches"])
            and abs(second_gap) <= 1.0e-8
            and abs(float(second_score["heavy_rmsd"]) - float(best_score["heavy_rmsd"]))
            <= 1.0e-6
        )
    metrics = MappingMetrics(
        method=(
            "trusted_experimental_labels"
            if reference.structure.path == experimental.path
            else "reference_graph_zmatrix_internal_rmsd"
        ),
        heavy_mapping_count=min(mapping_count, max_heavy_mappings),
        candidates_truncated=truncated,
        gross_bond_angle_mismatches=int(best_score["gross_bond_angle_mismatches"]),
        internal_coordinate_score=float(best_score["internal_coordinate_score"]),
        heavy_kabsch_rmsd_angstrom=float(best_score["heavy_rmsd"]),
        second_best_score_gap=second_gap,
        torsion_orientation=str(best_score["torsion_orientation"]),
        mapping_ambiguous=mapping_ambiguous,
        reference_order_assumed=reference.order_assumed,
    )
    return best_mapping, metrics


def _assign_hydrogens(
    reference: _ReferenceContext,
    experimental: StructureData,
    experimental_connectivity: _Connectivity,
    heavy_mapping: Mapping[int, int],
) -> tuple[dict[int, int], float]:
    reference_order = sorted(heavy_mapping)
    experimental_order = [heavy_mapping[index] for index in reference_order]
    reference_coords = np.asarray(
        [reference.unwrapped[index] for index in reference_order], dtype=float
    )
    experimental_coords = np.asarray(
        [experimental_connectivity.unwrapped[index] for index in experimental_order],
        dtype=float,
    )
    rotation, mobile_center, target_center, heavy_rmsd = _kabsch_transform(
        experimental_coords, reference_coords,
    )

    result = dict(heavy_mapping)
    for reference_heavy_atom, experimental_heavy_atom in heavy_mapping.items():
        reference_hydrogens = sorted(
            node
            for node in reference.graph.neighbors(reference_heavy_atom)
            if reference.structure.atoms[node].element == "H"
        )
        experimental_hydrogens = sorted(
            node
            for node in experimental_connectivity.graph.neighbors(
                experimental_heavy_atom
            )
            if experimental.atoms[node].element == "H"
        )
        if len(reference_hydrogens) != len(experimental_hydrogens):
            raise ValueError("Hydrogen count differs for mapped heavy atoms")
        if not reference_hydrogens:
            continue
        best_permutation = min(
            itertools.permutations(experimental_hydrogens),
            key=lambda permutation: (
                sum(
                    float(
                        np.linalg.norm(
                            _apply_kabsch(
                                experimental_connectivity.unwrapped[experimental_index],
                                rotation,
                                mobile_center,
                                target_center,
                            )
                            - reference.unwrapped[reference_index]
                        )
                    )
                    ** 2
                    for reference_index, experimental_index in zip(
                        reference_hydrogens, permutation
                    )
                ),
                permutation,
            ),
        )
        result.update(zip(reference_hydrogens, best_permutation))
    if len(result) != len(reference.component) or len(set(result.values())) != len(
        result
    ):
        raise ValueError("Incomplete or non-bijective full-atom mapping")
    return result, heavy_rmsd


def _score_zmatrix_mapping(
    reference: _ReferenceContext,
    experimental: StructureData,
    experimental_connectivity: _Connectivity,
    mapping: Mapping[str, int],
) -> dict[str, object]:
    reference_values = _internal_values(
        reference.topology, reference.canonical_to_reference, reference.unwrapped,
    )
    experimental_values = _internal_values(
        reference.topology, mapping, experimental_connectivity.unwrapped,
    )
    same_cost = _torsion_cost(reference_values, experimental_values, 1)
    inverted_cost = _torsion_cost(reference_values, experimental_values, -1)
    orientation = -1 if inverted_cost + 1.0e-10 < same_cost else 1

    weighted: list[float] = []
    bond_deltas: list[float] = []
    angle_deltas: list[float] = []
    for reference_value, experimental_value in zip(
        reference_values, experimental_values
    ):
        if reference_value[0] is not None:
            delta = float(experimental_value[0]) - float(reference_value[0])
            bond_deltas.append(delta)
            weighted.append((delta / 0.05) ** 2)
        if reference_value[1] is not None:
            delta = float(experimental_value[1]) - float(reference_value[1])
            angle_deltas.append(delta)
            weighted.append((delta / 5.0) ** 2)
        if reference_value[2] is not None:
            delta = _circular_difference(
                float(experimental_value[2]), orientation * float(reference_value[2]),
            )
            weighted.append((delta / 15.0) ** 2)
    gross = sum(abs(value) > 0.25 for value in bond_deltas)
    gross += sum(abs(value) > 15.0 for value in angle_deltas)
    return {
        "gross_bond_angle_mismatches": gross,
        "internal_coordinate_score": math.sqrt(sum(weighted) / max(1, len(weighted))),
        "torsion_orientation": "same" if orientation == 1 else "globally_inverted",
    }


def _internal_values(
    topology: CP2LamTopology,
    mapping: Mapping[str, int],
    coordinates: Mapping[int, np.ndarray],
) -> list[tuple[float | None, float | None, float | None]]:
    values = []
    for site in topology.sites:
        atom = coordinates[mapping[site.label]]
        bond = angle = dihedral = None
        if site.bond_to is not None:
            bond_atom = coordinates[mapping[topology.sites[site.bond_to - 1].label]]
            bond = float(np.linalg.norm(atom - bond_atom))
        if site.angle_to is not None and site.bond_to is not None:
            angle_atom = coordinates[mapping[topology.sites[site.angle_to - 1].label]]
            angle = _angle_degrees(atom, bond_atom, angle_atom)
        if (
            site.dihedral_to is not None
            and site.angle_to is not None
            and site.bond_to is not None
        ):
            dihedral_atom = coordinates[
                mapping[topology.sites[site.dihedral_to - 1].label]
            ]
            dihedral = _dihedral_degrees(atom, bond_atom, angle_atom, dihedral_atom)
        values.append((bond, angle, dihedral))
    return values


def _assign_components(
    input_definition: CP2InputDefinition,
    components: Sequence[tuple[int, ...]],
    candidates: Mapping[tuple[int, int], tuple[dict[str, int], MappingMetrics]],
) -> dict[int, tuple[int, ...]]:
    best: tuple[tuple[object, ...], dict[int, tuple[int, ...]]] | None = None

    def search(
        type_position: int, available: set[int], selected: dict[int, tuple[int, ...]]
    ) -> None:
        nonlocal best
        if type_position == len(input_definition.molecular_types):
            if available:
                return
            metrics = [
                candidates[(type_index, component_index)][1]
                for type_index, component_indices in selected.items()
                for component_index in component_indices
            ]
            rank = (
                sum(item.gross_bond_angle_mismatches for item in metrics),
                sum(item.internal_coordinate_score for item in metrics),
                sum(item.heavy_kabsch_rmsd_angstrom for item in metrics),
                tuple((key, selected[key]) for key in sorted(selected)),
            )
            if best is None or rank < best[0]:
                best = rank, dict(selected)
            return
        molecular_type = input_definition.molecular_types[type_position]
        compatible = sorted(
            component_index
            for component_index in available
            if (molecular_type.index, component_index) in candidates
        )
        for chosen in itertools.combinations(compatible, molecular_type.occurrences):
            selected[molecular_type.index] = tuple(chosen)
            search(type_position + 1, available - set(chosen), selected)
            del selected[molecular_type.index]

    search(0, set(range(1, len(components) + 1)), {})
    if best is None:
        raise ValueError(
            "Could not assign every experimental component to the CP2 molecular TYPE blocks"
        )
    return best[1]


def _write_expcrys_pdb(
    path: Path,
    experimental: StructureData,
    input_definition: CP2InputDefinition,
    topologies: Mapping[int, CP2LamTopology],
    mappings: Sequence[ComponentMapping],
    connectivity: _Connectivity,
    *,
    cp2_space_group: str,
) -> list[dict[str, object]]:
    a, b, c, alpha, beta, gamma = experimental.cell_parameters
    if any(value >= 1000.0 for value in (a, b, c)):
        raise ValueError(
            "CP2 expcrys.pdb CRYST1 field supports cell lengths below 1000 A"
        )
    if len(cp2_space_group) > 11:
        raise ValueError(
            f"Space-group symbol is too long for CP2 PDB field: {cp2_space_group!r}"
        )
    lines = [
        f"HEADER    CP2 EXPERIMENTAL LOCAL MIN {path.parent.name}",
        f"CRYST1{a:9.4f}{b:9.4f}{c:9.4f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} "
        f"{cp2_space_group:<11}",
    ]
    rows: list[dict[str, object]] = []
    serial = 0
    mapping_by_slot = {
        (item.molecular_type_index, item.occurrence_index): item for item in mappings
    }
    seen_experimental: set[int] = set()
    for molecular_type in input_definition.molecular_types:
        topology = topologies[molecular_type.index]
        for occurrence_index in range(1, molecular_type.occurrences + 1):
            component_mapping = mapping_by_slot[
                (molecular_type.index, occurrence_index)
            ]
            canonical_to_experimental = component_mapping.as_index_mapping()
            for site in topology.sites:
                serial += 1
                experimental_index = canonical_to_experimental[site.label]
                if experimental_index in seen_experimental:
                    raise ValueError("Experimental atom mapping is not bijective")
                seen_experimental.add(experimental_index)
                source_atom = experimental.atoms[experimental_index]
                if source_atom.element != site.element:
                    raise ValueError(
                        f"Element mismatch for CP2 site {site.label}: {site.element} vs {source_atom.element}"
                    )
                coordinate = connectivity.unwrapped[experimental_index]
                if any(abs(float(value)) >= 1000.0 for value in coordinate):
                    raise ValueError(
                        "CP2 expcrys.pdb coordinate exceeds the fixed f8.3 field"
                    )
                x, y, z = (round(float(value), 3) for value in coordinate)
                atom_name = site.label[:4]
                chain = chr(ord("A") + ((molecular_type.index - 1) % 26))
                lines.append(
                    f"HETATM{serial:5d} {atom_name:<4s} MOL {chain}{occurrence_index:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {site.element:>2s}"
                )
                rows.append(
                    {
                        "molecular_type_index": molecular_type.index,
                        "molecular_type_name": molecular_type.name,
                        "occurrence_index": occurrence_index,
                        "experimental_component_index": component_mapping.experimental_component_index,
                        "cp2_site_index_1based": site.index,
                        "cp2_label": site.label,
                        "element": site.element,
                        "experimental_atom_index_0based": experimental_index,
                        "experimental_label": source_atom.label,
                        "original_fractional_x": f"{source_atom.fractional[0]:.10f}",
                        "original_fractional_y": f"{source_atom.fractional[1]:.10f}",
                        "original_fractional_z": f"{source_atom.fractional[2]:.10f}",
                        "pdb_serial": serial,
                        "pdb_x_angstrom": f"{x:.3f}",
                        "pdb_y_angstrom": f"{y:.3f}",
                        "pdb_z_angstrom": f"{z:.3f}",
                        "mapping_method": component_mapping.metrics.method,
                    }
                )
    if len(seen_experimental) != len(experimental.atoms):
        raise ValueError(
            f"Mapped {len(seen_experimental)} of {len(experimental.atoms)} experimental atoms"
        )
    lines.extend(["END", ""])
    path.write_text("\n".join(lines), encoding="utf-8")
    _validate_written_pdb(path, len(rows))
    return rows


def _validate_written_pdb(path: Path, expected_atoms: int) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    if len(atom_lines) != expected_atoms:
        raise AssertionError(
            f"PDB atom count mismatch: {len(atom_lines)} vs {expected_atoms}"
        )
    for serial, line in enumerate(atom_lines, start=1):
        if int(line[6:11]) != serial:
            raise AssertionError("PDB serial numbers are not contiguous")
        for field in (line[30:38], line[38:46], line[46:54]):
            float(field)


def _normalize_space_group_key(symbol: str) -> str:
    return re.sub(r"[\s()_]", "", symbol).upper()


def resolve_cp2_space_group(experimental_symbol: str) -> str:
    key = _normalize_space_group_key(experimental_symbol)
    matches = [
        symbol
        for symbol in CP2_SUPPORTED_SPACE_GROUPS
        if _normalize_space_group_key(symbol) == key
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Experimental space group {experimental_symbol!r} does not resolve uniquely "
            "to CrystalPredictor2 SpaceSupported spelling"
        )
    return matches[0]


def _is_parseable_cp2_pdb(path: Path) -> bool:
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        cryst1 = next(line for line in lines if line.startswith("CRYST1"))
        for field in (
            cryst1[6:15],
            cryst1[15:24],
            cryst1[24:33],
            cryst1[33:40],
            cryst1[40:47],
            cryst1[47:54],
        ):
            float(field)
        atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
        if not atom_lines:
            return False
        for line in atom_lines:
            float(line[30:38])
            float(line[38:46])
            float(line[46:54])
    except (OSError, StopIteration, ValueError):
        return False
    return True


def _find_log_candidates(root: Path) -> list[Path]:
    ordered: list[Path] = []
    for name in (
        "Minimisation_log.out",
        "minimisation.output",
        "cp2_local_min.stdout",
        "minimise.out",
    ):
        path = root / name
        if path.is_file():
            ordered.append(path)
    for pattern in ("*_minimisation.output", "*.minimisation.output"):
        for path in sorted(root.glob(pattern)):
            if path not in ordered and path.is_file():
                ordered.append(path)
    return ordered


def _parse_energy_report(text: str) -> dict[str, dict[str, float | None]]:
    result: dict[str, dict[str, float | None]] = {}
    for line in text.splitlines():
        match = re.match(
            r"^\s*(Utot|Uvdw|Uelec|Ureal|Urec|Umol_cor|Usurf|Upre|Uintra)\s*=\s*(.*)$",
            line,
            flags=re.IGNORECASE,
        )
        if not match:
            continue
        values = [
            _fortran_float(value)
            for value in re.findall(
                r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?", match.group(2)
            )
        ]
        if not values:
            continue
        key = match.group(1)
        if len(values) >= 2:
            result[key] = {"initial": values[0], "final": values[1]}
        elif key in result:
            if result[key]["initial"] is None:
                result[key]["initial"] = result[key]["final"]
            result[key]["final"] = values[0]
        else:
            result[key] = {"initial": None, "final": values[0]}
    return result


def _parse_result_marker(
    text: str,
) -> tuple[str | None, str | None, int | None, dict[str, float]]:
    """Parse the final structured result emitted by the CP2 local-min executable."""

    marker_lines = [
        line.strip()
        for line in text.splitlines()
        if line.lstrip().startswith("CP2_LOCAL_MIN_RESULT_V1")
    ]
    if not marker_lines:
        return None, None, None, {}
    fields = dict(re.findall(r"([A-Za-z0-9_]+)=\s*([^\s]+)", marker_lines[-1]))
    optimizer_status = fields.pop("status", None)
    raw_info = fields.pop("info", None)
    optimizer_info = None
    if raw_info is not None:
        try:
            optimizer_info = int(raw_info)
        except ValueError:
            optimizer_info = None
    values: dict[str, float] = {}
    for key, value in fields.items():
        try:
            values[key] = _fortran_float(value)
        except ValueError:
            continue
    return "CP2_LOCAL_MIN_RESULT_V1", optimizer_status, optimizer_info, values


def _resolve_lam_source(global_root: Path, reference: str) -> Path:
    path = Path(reference)
    return path if path.is_absolute() else global_root / path


def build_cp2_pbs_script(
    output_path: str | Path,
    job_dir: str | Path,
    cp2_executable: str | Path,
    runtime_paths: Sequence[str | Path],
    *,
    settings: CP2PBSSettings | None = None,
    job_name: str = "cp2_local_min",
) -> Path:
    """Write a portable single-core CX3 PBS script without submitting it."""

    destination = Path(output_path)
    root = Path(job_dir).absolute()
    executable = Path(cp2_executable)
    effective = settings or CP2PBSSettings()
    if not executable.is_file():
        raise FileNotFoundError(f"CP2 Minimise executable not found: {executable}")
    try:
        executable_relative = executable.absolute().relative_to(root)
    except ValueError as error:
        raise ValueError(
            "CP2 executable must be staged inside the prepared job directory"
        ) from error
    if effective.memory_gb < 1:
        raise ValueError("PBS memory_gb must be positive")
    if not re.fullmatch(r"\d+:[0-5]\d:[0-5]\d", effective.walltime):
        raise ValueError(f"Invalid PBS walltime: {effective.walltime!r}")
    for value in (*effective.modules, effective.nagvars_command, effective.queue or ""):
        if "\n" in value or "\r" in value:
            raise ValueError("PBS settings may not contain newlines")

    relative_paths: list[Path] = []
    for raw_path in runtime_paths:
        path = Path(raw_path).absolute()
        try:
            relative = path.relative_to(root)
        except ValueError as error:
            raise ValueError(
                f"Runtime input is outside prepared job directory: {path}"
            ) from error
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe runtime input path: {relative}")
        relative_paths.append(relative)
    relative_paths = sorted(set(relative_paths), key=str)

    safe_job_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", job_name).strip("_.-")
    safe_job_name = (safe_job_name or "cp2_local_min")[:50]
    directives = [
        "#!/bin/bash",
        f"#PBS -N {safe_job_name}",
        f"#PBS -l walltime={effective.walltime}",
        f"#PBS -l select=1:ncpus=1:mem={effective.memory_gb}gb",
        "#PBS -j oe",
    ]
    if effective.queue:
        directives.append(f"#PBS -q {effective.queue}")

    module_lines = ["module purge"]
    # CX3 compute nodes can expose a stale Lmod cache even when the same
    # production module is visible on the login node.  Lmod itself recommends
    # bypassing that cache in this situation; doing so is deterministic and
    # avoids a setup-only job failure before CP2 starts.
    module_lines.extend(
        f"module --ignore_cache load {shlex.quote(module)}"
        for module in effective.modules
    )
    if effective.nagvars_command.strip():
        module_lines.append(effective.nagvars_command.strip())
    module_lines.extend(
        (
            "export OMP_NUM_THREADS=1",
            "export OMP_DYNAMIC=FALSE",
            "export MKL_NUM_THREADS=1",
            "export MKL_DYNAMIC=FALSE",
        )
    )

    copy_lines: list[str] = []
    for relative in relative_paths:
        quoted_relative = shlex.quote(str(relative))
        if relative.parent != Path("."):
            copy_lines.append(
                f'mkdir -p "$RUN_DIR"/{shlex.quote(str(relative.parent))}'
            )
        copy_lines.append(
            f'cp -L "$JOB_SOURCE"/{quoted_relative} "$RUN_DIR"/{quoted_relative}'
        )
    input_roots = sorted({path.parts[0] for path in relative_paths})
    input_root_array = " ".join(shlex.quote(value) for value in input_roots)

    lines = [
        *directives,
        "",
        "set -euo pipefail",
        "umask 077",
        *module_lines,
        "",
        'JOB_SOURCE="${PBS_O_WORKDIR:?Submit qsub from the prepared job directory}"',
        'if [[ ! -f "$JOB_SOURCE/cp2_local_min_manifest.json" ]]; then',
        '  echo "PBS_O_WORKDIR is not a prepared CP2 job; cd into the bundle before qsub" >&2',
        "  exit 2",
        "fi",
        'RUN_PARENT="${TMPDIR:-/tmp/${USER}}"',
        'RUN_DIR="${RUN_PARENT}/cp2_local_min_${PBS_JOBID:?PBS_JOBID is not set}"',
        'if [[ -e "$RUN_DIR" ]]; then',
        '  echo "Refusing to reuse existing run directory: $RUN_DIR" >&2',
        "  exit 2",
        "fi",
        'mkdir -p "$RUN_DIR"',
        *copy_lines,
        f'cp -L "$JOB_SOURCE"/{shlex.quote(str(executable_relative))} "$RUN_DIR/Minimise"',
        'chmod u+x "$RUN_DIR/Minimise"',
        'cd "$RUN_DIR"',
        "set +e",
        "./Minimise > cp2_local_min.stdout 2>&1",
        "cp2_status=$?",
        "set -e",
        "printf '%s\\n' \"$cp2_status\" > exit_status.txt",
        "shopt -s nullglob dotglob",
        f"input_roots=({input_root_array} Minimise)",
        'for item in "$RUN_DIR"/*; do',
        '  base="${item##*/}"',
        "  skip=0",
        '  for input_root in "${input_roots[@]}"; do',
        '    if [[ "$base" == "$input_root" ]]; then skip=1; break; fi',
        "  done",
        '  if [[ "$skip" -eq 1 ]]; then continue; fi',
        '  cp -a "$item" "$JOB_SOURCE"/',
        "done",
        'cd "$JOB_SOURCE"',
        'case "$RUN_DIR" in "$RUN_PARENT"/cp2_local_min_*) rm -rf -- "$RUN_DIR" ;; *) exit 3 ;; esac',
        'exit "$cp2_status"',
        "",
    ]
    destination.write_text("\n".join(lines), encoding="utf-8")
    destination.chmod(0o750)
    return destination


def _stage_file(source: Path, destination: Path, mode: str) -> None:
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to overwrite staged file {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if mode == "copy":
        shutil.copyfile(source, destination)
    else:
        destination.symlink_to(source.resolve())


def _file_provenance(source: Path, staged: Path) -> dict[str, object]:
    return {
        "source_path": str(source.resolve()),
        "staged_path": str(staged.absolute()),
        "staged_is_symlink": staged.is_symlink(),
        "sha256": _sha256(source),
        "size_bytes": source.stat().st_size,
    }


def _path_provenance(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": _sha256(path),
        "size_bytes": path.stat().st_size,
    }


def _write_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty CP2 atom-mapping table")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _next_matching_line(lines: Sequence[str], start: int, predicate) -> int:
    for index in range(start, len(lines)):
        if predicate(lines[index].strip()):
            return index
    raise ValueError("Unexpected end of CP2 input while looking for a TYPE block")


def _read_next_data_line(
    lines: Sequence[str], start: int, source: Path, description: str,
) -> tuple[str, int]:
    for index in range(start, len(lines)):
        stripped = lines[index].strip()
        if not stripped or set(stripped) <= {"-"} or stripped.startswith("#"):
            continue
        return stripped, index + 1
    raise ValueError(f"Unexpected end of {source} while reading {description}")


def _read_next_integer(
    lines: Sequence[str], start: int, source: Path, description: str,
) -> tuple[int, int]:
    line, cursor = _read_next_data_line(lines, start, source, description)
    try:
        return int(line.split()[0]), cursor
    except ValueError as error:
        raise ValueError(f"Invalid {description} in {source}: {line!r}") from error


def _element_from_site_label(label: str) -> str:
    match = re.match(r"^([A-Za-z]+)", label.strip())
    if match is None:
        raise ValueError(f"Could not infer element from CP2 site label {label!r}")
    return _normalize_element(match.group(1))


def _normalize_element(value: str) -> str:
    raw = re.sub(r"[^A-Za-z]", "", value)
    if raw.upper() in {"D", "T"}:
        return "H"
    candidates = (raw[:1].upper() + raw[1:].lower(), raw.upper(), raw[:1].upper())
    for candidate in candidates:
        if candidate in atomic_numbers:
            return candidate
    raise ValueError(f"Unsupported element symbol {value!r}")


def _fortran_float(value: str) -> float:
    return float(str(value).replace("D", "E").replace("d", "e"))


def _component_formula(
    structure: StructureData, component: Iterable[int]
) -> tuple[tuple[str, int], ...]:
    return tuple(
        sorted(Counter(structure.atoms[index].element for index in component).items())
    )


def _topology_formula(topology: CP2LamTopology) -> tuple[tuple[str, int], ...]:
    return tuple(sorted(Counter(topology.elements).items()))


def _kabsch_transform(
    mobile: np.ndarray, target: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    centered_mobile = mobile - mobile_center
    centered_target = target - target_center
    covariance = centered_mobile.T @ centered_target
    left, _singular, right_transpose = np.linalg.svd(covariance)
    rotation = left @ right_transpose
    if np.linalg.det(rotation) < 0.0:
        left[:, -1] *= -1.0
        rotation = left @ right_transpose
    aligned = centered_mobile @ rotation + target_center
    rmsd = math.sqrt(float(np.mean(np.sum((aligned - target) ** 2, axis=1))))
    return rotation, mobile_center, target_center, rmsd


def _apply_kabsch(
    coordinate: np.ndarray,
    rotation: np.ndarray,
    mobile_center: np.ndarray,
    target_center: np.ndarray,
) -> np.ndarray:
    return (np.asarray(coordinate) - mobile_center) @ rotation + target_center


def _angle_degrees(first: np.ndarray, center: np.ndarray, last: np.ndarray) -> float:
    left = first - center
    right = last - center
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator <= 1.0e-12:
        raise ValueError("Undefined bond angle in atom mapping")
    cosine = max(-1.0, min(1.0, float(np.dot(left, right) / denominator)))
    return math.degrees(math.acos(cosine))


def _dihedral_degrees(
    first: np.ndarray, second: np.ndarray, third: np.ndarray, fourth: np.ndarray,
) -> float:
    b0 = -(second - first)
    b1 = third - second
    b2 = fourth - third
    norm = float(np.linalg.norm(b1))
    if norm <= 1.0e-12:
        raise ValueError("Undefined dihedral in atom mapping")
    b1 = b1 / norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(
        math.atan2(float(np.dot(np.cross(b1, v), w)), float(np.dot(v, w)))
    )


def _circular_difference(left: float, right: float) -> float:
    return ((left - right + 180.0) % 360.0) - 180.0


def _torsion_cost(
    reference: Sequence[tuple[float | None, float | None, float | None]],
    experimental: Sequence[tuple[float | None, float | None, float | None]],
    sign: int,
) -> float:
    return sum(
        _circular_difference(float(exp[2]), sign * float(ref[2])) ** 2
        for ref, exp in zip(reference, experimental)
        if ref[2] is not None and exp[2] is not None
    )


__all__ = [
    "CP2_SUPPORTED_SPACE_GROUPS",
    "CP2CanonicalSite",
    "CP2CanonicalZMatrix",
    "CP2InputDefinition",
    "CP2LamSite",
    "CP2LamTopology",
    "CP2LocalMinArtifacts",
    "CP2LocalMinStatus",
    "CP2MolecularType",
    "CP2PBSSettings",
    "CP2Torsion",
    "ComponentMapping",
    "MappingMetrics",
    "StructureAtom",
    "StructureData",
    "build_cp2_pbs_script",
    "collect_cp2_local_min_status",
    "discover_global_search_dir",
    "match_experimental_to_cp2",
    "parse_cp2_canonical_zmatrix",
    "parse_cp2_input",
    "parse_cp2_lam_topology",
    "prepare_cp2_local_min_inputs",
    "read_structure",
    "resolve_cp2_space_group",
]
