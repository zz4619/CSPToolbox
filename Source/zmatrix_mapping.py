"""Engine-neutral atom mapping against an authoritative Z-matrix template.

The reference coordinates may come from any structure source, but the local-
minimisation workflow uses a CrystalPredictor2 global-search structure.  The
matcher deliberately knows nothing about CP2 or CSO-FM LAM file formats.  It
selects one graph-preserving atom permutation from all available Z-matrix
coordinates, including independent torsions, and uses fixed-correspondence,
all-atom RMSD only to break numerical ties.

LAM-domain checks still belong to the consuming energy model.  The shared
template-coordinate score keeps the atom permutation identical when the
resulting mapping is later used to prepare either a CP2 PDB or a CSO-FM RES
file, while preventing symmetry-related independent torsions from being
ignored during atom assignment.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import networkx as nx
import numpy as np
from networkx.algorithms import isomorphism as nx_isomorphism


MAPPING_METHOD_VERSION = (
    "cp2_global_reference_zmatrix_all_torsion_fixed_pair_rmsd_v4"
)
DEFAULT_MAX_HEAVY_MAPPINGS = 100_000
DEFAULT_MAX_FULL_MAPPINGS = 500_000


@dataclass(frozen=True, order=True)
class InternalCoordinateKey:
    """One coordinate in a Z-matrix row."""

    kind: str
    site_index: int

    def __post_init__(self) -> None:
        if self.kind not in {"bond", "angle", "dihedral"}:
            raise ValueError(f"Unknown internal-coordinate kind {self.kind!r}")
        if self.site_index < 1:
            raise ValueError("Z-matrix site indices are one-based")


@dataclass(frozen=True)
class ZMatrixSite:
    """Engine-neutral definition of one authoritative Z-matrix row."""

    index: int
    label: str
    element: str
    bond_to: int | None = None
    angle_to: int | None = None
    dihedral_to: int | None = None


@dataclass(frozen=True)
class InternalCoordinateValues:
    """Values measured for one Z-matrix row."""

    site_index: int
    label: str
    bond_angstrom: float | None
    angle_degrees: float | None
    dihedral_degrees: float | None


@dataclass(frozen=True)
class MappingScoreSettings:
    """Numerical scales and tie tolerances for the common matcher."""

    bond_scale_angstrom: float = 0.05
    angle_scale_degrees: float = 5.0
    torsion_scale_degrees: float = 15.0
    gross_bond_delta_angstrom: float = 0.25
    gross_angle_delta_degrees: float = 15.0
    primary_score_atol: float = 1.0e-8
    fixed_pair_rmsd_atol_angstrom: float = 1.0e-6


@dataclass(frozen=True)
class MappingCandidate:
    """A complete canonical-to-experimental atom permutation and its metrics."""

    canonical_to_experimental: tuple[tuple[str, int], ...]
    gross_bond_angle_mismatches: int
    internal_coordinate_score: float
    fixed_pair_all_atom_rmsd_angstrom: float
    heavy_atom_rmsd_angstrom: float
    bond_rms_delta_angstrom: float
    bond_max_abs_delta_angstrom: float
    angle_rms_delta_degrees: float
    angle_max_abs_delta_degrees: float
    rigid_torsion_rms_delta_degrees: float
    rigid_torsion_max_abs_delta_degrees: float
    independent_torsion_rms_delta_degrees: float
    independent_torsion_max_abs_delta_degrees: float
    reference_values: tuple[InternalCoordinateValues, ...]
    experimental_values: tuple[InternalCoordinateValues, ...]

    @property
    def deterministic_key(self) -> tuple[int, ...]:
        return tuple(index for _label, index in self.canonical_to_experimental)

    def as_index_mapping(self) -> dict[str, int]:
        return dict(self.canonical_to_experimental)


@dataclass(frozen=True)
class MappingDecision:
    """Auditable result of the staged candidate-selection process."""

    method_version: str
    selected: MappingCandidate
    heavy_mapping_count: int
    full_mapping_count: int
    candidates_truncated: bool
    primary_tie_count: int
    final_tie_count: int
    second_best_score_gap: float | None
    fixed_pair_rmsd_gap_angstrom: float | None
    mapping_ambiguous: bool
    selection_reason: str
    primary_score_atol: float
    fixed_pair_rmsd_atol_angstrom: float
    reflection_allowed: bool = False


def match_zmatrix_atoms(
    *,
    sites: Sequence[ZMatrixSite],
    reference_graph: nx.Graph,
    reference_coordinates: Mapping[int, Sequence[float]],
    canonical_to_reference: Mapping[str, int],
    experimental_graph: nx.Graph,
    experimental_coordinates: Mapping[int, Sequence[float]],
    flexible_coordinates: frozenset[InternalCoordinateKey] = frozenset(),
    settings: MappingScoreSettings = MappingScoreSettings(),
    max_heavy_mappings: int = DEFAULT_MAX_HEAVY_MAPPINGS,
    max_full_mappings: int = DEFAULT_MAX_FULL_MAPPINGS,
) -> MappingDecision:
    """Select a graph-preserving mapping using the shared Z-matrix contract.

    ``reference_graph`` and ``experimental_graph`` must contain one molecular
    component each.  Every node must carry an ``element`` attribute and every
    coordinate must be a molecule-contiguous Cartesian position in angstrom.
    Dihedrals listed in ``flexible_coordinates`` are classified as independent
    torsions for diagnostics, but they still contribute to the primary score.
    Flexible bonds and angles remain excluded when explicitly listed.
    """

    _validate_inputs(
        sites,
        reference_graph,
        reference_coordinates,
        canonical_to_reference,
        experimental_graph,
        experimental_coordinates,
        flexible_coordinates,
        settings,
        max_heavy_mappings,
        max_full_mappings,
    )
    reference_heavy = reference_graph.subgraph(
        [node for node in reference_graph if reference_graph.nodes[node]["element"] != "H"]
    ).copy()
    experimental_heavy = experimental_graph.subgraph(
        [
            node
            for node in experimental_graph
            if experimental_graph.nodes[node]["element"] != "H"
        ]
    ).copy()
    for node in reference_heavy:
        reference_heavy.nodes[node]["attached_hydrogens"] = sum(
            reference_graph.nodes[neighbor]["element"] == "H"
            for neighbor in reference_graph.neighbors(node)
        )
    for node in experimental_heavy:
        experimental_heavy.nodes[node]["attached_hydrogens"] = sum(
            experimental_graph.nodes[neighbor]["element"] == "H"
            for neighbor in experimental_graph.neighbors(node)
        )
    matcher = nx_isomorphism.GraphMatcher(
        reference_heavy,
        experimental_heavy,
        node_match=lambda left, right: (
            left["element"] == right["element"]
            and left["attached_hydrogens"] == right["attached_hydrogens"]
        ),
    )
    if not matcher.is_isomorphic():
        raise ValueError("Heavy-atom graphs are not isomorphic")

    candidates: list[MappingCandidate] = []
    heavy_mapping_count = 0
    full_mapping_count = 0
    candidates_truncated = False
    for heavy_mapping in matcher.isomorphisms_iter():
        heavy_mapping_count += 1
        if heavy_mapping_count > max_heavy_mappings:
            candidates_truncated = True
            break
        for reference_to_experimental in _enumerate_hydrogen_assignments(
            reference_graph, experimental_graph, heavy_mapping
        ):
            full_mapping_count += 1
            if full_mapping_count > max_full_mappings:
                candidates_truncated = True
                break
            canonical_mapping = {
                site.label: reference_to_experimental[
                    canonical_to_reference[site.label]
                ]
                for site in sites
            }
            candidates.append(
                score_mapping_candidate(
                    sites=sites,
                    reference_coordinates=reference_coordinates,
                    canonical_to_reference=canonical_to_reference,
                    experimental_coordinates=experimental_coordinates,
                    canonical_to_experimental=canonical_mapping,
                    flexible_coordinates=flexible_coordinates,
                    settings=settings,
                    include_internal_coordinates=False,
                )
            )
        if candidates_truncated:
            break

    if not candidates:
        raise ValueError("No complete graph-preserving atom mapping was generated")
    decision = select_mapping_candidate(
        candidates,
        settings=settings,
        heavy_mapping_count=min(heavy_mapping_count, max_heavy_mappings),
        full_mapping_count=min(full_mapping_count, max_full_mappings),
        candidates_truncated=candidates_truncated,
    )
    selected_mapping = decision.selected.as_index_mapping()
    selected_with_coordinates = score_mapping_candidate(
        sites=sites,
        reference_coordinates=reference_coordinates,
        canonical_to_reference=canonical_to_reference,
        experimental_coordinates=experimental_coordinates,
        canonical_to_experimental=selected_mapping,
        flexible_coordinates=flexible_coordinates,
        settings=settings,
        include_internal_coordinates=True,
    )
    return replace(decision, selected=selected_with_coordinates)


def score_mapping_candidate(
    *,
    sites: Sequence[ZMatrixSite],
    reference_coordinates: Mapping[int, Sequence[float]],
    canonical_to_reference: Mapping[str, int],
    experimental_coordinates: Mapping[int, Sequence[float]],
    canonical_to_experimental: Mapping[str, int],
    flexible_coordinates: frozenset[InternalCoordinateKey] = frozenset(),
    settings: MappingScoreSettings = MappingScoreSettings(),
    include_internal_coordinates: bool = True,
) -> MappingCandidate:
    """Measure one fixed atom permutation without rematching or reflection."""

    reference_values = internal_values(
        sites, canonical_to_reference, reference_coordinates
    )
    experimental_values = internal_values(
        sites, canonical_to_experimental, experimental_coordinates
    )
    weighted: list[float] = []
    bond_deltas: list[float] = []
    angle_deltas: list[float] = []
    rigid_torsion_deltas: list[float] = []
    independent_torsion_deltas: list[float] = []
    for site, reference_value, experimental_value in zip(
        sites, reference_values, experimental_values
    ):
        bond_key = InternalCoordinateKey("bond", site.index)
        if (
            reference_value.bond_angstrom is not None
            and bond_key not in flexible_coordinates
        ):
            delta = (
                float(experimental_value.bond_angstrom)
                - float(reference_value.bond_angstrom)
            )
            bond_deltas.append(delta)
            weighted.append((delta / settings.bond_scale_angstrom) ** 2)
        angle_key = InternalCoordinateKey("angle", site.index)
        if (
            reference_value.angle_degrees is not None
            and angle_key not in flexible_coordinates
        ):
            delta = (
                float(experimental_value.angle_degrees)
                - float(reference_value.angle_degrees)
            )
            angle_deltas.append(delta)
            weighted.append((delta / settings.angle_scale_degrees) ** 2)
        torsion_key = InternalCoordinateKey("dihedral", site.index)
        if reference_value.dihedral_degrees is not None:
            delta = circular_difference_degrees(
                float(experimental_value.dihedral_degrees),
                float(reference_value.dihedral_degrees),
            )
            if torsion_key in flexible_coordinates:
                independent_torsion_deltas.append(delta)
            else:
                rigid_torsion_deltas.append(delta)
            weighted.append((delta / settings.torsion_scale_degrees) ** 2)

    canonical_mapping = tuple(
        (site.label, int(canonical_to_experimental[site.label])) for site in sites
    )
    reference_all = np.asarray(
        [reference_coordinates[canonical_to_reference[site.label]] for site in sites],
        dtype=float,
    )
    experimental_all = np.asarray(
        [experimental_coordinates[canonical_to_experimental[site.label]] for site in sites],
        dtype=float,
    )
    heavy_positions = [
        position for position, site in enumerate(sites) if site.element != "H"
    ]
    return MappingCandidate(
        canonical_to_experimental=canonical_mapping,
        gross_bond_angle_mismatches=(
            sum(abs(value) > settings.gross_bond_delta_angstrom for value in bond_deltas)
            + sum(
                abs(value) > settings.gross_angle_delta_degrees
                for value in angle_deltas
            )
        ),
        internal_coordinate_score=math.sqrt(
            sum(weighted) / max(1, len(weighted))
        ),
        fixed_pair_all_atom_rmsd_angstrom=proper_rotation_rmsd(
            experimental_all, reference_all
        ),
        heavy_atom_rmsd_angstrom=proper_rotation_rmsd(
            experimental_all[heavy_positions], reference_all[heavy_positions]
        ),
        bond_rms_delta_angstrom=_rms(bond_deltas),
        bond_max_abs_delta_angstrom=_max_abs(bond_deltas),
        angle_rms_delta_degrees=_rms(angle_deltas),
        angle_max_abs_delta_degrees=_max_abs(angle_deltas),
        rigid_torsion_rms_delta_degrees=_rms(rigid_torsion_deltas),
        rigid_torsion_max_abs_delta_degrees=_max_abs(rigid_torsion_deltas),
        independent_torsion_rms_delta_degrees=_rms(
            independent_torsion_deltas
        ),
        independent_torsion_max_abs_delta_degrees=_max_abs(
            independent_torsion_deltas
        ),
        reference_values=reference_values if include_internal_coordinates else (),
        experimental_values=(
            experimental_values if include_internal_coordinates else ()
        ),
    )


def select_mapping_candidate(
    candidates: Sequence[MappingCandidate],
    *,
    settings: MappingScoreSettings = MappingScoreSettings(),
    heavy_mapping_count: int = 0,
    full_mapping_count: int | None = None,
    candidates_truncated: bool = False,
) -> MappingDecision:
    """Apply primary-score, fixed-RMSD, then deterministic selection stages."""

    if not candidates:
        raise ValueError("At least one mapping candidate is required")
    best_gross = min(item.gross_bond_angle_mismatches for item in candidates)
    gross_class = [
        item for item in candidates if item.gross_bond_angle_mismatches == best_gross
    ]
    gross_class.sort(
        key=lambda item: (item.internal_coordinate_score, item.deterministic_key)
    )
    best_primary = gross_class[0].internal_coordinate_score
    primary_ties = [
        item
        for item in gross_class
        if abs(item.internal_coordinate_score - best_primary)
        <= settings.primary_score_atol
    ]
    primary_ties.sort(
        key=lambda item: (
            item.fixed_pair_all_atom_rmsd_angstrom,
            item.deterministic_key,
        )
    )
    best_rmsd = primary_ties[0].fixed_pair_all_atom_rmsd_angstrom
    final_ties = [
        item
        for item in primary_ties
        if abs(item.fixed_pair_all_atom_rmsd_angstrom - best_rmsd)
        <= settings.fixed_pair_rmsd_atol_angstrom
    ]
    selected = min(final_ties, key=lambda item: item.deterministic_key)

    second_score_gap = (
        gross_class[1].internal_coordinate_score - best_primary
        if len(gross_class) > 1
        else None
    )
    rmsd_values = sorted(
        item.fixed_pair_all_atom_rmsd_angstrom for item in primary_ties
    )
    fixed_rmsd_gap = (
        rmsd_values[1] - rmsd_values[0] if len(rmsd_values) > 1 else None
    )
    if len(primary_ties) == 1:
        reason = "unique_primary_score"
    elif len(final_ties) == 1:
        reason = "fixed_pair_all_atom_rmsd"
    else:
        reason = "deterministic_unresolved_tie"
    return MappingDecision(
        method_version=MAPPING_METHOD_VERSION,
        selected=selected,
        heavy_mapping_count=heavy_mapping_count,
        full_mapping_count=(
            len(candidates) if full_mapping_count is None else full_mapping_count
        ),
        candidates_truncated=candidates_truncated,
        primary_tie_count=len(primary_ties),
        final_tie_count=len(final_ties),
        second_best_score_gap=second_score_gap,
        fixed_pair_rmsd_gap_angstrom=fixed_rmsd_gap,
        mapping_ambiguous=len(final_ties) > 1,
        selection_reason=reason,
        primary_score_atol=settings.primary_score_atol,
        fixed_pair_rmsd_atol_angstrom=(
            settings.fixed_pair_rmsd_atol_angstrom
        ),
        reflection_allowed=False,
    )


def internal_values(
    sites: Sequence[ZMatrixSite],
    mapping: Mapping[str, int],
    coordinates: Mapping[int, Sequence[float]],
) -> tuple[InternalCoordinateValues, ...]:
    """Measure the complete vector defined by ``sites`` for one mapping."""

    values: list[InternalCoordinateValues] = []
    for site in sites:
        atom = _coordinate(coordinates, mapping[site.label])
        bond = angle = dihedral = None
        if site.bond_to is not None:
            bond_atom = _coordinate(
                coordinates, mapping[sites[site.bond_to - 1].label]
            )
            bond = float(np.linalg.norm(atom - bond_atom))
        if site.angle_to is not None and site.bond_to is not None:
            angle_atom = _coordinate(
                coordinates, mapping[sites[site.angle_to - 1].label]
            )
            angle = angle_degrees(atom, bond_atom, angle_atom)
        if (
            site.dihedral_to is not None
            and site.angle_to is not None
            and site.bond_to is not None
        ):
            dihedral_atom = _coordinate(
                coordinates, mapping[sites[site.dihedral_to - 1].label]
            )
            dihedral = dihedral_degrees(
                atom, bond_atom, angle_atom, dihedral_atom
            )
        values.append(
            InternalCoordinateValues(
                site_index=site.index,
                label=site.label,
                bond_angstrom=bond,
                angle_degrees=angle,
                dihedral_degrees=dihedral,
            )
        )
    return tuple(values)


def proper_rotation_rmsd(mobile: np.ndarray, target: np.ndarray) -> float:
    """Fixed-pair RMSD after translation and a determinant-positive rotation."""

    mobile_array = np.asarray(mobile, dtype=float)
    target_array = np.asarray(target, dtype=float)
    if (
        mobile_array.shape != target_array.shape
        or mobile_array.ndim != 2
        or mobile_array.shape[1] != 3
        or mobile_array.shape[0] < 1
    ):
        raise ValueError("Kabsch arrays must have matching non-empty (N, 3) shape")
    centered_mobile = mobile_array - mobile_array.mean(axis=0)
    centered_target = target_array - target_array.mean(axis=0)
    covariance = centered_mobile.T @ centered_target
    left, _singular, right_transpose = np.linalg.svd(covariance)
    rotation = left @ right_transpose
    if np.linalg.det(rotation) < 0.0:
        left[:, -1] *= -1.0
        rotation = left @ right_transpose
    if np.linalg.det(rotation) <= 0.0:
        raise ValueError("Proper-rotation Kabsch alignment has non-positive determinant")
    aligned = centered_mobile @ rotation
    return math.sqrt(
        float(np.mean(np.sum((aligned - centered_target) ** 2, axis=1)))
    )


def circular_difference_degrees(left: float, right: float) -> float:
    """Signed shortest-periodic difference without reflection/sign inversion."""

    return ((left - right + 180.0) % 360.0) - 180.0


def read_mapping_artifact(
    path: str | Path,
    *,
    experimental_structure: str | Path | None = None,
    reference_paths: Mapping[int, str | Path] | None = None,
    zmatrix_paths: Mapping[int, str | Path] | None = None,
    require_csofm_reusable: bool = False,
) -> dict[str, object]:
    """Load a mapping artifact and optionally verify its reusable inputs."""

    source = Path(path)
    payload = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Mapping artifact must contain a JSON object: {source}")
    if payload.get("schema_version") != 1:
        raise ValueError(f"Unsupported mapping artifact schema in {source}")
    if payload.get("mapping_method_version") != MAPPING_METHOD_VERSION:
        raise ValueError(
            f"Stale mapping method in {source}: "
            f"{payload.get('mapping_method_version')!r}"
        )
    if payload.get("reflection_allowed") is not False:
        raise ValueError(f"Mapping artifact does not prohibit reflection: {source}")
    if require_csofm_reusable and payload.get("reusable_for_csofm") is not True:
        raise ValueError(f"Mapping artifact is not reusable for CSO-FM: {source}")
    components = payload.get("components")
    if not isinstance(components, list) or not components:
        raise ValueError(f"Mapping artifact has no molecular components: {source}")
    seen_experimental_indices: set[object] = set()
    for component in components:
        if not isinstance(component, dict):
            raise ValueError(f"Invalid component record in {source}")
        rows = component.get("atom_mapping")
        if not isinstance(rows, list) or not rows:
            raise ValueError(f"Mapping artifact component has no atoms: {source}")
        labels = [
            row.get("zmatrix_label") for row in rows if isinstance(row, dict)
        ]
        indices = [
            row.get("experimental_atom_index_0based")
            for row in rows
            if isinstance(row, dict)
        ]
        site_indices = [
            row.get("zmatrix_site_index_1based")
            for row in rows
            if isinstance(row, dict)
        ]
        fractional_coordinates = [
            row.get("experimental_fractional_coordinates")
            for row in rows
            if isinstance(row, dict)
        ]
        if (
            len(labels) != len(rows)
            or any(not isinstance(label, str) or not label for label in labels)
            or len(labels) != len(set(labels))
            or len(indices) != len(rows)
            or any(not isinstance(index, int) or index < 0 for index in indices)
            or len(indices) != len(set(indices))
            or site_indices != list(range(1, len(rows) + 1))
            or any(
                not isinstance(coordinates, list)
                or len(coordinates) != 3
                or any(
                    not isinstance(value, (int, float))
                    or not math.isfinite(float(value))
                    for value in coordinates
                )
                for coordinates in fractional_coordinates
            )
        ):
            raise ValueError(f"Mapping artifact component is not bijective: {source}")
        overlap = seen_experimental_indices.intersection(indices)
        if overlap:
            raise ValueError(
                f"Mapping artifact reuses experimental atom indices {sorted(overlap)}"
            )
        seen_experimental_indices.update(indices)
    if experimental_structure is not None:
        expected = payload.get("experimental_structure")
        if not isinstance(expected, dict) or expected.get("sha256") != _sha256_file(
            Path(experimental_structure)
        ):
            raise ValueError("Experimental structure hash differs from mapping artifact")
    if reference_paths is not None:
        recorded = payload.get("reference_structures")
        if not isinstance(recorded, dict):
            raise ValueError("Mapping artifact has no reference-structure hashes")
        for index, reference in reference_paths.items():
            item = recorded.get(str(index))
            if not isinstance(item, dict) or item.get("sha256") != _sha256_file(
                Path(reference)
            ):
                raise ValueError(
                    f"Reference hash for TYPE {index} differs from mapping artifact"
                )
    if zmatrix_paths is not None:
        recorded = payload.get("canonical_zmatrices")
        if not isinstance(recorded, dict):
            raise ValueError("Mapping artifact has no canonical Zmatrix hashes")
        for index, zmatrix in zmatrix_paths.items():
            item = recorded.get(str(index))
            if not isinstance(item, dict) or item.get("sha256") != _sha256_file(
                Path(zmatrix)
            ):
                raise ValueError(
                    f"Zmatrix hash for TYPE {index} differs from mapping artifact"
                )
    return payload


def angle_degrees(first: np.ndarray, center: np.ndarray, last: np.ndarray) -> float:
    left = first - center
    right = last - center
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator <= 1.0e-12:
        raise ValueError("Undefined bond angle in atom mapping")
    cosine = max(-1.0, min(1.0, float(np.dot(left, right) / denominator)))
    return math.degrees(math.acos(cosine))


def dihedral_degrees(
    first: np.ndarray,
    second: np.ndarray,
    third: np.ndarray,
    fourth: np.ndarray,
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


def _enumerate_hydrogen_assignments(
    reference_graph: nx.Graph,
    experimental_graph: nx.Graph,
    heavy_mapping: Mapping[int, int],
):
    groups: list[tuple[tuple[tuple[int, int], ...], ...]] = []
    for reference_heavy, experimental_heavy in sorted(heavy_mapping.items()):
        reference_hydrogens = sorted(
            node
            for node in reference_graph.neighbors(reference_heavy)
            if reference_graph.nodes[node]["element"] == "H"
        )
        experimental_hydrogens = sorted(
            node
            for node in experimental_graph.neighbors(experimental_heavy)
            if experimental_graph.nodes[node]["element"] == "H"
        )
        if len(reference_hydrogens) != len(experimental_hydrogens):
            raise ValueError("Hydrogen count differs for mapped heavy atoms")
        if reference_hydrogens:
            groups.append(
                tuple(
                    tuple(zip(reference_hydrogens, permutation))
                    for permutation in itertools.permutations(experimental_hydrogens)
                )
            )
    combinations = itertools.product(*groups) if groups else [()]
    for combination in combinations:
        result = dict(heavy_mapping)
        for assignments in combination:
            result.update(assignments)
        if len(result) != len(reference_graph) or len(set(result.values())) != len(result):
            raise ValueError("Incomplete or non-bijective full-atom mapping")
        yield result


def _validate_inputs(
    sites: Sequence[ZMatrixSite],
    reference_graph: nx.Graph,
    reference_coordinates: Mapping[int, Sequence[float]],
    canonical_to_reference: Mapping[str, int],
    experimental_graph: nx.Graph,
    experimental_coordinates: Mapping[int, Sequence[float]],
    flexible_coordinates: frozenset[InternalCoordinateKey],
    settings: MappingScoreSettings,
    max_heavy_mappings: int,
    max_full_mappings: int,
) -> None:
    if not sites or tuple(site.index for site in sites) != tuple(
        range(1, len(sites) + 1)
    ):
        raise ValueError("Z-matrix site indices must be contiguous and one-based")
    labels = tuple(site.label for site in sites)
    if len(labels) != len(set(labels)) or set(labels) != set(canonical_to_reference):
        raise ValueError("Canonical reference mapping must match unique Z-matrix labels")
    if len(reference_graph) != len(experimental_graph) or len(reference_graph) != len(sites):
        raise ValueError("Reference, experimental, and Z-matrix atom counts differ")
    if set(reference_graph) - set(reference_coordinates):
        raise ValueError("Missing reference Cartesian coordinates")
    if set(experimental_graph) - set(experimental_coordinates):
        raise ValueError("Missing experimental Cartesian coordinates")
    for graph_name, graph in (
        ("reference", reference_graph),
        ("experimental", experimental_graph),
    ):
        if not nx.is_connected(graph):
            raise ValueError(f"The {graph_name} graph must contain one connected molecule")
        if any("element" not in graph.nodes[node] for node in graph):
            raise ValueError(f"The {graph_name} graph lacks atom elements")
        if not any(graph.nodes[node]["element"] != "H" for node in graph):
            raise ValueError(f"The {graph_name} graph has no heavy atom")
    for site in sites:
        reference_index = canonical_to_reference[site.label]
        if reference_graph.nodes[reference_index]["element"] != site.element:
            raise ValueError(f"Reference element mismatch for {site.label}")
        for relation in (site.bond_to, site.angle_to, site.dihedral_to):
            if relation is not None and not 1 <= relation < site.index:
                raise ValueError(f"Invalid Z-matrix dependency in row {site.index}")
    valid_keys: set[InternalCoordinateKey] = set()
    for site in sites:
        if site.bond_to is not None:
            valid_keys.add(InternalCoordinateKey("bond", site.index))
        if site.angle_to is not None:
            valid_keys.add(InternalCoordinateKey("angle", site.index))
        if site.dihedral_to is not None:
            valid_keys.add(InternalCoordinateKey("dihedral", site.index))
    if not flexible_coordinates.issubset(valid_keys):
        raise ValueError("Flexible-coordinate set refers to a missing Z-matrix row")
    positive = (
        settings.bond_scale_angstrom,
        settings.angle_scale_degrees,
        settings.torsion_scale_degrees,
        settings.gross_bond_delta_angstrom,
        settings.gross_angle_delta_degrees,
        settings.primary_score_atol,
        settings.fixed_pair_rmsd_atol_angstrom,
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in positive):
        raise ValueError("Mapping scales and tolerances must be positive and finite")
    if max_heavy_mappings < 1 or max_full_mappings < 1:
        raise ValueError("Mapping enumeration limits must be positive")


def _coordinate(
    coordinates: Mapping[int, Sequence[float]], index: int
) -> np.ndarray:
    value = np.asarray(coordinates[index], dtype=float)
    if value.shape != (3,) or not np.all(np.isfinite(value)):
        raise ValueError(f"Invalid Cartesian coordinate for atom index {index}")
    return value


def _rms(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    return math.sqrt(sum(value * value for value in values) / len(values))


def _max_abs(values: Sequence[float]) -> float:
    return max((abs(value) for value in values), default=0.0)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
