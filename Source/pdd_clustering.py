"""Exact AMD radius lookup followed by PDD/EMD representative clustering.

The index is a candidate filter, never an AMD-only merging decision. Results
match the exhaustive implementation for the same order, gates and PDD settings.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
import math
import time
from typing import Callable, Sequence

import numpy as np
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

from .pdd_descriptor import PDDDescriptor, pdd_to_amd, _earth_movers_distance

__all__ = ['PDDClusterRecord', 'PDDClusterSettings', 'PDDClusterResult',
           'cluster_pdd', 'threshold_agreement', 'make_compack_confirmation']


@dataclass(frozen=True)
class PDDClusterRecord:
    identifier: str
    energy: float
    density: float
    descriptor: PDDDescriptor


@dataclass(frozen=True)
class PDDClusterSettings:
    # No unvalidated default: 0.20 Å COMPACK RMSD is NOT 0.20 Å PDD EMD.
    threshold: float
    energy_tolerance: float | None = 1.0
    density_tolerance: float | None = 20.0
    index_dimensions: int = 8
    block_size: int = 64
    numeric_slack: float = 1e-10

    def __post_init__(self):
        if not math.isfinite(self.threshold) or self.threshold < 0:
            raise ValueError("PDD threshold must be finite and nonnegative.")
        for value in (self.energy_tolerance, self.density_tolerance):
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError("Energy/density tolerances must be positive, or None.")
        for value in (self.index_dimensions, self.block_size):
            if isinstance(value, bool) or int(value) != value or value < 1:
                raise ValueError("Index dimensions and block size must be positive integers.")
        if not math.isfinite(self.numeric_slack) or self.numeric_slack < 0:
            raise ValueError("Numeric slack must be finite and nonnegative.")


@dataclass
class PDDClusterResult:
    assignments: dict[str, str] = field(default_factory=dict)
    representatives: list[str] = field(default_factory=list)
    populations: dict[str, int] = field(default_factory=dict)
    statistics: dict[str, float | int] = field(default_factory=dict)


class _RadiusForest:
    """Immutable KDTree blocks plus a small always-searched pending buffer.

    Equal-size blocks merge like a binary counter. This avoids rebuilding one
    growing tree every few insertions. All blocks and pending entries are queried;
    newly created representatives cannot disappear between index rebuilds.
    """
    def __init__(self, vectors, dimensions, block_size):
        self.vectors = vectors
        self.dimensions = dimensions
        self.block_size = block_size
        self.pending = []
        self.blocks = {}
        self.builds = 0

    def add(self, index):
        self.pending.append(index)
        if len(self.pending) < self.block_size:
            return
        indexes = np.asarray(self.pending, dtype=int)
        self.pending = []
        level = 0
        while level in self.blocks:
            previous, _ = self.blocks.pop(level)
            indexes = np.concatenate((previous, indexes))
            level += 1
        self.blocks[level] = (indexes, KDTree(self.vectors[indexes][:, self.dimensions], copy_data=True))
        self.builds += 1

    def query(self, vector, radius):
        found = []
        for indexes, tree in self.blocks.values():
            positions = tree.query_ball_point(vector[self.dimensions], radius, p=np.inf, eps=0, workers=1)
            found.extend(indexes[positions].tolist())
        if self.pending:
            indexes = np.asarray(self.pending)
            distances = np.max(np.abs(self.vectors[indexes][:, self.dimensions] - vector[self.dimensions]), axis=1)
            found.extend(indexes[distances <= radius].tolist())
        return sorted(found)  # record indexes follow the deterministic energy order


def _composition_key(descriptor):
    if not descriptor.typed:
        return (False,)
    # Do not hash rounded fractions: roundoff near a bin boundary could separate
    # otherwise compatible supercells. Check fractions explicitly before EMD.
    return (True, tuple(sorted(set(descriptor.center_elements))))


def _prepare_comparison(descriptor):
    """Prepare an already validated descriptor once per record/representative.

    This private path is used only after pdd_to_amd validates every input. Keep
    just the current record and accepted representatives, not another full copy
    of a potentially memory-mapped input dataset.
    """
    groups = sorted(set(descriptor.center_elements)) if descriptor.typed else [None]
    prepared = []
    for element in groups:
        selected = (np.asarray(descriptor.center_elements) == element if element is not None
                    else np.ones(len(descriptor.weights), dtype=bool))
        weights = np.asarray(descriptor.weights)[selected]
        distances = np.asarray(descriptor.distances)[selected]
        uniform = float(weights[0]) if np.all(weights == weights[0]) else None
        prepared.append((weights, distances, uniform))
    return prepared


def _prepared_distance(a, b):
    """Same typed Chebyshev EMD, with membership/fractions checked by caller."""
    total = 0.
    for (wa, da, ua), (wb, db, ub) in zip(a, b):
        cost = cdist(da, db, metric='chebyshev')
        if ua is not None and ub == ua and len(wa) == len(wb):
            rows, cols = linear_sum_assignment(cost)
            total += float(ua * cost[rows, cols].sum())
        else:
            total += _earth_movers_distance(wa, wb, cost)
    return total


def cluster_pdd(
    records: Sequence[PDDClusterRecord], settings: PDDClusterSettings, *,
    indexed: bool = True,
    confirm: Callable[[PDDClusterRecord, PDDClusterRecord], bool] | None = None,
) -> PDDClusterResult:
    """Cluster by fixed representatives in (energy, identifier) order.

    Energy is kJ/mol and density kg/m3 when reproducing AXOSOW eligibility gates;
    callers must use a consistent normalization. Gates use strict inequalities,
    as in AXOSOW; EMD uses <= threshold. SG and cell atom count do not partition
    the index. Element sets partition it; fractions are checked before EMD.
    Optional confirmation is called
    only after EMD passes (e.g. COMPACK); AMD/EMD pruning is not a proven bound
    on a different confirmation metric.
    """
    start = time.perf_counter()
    ordered = sorted(records, key=lambda record: (record.energy, record.identifier))
    if len({r.identifier for r in ordered}) != len(ordered):
        raise ValueError("Structure identifiers must be unique.")
    if any(not isinstance(r.identifier, str) or not r.identifier for r in ordered):
        raise ValueError("Structure identifiers must be nonempty strings.")
    if any(not math.isfinite(r.energy) or not math.isfinite(r.density) for r in ordered):
        raise ValueError("Energy and density must be finite.")
    result = PDDClusterResult()
    stats = dict(structures=len(ordered), index_candidates=0, amd_checks=0,
                 eligibility_passes=0, emd_calls=0, confirmation_calls=0,
                 exhaustive_representatives=0)
    if not ordered:
        result.statistics = dict(stats, elapsed_s=time.perf_counter() - start, index_builds=0)
        return result
    if len({(r.descriptor.k, r.descriptor.typed) for r in ordered}) != 1:
        raise ValueError("All records must use the same k and typing mode.")
    if len({(r.descriptor.collapse, r.descriptor.collapse_tol) for r in ordered}) != 1:
        raise ValueError("All records must use consistent PDD collapse settings.")
    amds = np.asarray([pdd_to_amd(r.descriptor) for r in ordered])
    fractions = [np.asarray([r.descriptor.rows_for_element(e)[0].sum()
                            for e in sorted(set(r.descriptor.center_elements))]) for r in ordered]
    dimensions = np.unique(np.linspace(0, amds.shape[1] - 1,
                          min(settings.index_dimensions, amds.shape[1]), dtype=int))
    forests, representatives, prepared_representatives = {}, {}, {}
    base_radius = settings.threshold + settings.numeric_slack
    radius = np.nextafter(base_radius, np.inf) if base_radius > 0 else 0.
    # Search metadata ranges in the same tree, rather than testing energy and
    # density only after retrieving a large AMD shortlist. Final gates below
    # retain their original strict inequalities. Scaling does not redefine EMD.
    if radius > 0:
        parts = [amds[:, dimensions] / radius]
        for values, tolerance in (([r.energy for r in ordered], settings.energy_tolerance),
                                  ([r.density for r in ordered], settings.density_tolerance)):
            if tolerance is not None:
                values = np.asarray(values)
                parts.append(((values - values.min()) / tolerance)[:, None])
        index_vectors = np.column_stack(parts)
        if not np.all(np.isfinite(index_vectors)):
            raise ValueError("Metadata/threshold scaling overflows the index coordinates.")
        # Include a conservative arithmetic margin when searching normalized
        # coordinates, then recheck unscaled values for final membership.
        query_radius = 1 + 64 * np.finfo(float).eps * max(1., float(np.max(np.abs(index_vectors))))
    else:
        index_vectors = amds[:, dimensions]
        query_radius = 0.
    index_dimensions = np.arange(index_vectors.shape[1])
    for i, record in enumerate(ordered):
        key = _composition_key(record.descriptor)
        if key not in forests:
            forests[key] = _RadiusForest(index_vectors, index_dimensions, settings.block_size)
            representatives[key] = []
        forest = forests[key]
        previous = representatives[key]
        stats['exhaustive_representatives'] += len(previous)
        candidates = forest.query(index_vectors[i], query_radius) if indexed else previous
        stats['index_candidates'] += len(candidates)
        # Both paths use the same complete-vector test. The exhaustive path is
        # intentionally vectorized, rather than an artificially slow Python scan.
        if len(candidates):
            stats['amd_checks'] += len(candidates)
            differences = np.max(np.abs(amds[candidates] - amds[i]), axis=1)
            candidates = np.asarray(candidates)[differences <= radius].tolist()
        representative = None
        prepared = None
        for j in candidates:
            other = ordered[j]
            if settings.energy_tolerance is not None and abs(record.energy - other.energy) >= settings.energy_tolerance:
                continue
            if settings.density_tolerance is not None and abs(record.density - other.density) >= settings.density_tolerance:
                continue
            if record.descriptor.typed and not np.allclose(fractions[i], fractions[j], atol=1e-12, rtol=0):
                continue
            stats['eligibility_passes'] += 1
            stats['emd_calls'] += 1
            if prepared is None:
                prepared = _prepare_comparison(record.descriptor)
            distance = _prepared_distance(prepared, prepared_representatives[j])
            if distance > settings.threshold:
                continue
            if confirm is not None:
                stats['confirmation_calls'] += 1
                if not confirm(record, other):
                    continue
            representative = other.identifier
            break
        if representative is None:
            representative = record.identifier
            result.representatives.append(representative)
            result.populations[representative] = 0
            previous.append(i)
            forest.add(i)
            prepared_representatives[i] = prepared if prepared is not None else _prepare_comparison(record.descriptor)
        result.assignments[record.identifier] = representative
        result.populations[representative] += 1
    result.statistics = dict(stats, elapsed_s=time.perf_counter() - start,
                             index_builds=sum(f.builds for f in forests.values()))
    return result


def threshold_agreement(distances, reference_matches, thresholds):
    """Evaluate labelled pairs; do not mistake a mostly-unique count for accuracy.

    Call separately for calibration and held-out pairs. Splits should be made by
    structure/family before constructing pairs, not by randomly splitting pairs.
    """
    distances = np.asarray(distances, dtype=float)
    matches = np.asarray(reference_matches, dtype=bool)
    if distances.ndim != 1 or distances.shape != matches.shape or not np.all(np.isfinite(distances)) or np.any(distances < 0):
        raise ValueError("Provide aligned finite nonnegative distances and labels.")
    output = []
    for threshold in thresholds:
        if not math.isfinite(threshold) or threshold < 0:
            raise ValueError("Thresholds must be finite and nonnegative.")
        predicted = distances <= threshold
        tp = int(np.count_nonzero(predicted & matches))
        fp = int(np.count_nonzero(predicted & ~matches))
        fn = int(np.count_nonzero(~predicted & matches))
        tn = int(np.count_nonzero(~predicted & ~matches))
        output.append(dict(threshold=float(threshold), true_positive=tp, false_merge=fp,
                           missed_duplicate=fn, true_negative=tn,
                           precision=tp / (tp + fp) if tp + fp else None,
                           recall=tp / (tp + fn) if tp + fn else None))
    return output


def make_compack_confirmation(paths, *, shell_size=15, rmsd_tolerance=.20, cache_size=256):
    """Optional licensed CCDC confirmation with the AXOSOW packing settings.

    ``paths`` maps record identifiers to structure files. CCDC is imported only
    when requested; this function never installs or configures that dependency.
    """
    from ccdc.io import CrystalReader
    from ccdc.crystal import PackingSimilarity
    if shell_size < 1 or rmsd_tolerance <= 0 or not math.isfinite(rmsd_tolerance) or cache_size < 1:
        raise ValueError("Packing shell, RMSD tolerance and cache size must be positive.")
    similarity = PackingSimilarity()
    similarity.settings.ignore_hydrogen_counts = False
    similarity.settings.ignore_bond_counts = False
    similarity.settings.packing_shell_size = shell_size

    @lru_cache(maxsize=cache_size)
    def read(identifier):
        return CrystalReader(str(paths[identifier]))[0]

    def confirm(a, b):
        match = similarity.compare(read(a.identifier), read(b.identifier))
        return (match is not None and match.nmatched_molecules == shell_size
                and match.rmsd < rmsd_tolerance)

    return confirm
