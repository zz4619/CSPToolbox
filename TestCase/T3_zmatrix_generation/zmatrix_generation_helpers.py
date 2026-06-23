"""Helpers for T3 Z-matrix generation tests."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import math
import os
from pathlib import Path
import sys


CASE_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = CASE_ROOT.parents[1]
CACHE_ROOT = Path("/private/tmp/csptoolbox_t3_zmatrix_cache")
os.environ.setdefault("MPLCONFIGDIR", str(CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_ROOT / "xdg"))
for cache_dir in (Path(os.environ["MPLCONFIGDIR"]), Path(os.environ["XDG_CACHE_HOME"])):
    cache_dir.mkdir(parents=True, exist_ok=True)

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import (  # noqa: E402
    DEFAULT_COVALENT_SCALE,
    AtomRecord,
    CrystalStructure,
    ZMatrixEntry,
    ZMatrixRepresentation,
)


@dataclass(frozen=True)
class ZMatrixGenerationCase:
    name: str
    structure: CrystalStructure
    expected_improper_labels: tuple[str, ...] = ()
    expected_proper_label_sets: tuple[frozenset[str], ...] = ()


@dataclass(frozen=True)
class DihedralReferenceClassification:
    row_index: int
    atom_label: str
    bond_label: str
    angle_label: str
    dihedral_label: str
    kind: str

    @property
    def label_set(self) -> frozenset[str]:
        return frozenset(
            (self.atom_label, self.bond_label, self.angle_label, self.dihedral_label)
        )

    @property
    def quartet(self) -> str:
        return (
            f"{self.atom_label}-{self.bond_label}-"
            f"{self.angle_label}-{self.dihedral_label}"
        )


@dataclass(frozen=True)
class ZMatrixCaseReport:
    case_name: str
    atom_count: int
    zmat_count: int
    proper_count: int
    improper_count: int
    fallback_count: int
    ambiguous_count: int
    warning_count: int
    validation_error_count: int


def sample_cases() -> tuple[ZMatrixGenerationCase, ...]:
    return (butane_chain_case(), methanol_case())


def butane_chain_case() -> ZMatrixGenerationCase:
    atoms = [
        AtomRecord("C1", "C", (10.00, 10.00, 10.00)),
        AtomRecord("C2", "C", (11.54, 10.00, 10.00)),
        AtomRecord("C3", "C", (12.10, 11.43, 10.00)),
        AtomRecord("C4", "C", (13.64, 11.43, 10.60)),
    ]
    structure = _gas_phase_structure("butane_chain", atoms)
    return ZMatrixGenerationCase(
        name="butane_chain",
        structure=structure,
        expected_proper_label_sets=(frozenset(("C1", "C2", "C3", "C4")),),
    )


def methanol_case() -> ZMatrixGenerationCase:
    atoms = [
        AtomRecord("C1", "C", (10.00, 10.00, 10.00)),
        AtomRecord("O1", "O", (11.43, 10.00, 10.00)),
        AtomRecord("H1", "H", (9.45, 10.90, 10.00)),
        AtomRecord("H2", "H", (9.45, 9.10, 10.00)),
        AtomRecord("H3", "H", (10.00, 10.00, 11.09)),
        AtomRecord("H4", "H", (11.90, 10.75, 10.00)),
    ]
    structure = _gas_phase_structure("methanol", atoms)
    return ZMatrixGenerationCase(
        name="methanol",
        structure=structure,
        expected_improper_labels=("H2", "H3"),
    )


def _gas_phase_structure(name: str, atoms: list[AtomRecord]) -> CrystalStructure:
    return CrystalStructure(
        atoms=atoms,
        cell_parameters=(30.0, 30.0, 30.0, 90.0, 90.0, 90.0),
        space_group="P 1",
        name=name,
        explict_unit_cell=True,
    )


def first_zmatrix(structure: CrystalStructure) -> ZMatrixRepresentation:
    zmatrices = structure.generate_zmatrices()
    if len(zmatrices) != 1:
        raise AssertionError(f"{structure.name} generated {len(zmatrices)} Z-matrices, expected 1.")
    return zmatrices[0]


def validate_zmatrix(
    structure: CrystalStructure,
    zmatrix: ZMatrixRepresentation,
) -> list[str]:
    errors: list[str] = []
    structure_labels = {atom.label for atom in structure.atoms}

    if len(zmatrix.entries) != len(zmatrix.ordered_atom_labels):
        errors.append("entry count does not match ordered label count")
    if len(zmatrix.ordered_atom_labels) != len(set(zmatrix.ordered_atom_labels)):
        errors.append("ordered atom labels are not unique")
    missing_labels = sorted(set(zmatrix.ordered_atom_labels) - structure_labels)
    if missing_labels:
        errors.append(f"ordered atom labels not present in structure: {missing_labels}")

    for row_index, entry in enumerate(zmatrix.entries, start=1):
        if row_index <= len(zmatrix.ordered_atom_labels):
            expected_label = zmatrix.ordered_atom_labels[row_index - 1]
            if entry.label != expected_label:
                errors.append(
                    f"row {row_index} entry label {entry.label} does not match order label {expected_label}"
                )
        errors.extend(_validate_entry_references(row_index, entry))

    return errors


def _validate_entry_references(row_index: int, entry: ZMatrixEntry) -> list[str]:
    errors: list[str] = []

    if row_index == 1:
        if any(
            value is not None
            for value in (
                entry.bond_to,
                entry.bond_length,
                entry.angle_to,
                entry.angle_degrees,
                entry.dihedral_to,
                entry.dihedral_degrees,
            )
        ):
            errors.append("row 1 should not have internal-coordinate references")
        return errors

    if entry.bond_to is None or entry.bond_length is None:
        errors.append(f"row {row_index} is missing a bond reference/value")
    elif not 1 <= entry.bond_to < row_index:
        errors.append(f"row {row_index} bond reference does not point to an earlier row")
    elif not _is_positive_finite(entry.bond_length):
        errors.append(f"row {row_index} bond length is not positive finite")

    if row_index >= 3:
        if entry.angle_to is None or entry.angle_degrees is None:
            errors.append(f"row {row_index} is missing an angle reference/value")
        elif not 1 <= entry.angle_to < row_index:
            errors.append(f"row {row_index} angle reference does not point to an earlier row")
        elif not _is_finite(entry.angle_degrees):
            errors.append(f"row {row_index} angle is not finite")

    if row_index >= 4:
        if entry.dihedral_to is None or entry.dihedral_degrees is None:
            errors.append(f"row {row_index} is missing a dihedral reference/value")
        elif not 1 <= entry.dihedral_to < row_index:
            errors.append(f"row {row_index} dihedral reference does not point to an earlier row")
        elif not _is_finite(entry.dihedral_degrees):
            errors.append(f"row {row_index} dihedral is not finite")

    return errors


def _is_positive_finite(value: float) -> bool:
    return _is_finite(value) and value > 0.0


def _is_finite(value: float) -> bool:
    return math.isfinite(float(value))


def classify_dihedral_references(
    structure: CrystalStructure,
    zmatrix: ZMatrixRepresentation,
    *,
    covalent_scale: float = DEFAULT_COVALENT_SCALE,
) -> list[DihedralReferenceClassification]:
    bond_pairs = bond_pairs_by_label(structure, covalent_scale=covalent_scale)
    labels = zmatrix.ordered_atom_labels
    classifications: list[DihedralReferenceClassification] = []

    for row_index, entry in enumerate(zmatrix.entries, start=1):
        if entry.bond_to is None or entry.angle_to is None or entry.dihedral_to is None:
            continue

        atom_label = labels[row_index - 1]
        bond_label = labels[entry.bond_to - 1]
        angle_label = labels[entry.angle_to - 1]
        dihedral_label = labels[entry.dihedral_to - 1]

        has_atom_bond = _bonded(atom_label, bond_label, bond_pairs)
        has_bond_angle = _bonded(bond_label, angle_label, bond_pairs)
        is_proper = (
            has_atom_bond
            and has_bond_angle
            and _bonded(angle_label, dihedral_label, bond_pairs)
        )
        is_improper = (
            has_atom_bond
            and has_bond_angle
            and _bonded(bond_label, dihedral_label, bond_pairs)
        )

        if is_proper and is_improper:
            kind = "ambiguous"
        elif is_proper:
            kind = "proper"
        elif is_improper:
            kind = "improper"
        else:
            kind = "fallback"

        classifications.append(
            DihedralReferenceClassification(
                row_index=row_index,
                atom_label=atom_label,
                bond_label=bond_label,
                angle_label=angle_label,
                dihedral_label=dihedral_label,
                kind=kind,
            )
        )

    return classifications


def bond_pairs_by_label(
    structure: CrystalStructure,
    *,
    covalent_scale: float = DEFAULT_COVALENT_SCALE,
) -> set[frozenset[str]]:
    adjacency = structure._build_connectivity(covalent_scale)
    pairs: set[frozenset[str]] = set()
    for atom_index, edges in adjacency.items():
        left_label = structure.atoms[atom_index].label
        for edge in edges:
            right_label = structure.atoms[edge.neighbor].label
            pairs.add(frozenset((left_label, right_label)))
    return pairs


def _bonded(left: str, right: str, bond_pairs: set[frozenset[str]]) -> bool:
    return frozenset((left, right)) in bond_pairs


def render_numeric_zmat_text(title: str, zmatrix: ZMatrixRepresentation) -> str:
    lines = [f"# ZMAT v1", f"# title: {title}"]
    for entry in zmatrix.entries:
        row = [entry.element]
        if entry.bond_to is not None and entry.bond_length is not None:
            row.extend([str(entry.bond_to), f"{entry.bond_length:.10f}"])
        if entry.angle_to is not None and entry.angle_degrees is not None:
            row.extend([str(entry.angle_to), f"{entry.angle_degrees:.10f}"])
        if entry.dihedral_to is not None and entry.dihedral_degrees is not None:
            row.extend([str(entry.dihedral_to), f"{entry.dihedral_degrees:.10f}"])
        lines.append(" ".join(row))
    return "\n".join(lines) + "\n"


def case_report(case: ZMatrixGenerationCase) -> ZMatrixCaseReport:
    zmatrices = case.structure.generate_zmatrices()
    if zmatrices:
        zmatrix = zmatrices[0]
        classifications = classify_dihedral_references(case.structure, zmatrix)
        validation_error_count = len(validate_zmatrix(case.structure, zmatrix))
        warning_count = len(zmatrix.warnings)
    else:
        classifications = []
        validation_error_count = 1
        warning_count = 0

    return ZMatrixCaseReport(
        case_name=case.name,
        atom_count=len(case.structure.atoms),
        zmat_count=len(zmatrices),
        proper_count=sum(item.kind == "proper" for item in classifications),
        improper_count=sum(item.kind == "improper" for item in classifications),
        fallback_count=sum(item.kind == "fallback" for item in classifications),
        ambiguous_count=sum(item.kind == "ambiguous" for item in classifications),
        warning_count=warning_count,
        validation_error_count=validation_error_count,
    )


def write_case_artifacts(output_dir: Path) -> list[ZMatrixCaseReport]:
    output_dir.mkdir(parents=True, exist_ok=True)
    reports: list[ZMatrixCaseReport] = []
    for case in sample_cases():
        zmatrix = first_zmatrix(case.structure)
        (output_dir / f"{case.name}.zmat").write_text(
            render_numeric_zmat_text(case.name, zmatrix),
            encoding="utf-8",
        )
        reports.append(case_report(case))
    return reports


def write_csv_report(path: Path, reports: list[ZMatrixCaseReport]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_name",
                "atom_count",
                "zmat_count",
                "proper_count",
                "improper_count",
                "fallback_count",
                "ambiguous_count",
                "warning_count",
                "validation_error_count",
            ],
        )
        writer.writeheader()
        for report in reports:
            writer.writerow(report.__dict__)


def markdown_report(reports: list[ZMatrixCaseReport]) -> str:
    lines = [
        "# T3 Z-Matrix Generation Report",
        "",
        "| case | atoms | zmatrices | proper | improper | fallback | ambiguous | warnings | validation errors |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for report in reports:
        lines.append(
            "| "
            f"{report.case_name} | "
            f"{report.atom_count} | "
            f"{report.zmat_count} | "
            f"{report.proper_count} | "
            f"{report.improper_count} | "
            f"{report.fallback_count} | "
            f"{report.ambiguous_count} | "
            f"{report.warning_count} | "
            f"{report.validation_error_count} |"
        )
    lines.append("")
    return "\n".join(lines)
