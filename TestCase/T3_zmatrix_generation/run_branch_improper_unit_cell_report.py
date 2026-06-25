"""Generate branch-improper Z-matrices from exported unit-cell CIFs."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import os
from pathlib import Path
import sys


CASE_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = CASE_ROOT.parents[1]
DEFAULT_INPUT_DIR = (
    CASE_ROOT
    / "HA_pair_TPSS_CONTCAR_labelled_diagrams"
    / "unit_cells"
)
DEFAULT_OUTPUT_DIR = CASE_ROOT / "BranchImproper_API_Zmatrices"
DEFAULT_SUMMARY_CSV = CASE_ROOT / "branch_improper_unit_cell_report.csv"
DEFAULT_SUMMARY_MD = CASE_ROOT / "branch_improper_unit_cell_report.md"
CACHE_ROOT = Path("/private/tmp/csptoolbox_t3_branch_improper_cache")

os.environ.setdefault("MPLCONFIGDIR", str(CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_ROOT / "xdg"))
for cache_dir in (Path(os.environ["MPLCONFIGDIR"]), Path(os.environ["XDG_CACHE_HOME"])):
    cache_dir.mkdir(parents=True, exist_ok=True)

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import AtomRecord, CrystalStructure  # noqa: E402
from zmatrix_generation_helpers import (  # noqa: E402
    classify_dihedral_references,
    render_numeric_zmat_text,
    validate_zmatrix,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read exported unit-cell CIFs, select the organic non-water API molecule, "
            "and generate branch-improper Z-matrices."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Directory containing *_unit_cell.cif files. Default: {DEFAULT_INPUT_DIR}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for generated .zmat files. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=DEFAULT_SUMMARY_CSV,
        help=f"CSV report path. Default: {DEFAULT_SUMMARY_CSV}",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=DEFAULT_SUMMARY_MD,
        help=f"Markdown report path. Default: {DEFAULT_SUMMARY_MD}",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional maximum number of CIF files to process.",
    )
    return parser.parse_args()


def is_water_molecule(molecule: list[AtomRecord]) -> bool:
    return Counter(atom.element for atom in molecule) == Counter({"H": 2, "O": 1})


def is_organic_api_molecule(molecule: list[AtomRecord]) -> bool:
    elements = {atom.element for atom in molecule}
    return "C" in elements and not is_water_molecule(molecule)


def select_api_molecule_group(structure: CrystalStructure):
    groups = structure.deduplicate_molecules()
    organic_groups = [
        group
        for group in groups
        if is_organic_api_molecule(group.representative_molecule)
    ]
    if organic_groups:
        candidates = organic_groups
    else:
        candidates = [
            group
            for group in groups
            if not is_water_molecule(group.representative_molecule)
        ]
    if not candidates:
        raise ValueError("No non-water molecule group was found.")

    return max(
        candidates,
        key=lambda group: (
            len(group.representative_molecule),
            sum(atom.element == "C" for atom in group.representative_molecule),
            len(group.duplicate_molecules),
        ),
    )


def clear_zmat_files(directory: Path) -> None:
    if not directory.exists():
        return
    for path in directory.glob("*.zmat"):
        path.unlink()


def system_name_from_cif(path: Path) -> str:
    return path.stem.removesuffix("_unit_cell")


def molecule_formula(molecule: list[AtomRecord]) -> str:
    counts = Counter(atom.element for atom in molecule)
    return "".join(f"{element}{counts[element]}" for element in sorted(counts))


def run_report(
    *,
    input_dir: Path,
    output_dir: Path,
    limit: int | None = None,
) -> list[dict[str, object]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    clear_zmat_files(output_dir)

    cif_paths = sorted(input_dir.glob("*_unit_cell.cif"))
    if limit is not None:
        cif_paths = cif_paths[:limit]

    rows: list[dict[str, object]] = []
    for index, cif_path in enumerate(cif_paths, start=1):
        system_name = system_name_from_cif(cif_path)
        row: dict[str, object] = {
            "system_name": system_name,
            "status": "ok",
            "cif_path": str(cif_path),
            "selected_formula": "",
            "selected_atom_count": "",
            "duplicate_count": "",
            "zmatrix_path": "",
            "entry_count": "",
            "proper_count": "",
            "improper_count": "",
            "fallback_count": "",
            "ambiguous_count": "",
            "warning_count": "",
            "validation_error_count": "",
            "warnings": "",
            "validation_errors": "",
            "error": "",
        }

        try:
            structure = CrystalStructure.from_file(cif_path, fmt="cif")
            group = select_api_molecule_group(structure)
            molecule = group.representative_molecule
            gas_structure = structure.generate_gas_phase_vasp_structure(
                molecule,
                name=f"{system_name}_api",
            )
            zmatrices = gas_structure.generate_zmatrices(zmatrix_mode="branch_improper")
            if len(zmatrices) != 1:
                raise ValueError(f"Expected 1 API Z-matrix, generated {len(zmatrices)}.")

            zmatrix = zmatrices[0]
            validation_errors = validate_zmatrix(gas_structure, zmatrix)
            classifications = classify_dihedral_references(gas_structure, zmatrix)
            zmat_path = output_dir / f"{system_name}_branch_improper.zmat"
            zmat_path.write_text(
                render_numeric_zmat_text(f"{system_name} branch improper API", zmatrix),
                encoding="utf-8",
            )

            row.update(
                {
                    "selected_formula": molecule_formula(molecule),
                    "selected_atom_count": len(molecule),
                    "duplicate_count": len(group.duplicate_molecules),
                    "zmatrix_path": str(zmat_path),
                    "entry_count": len(zmatrix.entries),
                    "proper_count": sum(item.kind == "proper" for item in classifications),
                    "improper_count": sum(item.kind == "improper" for item in classifications),
                    "fallback_count": sum(item.kind == "fallback" for item in classifications),
                    "ambiguous_count": sum(item.kind == "ambiguous" for item in classifications),
                    "warning_count": len(zmatrix.warnings),
                    "validation_error_count": len(validation_errors),
                    "warnings": "; ".join(zmatrix.warnings),
                    "validation_errors": "; ".join(validation_errors),
                }
            )
        except Exception as exc:
            row["status"] = "error"
            row["error"] = str(exc)

        rows.append(row)
        if index % 10 == 0 or index == len(cif_paths):
            print(f"processed={index}/{len(cif_paths)}", flush=True)

    return rows


def write_csv_report(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else ["system_name", "status", "error"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def markdown_report(rows: list[dict[str, object]]) -> str:
    success_rows = [row for row in rows if row["status"] == "ok"]
    error_rows = [row for row in rows if row["status"] != "ok"]
    validation_error_rows = [
        row
        for row in success_rows
        if int(row["validation_error_count"]) > 0
    ]
    total_improper = sum(int(row["improper_count"]) for row in success_rows)
    total_proper = sum(int(row["proper_count"]) for row in success_rows)
    total_fallback = sum(int(row["fallback_count"]) for row in success_rows)
    total_ambiguous = sum(int(row["ambiguous_count"]) for row in success_rows)

    lines = [
        "# Branch-Improper Unit-Cell Z-Matrix Report",
        "",
        f"- CIF files processed: {len(rows)}",
        f"- Successful Z-matrices: {len(success_rows)}",
        f"- Errors: {len(error_rows)}",
        f"- Validation-error rows: {len(validation_error_rows)}",
        f"- Total proper references: {total_proper}",
        f"- Total improper references: {total_improper}",
        f"- Total fallback references: {total_fallback}",
        f"- Total ambiguous references: {total_ambiguous}",
        "",
        "| system | atoms | duplicates | proper | improper | fallback | ambiguous | warnings | validation errors |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]

    for row in success_rows:
        lines.append(
            "| "
            f"{row['system_name']} | "
            f"{row['selected_atom_count']} | "
            f"{row['duplicate_count']} | "
            f"{row['proper_count']} | "
            f"{row['improper_count']} | "
            f"{row['fallback_count']} | "
            f"{row['ambiguous_count']} | "
            f"{row['warning_count']} | "
            f"{row['validation_error_count']} |"
        )

    if error_rows:
        lines.extend(["", "## Errors", ""])
        for row in error_rows:
            lines.append(f"- {row['system_name']}: {row['error']}")

    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    rows = run_report(
        input_dir=args.input_dir.resolve(),
        output_dir=args.output_dir.resolve(),
        limit=args.limit,
    )
    if not rows:
        raise FileNotFoundError(f"No *_unit_cell.cif files found under {args.input_dir.resolve()}")

    write_csv_report(args.summary_csv.resolve(), rows)
    args.summary_md.resolve().write_text(markdown_report(rows), encoding="utf-8")

    success_count = sum(row["status"] == "ok" for row in rows)
    error_count = len(rows) - success_count
    validation_error_count = sum(
        int(row["validation_error_count"])
        for row in rows
        if row["status"] == "ok"
    )

    print(f"systems={len(rows)}")
    print(f"successes={success_count}")
    print(f"errors={error_count}")
    print(f"validation_errors={validation_error_count}")
    print(f"output_dir={args.output_dir.resolve()}")
    print(f"summary_csv={args.summary_csv.resolve()}")
    print(f"summary_md={args.summary_md.resolve()}")

    return 1 if error_count or validation_error_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
