from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.crystal_structure import CrystalStructure


DEFAULT_INPUT_ROOT = Path("/Users/zianzhan/Desktop/CSP_sandbox/CSP-personal/2_VASP/01_HA/Optimised_CIF")
DEFAULT_OUTPUT_CSV = PROJECT_ROOT / "TestCase" / "VASP" / "optimised_cif_contcar_chain_asymmetric_unit_molecule_counts.csv"
FILE_PATTERN = re.compile(r"^(?P<system>.+)_(?P<label>.+)_CONTCAR\.cif$")
POSCAR_PATTERN = re.compile(r"^(?P<system>.+)_(?P=system)_POSCAR\.cif$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reduce Optimised_CIF CONTCAR CIFs to asymmetric unit and count molecules per optimisation stage."
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help=f"Directory containing SYSTEM_X_CONTCAR.cif files. Default: {DEFAULT_INPUT_ROOT}",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=DEFAULT_OUTPUT_CSV,
        help=f"Destination CSV path. Default: {DEFAULT_OUTPUT_CSV}",
    )
    return parser.parse_args()


def _label_sort_key(system_name: str, label: str) -> tuple[int, int | str]:
    if label == system_name:
        return (0, label)
    if label.isdigit():
        return (1, int(label))
    return (2, label)


def _molecule_size_signature(path: Path) -> tuple[int, list[int]]:
    structure = CrystalStructure.from_file(path, fmt="cif")
    reduced = structure.reduce_to_asymmetric_unit(symprec=0.05)
    sizes = sorted((len(molecule) for molecule in reduced.detect_molecules()), reverse=True)
    return len(reduced.atoms), sizes


def summarize(input_root: Path) -> tuple[list[dict[str, object]], Counter[int], Counter[int]]:
    grouped: dict[str, list[tuple[str, Path]]] = defaultdict(list)
    starting_poscars: dict[str, Path] = {}
    max_numeric_label = 0

    for path in sorted(input_root.iterdir()):
        if not path.is_file():
            continue
        poscar_match = POSCAR_PATTERN.match(path.name)
        if poscar_match is not None:
            starting_poscars[poscar_match.group("system")] = path

        match = FILE_PATTERN.match(path.name)
        if match is None:
            continue
        system_name = match.group("system")
        label = match.group("label")
        grouped[system_name].append((label, path))
        if label.isdigit():
            max_numeric_label = max(max_numeric_label, int(label))

    rows: list[dict[str, object]] = []
    distribution: Counter[int] = Counter()
    changed_distribution: Counter[int] = Counter()

    for system_name in sorted(grouped):
        entries = sorted(
            grouped[system_name],
            key=lambda item: _label_sort_key(system_name, item[0]),
        )
        starting_poscar_path = starting_poscars.get(system_name)
        row: dict[str, object] = {
            "system_name": system_name,
            "starting_poscar_path": str(starting_poscar_path) if starting_poscar_path is not None else "",
            "starting_poscar_asymmetric_unit_atom_count": "",
            "starting_poscar_molecule_count": "",
            "starting_poscar_molecule_sizes": "",
            "stage_count": len(entries),
            "stage_labels": ";".join(label for label, _ in entries),
            "stage_molecule_counts": "",
            "stage_molecule_sizes": "",
            "final_stage_label": entries[-1][0],
            "final_stage_molecule_count": "",
            "final_stage_molecule_sizes": "",
            "molecule_size_signature_changed": "",
            "changed_steps": "",
            "first_changed_stage": "",
            "error": "",
        }

        for column_label in [system_name, *[str(index) for index in range(1, max_numeric_label + 1)]]:
            column_name = "stage_SYSTEM_molecule_count" if column_label == system_name else f"stage_{column_label}_molecule_count"
            row[column_name] = ""
            sizes_column_name = "stage_SYSTEM_molecule_sizes" if column_label == system_name else f"stage_{column_label}_molecule_sizes"
            row[sizes_column_name] = ""

        stage_counts: list[str] = []
        stage_sizes: list[str] = []
        try:
            previous_signature: tuple[int, ...] | None = None
            previous_label = "POSCAR"
            changed_steps: list[str] = []
            first_changed_stage = ""

            if starting_poscar_path is not None:
                starting_atom_count, starting_sizes = _molecule_size_signature(starting_poscar_path)
                row["starting_poscar_asymmetric_unit_atom_count"] = starting_atom_count
                row["starting_poscar_molecule_count"] = len(starting_sizes)
                row["starting_poscar_molecule_sizes"] = ";".join(str(size) for size in starting_sizes)
                previous_signature = tuple(starting_sizes)

            for label, path in entries:
                _, sizes = _molecule_size_signature(path)
                count = len(sizes)
                signature = tuple(sizes)
                stage_counts.append(f"{label}:{count}")
                stage_sizes.append(f"{label}:{';'.join(str(size) for size in sizes)}")
                column_name = "stage_SYSTEM_molecule_count" if label == system_name else f"stage_{label}_molecule_count"
                row[column_name] = count
                sizes_column_name = "stage_SYSTEM_molecule_sizes" if label == system_name else f"stage_{label}_molecule_sizes"
                row[sizes_column_name] = ";".join(str(size) for size in sizes)

                if previous_signature is not None and signature != previous_signature:
                    changed_steps.append(f"{previous_label}->{label}")
                    if not first_changed_stage:
                        first_changed_stage = label
                previous_signature = signature
                previous_label = label

            final_count = int(stage_counts[-1].split(":")[1])
            row["stage_molecule_counts"] = ";".join(stage_counts)
            row["stage_molecule_sizes"] = "|".join(stage_sizes)
            row["final_stage_molecule_count"] = final_count
            row["final_stage_molecule_sizes"] = row[
                "stage_SYSTEM_molecule_sizes" if entries[-1][0] == system_name else f"stage_{entries[-1][0]}_molecule_sizes"
            ]
            row["molecule_size_signature_changed"] = bool(changed_steps)
            row["changed_steps"] = ";".join(changed_steps)
            row["first_changed_stage"] = first_changed_stage
            distribution[final_count] += 1
            changed_distribution[int(bool(changed_steps))] += 1
        except Exception as exc:
            row["error"] = str(exc)

        rows.append(row)

    return rows, distribution, changed_distribution


def write_csv(rows: list[dict[str, object]], output_csv: Path) -> None:
    if not rows:
        raise ValueError("No SYSTEM_X_CONTCAR.cif files were summarised.")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows, distribution, changed_distribution = summarize(args.input_root.resolve())
    write_csv(rows, args.output_csv.resolve())

    print(f"systems={len(rows)}")
    print(f"output_csv={args.output_csv.resolve()}")
    for key in sorted(distribution):
        print(f"final_stage_molecule_count_{key}={distribution[key]}")
    print(f"molecule_signature_unchanged={changed_distribution[0]}")
    print(f"molecule_signature_changed={changed_distribution[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
