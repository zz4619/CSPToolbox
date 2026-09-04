"""Prepare many experimental CP2 cases and assemble one sequential PBS job."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from Source.cp2_local_min import (
    CP2PBSSettings,
    DEFAULT_CX3_MODULES,
    build_cp2_local_min_batch,
    prepare_cp2_local_min_inputs,
)


STATUS_FILENAME = "cp2_local_min_preparation_status.tsv"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare every row in a CP2 experimental-structure inventory and "
            "assemble the runnable cases into one sequential PBS job."
        )
    )
    parser.add_argument("inventory", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--cp2-executable", type=Path, required=True)
    parser.add_argument("--stage-mode", choices=("copy", "symlink"), default="copy")
    parser.add_argument("--allow-unvalidated", action="store_true")
    parser.add_argument("--allow-unvalidated-mapping", action="store_true")
    parser.add_argument(
        "--auto-single-type-occurrences",
        action="store_true",
        help=(
            "For one-TYPE systems, adapt a case-local input.in occurrence "
            "count to the experimental asymmetric-unit component count."
        ),
    )
    parser.add_argument("--allow-generated-labels", action="store_true")
    parser.add_argument("--max-heavy-mappings", type=int, default=100_000)
    parser.add_argument("--pbs-walltime", default="24:00:00")
    parser.add_argument("--pbs-memory-gb", type=int, default=8)
    parser.add_argument("--pbs-queue")
    parser.add_argument("--per-case-timeout", default="15m")
    parser.add_argument("--job-name", default="cp2lm_all_exp")
    parser.add_argument(
        "--pbs-module",
        action="append",
        help="Override the default CX3 runtime modules; repeat in load order.",
    )
    parser.add_argument(
        "--nag-license-file",
        default="$HOME/.nag/ChemEngDept-nag_keys-2026-linux.txt",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = _read_inventory(args.inventory)
    destination = args.output_dir
    if destination.exists() and not destination.is_dir():
        raise FileExistsError(f"Batch output path is not a directory: {destination}")
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError(
            f"Refusing to overwrite non-empty CP2 batch directory {destination}"
        )
    destination.mkdir(parents=True, exist_ok=True)

    settings = CP2PBSSettings(
        walltime=args.pbs_walltime,
        memory_gb=args.pbs_memory_gb,
        modules=tuple(args.pbs_module or DEFAULT_CX3_MODULES),
        nag_license_file=args.nag_license_file,
        queue=args.pbs_queue,
    )
    status_rows: list[dict[str, object]] = []
    runnable_cases: list[Path] = []
    for row in rows:
        system_name = row["system_name"]
        refcode = row["refcode"]
        case_dir = destination / "cases" / system_name / refcode
        try:
            references = _indexed_paths(row, "reference_")
            zmatrices = _indexed_paths(row, "zmatrix_")
            compack_metadata = _optional_json(row.get("compack_metadata_json", ""))
            artifacts = prepare_cp2_local_min_inputs(
                row["system_dir"],
                row["experimental_structure"],
                case_dir,
                reference_paths=references,
                zmatrix_paths=zmatrices,
                global_search_dir=row.get("global_search_dir") or None,
                stage_mode=args.stage_mode,
                space_group=row.get("space_group") or None,
                system_name=system_name,
                refcode=refcode,
                compack_metadata=compack_metadata,
                max_heavy_mappings=args.max_heavy_mappings,
                allow_generated_labels=args.allow_generated_labels,
                allow_unvalidated=args.allow_unvalidated,
                allow_unvalidated_mapping=args.allow_unvalidated_mapping,
                auto_single_type_occurrences=args.auto_single_type_occurrences,
                trust_experimental_labels=_optional_bool(
                    row.get("trust_experimental_labels", "")
                ),
                cp2_executable=args.cp2_executable,
                pbs_settings=settings,
            )
            manifest = json.loads(
                artifacts.manifest_path.read_text(encoding="utf-8")
            )
            runnable_cases.append(case_dir)
            status_rows.append(
                {
                    "system_name": system_name,
                    "refcode": refcode,
                    "preparation_status": "prepared_runnable",
                    "supported_scope": manifest["supported_scope"],
                    "structurally_supported": manifest["structurally_supported"],
                    "mapping_validated": manifest["mapping_validated"],
                    "case_dir": str(case_dir.resolve()),
                    "error": "",
                }
            )
        except Exception as error:  # preserve every failed inventory row for audit
            status_rows.append(
                {
                    "system_name": system_name,
                    "refcode": refcode,
                    "preparation_status": "preparation_failed",
                    "supported_scope": "",
                    "structurally_supported": "",
                    "mapping_validated": "",
                    "case_dir": str(case_dir.absolute()),
                    "error": f"{type(error).__name__}: {error}",
                }
            )

    status_path = destination / STATUS_FILENAME
    _write_tsv(status_path, status_rows)
    if not runnable_cases:
        raise RuntimeError(
            f"No runnable CP2 cases were prepared; inspect {status_path}"
        )
    batch = build_cp2_local_min_batch(
        destination,
        runnable_cases,
        settings=settings,
        job_name=args.job_name,
        per_case_timeout=args.per_case_timeout,
    )
    print(f"inventory_rows={len(rows)}")
    print(f"prepared_runnable={len(runnable_cases)}")
    print(f"preparation_failed={len(rows) - len(runnable_cases)}")
    print(f"preparation_status={status_path}")
    print(f"cases_manifest={batch.cases_manifest_path}")
    print(f"pbs_script={batch.pbs_script_path}")
    return 0 if len(runnable_cases) == len(rows) else 2


def _read_inventory(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "system_name",
            "refcode",
            "system_dir",
            "experimental_structure",
            "reference_1",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "CP2 batch inventory is missing columns: " + ", ".join(sorted(missing))
            )
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError(f"CP2 batch inventory is empty: {path}")
    identities = [(row["system_name"], row["refcode"]) for row in rows]
    if len(identities) != len(set(identities)):
        raise ValueError("CP2 batch inventory contains duplicate system/refcode rows")
    for row in rows:
        for key in required:
            if not row.get(key):
                raise ValueError(
                    f"CP2 batch inventory has an empty {key} for "
                    f"{row.get('system_name')}/{row.get('refcode')}"
                )
    return rows


def _indexed_paths(row: dict[str, str], prefix: str) -> dict[int, Path]:
    paths: dict[int, Path] = {}
    for key, value in row.items():
        if not key.startswith(prefix) or not value:
            continue
        suffix = key[len(prefix) :]
        if not suffix.isdigit() or int(suffix) < 1:
            continue
        paths[int(suffix)] = Path(value)
    return paths


def _optional_json(value: str) -> dict[str, object] | None:
    if not value:
        return None
    parsed = json.loads(value)
    if not isinstance(parsed, dict):
        raise ValueError("compack_metadata_json must contain a JSON object")
    return parsed


def _optional_bool(value: str) -> bool:
    if not value:
        return False
    normalized = value.strip().lower()
    if normalized in {"true", "yes", "1"}:
        return True
    if normalized in {"false", "no", "0"}:
        return False
    raise ValueError(f"Invalid boolean value {value!r}")


def _write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
