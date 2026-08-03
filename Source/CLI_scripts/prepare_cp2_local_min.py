"""Command-line entry point for CP2 experimental local minimisation."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

from Source.cp2_local_min import (
    CP2PBSSettings,
    DEFAULT_CX3_MODULES,
    collect_cp2_local_min_status,
    prepare_cp2_local_min_inputs,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare or inspect a CrystalPredictor2 experimental local-minimisation job."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser(
        "prepare",
        help="Map an experimental structure into CP2 LAM order and stage a job.",
    )
    prepare.add_argument("system_dir", type=Path)
    prepare.add_argument("experimental_structure", type=Path)
    prepare.add_argument("output_dir", type=Path)
    prepare.add_argument(
        "--global-search-dir",
        type=Path,
        help="Explicit global-search directory; otherwise discover 5_Globalsearch/5_GlobSrch.",
    )
    prepare.add_argument(
        "--reference",
        action="append",
        default=[],
        metavar="TYPE_INDEX=PATH",
        help="Canonical landscape/reference structure for a CP2 molecular type; repeat as needed.",
    )
    prepare.add_argument(
        "--zmatrix",
        action="append",
        default=[],
        metavar="TYPE_INDEX=PATH",
        help="Canonical system Zmatrix for a molecular type; TYPE 1 is discovered automatically.",
    )
    prepare.add_argument("--stage-mode", choices=("symlink", "copy"), default="symlink")
    prepare.add_argument(
        "--space-group", help="Override the experimental space-group symbol."
    )
    prepare.add_argument("--system-name")
    prepare.add_argument("--refcode")
    prepare.add_argument(
        "--compack-metadata",
        type=Path,
        help="Optional JSON metadata, such as the selected COMPACK RMSD_1 reference.",
    )
    prepare.add_argument("--max-heavy-mappings", type=int, default=100_000)
    prepare.add_argument("--allow-generated-labels", action="store_true")
    prepare.add_argument(
        "--allow-unvalidated",
        action="store_true",
        help="Allow experimental Z'>1 or multicomponent preparation.",
    )
    prepare.add_argument(
        "--allow-unvalidated-mapping",
        action="store_true",
        help=(
            "Allow a runnable job for other explicitly unvalidated mappings. "
            "This never waives a globally inverted mapping because no automatic "
            "achirality or symmetry-equivalence proof is implemented."
        ),
    )
    prepare.add_argument(
        "--cp2-executable",
        type=Path,
        help="Generate a runnable CX3 PBS bundle using this Minimise executable.",
    )
    prepare.add_argument("--pbs-walltime", default="24:00:00")
    prepare.add_argument("--pbs-memory-gb", type=int, default=8)
    prepare.add_argument(
        "--pbs-module",
        action="append",
        help="Override the default verified CX3 module stack; repeat in load order.",
    )
    prepare.add_argument("--pbs-queue")
    prepare.add_argument(
        "--nag-license-file",
        default="$HOME/.nag/license.dat",
        help=(
            "NAG Kusari licence path on the execution host; use an absolute path "
            "or a path beginning with '$HOME/'. Licence contents are never staged."
        ),
    )

    status = subparsers.add_parser(
        "status",
        help="Collect CP2 output, IFAIL, energy, and completion status as JSON.",
    )
    status.add_argument("job_dir", type=Path)
    status.add_argument("--log", type=Path, help="Explicit minimisation log path.")
    status.add_argument(
        "--no-write",
        action="store_true",
        help="Print status without writing cp2_local_min_status.json.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.command == "status":
        result = collect_cp2_local_min_status(
            args.job_dir, log_path=args.log, write_json=not args.no_write,
        )
        print(json.dumps(asdict(result), indent=2, sort_keys=True))
        return 0 if result.status == "optimizer_converged" else 1

    references = _parse_indexed_paths(args.reference, option_name="--reference")
    zmatrices = _parse_indexed_paths(args.zmatrix, option_name="--zmatrix")
    compack_metadata = None
    if args.compack_metadata is not None:
        compack_metadata = json.loads(args.compack_metadata.read_text(encoding="utf-8"))
        if not isinstance(compack_metadata, dict):
            raise ValueError("--compack-metadata must contain a JSON object")
    artifacts = prepare_cp2_local_min_inputs(
        args.system_dir,
        args.experimental_structure,
        args.output_dir,
        reference_paths=references,
        zmatrix_paths=zmatrices,
        global_search_dir=args.global_search_dir,
        stage_mode=args.stage_mode,
        space_group=args.space_group,
        system_name=args.system_name,
        refcode=args.refcode,
        compack_metadata=compack_metadata,
        max_heavy_mappings=args.max_heavy_mappings,
        allow_generated_labels=args.allow_generated_labels,
        allow_unvalidated=args.allow_unvalidated,
        allow_unvalidated_mapping=args.allow_unvalidated_mapping,
        cp2_executable=args.cp2_executable,
        pbs_settings=CP2PBSSettings(
            walltime=args.pbs_walltime,
            memory_gb=args.pbs_memory_gb,
            modules=tuple(args.pbs_module or DEFAULT_CX3_MODULES),
            nag_license_file=args.nag_license_file,
            queue=args.pbs_queue,
        ),
    )
    print(f"output_dir={artifacts.output_dir}")
    print(f"expcrys_pdb={artifacts.expcrys_pdb_path}")
    print(f"manifest={artifacts.manifest_path}")
    print(f"atom_mapping={artifacts.mapping_tsv_path}")
    if artifacts.pbs_script_path is not None:
        print(f"pbs_script={artifacts.pbs_script_path}")
    return 0


def _parse_indexed_paths(values: list[str], *, option_name: str) -> dict[int, Path]:
    paths: dict[int, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(
                f"Invalid {option_name} {value!r}; expected TYPE_INDEX=PATH"
            )
        raw_index, raw_path = value.split("=", 1)
        type_index = int(raw_index)
        if type_index < 1 or not raw_path:
            raise ValueError(
                f"Invalid {option_name} {value!r}; expected TYPE_INDEX=PATH"
            )
        if type_index in paths:
            raise ValueError(f"Duplicate {option_name} for molecular TYPE {type_index}")
        paths[type_index] = Path(raw_path)
    return paths


if __name__ == "__main__":
    raise SystemExit(main())
