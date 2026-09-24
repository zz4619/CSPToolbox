"""Cluster a CSV manifest of crystal files, reusing a provenance-checked cache."""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path
import time

from Source.pdd_clustering_io import load_pdd_cluster_manifest
from Source.pdd_clustering import PDDClusterSettings, cluster_pdd, make_compack_confirmation


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--threshold', type=float, required=True,
                        help='PDD EMD threshold in Å; must be calibrated, not copied from RMSD15.')
    parser.add_argument('--k', type=int, default=100)
    parser.add_argument('--energy-tolerance', type=float, default=1.0, help='kJ/mol; 0 disables this gate.')
    parser.add_argument('--density-tolerance', type=float, default=20.0, help='kg/m3; 0 disables this gate.')
    parser.add_argument('--index-dimensions', type=int, default=8)
    parser.add_argument('--block-size', type=int, default=64)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--exhaustive', action='store_true')
    parser.add_argument('--compack', action='store_true',
                        help='Confirm surviving pairs with licensed CCDC RMSD15 <0.20 Å.')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error('Output already exists; choose a new path.')
    start = time.perf_counter()
    records, identity = load_pdd_cluster_manifest(args.manifest, k=args.k, cache=args.cache)
    descriptor_s = time.perf_counter() - start
    settings = PDDClusterSettings(args.threshold, args.energy_tolerance or None,
                                  args.density_tolerance or None, args.index_dimensions, args.block_size)
    confirm = None
    if args.compack:
        with args.manifest.open(newline='') as stream:
            paths = {row['id']: (args.manifest.resolve().parent / row['path']).resolve()
                     for row in csv.DictReader(stream)}
        confirm = make_compack_confirmation(paths)
    result = cluster_pdd(records, settings, indexed=not args.exhaustive, confirm=confirm)
    payload = dict(asdict(result), settings=asdict(settings), k=args.k,
                   descriptor_load_or_build_s=descriptor_s, input_identity=identity,
                   elapsed_total_s=time.perf_counter() - start, compack_confirmation=args.compack)
    with args.output.open('x') as stream:
        json.dump(payload, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(structures=len(records), clusters=len(result.representatives),
                          statistics=result.statistics, elapsed_total_s=payload['elapsed_total_s']), indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
