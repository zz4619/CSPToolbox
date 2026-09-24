"""Reproducible comparison with Ben's saved Nicotinamide Z'=1 clustering.

Consumes every native first-representative RES in all four saved analysis
batches. The archived reference logs are read-only; no CCDC is required.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import re
import resource
import sys
import time

import numpy as np

from Source.pdd_clustering import PDDClusterSettings, cluster_pdd
from Source.pdd_clustering_io import load_pdd_cluster_manifest

ROOT = Path(__file__).resolve().parent


def write_json(path, payload):
    with path.open('x') as stream:
        json.dump(payload, stream, indent=2)
        stream.write('\n')


def identifier(batch, rank):
    return f'A{batch}_{int(rank):05d}'


def prepare():
    rows, batches = [], {}
    for batch in range(1, 5):
        folder = ROOT / f'Analyse_{batch}'
        native = (folder / 'Analyse_log.out').read_text()
        log = (folder / 'Clustering_log.out').read_text()
        assert log.rstrip().endswith('Clustering terminated successfully')
        blocks = re.split(r'Cluster:\s+(\d+) has\s+(\d+) structures', native)[1:]
        metadata = {}
        for offset in range(0, len(blocks), 3):
            rank, population, body = blocks[offset:offset + 3]
            metadata[int(rank)] = dict(
                energy=float(re.search(r'Utot:\s+([-\d.]+)', body).group(1)),
                density=float(re.search(r'Density:\s+([-\d.]+)', body).group(1)),
                native_population=int(population),
            )
        cluster_to_source, assignment = {}, {}
        comparisons = 0
        for line in log.splitlines():
            match = re.match(r'(\d+) is unique, making it clustered/(\d+)', line)
            if match:
                rank, cluster = map(int, match.groups())
                cluster_to_source[cluster] = rank
                assignment[rank] = rank
            match = re.match(r'(\d+) matches (\d+) ! not unique', line)
            if match:
                rank, cluster = map(int, match.groups())
                assignment[rank] = cluster_to_source[cluster]
            comparisons += bool(re.match(r'\d+ and \d+ results:', line))
        assert set(assignment) == set(metadata)
        assert len(assignment) == int(re.search(r'there are (\d+) structures', log).group(1))
        native_count = int(re.search(r'NUMBER OF MINIMA FOUND:\s+(\d+)', native).group(1))
        assert sum(r['native_population'] for r in metadata.values()) == native_count
        for rank in sorted(metadata):
            relpath = f'Analyse_{batch}/unique_pool/{rank}-1/{rank}.res'
            assert (ROOT / relpath).is_file(), relpath
            rows.append(dict(id=identifier(batch, rank), path=relpath,
                             **metadata[rank], batch=batch,
                             reference=identifier(batch, assignment[rank])))
        cpu_progress = re.findall(r'processed (\d+) structures\s+([\d.]+) mins', log)
        batches[str(batch)] = dict(
            successful_raw_minima=native_count, native_exports=len(metadata),
            reference_clusters=len(cluster_to_source),
            reference_duplicate_assignments=len(assignment)-len(cluster_to_source),
            logged_compack_comparisons=comparisons,
            reference_log_sha256=hashlib.sha256(log.encode()).hexdigest(),
            historical_last_progress_cpu_s=float(cpu_progress[-1][1])*60,
            historical_last_progress_structure=int(cpu_progress[-1][0]),
            energy_min=min(m['energy'] for m in metadata.values()),
            energy_max=max(m['energy'] for m in metadata.values()),
        )
    with (ROOT / 'manifest.csv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    audit = dict(source='/rds/general/user/zz4619/home/03_Ben_CSP/Nicotinamide/5_GlobSrch',
                 batches=batches, exported_inputs=len(rows),
                 raw_minima=sum(b['successful_raw_minima'] for b in batches.values()),
                 reference_clusters=sum(b['reference_clusters'] for b in batches.values()),
                 note='All first-representative exports, not all raw minima; original four batch boundaries retained.')
    write_json(ROOT / 'reference_audit.json', audit)
    print(json.dumps(audit, indent=2), flush=True)


def partition_metrics(reference, predicted):
    """Label-invariant agreement, emphasizing same-cluster pairs and splits."""
    assert set(reference) == set(predicted)
    truth = Counter(reference.values())
    guess = Counter(predicted.values())
    cells = Counter((reference[key], predicted[key]) for key in reference)
    choose2 = lambda n: n*(n-1)//2
    true_pairs = sum(choose2(n) for n in truth.values())
    predicted_pairs = sum(choose2(n) for n in guess.values())
    tp = sum(choose2(n) for n in cells.values())
    total_pairs = choose2(len(reference))
    expected = true_pairs * predicted_pairs / total_pairs if total_pairs else 0.
    denominator = .5*(true_pairs+predicted_pairs)-expected
    ref_splits = Counter(r for r, p in cells)
    pred_merges = Counter(p for r, p in cells)
    b3_precision = sum(n*n/guess[p] for (r, p), n in cells.items()) / len(reference)
    b3_recall = sum(n*n/truth[r] for (r, p), n in cells.items()) / len(reference)
    return dict(structures=len(reference), reference_clusters=len(truth), predicted_clusters=len(guess),
                reference_same_cluster_pairs=true_pairs, predicted_same_cluster_pairs=predicted_pairs,
                shared_same_cluster_pairs=tp,
                pair_precision=tp/predicted_pairs if predicted_pairs else None,
                pair_recall=tp/true_pairs if true_pairs else None,
                adjusted_rand=(tp-expected)/denominator if denominator else 1.,
                b_cubed_precision=b3_precision, b_cubed_recall=b3_recall,
                reference_clusters_split=sum(n>1 for n in ref_splits.values()),
                predicted_clusters_merging_reference=sum(n>1 for n in pred_merges.values()),
                identical_representative_assignments=sum(reference[k]==predicted[k] for k in reference))


def run(thresholds):
    start = time.perf_counter()
    rows = list(csv.DictReader((ROOT / 'manifest.csv').open()))
    cache = ROOT / 'pdd_k100.npz'
    cache_existed = cache.exists()
    print(f'Descriptors: {len(rows)} structures; existing cache={cache_existed}', flush=True)
    records, identity = load_pdd_cluster_manifest(ROOT / 'manifest.csv', cache=cache)
    preparation_s = time.perf_counter()-start
    nrows = [len(r.descriptor.weights) for r in records]
    prep = dict(input_identity=identity, structures=len(records), descriptor_k=100,
                descriptor_load_or_build_s=preparation_s, loaded_existing_cache=cache_existed,
                descriptor_rows_min=min(nrows), descriptor_rows_max=max(nrows),
                descriptor_rows_median=float(np.median(nrows)),
                distance_bytes=sum(r.descriptor.distances.nbytes for r in records))
    prep_name = 'cache_load.json' if cache_existed else 'descriptor_build.json'
    if not (ROOT / prep_name).exists():
        write_json(ROOT / prep_name, prep)
    print(json.dumps(prep), flush=True)
    by_batch = {batch: [r for r in records if r.identifier.startswith(f'A{batch}_')]
                for batch in range(1,5)}
    reference = {row['id']: row['reference'] for row in rows}
    for threshold in thresholds:
        output = ROOT / f'pdd_only_{threshold:.3f}.json'
        if output.exists():
            raise FileExistsError(output)
        settings = PDDClusterSettings(threshold=threshold)
        batches, assignments = {}, {}
        timer = time.perf_counter()
        for batch, subset in by_batch.items():
            result = cluster_pdd(subset, settings)
            assignments.update(result.assignments)
            subset_ref = {r.identifier:reference[r.identifier] for r in subset}
            batches[str(batch)] = dict(statistics=result.statistics,
                                      metrics=partition_metrics(subset_ref, result.assignments))
            print(f'threshold={threshold:g} batch={batch} structures={len(subset)} '
                  f'clusters={len(result.representatives)} seconds={result.statistics["elapsed_s"]:.3f}', flush=True)
        report = dict(settings=asdict(settings), input_identity=identity,
                      compack_confirmation=False, independent_original_batches=True,
                      clustering_wall_s=time.perf_counter()-timer,
                      metrics=partition_metrics(reference, assignments), batches=batches,
                      assignments=assignments,
                      process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                                            *(1 if sys.platform=='darwin' else 1024))
        write_json(output, report)
        print(json.dumps(dict(threshold=threshold,metrics=report['metrics'],seconds=report['clustering_wall_s'])),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','run'])
    parser.add_argument('--thresholds',type=float,nargs='+',default=[.05,.10,.15,.20,.30])
    args=parser.parse_args()
    if args.action=='prepare':
        prepare()
    else:
        run(args.thresholds)
