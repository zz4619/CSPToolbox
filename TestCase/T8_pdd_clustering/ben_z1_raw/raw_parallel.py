"""Parallelize independent historical batches, retaining serial cluster order."""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import argparse
import csv
import json
import multiprocessing
import resource
import sys
import time
import numpy as np
from cp2_raw_reconstruction import ROOT,descriptor
from raw_benchmark import anchors,save
from benchmark import partition_metrics
from Source.pdd_clustering import PDDClusterRecord,PDDClusterSettings,cluster_pdd


def one_batch(task):
    b,threshold=task
    a=np.loadtxt(ROOT/f'Analyse_{b}/crystals_stable.out')
    distances=np.load(ROOT/f'raw_pdd_A{b}.npy',mmap_mode='r')
    densities=np.load(ROOT/f'raw_density_A{b}.npy',mmap_mode='r')
    records=[PDDClusterRecord(f'R{b}_{i+1:06d}',r[4],densities[i],descriptor(distances[i])) for i,r in enumerate(a)]
    print(f'start batch {b}, threshold {threshold}, structures {len(records)}',flush=True)
    result=cluster_pdd(records,PDDClusterSettings(threshold))
    labels=np.array([int(result.assignments[r.identifier].split('_')[1])-1 for r in records],dtype=np.int32)
    with (ROOT/f'raw_labels_{threshold:.3f}_A{b}.npy').open('xb') as f:np.save(f,labels)
    predicted={name:result.assignments[records[index].identifier] for name,(batch,index) in anchors().items() if batch==b}
    report=dict(clusters=len(result.representatives),statistics=result.statistics,
                representative_lines=[int(s.split('_')[1]) for s in result.representatives],
                populations=result.populations,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
    save(ROOT/f'raw_batch_{threshold:.3f}_A{b}.json',dict(batch=b,report=report,anchors=predicted))
    print(f'complete batch {b}, clusters {len(result.representatives)}, seconds {result.statistics["elapsed_s"]:.1f}',flush=True)
    return b,report,predicted


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threshold',type=float,required=True)
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    if not 1<=args.workers<=4:parser.error('workers must be between 1 and 4')
    target=ROOT/f'raw_clusters_{args.threshold:.3f}.json'
    if target.exists():raise FileExistsError(target)
    # Refuse collisions before starting workers, not after long calculations.
    for b in range(1,5):
        for name in [f'raw_labels_{args.threshold:.3f}_A{b}.npy',f'raw_batch_{args.threshold:.3f}_A{b}.json']:
            if (ROOT/name).exists():raise FileExistsError(ROOT/name)
    assert json.loads((ROOT/'raw_descriptor_build.json').read_text())['structures']==548435
    manifest=list(csv.DictReader((ROOT/'manifest.csv').open()))
    reference={row['id']:row['reference'] for row in manifest}
    batches={};predicted={};start=time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        for b,report,labels in pool.map(one_batch,[(b,args.threshold) for b in range(1,5)]):
            batches[str(b)]=report;predicted.update(labels)
    metrics=partition_metrics(reference,predicted)
    metrics.pop('identical_representative_assignments',None)
    report=dict(structures=548435,clusters=sum(b['clusters'] for b in batches.values()),
                settings=asdict(PDDClusterSettings(args.threshold)),independent_original_batches=True,
                native_reduction=False,compack=False,workers=args.workers,
                wall_s=time.perf_counter()-start,batches=batches,anchor_metrics=metrics,
                predicted_anchor_labels=predicted,
                peak_rss_bytes=max(b['peak_rss_bytes'] for b in batches.values()),
                peak_rss_scope='largest_worker_not_aggregate',
                parent_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
    save(target,report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('batches','predicted_anchor_labels')}),flush=True)


if __name__=='__main__':main()
