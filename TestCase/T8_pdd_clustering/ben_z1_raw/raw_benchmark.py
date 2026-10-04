"""Full raw-minimum benchmark: no native fingerprint reduction or COMPACK."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import asdict
import hashlib
import json
import multiprocessing
import re
import resource
import sys
import time
import numpy as np
from benchmark import partition_metrics
from cp2_raw_reconstruction import ROOT, HistoricalLAM, ELEMENTS, asu_pdd, descriptor
from Source.pdd_descriptor import PDDDescriptor, calculate_pdd_from_arrays, pdd_descriptor_distance
from Source.pdd_clustering import PDDClusterRecord, PDDClusterSettings, cluster_pdd


def save(path, payload):
    with path.open('x') as stream:
        json.dump(payload,stream,indent=2); stream.write('\n')


def anchors():
    out = {}
    for batch in range(1,5):
        text = (ROOT/f'Analyse_{batch}/Analyse_log.out').read_text()
        blocks = re.split(r'Cluster:\s+(\d+) has\s+(\d+) structures',text)[1:]
        for i in range(0,len(blocks),3):
            rank, population, body = blocks[i:i+3]
            points = re.search(r'Participating minima:\s*\n([^\n]+)',body).group(1).split()
            out[f'A{batch}_{int(rank):05d}'] = (batch,int(points[0])-1)
    assert len(out)==4953
    return out


def raw_arrays():
    arrays={b:np.loadtxt(ROOT/f'Analyse_{b}/crystals_stable.out') for b in range(1,5)}
    assert [len(arrays[b]) for b in range(1,5)] == [21588,153097,217818,155932]
    assert all(a.shape[1]==21 and np.isfinite(a).all() for a in arrays.values())
    return arrays


def validate():
    start=time.perf_counter(); model=HistoricalLAM(); arrays=raw_arrays(); mapping=anchors()
    raw=np.loadtxt(ROOT/'crystals_stable.out')
    combined=np.concatenate(list(arrays.values()))
    # Original generation IDs are unique here, so compare every numeric record.
    assert len(set(raw[:,0]))==len(raw)==548435
    assert np.array_equal(raw[np.argsort(raw[:,0])],combined[np.argsort(combined[:,0])])
    cache=np.load(ROOT/'pdd_k100.npz',allow_pickle=False)
    manifest=list(csv.DictReader((ROOT/'manifest.csv').open()))
    errors=[]; symmetry_checks={}; worst=[]
    for i,row in enumerate(manifest):
        batch,index=mapping[row['id']]
        record=arrays[batch][index]
        motif,cell,density=model.reconstruct(record)
        desc=descriptor(asu_pdd(motif,cell))
        ref=PDDDescriptor(row['id'],100,cache[f'w{i}'],cache[f'd{i}'],tuple(cache[f'e{i}']),True,True,1e-4)
        error=pdd_descriptor_distance(desc,ref)
        errors.append(error);worst.append((error,row['id'],int(record[0]),int(record[1])))
        sg=int(record[1])
        if sg not in symmetry_checks:
            full=calculate_pdd_from_arrays(motif,cell,ELEMENTS*(len(motif)//15),collapse=False)
            discrepancy=pdd_descriptor_distance(desc,full)
            assert discrepancy < 1e-10,(sg,discrepancy)
            symmetry_checks[sg]=discrepancy
        if (i+1)%500==0: print(f'validated {i+1}/4953; max PDD error={max(errors):.6g} A',flush=True)
    cache.close()
    report=dict(passed=max(errors)<.002,structures=len(errors),lam_points=426,
                raw_record_union_exact=True,raw_count=len(raw),raw_sha256=hashlib.sha256((ROOT/'crystals_stable.out').read_bytes()).hexdigest(),
                historical_smooth=.4,historical_cutoff=1000.,
                archived_pdd_errors_angstrom=dict(max=max(errors),median=float(np.median(errors)),p99=float(np.percentile(errors,99))),
                worst=sorted(worst,reverse=True)[:20],asu_full_cell_equivalence_by_sg=symmetry_checks,
                seconds=time.perf_counter()-start)
    save(ROOT/'raw_reconstruction_validation.json',report)
    print(json.dumps(report,indent=2),flush=True)
    assert report['passed']


_MODEL=None
def init_worker():
    global _MODEL
    _MODEL=HistoricalLAM()


def make_chunk(payload):
    batch,start,records=payload
    distances=np.empty((len(records),15,100)); densities=np.empty(len(records))
    for i,record in enumerate(records):
        motif,cell,densities[i]=_MODEL.reconstruct(record)
        distances[i]=asu_pdd(motif,cell)
    return batch,start,distances,densities,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)


def generate(workers):
    assert json.loads((ROOT/'raw_reconstruction_validation.json').read_text())['passed']
    start=time.perf_counter(); arrays=raw_arrays(); distance_maps={};density_maps={}
    for b,a in arrays.items():
        for kind,shape,target in [('pdd',(len(a),15,100),distance_maps),('density',(len(a),),density_maps)]:
            path=ROOT/f'raw_{kind}_A{b}.npy'
            if path.exists(): raise FileExistsError(path)
            target[b]=np.lib.format.open_memmap(path,mode='w+',dtype=np.float64,shape=shape)
    tasks=((b,i,a[i:i+250]) for b,a in arrays.items() for i in range(0,len(a),250))
    done=0;peak=0;last=time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers,initializer=init_worker,mp_context=multiprocessing.get_context('spawn')) as pool:
        for b,i,distances,density,rss in pool.map(make_chunk,tasks,chunksize=1):
            distance_maps[b][i:i+len(distances)]=distances
            density_maps[b][i:i+len(density)]=density
            done+=len(distances);peak=max(peak,rss)
            if time.perf_counter()-last>15 or done==548435:
                print(f'descriptors {done}/548435; elapsed={time.perf_counter()-start:.1f}s; rate={done/(time.perf_counter()-start):.1f}/s',flush=True)
                last=time.perf_counter()
    for item in list(distance_maps.values())+list(density_maps.values()): item.flush()
    save(ROOT/'raw_descriptor_build.json',dict(structures=done,workers=workers,elapsed_s=time.perf_counter()-start,
         largest_worker_peak_rss_bytes=peak,parent_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
         k=100,typed=True,asu_orbit_rows=15,collapse=False,distance_bytes=548435*15*100*8))


def run(thresholds):
    assert json.loads((ROOT/'raw_descriptor_build.json').read_text())['structures']==548435
    arrays=raw_arrays(); mapping=anchors()
    manifest=list(csv.DictReader((ROOT/'manifest.csv').open()))
    reference={row['id']:row['reference'] for row in manifest}
    for threshold in thresholds:
        output=ROOT/f'raw_clusters_{threshold:.3f}.json'
        if output.exists(): raise FileExistsError(output)
        started=time.perf_counter(); batches={};predicted_anchors={};total_clusters=0
        for b,a in arrays.items():
            distances=np.load(ROOT/f'raw_pdd_A{b}.npy',mmap_mode='r')
            densities=np.load(ROOT/f'raw_density_A{b}.npy',mmap_mode='r')
            records=[PDDClusterRecord(f'R{b}_{i+1:06d}',r[4],densities[i],descriptor(distances[i])) for i,r in enumerate(a)]
            settings=PDDClusterSettings(threshold)
            print(f'clustering threshold={threshold:g}, batch={b}, minima={len(records)}',flush=True)
            result=cluster_pdd(records,settings)
            total_clusters+=len(result.representatives)
            labels=np.array([int(result.assignments[r.identifier].split('_')[1])-1 for r in records],dtype=np.int32)
            with (ROOT/f'raw_labels_{threshold:.3f}_A{b}.npy').open('xb') as f: np.save(f,labels)
            for name,(batch,index) in mapping.items():
                if batch==b: predicted_anchors[name]=result.assignments[records[index].identifier]
            batches[str(b)]=dict(clusters=len(result.representatives),statistics=result.statistics,
                representative_lines=[int(s.split('_')[1]) for s in result.representatives],
                populations=result.populations)
            print(f'completed threshold={threshold:g}, batch={b}, clusters={len(result.representatives)}, seconds={result.statistics["elapsed_s"]:.1f}',flush=True)
            del records,result,distances,densities
        report=dict(structures=548435,clusters=total_clusters,settings=asdict(settings),
                    independent_original_batches=True,native_reduction=False,compack=False,
                    workers=1,peak_rss_scope='single_process',
                    wall_s=time.perf_counter()-started,batches=batches,
                    anchor_metrics=partition_metrics(reference,predicted_anchors),
                    predicted_anchor_labels=predicted_anchors,
                    peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
        # Native anchor names and raw representative line IDs are different
        # namespaces; direct string equality of representative IDs is undefined.
        report['anchor_metrics'].pop('identical_representative_assignments',None)
        save(output,report); print(json.dumps({k:v for k,v in report.items() if k not in ('batches','predicted_anchor_labels')}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['validate','generate','run'])
    parser.add_argument('--workers',type=int,default=6)
    parser.add_argument('--thresholds',type=float,nargs='+',default=[.10,.15])
    args=parser.parse_args()
    if args.action=='validate': validate()
    elif args.action=='generate': generate(args.workers)
    else: run(args.thresholds)
