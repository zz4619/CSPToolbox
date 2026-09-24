"""Ordered batch clustering for immutable, uniform-row PDD array datasets.

Workers search a frozen representative snapshot. Only unmatched rows need ordered
checks against representatives created inside the batch. This preserves the
serial energy/identifier order; it is not independent clustering of shards.
The array path requires identical central-element order and uniform row weights.
Use cluster_pdd for mixed/collapsed descriptors or COMPACK confirmation.
"""
from __future__ import annotations

from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import time

import numpy as np

from .pdd_clustering import PDDClusterSettings, _RadiusForest, _prepare_comparison, _prepared_distance
from .pdd_descriptor import PDDDescriptor

__all__ = ['prepare_pdd_batch_store', 'cluster_pdd_batched']
_AMD_SOURCE = None
_WORKER = None


def _save_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def _file_identity(path):
    p = Path(path).resolve(); stat = p.stat()
    return dict(path=str(p), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def _digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def _amd_init(path):
    global _AMD_SOURCE
    _AMD_SOURCE = np.load(path, mmap_mode='r')


def _amd_chunk(bounds):
    start, end = bounds
    data = _AMD_SOURCE[start:end]
    if (not np.all(np.isfinite(data)) or np.any(data < 0)
            or np.any(np.diff(data, axis=2) < -1e-12)):
        raise ValueError(f'Invalid PDD distances at rows {start}:{end}')
    # Same uniform weighted column means as the serial PDD descriptor.
    weights = np.full(data.shape[1], 1. / data.shape[1])
    return start, weights @ data


def prepare_pdd_batch_store(pdd_path, energies, densities, identifiers, elements,
                            directory, *, workers=1, source_sha256=None):
    """Validate/cache AMD and ordering once for both thresholds.

PDD is a float64 NPY array (structures, ASU rows, neighbours), with equal positive
row weights and the same element sequence in every record. Numerical metadata
are in the same raw order. Directory must be new; completed stores are immutable.
"""
    root = Path(directory); root.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter(); identity = _file_identity(pdd_path)
    data = np.load(pdd_path, mmap_mode='r')
    if data.ndim != 3 or data.dtype != np.float64 or min(data.shape) < 1:
        raise ValueError('Expected a nonempty float64 (N, rows, k) NPY PDD array')
    n, rows, k = data.shape
    energies = np.asarray(energies, dtype=np.float64)
    densities = np.asarray(densities, dtype=np.float64)
    identifiers = np.asarray(identifiers)
    if energies.shape != (n,) or densities.shape != (n,) or identifiers.shape != (n,):
        raise ValueError('Metadata must have one entry for every PDD')
    if identifiers.dtype.kind not in 'US' or len(set(identifiers.tolist())) != n or np.any(identifiers == ''):
        raise ValueError('Identifiers must be unique nonempty strings')
    if not np.all(np.isfinite(energies)) or not np.all(np.isfinite(densities)):
        raise ValueError('Nonfinite energy/density metadata')
    if len(elements) != rows or any(not isinstance(e, str) or not e for e in elements):
        raise ValueError('Every row needs a central element')
    if workers < 1:
        raise ValueError('workers must be positive')
    order = np.lexsort((identifiers, energies))
    np.save(root/'order.npy', order)
    np.save(root/'energies.npy', energies[order])
    np.save(root/'densities.npy', densities[order])
    np.save(root/'identifiers.npy', identifiers)
    raw_amd = np.lib.format.open_memmap(root/'raw_amd.npy', mode='w+', dtype=np.float64, shape=(n, k))
    tasks = [(i, min(n, i+1024)) for i in range(0, n, 1024)]
    if workers == 1:
        _amd_init(str(Path(pdd_path).resolve()))
        for start, value in map(_amd_chunk, tasks): raw_amd[start:start+len(value)] = value
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=_amd_init,
                initargs=(str(Path(pdd_path).resolve()),), mp_context=multiprocessing.get_context('spawn')) as pool:
            for start, value in pool.map(_amd_chunk, tasks, chunksize=1):
                raw_amd[start:start+len(value)] = value
    raw_amd.flush()
    amd = np.lib.format.open_memmap(root/'amd.npy', mode='w+', dtype=np.float64, shape=(n,k))
    for start in range(0,n,4096): amd[start:start+4096] = raw_amd[order[start:start+4096]]
    amd.flush(); del amd, raw_amd
    (root/'raw_amd.npy').unlink()  # Task-owned temporary; sorted AMD is now durable.
    if identity != _file_identity(pdd_path): raise ValueError('PDD source changed during preparation')
    meta = dict(version=1, structures=n, rows=rows, k=k, elements=list(elements),
                uniform_weights=True, collapse=False, pdd_identity=identity,
                source_sha256=source_sha256, seconds=time.perf_counter()-started,
                array_hashes={p.name:_digest(p) for p in root.glob('*.npy')})
    _save_json(root/'store.json', meta)
    return meta


def _new_stats():
    return dict(index_candidates=0, amd_checks=0, emd_calls=0, index_seconds=0.,
                amd_filter_seconds=0., emd_seconds=0., index_sync_seconds=0.)


def _add_stats(left, right):
    for key, value in right.items(): left[key] += value


class _State:
    def __init__(self, store, output, settings, cache_records):
        self.root = Path(store); self.output = Path(output); self.settings = settings
        self.meta = json.loads((self.root/'store.json').read_text())
        self.data = np.load(self.meta['pdd_identity']['path'], mmap_mode='r')
        self.order = np.load(self.root/'order.npy', mmap_mode='r')
        self.amd = np.load(self.root/'amd.npy', mmap_mode='r')
        self.energy = np.load(self.root/'energies.npy', mmap_mode='r')
        self.density = np.load(self.root/'densities.npy', mmap_mode='r')
        self.vectors = np.load(self.output/'index_vectors.npy', mmap_mode='r')
        self.reps = np.load(self.output/'published_representatives.npy', mmap_mode='r')
        self.elements = tuple(self.meta['elements'])
        self.weights = np.full(self.meta['rows'], 1. / self.meta['rows'])
        config = json.loads((self.output/'index.json').read_text())
        self.radius, self.query_radius = config['amd_radius'], config['query_radius']
        self.forest = self.new_forest(); self.published = 0
        self.cache = OrderedDict(); self.cache_records = cache_records

    def new_forest(self):
        return _RadiusForest(self.vectors, np.arange(self.vectors.shape[1]), self.settings.block_size)

    def prepared(self, i):
        if i in self.cache:
            self.cache.move_to_end(i); return self.cache[i]
        desc = PDDDescriptor(str(i), self.meta['k'], self.weights,
                self.data[self.order[i]], self.elements, True, False, 0.)
        value = _prepare_comparison(desc)
        self.cache[i] = value
        if len(self.cache) > self.cache_records: self.cache.popitem(last=False)
        return value

    def find(self, i, forest, stats):
        timer = time.perf_counter()
        candidates = forest.query(self.vectors[i], self.query_radius)
        stats['index_seconds'] += time.perf_counter()-timer
        stats['index_candidates'] += len(candidates)
        if not candidates: return -1
        timer = time.perf_counter()
        stats['amd_checks'] += len(candidates)
        indexes = np.asarray(candidates)
        keep = np.max(np.abs(self.amd[indexes]-self.amd[i]), axis=1) <= self.radius
        if self.settings.energy_tolerance is not None:
            keep &= np.abs(self.energy[indexes]-self.energy[i]) < self.settings.energy_tolerance
        if self.settings.density_tolerance is not None:
            keep &= np.abs(self.density[indexes]-self.density[i]) < self.settings.density_tolerance
        candidates = indexes[keep]
        stats['amd_filter_seconds'] += time.perf_counter()-timer
        if not len(candidates): return -1
        timer = time.perf_counter(); current = self.prepared(i)
        for j in candidates:
            stats['emd_calls'] += 1
            # Same scalar EMD and acceptance threshold as cluster_pdd. Element
            # fractions are equal by this fixed-layout dataset's input contract.
            if _prepared_distance(current, self.prepared(int(j))) <= self.settings.threshold:
                stats['emd_seconds'] += time.perf_counter()-timer
                return int(j)
        stats['emd_seconds'] += time.perf_counter()-timer
        return -1

    def task(self, task):
        start, end, count = task
        stats = _new_stats(); timer = time.perf_counter()
        if count < self.published: raise RuntimeError('Representative snapshot moved backwards')
        for j in self.reps[self.published:count]: self.forest.add(int(j))
        self.published = count; stats['index_sync_seconds'] = time.perf_counter()-timer
        matches = np.asarray([self.find(i, self.forest, stats) for i in range(start,end)], dtype=np.int64)
        return start, matches, stats


def _worker_init(store, output, settings, cache_records):
    global _WORKER
    _WORKER = _State(store, output, settings, cache_records)


def _worker_task(task): return _WORKER.task(task)


def _make_index(root, output, settings):
    amd = np.load(root/'amd.npy', mmap_mode='r'); n, k = amd.shape
    dimensions = np.unique(np.linspace(0,k-1,min(settings.index_dimensions,k),dtype=int))
    radius = settings.threshold+settings.numeric_slack
    radius = np.nextafter(radius,np.inf) if radius>0 else 0.
    meta = []
    if radius>0:
        for name, tolerance in [('energies',settings.energy_tolerance),('densities',settings.density_tolerance)]:
            if tolerance is not None:
                values=np.load(root/f'{name}.npy',mmap_mode='r')
                meta.append((values,tolerance,float(values.min())))
    vectors=np.lib.format.open_memmap(output/'index_vectors.npy',mode='w+',dtype=np.float64,shape=(n,len(dimensions)+len(meta)))
    maximum=1.
    for start in range(0,n,8192):
        block=amd[start:start+8192,dimensions]
        parts=[block/radius if radius>0 else block]
        parts.extend(((v[start:start+8192]-lo)/tol)[:,None] for v,tol,lo in meta)
        value=np.column_stack(parts)
        if not np.all(np.isfinite(value)):raise ValueError('Index coordinates overflow')
        vectors[start:start+len(value)]=value
        maximum=max(maximum,float(np.abs(value).max()))
    vectors.flush()
    query_radius=1+64*np.finfo(float).eps*maximum if radius>0 else 0.
    _save_json(output/'index.json',dict(amd_radius=radius,query_radius=query_radius))
    reps=np.lib.format.open_memmap(output/'published_representatives.npy',mode='w+',dtype=np.int64,shape=(n,))
    reps.flush()


def cluster_pdd_batched(store, settings: PDDClusterSettings, output, *, workers=15,
                        batch_size=1024, task_size=32, cache_records=256,
                        checkpoint_seconds=120., progress=None, max_batches=None):
    """Return raw-row representative labels; resume only completed batch checkpoints.

workers is the number of comparison processes, plus one lightweight coordinator.
workers=1 runs in-process. Every batch is finalized in the serial reference order.
Changing worker/batch counts on resume is allowed; settings/data must stay fixed.
The immutable input files must not be modified while any worker is running.
"""
    if any(isinstance(v,bool) or int(v)!=v or v<1 for v in [workers,batch_size,task_size,cache_records]):
        raise ValueError('Worker, batch, task and cache sizes must be positive integers')
    root=Path(store).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=True)
    meta=json.loads((root/'store.json').read_text())
    if _file_identity(meta['pdd_identity']['path'])!=meta['pdd_identity']:
        raise ValueError('PDD source identity changed')
    for filename,expected in meta['array_hashes'].items():
        if _digest(root/filename)!=expected:raise ValueError(f'Changed array store: {filename}')
    identity=dict(version=1,store_sha256=_digest(root/'store.json'),settings=asdict(settings))
    identity_path=out/'identity.json'
    if identity_path.exists():
        if json.loads(identity_path.read_text())!=identity:raise ValueError('Checkpoint settings/data mismatch')
    else:
        _save_json(identity_path,identity)
    # Reject concurrent writers. An unclean kill leaves this small lock for an
    # operator to remove only after confirming the previous process has stopped.
    lock=out/'writer.lock';handle=os.open(lock,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    os.write(handle,f'{os.getpid()}\n'.encode());os.close(handle)
    pool=None
    try:
        if (out/'result.json').exists():
            result=json.loads((out/'result.json').read_text())
            if not result['complete'] or result['labels_sha256']!=_digest(out/'labels.npy'):
                raise ValueError('Completed labels are inconsistent')
            return result
        n=meta['structures'];labels=np.full(n,-1,dtype=np.int64);representatives=[];processed=0
        snapshot_stats=_new_stats();ordered_stats=_new_stats();prior_seconds=0.;batch_count=0
        if (out/'checkpoint.npz').exists():
            with np.load(out/'checkpoint.npz',allow_pickle=False) as cp:
                processed=int(cp['processed']);labels[:processed]=cp['labels'];representatives=cp['representatives'].tolist()
                saved=json.loads(str(cp['metadata']))
                snapshot_stats=saved['snapshot'];ordered_stats=saved['ordered'];prior_seconds=saved['elapsed_s'];batch_count=saved['batches']
            if any(r>=processed or labels[r]!=r for r in representatives):raise ValueError('Invalid checkpoint representatives')
        _make_index(root,out,settings)
        published=np.load(out/'published_representatives.npy',mmap_mode='r+')
        published[:len(representatives)]=representatives;published.flush()
        local=_State(root,out,settings,cache_records)
        if workers>1:
            pool=ProcessPoolExecutor(max_workers=workers,initializer=_worker_init,
                initargs=(str(root),str(out),settings,cache_records),mp_context=multiprocessing.get_context('spawn'))
        started=time.perf_counter();last_checkpoint=started;run_batches=0;reconcile_seconds=0.;parallel_seconds=0.
        def checkpoint():
            value=dict(snapshot=snapshot_stats,ordered=ordered_stats,elapsed_s=prior_seconds+time.perf_counter()-started,batches=batch_count)
            tmp=out/'checkpoint.npz.tmp'
            with tmp.open('wb') as stream:
                np.savez(stream,processed=processed,labels=labels[:processed],representatives=np.asarray(representatives,dtype=np.int64),metadata=np.asarray(json.dumps(value)))
            tmp.replace(out/'checkpoint.npz')
        for start in range(processed,n,batch_size):
            end=min(n,start+batch_size);count=len(representatives);timer=time.perf_counter()
            tasks=[(i,min(end,i+task_size),count) for i in range(start,end,task_size)]
            results=pool.map(_worker_task,tasks,chunksize=1) if pool else map(local.task,tasks)
            matches=np.full(end-start,-1,dtype=np.int64)
            for first,value,stats in results:
                matches[first-start:first-start+len(value)]=value;_add_stats(snapshot_stats,stats)
            parallel_seconds+=time.perf_counter()-timer
            # Workers have finished reading the snapshot before it is extended.
            timer=time.perf_counter();new=local.new_forest();new_count=0;matched_snapshot=int(np.sum(matches>=0))
            for offset,i in enumerate(range(start,end)):
                match=int(matches[offset])
                if match<0:match=local.find(i,new,ordered_stats)
                if match<0:
                    match=i;representatives.append(i);new.add(i);new_count+=1
                labels[i]=match
            published[count:len(representatives)]=representatives[count:]
            reconcile_seconds+=time.perf_counter()-timer
            processed=end;batch_count+=1;run_batches+=1
            event=dict(processed=processed,total=n,representatives=len(representatives),
                batch_new_representatives=new_count,batch_snapshot_matches=matched_snapshot,
                batch_size=end-start,workers=workers,elapsed_s=prior_seconds+time.perf_counter()-started,
                current_run_parallel_wall_s=parallel_seconds,current_run_ordered_wall_s=reconcile_seconds,
                snapshot=snapshot_stats.copy(),ordered=ordered_stats.copy(),batches=batch_count)
            _save_json(out/'progress.json',event)
            if progress is not None:progress(event)
            if time.perf_counter()-last_checkpoint>=checkpoint_seconds or processed==n or (max_batches is not None and run_batches>=max_batches):
                checkpoint();last_checkpoint=time.perf_counter()
            if max_batches is not None and run_batches>=max_batches:break
        if processed<n:return dict(complete=False,processed=processed,representatives=len(representatives))
        order=np.load(root/'order.npy',mmap_mode='r');raw_labels=np.empty(n,dtype=np.int64)
        raw_labels[order]=order[labels];rep_raw=order[representatives]
        if not np.array_equal(raw_labels[rep_raw],rep_raw):raise RuntimeError('Representative self-assignment failed')
        if set(raw_labels)!=set(rep_raw):raise RuntimeError('Unrepresented label')
        if np.any(local.energy[labels]>local.energy):raise RuntimeError('Higher-energy representative')
        np.save(out/'labels.npy',raw_labels);np.save(out/'representatives.npy',rep_raw)
        result=dict(complete=True,structures=n,clusters=len(representatives),settings=asdict(settings),
                    workers=workers,batch_size=batch_size,task_size=task_size,cache_records=cache_records,
                    elapsed_s=prior_seconds+time.perf_counter()-started,snapshot=snapshot_stats,ordered=ordered_stats,
                    current_run_parallel_wall_s=parallel_seconds,current_run_ordered_wall_s=reconcile_seconds,
                    labels_sha256=_digest(out/'labels.npy'),populations_sum=int(np.bincount(raw_labels,minlength=n).sum()),
                    representative_order='energy, identifier',identity=identity)
        _save_json(out/'result.json',result)
        return result
    finally:
        if pool is not None:pool.shutdown(wait=True,cancel_futures=True)
        lock.unlink(missing_ok=True)
