"""Separate reconstruction effects from metadata/representative changes."""
import csv
import json
import time
import numpy as np
from cp2_raw_reconstruction import ROOT,descriptor
from raw_benchmark import anchors,raw_arrays,save
from benchmark import partition_metrics
from Source.pdd_clustering import PDDClusterRecord,PDDClusterSettings,cluster_pdd

mapping=anchors();arrays=raw_arrays()
rows=list(csv.DictReader((ROOT/'manifest.csv').open()))
reference={r['id']:r['reference'] for r in rows}
expected=json.loads((ROOT/'pdd_only_0.100.json').read_text())['assignments']
reports=[]
for use_native_averages in [True,False]:
    predictions={};start=time.perf_counter()
    for b,a in arrays.items():
        d=np.load(ROOT/f'raw_pdd_A{b}.npy',mmap_mode='r')
        rho=np.load(ROOT/f'raw_density_A{b}.npy',mmap_mode='r')
        records=[]
        for row in rows:
            batch,i=mapping[row['id']]
            if batch!=b:continue
            energy=float(row['energy']) if use_native_averages else a[i,4]
            density=float(row['density']) if use_native_averages else rho[i]
            records.append(PDDClusterRecord(row['id'],energy,density,descriptor(d[i])))
        predictions.update(cluster_pdd(records,PDDClusterSettings(.1)).assignments)
    if use_native_averages:
        assert predictions==expected, 'Reconstructed geometry changed anchor-only assignments'
    reports.append(dict(native_average_metadata=use_native_averages,
                        structures=len(predictions),threshold=.1,
                        exactly_matches_exported_geometry_assignments=predictions==expected,
                        metrics=partition_metrics(reference,predictions),
                        seconds=time.perf_counter()-start))
    print(json.dumps(reports[-1]),flush=True)
save(ROOT/'raw_metadata_controls.json',dict(controls=reports))
