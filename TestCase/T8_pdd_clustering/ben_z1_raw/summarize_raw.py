"""Audit all labels, summarize full-run counts and separately labelled anchors."""
import json
import csv
import numpy as np
from cp2_raw_reconstruction import ROOT
from raw_benchmark import raw_arrays, anchors, save

arrays=raw_arrays();mapping=anchors();summaries=[]
for threshold in [.1,.15]:
    report=json.loads((ROOT/f'raw_clusters_{threshold:.3f}.json').read_text())
    # An older runner included this generic metric, but native A... and raw
    # R... labels are disjoint namespaces. Keep only label-invariant metrics.
    report['anchor_metrics'].pop('identical_representative_assignments',None)
    anchor_representatives=set(report['predicted_anchor_labels'].values())
    without_anchor=without_anchor_population=0
    exports=[]
    for b,a in arrays.items():
        labels=np.load(ROOT/f'raw_labels_{threshold:.3f}_A{b}.npy')
        assert labels.shape==(len(a),) and labels.min()>=0 and labels.max()<len(a)
        reps,populations=np.unique(labels,return_counts=True)
        source_lines=(ROOT/f'Analyse_{b}/crystals_stable.out').read_text().splitlines(keepends=True)
        density=np.load(ROOT/f'raw_density_A{b}.npy',mmap_mode='r')
        assert np.all(labels[reps]==reps)
        assert len(reps)==report['batches'][str(b)]['clusters']
        # Fixed representatives must have the lowest energy among their members.
        assert np.all(a[labels,4]<=a[:,4])
        for rep,population in zip(reps,populations):
            name=f'R{b}_{rep+1:06d}'
            assert int(population)==report['batches'][str(b)]['populations'][name]
            exports.append((a[rep,4],int(a[rep,0]),source_lines[rep],dict(
                representative=name,batch=b,original_batch_line=int(rep)+1,
                generation_id=int(a[rep,0]),energy_kj_mol=float(a[rep,4]),
                density_kg_m3=float(density[rep]),population=int(population))))
            if name not in anchor_representatives:
                without_anchor+=1;without_anchor_population+=int(population)
        for name,(batch,index) in mapping.items():
            if batch==b:
                assert report['predicted_anchor_labels'][name]==f'R{b}_{labels[index]+1:06d}'
    stats={key:sum(batch['statistics'][key] for batch in report['batches'].values())
           for key in ['index_candidates','amd_checks','emd_calls','exhaustive_representatives']}
    exports.sort(key=lambda item:(item[0],item[1]))
    # Preserve native record text exactly; do not round reserialized geometry.
    with (ROOT/f'raw_representatives_{threshold:.3f}.out').open('x') as f:
        f.writelines(item[2] for item in exports)
    with (ROOT/f'raw_representatives_{threshold:.3f}.csv').open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(exports[0][3]))
        writer.writeheader();writer.writerows(item[3] for item in exports)
    summaries.append(dict(threshold=threshold,structures=548435,clusters=report['clusters'],
        historical_reference_clusters=2777,clustering_wall_s=report['wall_s'],
        clustering_workers=report.get('workers',1),
        peak_rss_scope=report.get('peak_rss_scope','single_process'),
        peak_rss_bytes=report['peak_rss_bytes'],anchor_metrics=report['anchor_metrics'],
        clusters_without_first_representative_anchor=without_anchor,
        raw_members_in_those_clusters=without_anchor_population,
        per_batch_clusters={b:v['clusters'] for b,v in report['batches'].items()},
        statistics=stats,labels_validated=True))
result=dict(scope='All raw successful minima; replaces both native reduction and COMPACK',
            descriptor_build=json.loads((ROOT/'raw_descriptor_build.json').read_text()),
            reconstruction_validation=json.loads((ROOT/'raw_reconstruction_validation.json').read_text()),
            results=summaries,
            limitations=['Reference pair metrics cover 4,953 saved anchors, not all raw memberships.',
                         'Original four analysis partitions retained; no cross-batch deduplication.',
                         'Historical native executable has not been hash-matched.',
                         'Local times are observed, not a matched CX3 speedup benchmark.',
                         'No production threshold selected.'])
save(ROOT/'raw_results_summary.json',result)
print(json.dumps(summaries,indent=2))
