#!/usr/bin/env python3
"""Collect every row/check from completed chunks and hash the evidence archives."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    args = parser.parse_args()
    result_dir = args.study/'results'
    result_dir.mkdir(exist_ok=True)
    summary = {}
    for case in ('Ice_Z4', 'Nicotinamide_Z4'):
        chunks = sorted((args.study/'runs'/case).glob('*.json.gz'))
        records, checks, archives, protocols = [], [], [], set()
        for path in chunks:
            with gzip.open(str(path), 'rt') as stream:
                chunk = json.load(stream)
            protocols.add(chunk['protocol_sha256'])
            records.extend(chunk['records'])
            checks.extend(chunk['checks'])
            archive = path.with_name(path.name.replace('.json.gz', '.tar.gz'))
            assert digest(archive) == chunk['archive_sha256']
            archives.append(dict(path=str(archive.relative_to(args.study)),
                sha256=chunk['archive_sha256'], bytes=archive.stat().st_size))
        assert len(protocols) == 1
        records.sort(key=lambda row: row['row'])
        checks.sort(key=lambda row: (row['row'], row['variant'], row['tolerance_a']))
        assert [row['row'] for row in records] == list(range(1, 10001)), 'Not all native returns are accounted for'
        keys = [(c['row'], c['variant'], c['tolerance_a']) for c in checks]
        assert len(keys) == len(set(keys))
        assert set(keys) == {(r['row'], v, t) for r in records if r['export_ok']
                            for v in ('ordinary', 'hydrogen_points') for t in (.01, .05, .1)}
        with gzip.open(str(result_dir/(case+'_records.json.gz')), 'wt') as stream:
            json.dump(records, stream)
        with gzip.open(str(result_dir/(case+'_checks.json.gz')), 'wt') as stream:
            json.dump(checks, stream)
        fields = ['row', 'status', 'sg', 'energy', 'variant', 'tolerance_a', 'checked',
                  'space_group', 'space_group_number', 'hall_symbol', 'reported_multiplicity',
                  'input_to_conventional_volume_ratio', 'included_atoms', 'expected_included_atoms',
                  'unitcell_atoms', 'expected_unitcell_atoms', 'process_exit', 'wall_s']
        by_row = {r['row']: r for r in records}
        with gzip.open(str(result_dir/(case+'_checks.csv.gz')), 'wt') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
            writer.writeheader()
            for check in checks:
                writer.writerow(dict(by_row[check['row']], **check))
        summary[case] = dict(returned_rows=len(records), native_status_counts=dict(Counter(r['status'] for r in records)),
            exported_rows=sum(r['export_ok'] for r in records),
            export_failures=[r for r in records if not r['export_ok']],
            attempted_checks=len(checks), passed_checks=sum(c['checked'] for c in checks),
            failed_checks=[c for c in checks if not c['checked']],
            protocol_sha256=next(iter(protocols)),
            maximum_native_expansion_error_a=max(r.get('native_expansion_error_a', 0) for r in records),
            maximum_native_metric_relative_error=max(r.get('native_metric_relative_error', 0) for r in records),
            compressed_evidence_bytes=sum(a['bytes'] for a in archives), archives=archives,
            settings={str(t)+'_'+v:dict(Counter(c['space_group'] for c in checks
                        if c['checked'] and c['variant']==v and c['tolerance_a']==t))
                      for t in (.01, .05, .1) for v in ('ordinary', 'hydrogen_points')})
    (result_dir/'collection_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:{f:d[f] for f in ('returned_rows','exported_rows','attempted_checks','passed_checks','compressed_evidence_bytes')}
                      for k,d in summary.items()}, indent=2))


if __name__ == '__main__':
    main()
