#!/usr/bin/env python3
"""Summarize PLATON proposals without mistaking alternate settings for new SGs.

Use international numbers and operations per unit volume. A proposed Z-prime is
a geometric estimate based on general-position multiplicity, not phase identity
or an independently verified molecular symmetry assignment.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import re

import spglib


def normalized(symbol):
    return re.sub(r'[\s_]', '', symbol).split(':')[0]


def group_map():
    numbers, orders = defaultdict(set), defaultdict(set)
    for hall in range(1, 531):
        info = spglib.get_spacegroup_type(hall)
        order = len(spglib.get_symmetry_from_database(hall)['rotations'])
        for symbol in (info.international_short, info.international_full):
            numbers[normalized(symbol)].add(info.number)
        orders[info.number].add(order)
    return numbers, orders


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    args = parser.parse_args()
    numbers, orders = group_map()
    overview = {}
    for case in ('Ice_Z4', 'Nicotinamide_Z4'):
        with gzip.open(args.results/(case+'_records.json.gz'), 'rt') as stream:
            records = {r['row']:r for r in json.load(stream)}
        with gzip.open(args.results/(case+'_checks.json.gz'), 'rt') as stream:
            checks = json.load(stream)
        settings, proposals, uncertain = {}, [], []
        for check in checks:
            record = records[check['row']]
            selection = '%s_%s_%s' % (check['variant'], check['tolerance_a'],
                                       'status_zero' if record['status']==0 else 'nonzero_status')
            count = settings.setdefault(selection, Counter())
            count['attempted'] += 1
            if not check['checked']:
                count['check_failed'] += 1
                continue
            count['passed'] += 1
            found = check.get('space_group_number')
            if found is None:
                matches = numbers[normalized(check['space_group'])]
                if len(matches)==1:
                    found = next(iter(matches))
            if found is None:
                count['unresolved_group_name'] += 1
                uncertain.append(check)
                continue
            count['same_international_group' if found==record['sg'] else 'different_international_group'] += 1
            ratio = check.get('input_to_conventional_volume_ratio')
            mult = check.get('reported_multiplicity')
            if mult is None and len(orders[found])==1:
                mult = next(iter(orders[found]))
            if mult is None and found==record['sg'] and ratio is not None and abs(abs(ratio)-1)<0.005:
                mult = record['full_molecules']/4
            if mult is None or ratio is None or ratio==0:
                count['unresolved_symmetry_gain'] += 1
                uncertain.append(check)
                continue
            gain = mult/(record['full_molecules']/4*abs(ratio))
            if gain>1.01:
                count['higher_symmetry_proposed'] += 1
                proposals.append(dict(row=record['row'], status=record['status'], energy=record['energy'],
                    original_group=record['sg'], proposed_group=found, proposed_symbol=check['space_group'],
                    variant=check['variant'], tolerance_a=check['tolerance_a'],
                    operations_per_volume_ratio=gain, suggested_effective_zprime=4/gain))
            elif gain<0.99:
                count['lower_symmetry_reported_needs_review'] += 1
                uncertain.append(check)
            else:
                count['same_symmetry_density'] += 1
        result = dict(settings=settings, proposals=proposals, unresolved=uncertain)
        (args.results/(case+'_interpretation.json')).write_text(json.dumps(result, indent=2)+'\n')
        overview[case] = dict(settings=settings, proposal_rows_any_setting=len({p['row'] for p in proposals}),
            unresolved_checks=len(uncertain), full_atom_proposals_at_005=[p for p in proposals
                if p['variant']=='hydrogen_points' and p['tolerance_a']==0.05])
    (args.results/'interpretation_summary.json').write_text(json.dumps(overview, indent=2)+'\n')
    print(json.dumps({case:result['settings'] for case,result in overview.items()}, indent=2))


if __name__ == '__main__':
    main()
