#!/usr/bin/env python3
"""Run explicit, zero-nonfit ADDSYM checks; compatible with CX3 Python 3.6.

Ordinary O/H CIFs validate only non-H coordinates because PLATON omits H.
The H-to-C point-label proxy is an auxiliary geometric test, not a chemical CIF.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--platon', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    reports = []
    for source in sorted(args.inputs.glob('*/*.cif')):
        if source.name not in ('full_OH.cif','oxygen_only.cif','hydrogen_proxy.cif'):
            continue
        for angle in (.01,.1,1.):
            for tolerance in (.001,.005,.01,.02,.05,.1,.2):
                name = source.parent.name+'__'+source.stem+'__'+str(angle)+'__'+str(tolerance)
                work = args.output/name
                work.mkdir()
                shutil.copyfile(str(source),str(work/'input.cif'))
                commands = 'CALC ADDSYM EXACT {a} {d} {d} {d}\nEND\n'.format(a=angle,d=tolerance)
                (work/'commands.txt').write_text(commands)
                started = time.monotonic()
                p = subprocess.run([str(args.platon.resolve()),'-o','input.cif'],input=commands.encode(),
                                   stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=str(work),timeout=30)
                (work/'stdout.txt').write_bytes(p.stdout)
                (work/'stderr.txt').write_bytes(p.stderr)
                stdout = p.stdout.decode(errors='replace')
                group = re.findall(r'SpaceGroup\s*=\s*(\S+)',stdout)
                # A detected increase uses a different native PLATON report block.
                if not group:
                    group = [g.replace(' ','') for g in re.findall(
                        r'Space Group\s+H-M:\s+(.+?)\s+Laue:',stdout)]
                included = re.findall(r'Number of Input Atoms Included in Search\s+(\d+)',stdout)
                success = p.returncode==0 and 'NORMAL END of PLATON' in stdout and len(group)==1 and bool(included)
                report = dict(name=name,structure=source.parent.name,variant=source.stem,
                    angle_tolerance_deg=angle,distance_tolerances_a=[tolerance]*3,
                    process_exit=p.returncode,checked=success,wall_s=time.monotonic()-started,
                    space_group=group[0] if group else None,included_atoms=int(included[-1]) if included else None,
                    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
                reports.append(report)
                if not success:
                    raise RuntimeError(json.dumps(report))
    (args.output/'summary.json').write_text(json.dumps(dict(schema=1,platon_sha256=
        hashlib.sha256(args.platon.read_bytes()).hexdigest(),runs=reports),indent=2)+'\n')
    print(json.dumps(dict(runs=len(reports),all_checked=all(r['checked'] for r in reports),
                         groups=sorted(set(r['space_group'] for r in reports)))))


if __name__ == '__main__':
    main()
