#!/usr/bin/env python3
"""Native CP2 terminal export and exhaustive PLATON checks, without minimization.

Each input return is preserved, including failed minimizations. Native export
failures and PLATON failures are explicit results. Work and full logs are archived
by chunk; completed chunks are resumable. No clustering or input rewriting occurs.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile
import time

import numpy as np
from native_geometry import parse_native, validate_expansion, canonical_asu, symmetry_string

CASES = ('Ice_Z4', 'Nicotinamide_Z4')
TOLS = (0.01, 0.05, 0.1)
VARIANTS = ('ordinary', 'hydrogen_points')
CHECKED_GROUPS = set()


def validate_group(native):
    rotations = np.array([r for r, _ in native['symops']])
    translations = np.array([t for _, t in native['symops']])
    key = (rotations.tobytes(), translations.tobytes())
    if key not in CHECKED_GROUPS:
        for first, shift in zip(rotations, translations):
            for second, delta in zip(rotations, translations):
                product = first @ second
                moved = first @ delta+shift
                residual = moved-translations
                residual -= np.rint(residual)
                matching = np.all(rotations == product, axis=(1, 2)) & (np.max(np.abs(residual), axis=1) < 1e-9)
                if not matching.any():
                    raise ValueError('Native symmetry operations are not closed under composition')
        CHECKED_GROUPS.add(key)
    metric = native['cell'] @ native['cell'].T
    errors = [np.max(np.abs(r.T @ metric @ r-metric)) for r in rotations]
    relative = float(max(errors)/np.max(np.abs(metric)))
    if relative > 1e-10:
        raise ValueError('Native symmetry rotations do not preserve the cell metric')
    return relative


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def poses(records):
    lines = []
    for r in records:
        p, nt, n = r['final'], r['nt'], r['nasm']
        assert n == 4 and len(p) == 6+n*(6+nt)
        assert np.isfinite(p).all() and all(t == 1 for t in r['types'])
        lines.append('%d 1 %d %d %d' % (r['row'], r['status'], r['sg'], n))
        lines.append(' '.join('%.17g' % v for v in p[:6]))
        for m in range(n):
            lines.append('1 '+' '.join('%.17g' % v for v in p[6+m*(6+nt):6+(m+1)*(6+nt)]))
    return '\n'.join(lines)+'\n'


def export_records(study, case, records, work):
    """Recover complete prefix records and isolate failures without skipping rows."""
    todo = list(records)
    attempt = 0
    while todo:
        dest = work/('export_%03d' % attempt)
        dest.mkdir()
        inputs = study/'inputs'/case
        for name in ('input.in', 'potential.in', 'rigid_lam_intra', 'flexible_lam_intra', 'Zmatrix'):
            if (inputs/name).exists():
                shutil.copyfile(str(inputs/name), str(dest/name))
        (dest/'export_poses.in').write_text(poses(todo))
        started = time.monotonic()
        with (dest/'stdout.txt').open('wb') as log:
            try:
                p = subprocess.run([str(study/'software/export.x')], cwd=str(dest),
                                   stdout=log, stderr=subprocess.STDOUT, timeout=180)
                exit_code = p.returncode
            except subprocess.TimeoutExpired:
                exit_code = 'timeout'
        count = 0
        parse_error = None
        if (dest/'export_geometry.txt').exists():
            try:
                for native in parse_native(dest/'export_geometry.txt'):
                    original = todo[count]
                    assert (native['row'], native['status'], native['sg']) == (
                        original['row'], original['status'], original['sg'])
                    count += 1
                    yield original, native, None
            except (ValueError, AssertionError, IndexError) as exc:
                parse_error = str(exc) or type(exc).__name__
        (dest/'process.json').write_text(json.dumps(dict(exit_code=exit_code,
            parsed_records=count, parse_error=parse_error, wall_s=time.monotonic()-started))+'\n')
        todo = todo[count:]
        if todo:
            # Retry a failed batch head alone before attributing failure to it.
            if len(todo) > 1:
                first = todo.pop(0)
                isolated = dest/'isolated'
                isolated.mkdir()
                for result in export_records(study, case, [first], isolated):
                    yield result
            else:
                yield todo.pop(0), None, dict(reason='native_export_failed',
                    exit_code=exit_code, parse_error=parse_error, log=str(dest.relative_to(work)))
        attempt += 1


def cif(native, name, hydrogen_points):
    cell = native['cell']
    lengths = np.linalg.norm(cell, axis=1)
    angles = [np.degrees(np.arccos(np.clip(np.dot(cell[a], cell[b])/lengths[a]/lengths[b], -1, 1)))
              for a, b in ((1, 2), (0, 2), (0, 1))]
    lines = ['data_'+name, '_space_group_IT_number '+str(native['sg'])]
    for key, val in zip(('length_a', 'length_b', 'length_c', 'angle_alpha', 'angle_beta', 'angle_gamma'),
                        list(lengths)+angles):
        lines.append('_cell_'+key+' %.12f' % val)
    lines += ['loop_', '_space_group_symop_operation_xyz']
    lines += ["'"+symmetry_string(rot, tr)+"'" for rot, tr in native['symops']]
    lines += ['loop_',
              '_atom_site_label', '_atom_site_type_symbol', '_atom_site_fract_x',
              '_atom_site_fract_y', '_atom_site_fract_z', '_atom_site_occupancy']
    inv = np.linalg.inv(cell)
    heavy = 0
    atoms = canonical_asu(native)
    for i, (_, _, label, xyz) in enumerate(atoms, 1):
        element = re.match('[A-Za-z]+', label).group().capitalize()
        assert element in ('C', 'N', 'O', 'H')
        heavy += element != 'H'
        if hydrogen_points and element == 'H':
            element = 'F'  # F is absent from both real systems; keep C/N/O distinct.
        point = xyz @ inv
        lines.append('%s%d %s %.12f %.12f %.12f 1' % (element, i, element, *point))
    included = len(atoms) if hydrogen_points else heavy
    centering = sum(np.array_equal(rot, np.eye(3, dtype=int)) for rot, _ in native['symops'])
    assert centering in (1, 2, 3, 4) and len(native['symops']) % centering == 0
    # PLATON's "Unitcell" is the reduced primitive-cell atom set (PLA060).
    return '\n'.join(lines)+'\n', included, included*len(native['symops'])//centering


def platon_check(exe, input_text, expected, expected_cell, dest, tolerance, timeout):
    dest.mkdir()
    (dest/'input.cif').write_text(input_text)
    commands = 'CALC ADDSYM EXACT 0.1 {0} {0} {0}\nEND\n'.format(tolerance)
    (dest/'commands.txt').write_text(commands)
    started = time.monotonic()
    with (dest/'stdout.txt').open('wb') as out, (dest/'stderr.txt').open('wb') as err:
        try:
            p = subprocess.run([str(exe), '-o', 'input.cif'], input=commands.encode(),
                cwd=str(dest), stdout=out, stderr=err, timeout=timeout)
            exit_code = p.returncode
        except subprocess.TimeoutExpired:
            exit_code = 'timeout'
    text = (dest/'stdout.txt').read_text(errors='replace')
    groups = re.findall(r'SpaceGroup\s*=\s*(\S+)', text)
    if not groups:
        groups = [g.replace(' ', '') for g in re.findall(r'Space Group\s+H-M:\s+(.+?)\s+Laue:', text)]
    hall = re.findall(r'Space Group Hall:\s+(.+?)\s+\[Schoenflies:', text)
    included = re.findall(r'Number of Input Atoms Included in Search\s+(\d+)\s+\(Unitcell\s+(\d+)\)', text)
    n, ncell = map(int, included[-1]) if included else (None, None)
    numbers = re.findall(r'Multiplicity:\s+\d+\([^)]*\), No:\s*(\d+)', text)
    multiplicity = re.findall(r'Multiplicity:\s+(\d+)\([^)]*\), No:', text)
    determinant = re.findall(r'Det\(T\)\s*\n[^\n]*\n[^\n]*?\)\s*([-+0-9.]+)', text)
    deleted = [line.strip() for line in text.splitlines() if 'DELETED from INPUT STREAM' in line]
    return dict(process_exit=exit_code, wall_s=time.monotonic()-started,
        checked=exit_code == 0 and 'NORMAL END of PLATON' in text and len(groups) == 1 and n == expected and ncell == expected_cell,
        space_group=groups[0] if len(groups) == 1 else None, groups=groups,
        hall_symbol=hall[0].strip() if len(hall) == 1 else None,
        space_group_number=int(numbers[-1]) if numbers else None, deleted_atom_messages=deleted,
        reported_multiplicity=int(multiplicity[-1]) if multiplicity else None,
        input_to_conventional_volume_ratio=float(determinant[-1]) if determinant else None,
        included_atoms=n, expected_included_atoms=expected, unitcell_atoms=ncell,
        expected_unitcell_atoms=expected_cell, source_sha256=digest(dest/'input.cif'))


def worker(task):
    study_name, case, chunk, records, timeout, output_name, protocol_sha = task
    study = Path(study_name)
    output = study/output_name/case
    output.mkdir(parents=True, exist_ok=True)
    name = '%04d' % chunk
    result_path = output/(name+'.json.gz')
    if result_path.exists():
        with gzip.open(str(result_path), 'rt') as stream:
            data = json.load(stream)
        assert data['rows'] == [r['row'] for r in records]
        assert data['protocol_sha256'] == protocol_sha, 'Cannot resume a different scientific protocol'
        return dict(case=case, chunk=chunk, resumed=True, records=len(records))
    temp = Path(tempfile.mkdtemp(prefix='cp2-platon-'+case+'-'+name+'-'))
    data = dict(case=case, chunk=chunk, protocol_sha256=protocol_sha,
                rows=[r['row'] for r in records], records=[], checks=[])
    started = time.monotonic()
    try:
        for original, native, error in export_records(study, case, records, temp):
            record = {k: original[k] for k in ('row', 'status', 'sg', 'energy')}
            if error:
                record.update(export_ok=False, error=error)
            else:
                try:
                    if not np.isfinite(native['cell']).all() or not np.isfinite(
                            np.array([atom[3] for atom in native['full']])).all():
                        raise ValueError('Native export contains non-finite cell or atom coordinates')
                    record.update(export_ok=True, full_atoms=len(native['full']), full_molecules=native['nmol'],
                                  native_metric_relative_error=validate_group(native),
                                  native_expansion_error_a=validate_expansion(native))
                except (ValueError, AssertionError) as exc:
                    record.update(export_ok=False, error=dict(reason='native_expansion_invalid', detail=str(exc)))
            data['records'].append(record)
            if not record['export_ok']:
                continue
            for variant in VARIANTS:
                input_text, expected, expected_cell = cif(native, '%s_%05d' % (case, original['row']), variant == 'hydrogen_points')
                for tolerance in TOLS:
                    dest = temp/('%05d_%s_%s' % (original['row'], variant, tolerance))
                    checked = platon_check(study/'software/platon', input_text, expected, expected_cell, dest, tolerance, timeout)
                    checked.update(row=original['row'], variant=variant, tolerance_a=tolerance)
                    data['checks'].append(checked)
        assert [r['row'] for r in data['records']] == data['rows']
        data['wall_s'] = time.monotonic()-started
        archive = output/(name+'.tar.gz')
        with tarfile.open(str(archive)+'.tmp', 'w:gz', compresslevel=1) as tar:
            tar.add(str(temp), arcname=name)
        os.replace(str(archive)+'.tmp', str(archive))
        data['archive_sha256'] = digest(archive)
        with gzip.open(str(result_path)+'.tmp', 'wt') as stream:
            json.dump(data, stream)
        os.replace(str(result_path)+'.tmp', str(result_path))
        return dict(case=case, chunk=chunk, records=len(records), wall_s=data['wall_s'],
                    export_errors=sum(not r['export_ok'] for r in data['records']),
                    checks=len(data['checks']), check_errors=sum(not r['checked'] for r in data['checks']))
    finally:
        shutil.rmtree(str(temp))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study', type=Path, required=True)
    p.add_argument('--workers', type=int, default=32)
    p.add_argument('--chunk-size', type=int, default=100)
    p.add_argument('--timeout', type=float, default=30)
    p.add_argument('--pilot', action='store_true')
    p.add_argument('--output-name', default=None)
    args = p.parse_args()
    output_name = args.output_name or ('pilot' if args.pilot else 'runs')
    output = args.study/output_name
    output.mkdir(exist_ok=True)
    protocol = dict(schema=1, tolerances_a=TOLS, angle_tolerance_deg=0.1,
        variants=VARIANTS, timeout_s=args.timeout, numpy_version=np.__version__,
        platon_sha256=digest(args.study/'software/platon'), exporter_sha256=digest(args.study/'software/export.x'),
        script_sha256=digest(__file__), parser_sha256=digest(Path(__file__).with_name('native_geometry.py')),
        input_sha256={str(f.relative_to(args.study)):digest(f) for f in sorted((args.study/'inputs').glob('*/*')) if f.is_file()})
    protocol_text = json.dumps(protocol, sort_keys=True)
    protocol_sha = hashlib.sha256(protocol_text.encode()).hexdigest()
    protocol_path = output/'protocol.json'
    if protocol_path.exists():
        assert json.loads(protocol_path.read_text()) == json.loads(protocol_text)
    else:
        protocol_path.write_text(protocol_text+'\n')
        scripts = output/'scripts'
        scripts.mkdir()
        for source in (Path(__file__), Path(__file__).with_name('native_geometry.py')):
            shutil.copyfile(str(source), str(scripts/source.name))
    tasks = []
    for case in CASES:
        records = json.loads((args.study/'inputs'/case/'records.json').read_text())
        assert len(records) == 10000 and [r['row'] for r in records] == list(range(1, 10001))
        if args.pilot:
            # Exercise different groups plus the observed optimizer-status classes.
            selected = {}
            for r in records:
                key = ('sg', r['sg']) if r['status'] == 0 else ('status', r['status'])
                selected.setdefault(key, r)
            records = sorted(selected.values(), key=lambda r: r['row'])
        for start in range(0, len(records), args.chunk_size):
            tasks.append((str(args.study.resolve()), case, start//args.chunk_size,
                records[start:start+args.chunk_size], args.timeout, output_name, protocol_sha))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(worker, task) for task in tasks]
        for future in as_completed(futures):
            print(json.dumps(future.result()), flush=True)


if __name__ == '__main__':
    main()
