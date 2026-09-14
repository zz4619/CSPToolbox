#!/usr/bin/env python3
"""Freeze the already validated CX3 exporter, PLATON and native inputs."""
import hashlib
import json
from pathlib import Path
import shutil


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    studies = Path.home()/'CSP_studies'
    target = studies/'z4_platon_20260914'
    target.mkdir(exist_ok=False)
    (target/'software').mkdir()
    (target/'scripts').mkdir()
    source = studies/'cp2_high_z_scientific_20260914'
    manifest = []
    copies = [(source/'native_export/export.x', target/'software/export.x'),
              (source/'native_export/export_driver.f90', target/'software/export_driver.f90'),
              (source/'native_export/build_command.json', target/'software/build_command.json'),
              (studies/'water_z16_symmetry_20260914/software/platon/platon', target/'software/platon')]
    for case in ('Ice_Z4', 'Nicotinamide_Z4'):
        (target/'inputs'/case).mkdir(parents=True)
        for name in ('input.in', 'potential.in', 'rigid_lam_intra', 'flexible_lam_intra', 'Zmatrix', 'records.json'):
            if (source/case/name).exists():
                copies.append((source/case/name, target/'inputs'/case/name))
    for src, dst in copies:
        before = sha(src)
        shutil.copy2(str(src), str(dst))
        assert before == sha(dst) == sha(src), 'Source changed while freezing '+str(src)
        manifest.append(dict(source=str(src), path=str(dst.relative_to(target)),
                             sha256=before, bytes=dst.stat().st_size))
    assert b"'symop'" in (target/'software/export_driver.f90').read_bytes()
    (target/'input_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(dict(study=str(target), files=len(manifest))))


if __name__ == '__main__':
    main()
