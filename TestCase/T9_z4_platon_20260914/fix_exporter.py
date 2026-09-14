#!/usr/bin/env python3
"""Fix only the separate diagnostic exporter to use the run's named SG setting."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    dest = Path.home()/'CSP_studies/z4_platon_20260914/software'
    source = dest/'export_driver.f90'
    original = dest/'export_driver_original.f90'
    assert not original.exists()
    shutil.copyfile(str(source), str(original))
    old = '''    do i=1,NSpaceSupported
        if(SpaceSupported(i)%Number==export_group) ispace=i
    enddo
    if(ispace==0) error stop 'Unsupported export space group'
    group=SpaceSupported(ispace)'''
    new = '''    ! Original diagnostic lookup selected the last supported setting by number:
    ! do i=1,NSpaceSupported
    !     if(SpaceSupported(i)%Number==export_group) ispace=i
    ! enddo
    ! group=SpaceSupported(ispace)
    ! The run input selects a named setting. Multiple names can share a number.
    do i=1,NSpaceList
        if(SpaceList(i)%Number==export_group) then
            if(ispace/=0) error stop 'Ambiguous exported space group setting'
            ispace=i
        endif
    enddo
    if(ispace==0) error stop 'Export group absent from run SpaceList'
    group=SpaceList(ispace)'''
    text = source.read_text()
    assert text.count(old) == 1
    source.write_text(text.replace(old, new))
    command = json.loads((dest/'build_command.json').read_text())
    command = [str(source) if c.endswith('/export_driver.f90') else
               str(dest/'export_fixed.x') if c.endswith('/export.x') else c for c in command]
    assert str(source) in command and str(dest/'export_fixed.x') in command
    (dest/'fixed_build_command.json').write_text(json.dumps(command, indent=2)+'\n')
    with (dest/'fixed_build.log').open('wb') as stream:
        subprocess.run(command, cwd=str(dest), stdout=stream, stderr=subprocess.STDOUT, check=True)
    shutil.move(str(dest/'export.x'), str(dest/'export_legacy.x'))
    shutil.copy2(str(dest/'export_fixed.x'), str(dest/'export.x'))
    report = dict(original_driver_sha256=sha(original), fixed_driver_sha256=sha(source),
                  original_binary_sha256=sha(dest/'export_legacy.x'), fixed_binary_sha256=sha(dest/'export.x'),
                  scientific_source_modified=False,
                  correction='Use unique run-specific SpaceList setting, not last SpaceSupported numeric match')
    (dest/'exporter_fix.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
