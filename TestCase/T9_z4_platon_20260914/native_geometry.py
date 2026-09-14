#!/usr/bin/env python3
"""Lossless parsing and CIF conversion of the separate native CP2 diagnostic.

Arrays here use row-vector cells. Native printed poses are already rounded;
additional output precision does not recover the original minimizer trajectory.
No clustering, symmetry inference, fitting, or minimization is performed.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np


def parse_native(path):
    record = None
    with Path(path).open() as stream:
        for line in stream:
            f = line.split()
            if not f:
                continue
            if f[0] == 'record':
                if record is not None:
                    raise ValueError('Unterminated native record')
                record = dict(zip(('row','phase','status','sg','nasm','nmol'), map(int,f[1:7])))
                record.update(density=float(f[7]),volume=float(f[8]),asu=[],full=[],zmat=[],symops=[])
            elif f[0] == 'cell':
                record['cell'] = np.array(f[1:],float).reshape(3,3)
            elif f[0] == 'symop':
                record['symops'].append((np.array(f[1:10],int).reshape(3,3).T,np.array(f[10:],float)))
            elif f[0] in ('asu','full'):
                record[f[0]].append((int(f[1]),int(f[2]),f[3],np.array(f[4:],float)))
            elif f[0] == 'zmat':
                record['zmat'].append(tuple(map(int,f[1:6]))+tuple(map(float,f[6:])))
            elif f[0] == 'end_record':
                assert len(record['full']) == len(record['asu'])*record['nmol']//record['nasm']
                assert np.linalg.det(record['cell'])>0
                yield record
                record = None
            else:
                raise ValueError('Unknown export line: '+line[:80])
    if record is not None:
        raise ValueError('Truncated native export')


def symmetry_string(rotation, translation):
    fields=[]
    for row,delta in zip(rotation,translation):
        terms=[]
        for coef,var in zip(row,'xyz'):
            if coef:
                terms.append(('+' if coef>0 else '-')+('' if abs(coef)==1 else str(abs(coef))+'*')+var)
        frac=Fraction(float(delta)).limit_denominator(48)
        if abs(float(frac)-delta)>1e-10:
            raise ValueError('Nonrational symmetry translation')
        if frac:
            terms.append(('+' if frac>0 else '-')+str(abs(frac)))
        fields.append(''.join(terms).lstrip('+') or '0')
    return ','.join(fields)


def cif_text(record, identifier, *, oxygen_only=False, explicit_cell=False):
    cell=record['cell']; lengths=np.linalg.norm(cell,axis=1)
    angles=[np.degrees(np.arccos(np.clip(np.dot(cell[a],cell[b])/lengths[a]/lengths[b],-1,1)))
            for a,b in ((1,2),(0,2),(0,1))]
    operations=record['symops']
    if explicit_cell:
        operations=[(np.eye(3,dtype=int),np.zeros(3))]
    if not operations:
        if record['sg']!=1:
            raise ValueError('Non-P1 ASU export needs native symmetry operations')
        operations=[(np.eye(3,dtype=int),np.zeros(3))]
    lines=['data_'+identifier, "_symmetry_space_group_name_H-M 'P 1'" if explicit_cell or record['sg']==1 else
           '_space_group_IT_number '+str(record['sg'])]
    for name,value in zip(('length_a','length_b','length_c','angle_alpha','angle_beta','angle_gamma'),list(lengths)+angles):
        lines.append('_cell_'+name+' %.12f'%value)
    lines+=['loop_','_space_group_symop_operation_xyz']
    lines += ["'"+symmetry_string(rot,tr)+"'" for rot,tr in operations]
    lines+=['loop_','_atom_site_label','_atom_site_type_symbol','_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z','_atom_site_occupancy']
    inv=np.linalg.inv(cell)
    # build_unit_cell restores exact fixed cell angles/length relations and
    # rebuilds identity copies. Use those copies for a symmetry-consistent ASU,
    # rather than pre-build ASU positions made with rounded printed cell angles.
    atoms=record['full'] if explicit_cell else canonical_asu(record)
    for m,a,label,xyz in atoms:
        element=''.join(c for c in label if c.isalpha()).capitalize()
        if oxygen_only and element!='O':
            continue
        frac=xyz@inv
        lines.append('%s_m%d_a%d %s %.12f %.12f %.12f 1'%(label,m,a,element,*frac))
    return '\n'.join(lines)+'\n'


def canonical_asu(record):
    nops=record['nmol']//record['nasm'] if 'nmol' in record else len(record['symops']) or 1
    return [(1+(m-1)//nops,a,label,xyz) for m,a,label,xyz in record['full'] if (m-1)%nops==0]


def validate_expansion(record):
    nops=len(record['symops'])
    if not nops: raise ValueError('Native symmetry operations are required')
    nat=len(record['asu'])//record['nasm']; inv=np.linalg.inv(record['cell'])
    full=np.array([a[3] for a in record['full']]).reshape(record['nasm'],nops,nat,3)@inv
    expected=np.einsum('sij,maj->msai',np.array([s[0] for s in record['symops']]),full[:,0])
    expected+=np.array([s[1] for s in record['symops']])[None,:,None,:]
    residual=expected-full;residual-=np.rint(residual)
    error=float(np.max(np.linalg.norm(residual@record['cell'],axis=-1)))
    if error>1e-8: raise ValueError('CIF/native full-cell disagreement: %.8g A'%error)
    return error


def prepare(repo, output):
    """Retain native returned-row identity and initial/final pairing."""
    for case in ('Ice_Z4','Nicotinamide_Z4','Water_Z16_P1','Nicotinamide_Z16_P1'):
        source=repo/'test/CP2'/case/'speed/high_z_20260914'
        finals=(source/'crystals.out').read_text().splitlines()
        starts=(source/'starting_crystals.out').read_text().splitlines()
        assert len(finals)==len(starts)
        nt=0 if case in ('Ice_Z4','Water_Z16_P1') else 3
        records=[]; export=[]
        for row,(line,start) in enumerate(zip(finals,starts),1):
            f=line.split(); s=start.split(); nasm=int(f[4]); natt=int(f[5]); k=6+natt+nasm
            assert len(f)==k+12+6+nasm*(6+nt)
            assert len(s)==6+nasm*(6+nt)+3 and int(s[-2])==int(f[0])
            e=list(map(float,f[k:k+12])); final=list(map(float,f[k+12:])); initial=list(map(float,s[1:-2]))
            d=dict(row=row,status=int(f[1]),sg=int(f[0]),nasm=nasm,nt=nt,energy=e[0],uintra=e[1],
                   density=e[9],volume=e[10],cpu_s=e[11],types=list(map(int,f[6+natt:k])),initial=initial,final=final)
            assert np.isfinite(e+initial+final).all()
            records.append(d)
            if 'Z16' in case or d['status']==0:
                phases=[(0,initial),(1,final)] if 'Z16' in case else [(1,final)]
                for phase,pose in phases:
                    export.append('%d %d %d %d %d'%(row,phase,d['status'],d['sg'],nasm))
                    export.append(' '.join('%.17g'%v for v in pose[:6]))
                    for m in range(nasm):
                        export.append(str(d['types'][m])+' '+' '.join('%.17g'%v for v in pose[6+m*(6+nt):6+(m+1)*(6+nt)]))
        dest=output/case; dest.mkdir(parents=True,exist_ok=True)
        (dest/'records.json').write_text(json.dumps(records)+'\n')
        (dest/'export_poses.in').write_text('\n'.join(export)+'\n')
        (dest/'native_inputs_sha256.json').write_text(json.dumps({name:hashlib.sha256((source/name).read_bytes()).hexdigest()
                   for name in ('crystals.out','starting_crystals.out','input.in')},indent=2)+'\n')
        print(case,len(records),sum(d['status']==0 for d in records),flush=True)


def convert(study,case):
    dest=study/case/'cif'; dest.mkdir(exist_ok=True)
    records=[]; errors=[]
    for d in parse_native(study/case/'export_geometry.txt'):
        errors.append(validate_expansion(d))
        name='%s_%05d_%s'%(case,d['row'],'initial' if d['phase']==0 else 'final')
        path=dest/(name+'.cif'); path.write_text(cif_text(d,name,oxygen_only=case=='Ice_Z4'))
        records.append({'row':d['row'],'phase':d['phase'],'status':d['status'],'sg':d['sg'],'path':str(path.relative_to(study))})
    (study/case/'cif_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    (study/case/'cif_validation.json').write_text(json.dumps({'records':len(records),
        'maximum_native_expansion_error_A':max(errors),'clustering':False,
        'asu_source':'identity copies from native build_unit_cell after restoring fixed SG cell parameters'},indent=2)+'\n')
    print(case,'CIFs',len(records),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','convert'))
    p.add_argument('--repo',type=Path);p.add_argument('--study',type=Path,required=True);p.add_argument('--case')
    a=p.parse_args()
    if a.action=='prepare': prepare(a.repo,a.study)
    else: convert(a.study,a.case)
