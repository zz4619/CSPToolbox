"""Historical Nicotinamide Z'=1 reconstruction, with explicit narrow contracts.

Mirrors CP2 Analyse's LAM geometry, cell constraints and Euler convention.
This study adapter is not a general CP2 reader. No energies are recalculated.
The symmetry dump comes from the repository's unmodified Fortran group table.
"""
from pathlib import Path
from itertools import product
import hashlib
import re
import numpy as np
from scipy.spatial import KDTree
from ase.geometry import minkowski_reduce
from Source.pdd_descriptor import PDDDescriptor

ROOT = Path(__file__).resolve().parent
PI = 3.14159265358979
ELEMENTS = tuple('C C N C C C C N O H H H H H H'.split())
MASS = 6*12.0107 + 2*14.0067 + 15.9994 + 6*1.00794
EXPECTED_MODEL_INPUTS = {
    'input.in': 'c0fe5518860be46de66b279a0c6734484ec7109a909548bcae091a6f111eed7a',
    'flexible_lam_intra': '494821de3f401c08dff0ae99caee06a64e89259a02e16058883845670ad7f940',
    'cp2_groups.txt': 'fb91a774602cbcbe03eab1a07063bf4b125017d405077605505150acb4ee10b7',
}


def read_groups(path=ROOT/'cp2_groups.txt'):
    lines = iter(path.read_text().splitlines())
    groups = {}
    for line in lines:
        fields = line.split()
        assert fields[0] == 'GROUP', line
        number, count = int(fields[2]), int(fields[3])
        constraints = np.array(next(lines).split(), float)
        operations = np.array([next(lines).split() for _ in range(count)], float)
        groups.setdefault(number, (constraints, operations[:, :9].reshape(-1,3,3).transpose(0,2,1), operations[:, 9:]))
    return groups


class HistoricalLAM:
    def __init__(self, path=ROOT/'flexible_lam_intra'):
        # This adapter hard-codes the audited smoothing/torsion conventions.
        # Fail explicitly for any other model instead of silently reusing them.
        for filename,expected in EXPECTED_MODEL_INPUTS.items():
            source=path if filename=='flexible_lam_intra' else ROOT/filename
            if hashlib.sha256(source.read_bytes()).hexdigest()!=expected:
                raise ValueError(f'Unsupported historical model input: {filename}')
        header, *entries = re.split(r'^ENTRY ', path.read_text(), flags=re.M)
        atomlines = re.findall(r'^([CNHO]\d+\s+\d+\s+.*)$', header, flags=re.M)
        assert len(atomlines) == 15 and len(entries) == 426
        self.refs = np.full((15,3), -1, int)
        for atom, line in enumerate(atomlines):
            fields = line.split()
            for col in range(min(atom, 3)):
                self.refs[atom,col] = int(fields[3+col*2])-1
        starting = dict((n, float(v)) for n,v in re.findall(r'^((?:bnd|ang|dih)\d+)\s+([-\d.]+)\s*$',header,re.M))
        self.names = list(starting)
        assert len(self.names) == 39 and self.names[-3:] == ['dih8','dih10','dih11']
        self.locations = [(int(n[3:])-1, ['bnd','ang','dih'].index(n[:3])) for n in self.names]
        values, slopes = [], []
        scale = np.array([1/.529177249 if n.startswith('bnd') else 1 for n in self.names])
        for entry in entries:
            ztext, htext = entry.split('The Hessian matrix',1)
            rows = re.findall(r'^((?:bnd|ang|dih)\d+)\s+(\S+)\s+(\S+)',ztext,re.M)
            assert [r[0] for r in rows] == self.names
            value = np.array([float(r[1]) for r in rows])
            for i,n in enumerate(self.names[:36]):
                if n.startswith('dih'):
                    if starting[n]-value[i] > 180: value[i] += 360
                    elif starting[n]-value[i] < -180: value[i] -= 360
            hessian = np.full((39,39), np.nan)
            columns = None
            for line in htext.splitlines():
                fields = line.split()
                if not fields or fields[0] not in starting:
                    continue
                if len(fields) > 1 and fields[1] in starting:
                    columns = [self.names.index(n) for n in fields]
                else:
                    assert columns is not None
                    row = self.names.index(fields[0])
                    nums = [float(x) for x in fields[1:]]
                    for col,v in zip(columns,nums):
                        hessian[row,col] = hessian[col,row] = v
            assert np.isfinite(hessian).all()
            hessian *= scale[:,None]*scale[None,:]
            slope = -np.linalg.solve(hessian[:36,:36],hessian[:36,36:])
            slope *= np.array([1 if n.startswith('bnd') else 180/PI for n in self.names[:36]])[:,None]
            values.append(value)
            slopes.append(slope)
        self.values, self.slopes = np.array(values), np.array(slopes)
        self.centres = self.values[:,36:] * PI/180
        self.groups = read_groups()

    def zmatrix(self, torsions):
        torsions = np.array(torsions, copy=True)
        lo = np.array([-60.,-60.,90.])*PI/180
        hi = np.array([300.,300.,270.])*PI/180
        # Native Analyse's one-step guard for rounded output at input bounds.
        torsions[torsions > hi] -= .001
        torsions[torsions < lo] += .001
        diff = torsions[None,:] - self.centres
        diff[:,:2] = (diff[:,:2]+PI) % (2*PI)-PI
        weights = np.exp(-np.sum(diff**2,axis=1)/.4**2)
        assert weights.sum() > 0
        weights /= weights.sum()
        dependent = weights @ self.values[:,:36] + np.einsum('l,li,ldi->d',weights,diff,self.slopes)
        all_values = np.r_[dependent,torsions*180/PI]
        z = np.zeros((15,3))
        for (atom,col),value in zip(self.locations,all_values):
            z[atom,col] = value if col == 0 else value*PI/180
        return z

    def reconstruct(self, record):
        assert record.shape == (21,) and record[2] == record[3] == 1
        constraints, rotations, translations = self.groups[int(record[1])]
        a,b,c,alpha,beta,gamma = record[6:12]
        if constraints[1] == constraints[0]: b=a
        if constraints[2] == constraints[0]: c=a
        if constraints[2] == constraints[1]: c=b
        angles = np.array([alpha,beta,gamma])
        fixed = constraints[3:] > .001
        angles[fixed] = constraints[3:][fixed]
        ca,cb,cg = np.cos(angles)
        sg = np.sin(angles[2])
        determinant = 1+2*ca*cb*cg-ca*ca-cb*cb-cg*cg
        assert determinant > 0
        cell = np.array([[a,0,0],[b*cg,b*sg,0],[c*cb,c*(ca-cb*cg)/sg,c*np.sqrt(determinant)/sg]])
        z = self.zmatrix(record[18:])
        body = np.zeros((15,3))
        body[1,0] = z[1,0]
        assert list(self.refs[2,:2]) == [1,0]
        body[2,:2] = [z[1,0]-z[2,0]*np.cos(z[2,1]),-z[2,0]*np.sin(z[2,1])]
        for atom in range(3,15):
            bond,angle,torsion = z[atom]
            bidx,aidx,didx = self.refs[atom]
            u, v = body[bidx]-body[aidx], body[didx]-body[aidx]
            u /= np.linalg.norm(u)
            v /= np.linalg.norm(v)
            normal = np.cross(u,v)/np.sqrt(max(1-np.dot(u,v)**2,1e-7))
            inplane = np.cross(normal,u)
            body[atom] = body[bidx]+bond*(inplane*np.sin(angle)*np.cos(torsion)+normal*np.sin(angle)*np.sin(torsion)-u*np.cos(angle))
        cp,ct,cs = np.cos(record[15:18]); sp,st,ss = np.sin(record[15:18])
        rotation = np.array([[cp*cs-sp*ct*ss,sp*cs+cp*ct*ss,st*ss],
                             [-cp*ss-sp*ct*cs,-sp*ss+cp*ct*cs,st*cs],
                             [sp*st,-cp*st,ct]])
        asu = (body@rotation.T) @ np.linalg.inv(cell) + record[12:15]
        fractional = np.einsum('aj,sij->sai',asu,rotations) + translations[:,None,:]
        motif = fractional.reshape(-1,3)@cell
        density = MASS*len(rotations)/np.linalg.det(cell)*1e27/6.02214199e23
        return motif, cell, density


def asu_pdd(motif, cell, *, k=100):
    """Query all 15 ASU sites against the complete periodic cell.

    Exact symmetry means each ASU row appears NSymOp times in a full-cell PDD;
    retain equal 1/15 weights, without an approximate row-collapse operation.
    Caller must supply the complete general-position orbit, identity first.
    """
    cell, _ = minkowski_reduce(cell)
    inverse = np.linalg.inv(cell)
    wrapped = np.mod(motif@inverse,1.)@cell
    query = wrapped[:15]
    reciprocal_lengths = np.linalg.norm(inverse,axis=0)
    radius = max(1e-6,1.2*(3*(k+1)*abs(np.linalg.det(cell))/(4*np.pi*len(motif)))**(1/3))
    while True:
        bounds = np.ceil(radius*reciprocal_lengths+1e-12).astype(int)
        shifts = np.array(list(product(*(range(-b,b+1) for b in bounds))))@cell
        cloud = (shifts[:,None,:]+wrapped[None,:,:]).reshape(-1,3)
        distances = KDTree(cloud).query(query,k=k+1,workers=1)[0]
        kth = float(distances[:,-1].max())
        if np.isfinite(kth) and kth <= radius:
            return distances[:,1:]
        radius = max(radius*1.5,kth*(1+1e-12)) if np.isfinite(kth) else radius*2


def descriptor(distances, name='raw'):
    return PDDDescriptor(name,distances.shape[1],np.full(15,1/15),distances,ELEMENTS,True,False,0.)
