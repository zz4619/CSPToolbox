"""Parallel ordering, restart and uniform-array contracts against serial PDD."""
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
import numpy as np
from Source.pdd_clustering import PDDClusterRecord, PDDClusterSettings, cluster_pdd
from Source.pdd_descriptor import PDDDescriptor
from csptoolbox import prepare_pdd_batch_store, cluster_pdd_batched


class ParallelTests(unittest.TestCase):
    def fixture(self, root, values, energy=None, density=None, ids=None, elements=None):
        values=np.asarray(values,dtype=np.float64);n,rows,k=values.shape
        energy=np.arange(n,dtype=float)*.001 if energy is None else np.asarray(energy)
        density=np.full(n,1000.) if density is None else np.asarray(density)
        ids=[f'id_{i:05d}' for i in range(n)] if ids is None else ids
        elements=['C']*rows if elements is None else elements
        np.save(root/'pdd.npy',values)
        prepare_pdd_batch_store(root/'pdd.npy',energy,density,ids,elements,root/'store',workers=2)
        return [PDDClusterRecord(ids[i],energy[i],density[i],
                    PDDDescriptor(ids[i],k,np.full(rows,1/rows),values[i],tuple(elements),True,False,0.)) for i in range(n)]

    def compare(self, root, records, settings, name, **kwargs):
        expected=cluster_pdd(records,settings)
        result=cluster_pdd_batched(root/'store',settings,root/name,**kwargs)
        self.assertTrue(result['complete'])
        labels=np.load(root/name/'labels.npy')
        actual={r.identifier:records[labels[i]].identifier for i,r in enumerate(records)}
        self.assertEqual(actual,expected.assignments)
        self.assertEqual([records[i].identifier for i in np.load(root/name/'representatives.npy')],expected.representatives)
        self.assertEqual(result['snapshot']['emd_calls']+result['ordered']['emd_calls'],expected.statistics['emd_calls'])
        return result

    def test_chain_across_batches_and_earliest_old_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            # A-B and B-C pass, A-C fails; equal energies exercise ID ordering.
            records=self.fixture(root,np.array([1.,1.08,1.16,1.04,1.24,1.08])[:,None,None],
                                 energy=np.zeros(6),ids=['b','c','d','a','f','e'])
            settings=PDDClusterSettings(.10,block_size=1)
            for batch in [1,2,3,6]:
                self.compare(root,records,settings,f'b{batch}',workers=2,batch_size=batch,task_size=1,cache_records=1)

    def test_random_typed_amd_collisions_gates_and_tree_merges(self):
        rng=np.random.default_rng(102)
        base=np.sort(rng.random((70,6,8)),axis=2)
        values=np.concatenate((base,base+rng.normal(0,.001,base.shape)))
        values=np.sort(np.maximum(values,0),axis=2)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            records=self.fixture(root,values,energy=rng.integers(0,3,len(values))*.5,
                density=1000+rng.integers(0,3,len(values))*10,elements=['C','H','N','C','H','O'])
            for threshold in [.01,.10,.15]:
                settings=PDDClusterSettings(threshold,block_size=4)
                self.compare(root,records,settings,f'r{threshold}',workers=2,batch_size=23,task_size=4,cache_records=3)
            self.compare(root,records,replace(settings,energy_tolerance=None,density_tolerance=None),
                         'nogates',workers=1,batch_size=16)

    def test_checkpoint_resume_changed_workers_and_batch_size(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            records=self.fixture(root,np.repeat(np.arange(24.)[:,None,None],2,axis=0))
            settings=PDDClusterSettings(0.,numeric_slack=0.,energy_tolerance=None,block_size=2)
            partial=cluster_pdd_batched(root/'store',settings,root/'resume',workers=1,batch_size=7,max_batches=2)
            self.assertFalse(partial['complete']);self.assertEqual(partial['processed'],14)
            self.compare(root,records,settings,'resume',workers=2,batch_size=13,task_size=2)
            self.compare(root,records,settings,'resume',workers=1,batch_size=2)
            with self.assertRaises(ValueError):
                cluster_pdd_batched(root/'store',replace(settings,threshold=.1),root/'resume')

    def test_exact_distance_and_strict_metadata_boundaries(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            records=self.fixture(root,np.asarray([1.,1.125,1.25,1.,1.])[:,None,None],
                energy=[0.,.01,.02,1.,.5],density=[1000.,1000.,1000.,1000.,1020.],ids=list('ABCDE'))
            settings=PDDClusterSettings(.125,block_size=1)
            self.assertEqual(cluster_pdd(records,settings).assignments,dict(A='A',B='A',C='C',D='D',E='E'))
            self.compare(root,records,settings,'boundary',workers=2,batch_size=2,task_size=1)

    def test_reject_changed_source_and_invalid_distances(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            self.fixture(root,np.ones((3,2,4)))
            values=np.load(root/'pdd.npy');values[0,0,0]=0.;np.save(root/'pdd.npy',values)
            with self.assertRaises(ValueError):cluster_pdd_batched(root/'store',PDDClusterSettings(.1),root/'run')
            values[0,0,1]=np.nan;np.save(root/'bad.npy',values)
            with self.assertRaises(ValueError):
                prepare_pdd_batch_store(root/'bad.npy',[0,1,2],[1000]*3,['a','b','c'],['C','H'],root/'badstore')


if __name__=='__main__':unittest.main()
