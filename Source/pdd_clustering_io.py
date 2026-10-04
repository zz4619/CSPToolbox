"""Manifest ingestion and content-validated PDD caches for clustering."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import numpy as np
from .pdd_descriptor import PDDDescriptor, calculate_pdd
from .pdd_clustering import PDDClusterRecord

__all__ = ['load_pdd_cluster_manifest']


def load_pdd_cluster_manifest(manifest, *, k=100, cache=None):
    """CSV columns: id,path,energy,density. Paths are relative to the manifest.

    Cache identity includes file contents, manifest contents and descriptor
    settings/version. No pickle or installation of optional CCDC is needed.
    """
    manifest = Path(manifest).resolve()
    with manifest.open(newline='') as stream:
        entries = list(csv.DictReader(stream))
    required = {'id', 'path', 'energy', 'density'}
    if not entries or not required.issubset(entries[0]):
        raise ValueError("Manifest must have id,path,energy,density and at least one row.")
    paths = [(manifest.parent / entry['path']).resolve() for entry in entries]
    checksum = hashlib.sha256(manifest.read_bytes())
    checksum.update(f'certified-pdd-v1;k={k};typed=True;collapse_tol=0.0001'.encode())
    for path in paths:
        checksum.update(path.read_bytes())
    identity = checksum.hexdigest()
    cached = None
    if cache is not None and Path(cache).exists():
        cached = np.load(cache, allow_pickle=False)
        if str(cached['identity']) != identity:
            cached.close()
            raise ValueError("Descriptor cache does not match these inputs/settings; use a new cache path.")
    records, arrays = [], {}
    if cached is None:
        from Source.crystal_structure import CrystalStructure
    for i, (entry, path) in enumerate(zip(entries, paths)):
        if cached is None:
            structure = CrystalStructure.from_file(path).expand_to_explicit_unit_cell()
            descriptor = calculate_pdd(structure, k=k, workers=1)
            arrays[f'w{i}'] = descriptor.weights
            arrays[f'd{i}'] = descriptor.distances
            arrays[f'e{i}'] = np.asarray(descriptor.center_elements)
        else:
            descriptor = PDDDescriptor(entry['id'], k, cached[f'w{i}'], cached[f'd{i}'],
                                       tuple(cached[f'e{i}'].tolist()), True, True, 1e-4)
        records.append(PDDClusterRecord(entry['id'], float(entry['energy']), float(entry['density']), descriptor))
    if cached is not None:
        cached.close()
    elif cache is not None:
        with Path(cache).open('xb') as stream:
            np.savez_compressed(stream, identity=np.asarray(identity), **arrays)
    return records, identity
