"""Core source package for CSPToolbox.

Exports are loaded lazily so lightweight utilities can run without importing
the full scientific stack immediately.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any


_EXPORTS = {
    "AtomRecord": ("Source.crystal_structure", "AtomRecord"),
    "CifExpansionReport": ("Source.crystal_structure", "CifExpansionReport"),
    "CSORMSymmetrySanityCheck": ("Source.crystal_structure", "CSORMSymmetrySanityCheck"),
    "CrystalStructure": ("Source.crystal_structure", "CrystalStructure"),
    "MoleculeGroup": ("Source.crystal_structure", "MoleculeGroup"),
    "SpaceGroupDetection": ("Source.crystal_structure", "SpaceGroupDetection"),
    "ZMatrixEntry": ("Source.crystal_structure", "ZMatrixEntry"),
    "ZMatrixRepresentation": ("Source.crystal_structure", "ZMatrixRepresentation"),
    "GaussianInputBuilder": ("Source.gaussian_input", "GaussianInputBuilder"),
    "GaussianJobArtifacts": ("Source.gaussian_input", "GaussianJobArtifacts"),
    "GaussianSettings": ("Source.gaussian_input", "GaussianSettings"),
    "CSOFMInputBuilder": ("Source.csofm_input", "CSOFMInputBuilder"),
    "CSOFMJobArtifacts": ("Source.csofm_input", "CSOFMJobArtifacts"),
    "CSOFMMoleculeDefinition": ("Source.csofm_input", "CSOFMMoleculeDefinition"),
    "CSOFMSettings": ("Source.csofm_input", "CSOFMSettings"),
    "read_gaussian_final_energy": ("Source.csofm_input", "read_gaussian_final_energy"),
    "CSORMInputBuilder": ("Source.csorm_input", "CSORMInputBuilder"),
    "CSORMJobArtifacts": ("Source.csorm_input", "CSORMJobArtifacts"),
    "CSORMSettings": ("Source.csorm_input", "CSORMSettings"),
    "FITAtomType": ("Source.mie_typing", "FITAtomType"),
    "InterSpec": ("Source.mie_typing", "InterSpec"),
    "MieAtomType": ("Source.mie_typing", "MieAtomType"),
    "count_atom_types": ("Source.mie_typing", "count_atom_types"),
    "detect_fit_atom_types": ("Source.mie_typing", "detect_fit_atom_types"),
    "detect_mie_atom_types": ("Source.mie_typing", "detect_mie_atom_types"),
    "read_inter_spec": ("Source.mie_typing", "read_inter_spec"),
    "validate_fit_atom_types": ("Source.mie_typing", "validate_fit_atom_types"),
    "validate_inter_atom_types": ("Source.mie_typing", "validate_inter_atom_types"),
    "validate_mie_atom_types": ("Source.mie_typing", "validate_mie_atom_types"),
    "DEFAULT_PBE0_INCAR_TEMPLATE": ("Source.vasp_input", "DEFAULT_PBE0_INCAR_TEMPLATE"),
    "DEFAULT_TPSS_INCAR_TEMPLATE": ("Source.vasp_input", "DEFAULT_TPSS_INCAR_TEMPLATE"),
    "VASP_PRESETS": ("Source.vasp_input", "VASP_PRESETS"),
    "VaspInputBuilder": ("Source.vasp_input", "VaspInputBuilder"),
    "VaspJobArtifacts": ("Source.vasp_input", "VaspJobArtifacts"),
    "VaspSettings": ("Source.vasp_input", "VaspSettings"),
    "OutcarSummary": ("Source.vasp_results", "OutcarSummary"),
    "VaspCalculationHealth": ("Source.vasp_results", "VaspCalculationHealth"),
    "VaspCalculationResult": ("Source.vasp_results", "VaspCalculationResult"),
    "VaspOutSummary": ("Source.vasp_results", "VaspOutSummary"),
    "VaspSystemParser": ("Source.vasp_results", "VaspSystemParser"),
    "VaspSystemResult": ("Source.vasp_results", "VaspSystemResult"),
    "classify_vasp_calculation_dir": ("Source.vasp_results", "classify_vasp_calculation_dir"),
    "parse_outcar_summary": ("Source.vasp_results", "parse_outcar_summary"),
    "parse_vasp_out": ("Source.vasp_results", "parse_vasp_out"),
    "read_contcar_as_crystal": ("Source.vasp_results", "read_contcar_as_crystal"),
    "summarize_system_rows": ("Source.vasp_results", "summarize_system_rows"),
    "DEFAULT_VASP_FILENAMES": ("Source.vasp_file_manifest", "DEFAULT_VASP_FILENAMES"),
    "VaspManifestEntry": ("Source.vasp_file_manifest", "VaspManifestEntry"),
    "VaspManifestSummary": ("Source.vasp_file_manifest", "VaspManifestSummary"),
    "collect_vasp_manifest": ("Source.vasp_file_manifest", "collect_vasp_manifest"),
    "create_tar_from_manifest": ("Source.vasp_file_manifest", "create_tar_from_manifest"),
    "render_summary": ("Source.vasp_file_manifest", "render_summary"),
    "summarize_manifest": ("Source.vasp_file_manifest", "summarize_manifest"),
    "write_manifest_csv": ("Source.vasp_file_manifest", "write_manifest_csv"),
    "PDDDescriptor": ("Source.pdd_descriptor", "PDDDescriptor"),
    "calculate_pdd": ("Source.pdd_descriptor", "calculate_pdd"),
    "pdd_distance": ("Source.pdd_descriptor", "pdd_distance"),
    "pdd_distance_breakdown": ("Source.pdd_descriptor", "pdd_distance_breakdown"),
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str) -> Any:
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    value = getattr(import_module(module_name), attr_name)
    globals()[name] = value
    return value
