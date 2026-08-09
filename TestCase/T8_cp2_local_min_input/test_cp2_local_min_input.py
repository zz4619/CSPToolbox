"""Focused tests for the CrystalPredictor2 local-minimisation adapter."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.cp2_local_min import (  # noqa: E402
    CP2_SUPPORTED_SPACE_GROUPS,
    ComponentMapping,
    MappingMetrics,
    build_cp2_local_min_batch,
    collect_cp2_local_min_status,
    discover_global_search_dir,
    parse_cp2_input,
    parse_cp2_lam_topology,
    prepare_cp2_local_min_inputs,
    read_structure,
    resolve_cp2_space_group,
)
from Source.CLI_scripts.prepare_cp2_local_min_batch import (  # noqa: E402
    main as prepare_batch_main,
)


SINGLE_TYPE_INPUT = """MOLECULAR TYPES 1
-----------------------------------------------------------------------------
TYPE Synthetic
1
4
flexible_lam_intra
1
dih4 -180.0 180.0
-----------------------------------------------------------------------------
SPACE_GROUPS
P1
END_SPACE_GROUPS
"""


SYNTHETIC_LAM = """Intramolecular energy/gradients/hessian/charges for:

Synthetic

Across 1 dimensional grid:

dih4 -180.0 30.0 180.0

From starting Z-matrix:
C1  1 -0.10 0 bnd1 0 ang1 0 dih1
N1  2  0.20 1 bnd2 0 ang2 0 dih2
O1  3 -0.30 2 bnd3 1 ang3 0 dih3
H1  4  0.20 3 bnd4 2 ang4 1 dih4
"""


CANONICAL_ZMATRIX = """Z-matrix for molecule 1
C1
N1 C1
O1 N1 C1
H1 O1 N1 C1
"""


REFERENCE_RES = """TITL canonical reference
CELL 1.0 20.0 20.0 20.0 90.0 90.0 90.0
LATT -1
SFAC C N O H
C1 1 0.100000 0.100000 0.100000 11.0 0.03
N1 2 0.170000 0.100000 0.100000 11.0 0.03
O1 3 0.235000 0.115000 0.100000 11.0 0.03
H1 4 0.270000 0.135000 0.100000 11.0 0.03
END
"""


EXPERIMENTAL_RES = """TITL shuffled experimental
REM SPACE_GROUP P 1
CELL 1.0 20.0 20.0 20.0 90.0 90.0 90.0
LATT -1
SFAC C N O H
HX 4 0.265000 0.320000 0.100000 11.0 =
 0.03 0.03 0.03 0.00 0.00 0.00
OX 3 0.285000 0.285000 0.100000 11.0 0.03
CX 1 0.300000 0.150000 0.100000 11.0 0.03
NX 2 0.300000 0.220000 0.100000 11.0 0.03
END
"""


EXPERIMENTAL_ZPRIME2_RES = """TITL two shuffled experimental molecules
REM SPACE_GROUP P 1
CELL 1.0 20.0 20.0 20.0 90.0 90.0 90.0
LATT -1
SFAC C N O H
HXA 4 0.265000 0.320000 0.100000 11.0 0.03
OXA 3 0.285000 0.285000 0.100000 11.0 0.03
CXA 1 0.300000 0.150000 0.100000 11.0 0.03
NXA 2 0.300000 0.220000 0.100000 11.0 0.03
HXB 4 0.665000 0.720000 0.500000 11.0 0.03
OXB 3 0.685000 0.685000 0.500000 11.0 0.03
CXB 1 0.700000 0.550000 0.500000 11.0 0.03
NXB 2 0.700000 0.620000 0.500000 11.0 0.03
END
"""


MULTI_TYPE_INPUT = """MOLECULAR TYPES 2
-----------------------------------------------------------------------------
TYPE First
2
4
first_lam
0
-----------------------------------------------------------------------------
TYPE Second
1
2
second_lam
0
-----------------------------------------------------------------------------
SPACE_GROUPS
P1 P-1
END_SPACE_GROUPS
"""


class CP2LocalMinInputTests(unittest.TestCase):
    def test_discovers_current_and_legacy_global_search_names(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            current = root / "5_Globalsearch"
            current.mkdir()
            self.assertEqual(current, discover_global_search_dir(root))

            current.rmdir()
            legacy = root / "5_GlobSrch"
            legacy.mkdir()
            self.assertEqual(legacy, discover_global_search_dir(root))

            current.mkdir()
            with self.assertRaisesRegex(ValueError, "Ambiguous"):
                discover_global_search_dir(root)

    def test_cp2_input_parser_marks_unvalidated_scope(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "input.in"
            path.write_text(MULTI_TYPE_INPUT, encoding="utf-8")

            parsed = parse_cp2_input(path)

            self.assertEqual(2, parsed.molecular_type_count)
            self.assertEqual(3, parsed.asymmetric_unit_molecule_count)
            self.assertFalse(parsed.structurally_supported)
            self.assertEqual(
                "experimental_unvalidated_zprime_gt1_or_multicomponent",
                parsed.supported_scope,
            )

    def test_repeated_lam_site_types_use_unique_system_zmatrix_labels(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            lam = root / "flexible_lam_intra"
            zmat = root / "Zmatrix"
            lam.write_text(
                """From starting Z-matrix:
C1 1 -0.1
C1 2 -0.1 1 bnd2
O1 3 -0.3 2 bnd3 1 ang3
H1 4  0.5 3 bnd4 2 ang4 1 dih4
""",
                encoding="utf-8",
            )
            zmat.write_text(
                """Z-matrix for molecule 1
C1
C2 C1
O1 C2 C1
H1 O1 C2 C1
""",
                encoding="utf-8",
            )

            topology = parse_cp2_lam_topology(lam, 4, canonical_zmatrix_path=zmat,)

            self.assertEqual(("C1", "C2", "O1", "H1"), topology.labels)
            self.assertEqual(
                ["C1", "C1", "O1", "H1"], [site.site_type for site in topology.sites]
            )
            self.assertFalse(topology.labels_generated)

    def test_shelx_latt_symm_space_group_is_inferred(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "pna21.res"
            path.write_text(
                """TITL ordinary SHELX file
CELL 1.0 10.0 11.0 12.0 90.0 90.0 90.0
LATT -1
SYMM -x,-y,1/2+z
SYMM 1/2+x,1/2-y,z
SYMM 1/2-x,1/2+y,1/2+z
SFAC C
C1 1 0.1 0.2 0.3 11.0 0.03
END
""",
                encoding="utf-8",
            )

            self.assertEqual("Pna21", read_structure(path).space_group)

    def test_space_group_is_normalized_to_cp2_spelling(self) -> None:
        self.assertEqual(67, len(CP2_SUPPORTED_SPACE_GROUPS))
        self.assertNotIn("P2/m", CP2_SUPPORTED_SPACE_GROUPS)
        self.assertNotIn("Pcab", CP2_SUPPORTED_SPACE_GROUPS)
        self.assertEqual("P1", resolve_cp2_space_group("P 1"))
        self.assertEqual("P21/c", resolve_cp2_space_group("P 21/c"))

    def test_globally_inverted_mapping_is_not_automatically_validated(self) -> None:
        common = dict(
            method="reference_graph_zmatrix_internal_rmsd",
            heavy_mapping_count=2,
            candidates_truncated=False,
            gross_bond_angle_mismatches=0,
            internal_coordinate_score=0.2,
            heavy_kabsch_rmsd_angstrom=0.02,
            second_best_score_gap=1.0,
            mapping_ambiguous=False,
            reference_order_assumed=False,
        )
        self.assertTrue(MappingMetrics(torsion_orientation="same", **common).validation_safe)
        self.assertFalse(
            MappingMetrics(
                torsion_orientation="globally_inverted", **common
            ).validation_safe
        )
        self.assertTrue(
            MappingMetrics(
                torsion_orientation="planar_inversion_equivalent", **common
            ).validation_safe
        )
        self.assertEqual(
            ("torsion_orientation_not_same",),
            MappingMetrics(
                torsion_orientation="globally_inverted", **common
            ).validation_failures,
        )

    def test_globally_inverted_mapping_cannot_be_waived_for_runnable_job(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            system = root / "SyntheticSystem"
            global_search = system / "5_GlobSrch"
            global_search.mkdir(parents=True)
            (global_search / "input.in").write_text(
                SINGLE_TYPE_INPUT, encoding="utf-8"
            )
            (global_search / "flexible_lam_intra").write_text(
                SYNTHETIC_LAM, encoding="utf-8"
            )
            (global_search / "potential.in").write_text(
                "authoritative potential\n", encoding="utf-8"
            )
            (system / "Zmatrix").write_text(
                CANONICAL_ZMATRIX, encoding="utf-8"
            )
            reference = root / "reference.res"
            experimental = root / "experimental.res"
            reference.write_text(REFERENCE_RES, encoding="utf-8")
            experimental.write_text(EXPERIMENTAL_RES, encoding="utf-8")
            executable = root / "Minimise"
            executable.write_text("#!/bin/bash\nexit 0\n", encoding="utf-8")
            executable.chmod(0o750)
            inverted_metrics = MappingMetrics(
                method="reference_graph_zmatrix_internal_rmsd",
                heavy_mapping_count=1,
                candidates_truncated=False,
                gross_bond_angle_mismatches=0,
                internal_coordinate_score=0.1,
                heavy_kabsch_rmsd_angstrom=0.01,
                second_best_score_gap=None,
                torsion_orientation="globally_inverted",
                mapping_ambiguous=False,
                reference_order_assumed=False,
            )
            inverted_mapping = ComponentMapping(
                molecular_type_index=1,
                occurrence_index=1,
                experimental_component_index=1,
                canonical_to_experimental=(
                    ("C1", 2),
                    ("N1", 3),
                    ("O1", 1),
                    ("H1", 0),
                ),
                metrics=inverted_metrics,
            )

            with patch(
                "Source.cp2_local_min.match_experimental_to_cp2",
                return_value=(inverted_mapping,),
            ):
                with self.assertRaisesRegex(
                    ValueError, "globally inverted.*cannot waive"
                ):
                    prepare_cp2_local_min_inputs(
                        system,
                        experimental,
                        system / "Local_Min_CP2" / "INVERTED",
                        reference_paths={1: reference},
                        stage_mode="copy",
                        cp2_executable=executable,
                        allow_unvalidated_mapping=True,
                    )

    def test_centrosymmetric_primitive_res_is_not_silently_p1(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "pinv.res"
            path.write_text(
                """TITL inversion only
CELL 1.0 10.0 10.0 10.0 90.0 90.0 90.0
LATT 1
SFAC C
C1 1 0.1 0.2 0.3 11.0 0.03
END
""",
                encoding="utf-8",
            )

            self.assertEqual("P-1", read_structure(path).space_group)

    def test_preparation_reorders_coordinates_and_preserves_potential_bytes(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            system = root / "SyntheticSystem"
            global_search = system / "5_GlobSrch"
            global_search.mkdir(parents=True)
            (global_search / "input.in").write_text(SINGLE_TYPE_INPUT, encoding="utf-8")
            (global_search / "flexible_lam_intra").write_text(
                SYNTHETIC_LAM, encoding="utf-8",
            )
            potential_bytes = b"CP2 potential is authoritative\r\nsecond line\x00\n"
            (global_search / "potential.in").write_bytes(potential_bytes)
            (system / "Zmatrix").write_text(CANONICAL_ZMATRIX, encoding="utf-8")
            reference = root / "reference.res"
            experimental = root / "experimental.res"
            reference.write_text(REFERENCE_RES, encoding="utf-8")
            experimental.write_text(EXPERIMENTAL_RES, encoding="utf-8")
            output = system / "Local_Min_CP2" / "TEST01"
            executable = root / "Minimise"
            executable.write_text("#!/bin/bash\nexit 0\n", encoding="utf-8")
            executable.chmod(0o750)

            artifacts = prepare_cp2_local_min_inputs(
                system,
                experimental,
                output,
                reference_paths={1: reference},
                stage_mode="copy",
                compack_metadata={"rmsd_1_angstrom": 0.08},
                cp2_executable=executable,
            )

            self.assertEqual(potential_bytes, artifacts.potential_path.read_bytes())
            atom_lines = [
                line
                for line in artifacts.expcrys_pdb_path.read_text(
                    encoding="utf-8"
                ).splitlines()
                if line.startswith("HETATM")
            ]
            self.assertEqual(
                ["C1", "N1", "O1", "H1"], [line[12:16].strip() for line in atom_lines]
            )
            self.assertEqual(
                [6.0, 6.0, 5.7, 5.3], [float(line[30:38]) for line in atom_lines]
            )
            self.assertTrue(all(len(line) >= 78 for line in atom_lines))
            cryst1 = artifacts.expcrys_pdb_path.read_text(
                encoding="utf-8"
            ).splitlines()[1]
            self.assertEqual("P1", cryst1[55:66].strip())

            with artifacts.mapping_tsv_path.open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                ["CX", "NX", "OX", "HX"], [row["experimental_label"] for row in rows]
            )
            manifest = json.loads(artifacts.manifest_path.read_text(encoding="utf-8"))
            self.assertFalse(manifest["scientifically_validated"])
            self.assertTrue(manifest["mapping_validated"])
            self.assertTrue(manifest["potential"]["byte_identical"])
            self.assertEqual(0.08, manifest["compack_metadata"]["rmsd_1_angstrom"])
            self.assertEqual(
                ["C1", "N1", "O1", "H1"],
                manifest["molecular_types"][0]["canonical_atom_order"],
            )
            self.assertIsNotNone(artifacts.pbs_script_path)
            self.assertIsNotNone(artifacts.executable_path)
            self.assertEqual(
                executable.read_bytes(), artifacts.executable_path.read_bytes()
            )
            script = artifacts.pbs_script_path.read_text(encoding="utf-8")
            self.assertIn("select=1:ncpus=1", script)
            self.assertIn("cp -L", script)
            self.assertIn("PBS_JOBID", script)
            self.assertIn("imkl/2022.1.0", script)
            self.assertIn("module --ignore_cache load", script)
            self.assertNotIn("impi/", script)
            self.assertNotIn("$EBROOTNAGLIB/scripts/nagvars.sh", script)
            self.assertIn(
                'export NAG_KUSARI_FILE="${HOME}"/.nag/license.dat', script
            )
            self.assertIn('! -f "$NAG_KUSARI_FILE"', script)
            self.assertIn('! -r "$NAG_KUSARI_FILE"', script)
            self.assertNotIn("/sw-eb/software/NAGlib/license", script)
            self.assertNotIn("NAG key begin", script)
            self.assertEqual(
                "$HOME/.nag/license.dat",
                manifest["execution"]["pbs_settings"]["nag_license_file"],
            )
            self.assertNotIn(
                "nagvars_command", manifest["execution"]["pbs_settings"]
            )
            self.assertIn("OMP_NUM_THREADS=1", script)
            self.assertIn("MKL_NUM_THREADS=1", script)
            self.assertIn("MKL_DYNAMIC=FALSE", script)
            self.assertIn("exit_status.txt", script)
            self.assertIn("cp2_Minimise", script)
            self.assertNotIn(str(executable.resolve()), script)
            subprocess.run(["bash", "-n", str(artifacts.pbs_script_path)], check=True)

    def test_single_type_occurrence_count_can_follow_experimental_zprime(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            system = root / "SyntheticSystem"
            global_search = system / "5_GlobSrch"
            global_search.mkdir(parents=True)
            source_input = global_search / "input.in"
            source_input.write_text(SINGLE_TYPE_INPUT, encoding="utf-8")
            (global_search / "flexible_lam_intra").write_text(
                SYNTHETIC_LAM, encoding="utf-8"
            )
            potential = global_search / "potential.in"
            potential.write_bytes(b"authoritative potential\r\n")
            (system / "Zmatrix").write_text(
                CANONICAL_ZMATRIX, encoding="utf-8"
            )
            reference = root / "reference.res"
            reference.write_text(REFERENCE_RES, encoding="utf-8")
            experimental = root / "zprime2.res"
            experimental.write_text(EXPERIMENTAL_ZPRIME2_RES, encoding="utf-8")

            artifacts = prepare_cp2_local_min_inputs(
                system,
                experimental,
                root / "case",
                reference_paths={1: reference},
                stage_mode="copy",
                auto_single_type_occurrences=True,
                allow_unvalidated=True,
            )

            self.assertEqual(2, parse_cp2_input(artifacts.input_path).asymmetric_unit_molecule_count)
            self.assertEqual(1, parse_cp2_input(source_input).asymmetric_unit_molecule_count)
            self.assertEqual(b"authoritative potential\r\n", artifacts.potential_path.read_bytes())
            manifest = json.loads(artifacts.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(
                {"source_occurrences": 1, "experimental_occurrences": 2},
                manifest["input"]["single_type_occurrence_update"],
            )
            self.assertFalse(manifest["input"]["byte_identical"])
            self.assertEqual(
                "experimental_unvalidated_zprime_gt1_or_multicomponent",
                manifest["supported_scope"],
            )

    def test_audited_canonical_experimental_labels_can_be_trusted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            system = root / "SyntheticSystem"
            global_search = system / "5_GlobSrch"
            global_search.mkdir(parents=True)
            (global_search / "input.in").write_text(
                SINGLE_TYPE_INPUT, encoding="utf-8"
            )
            (global_search / "flexible_lam_intra").write_text(
                SYNTHETIC_LAM, encoding="utf-8"
            )
            (global_search / "potential.in").write_text(
                "authoritative potential\n", encoding="utf-8"
            )
            (system / "Zmatrix").write_text(
                CANONICAL_ZMATRIX, encoding="utf-8"
            )
            experimental = root / "audited.res"
            experimental.write_text(REFERENCE_RES, encoding="utf-8")

            artifacts = prepare_cp2_local_min_inputs(
                system,
                experimental,
                root / "case",
                reference_paths={1: experimental},
                stage_mode="copy",
                trust_experimental_labels=True,
            )

            manifest = json.loads(artifacts.manifest_path.read_text(encoding="utf-8"))
            self.assertTrue(manifest["trusted_experimental_labels"])
            self.assertEqual(
                "trusted_experimental_labels",
                manifest["mappings"][0]["metrics"]["method"],
            )
            self.assertTrue(manifest["mapping_validated"])

    def test_structured_result_marker_is_collected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            job = Path(temporary_directory)
            (job / "Minimisation_log.out").write_text(
                """VOLUME_PML_MAX: 1500. Consider increasing it if you get ifail=-4 errors.
IFAIL 1 0
Utot = -120.0 -110.0 10.0
energy failed at a recoverable trial point
CP2_LOCAL_MIN_RESULT_V1 status=OPTIMIZER_CONVERGED info=0 utot_kj_mol_entity=-1.0125D+02 uvdw_kj_mol_entity=-8.0D+01 uelec_kj_mol_entity=-2.3D+01 uintra_kj_mol_entity=7.5D-01 volume_per_asu_a3=2.5D+02 density_kg_m3=1.2D+03
""",
                encoding="utf-8",
            )
            (job / "minimise_final_structure.pdb").write_text(
                """CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P1
HETATM    1 C1   UNK     0       1.000   2.000   3.000  1.00  0.00           C
END
""",
                encoding="utf-8",
            )
            (job / "exit_status.txt").write_text("0\n", encoding="utf-8")

            status = collect_cp2_local_min_status(job)

            self.assertEqual("optimizer_converged", status.status)
            self.assertTrue(status.optimizer_converged)
            self.assertFalse(status.local_minimum_confirmed)
            self.assertEqual("CP2_LOCAL_MIN_RESULT_V1", status.result_schema)
            self.assertEqual("OPTIMIZER_CONVERGED", status.optimizer_status)
            self.assertEqual(0, status.optimizer_info)
            self.assertEqual((0,), status.ifail_history)
            self.assertAlmostEqual(-101.25, status.energies_kj_mol["Utot"]["final"])
            self.assertAlmostEqual(-120.0, status.energies_kj_mol["Utot"]["initial"])
            self.assertAlmostEqual(1200.0, status.result_values["density_kg_m3"])
            self.assertTrue((job / "cp2_local_min_status.json").is_file())

    def test_prepared_cases_are_assembled_into_one_resumable_batch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            cases = []
            executable_bytes = b"#!/bin/bash\nexit 0\n"
            executable_hash = hashlib.sha256(executable_bytes).hexdigest()
            pbs_settings = {
                "walltime": "24:00:00",
                "memory_gb": 8,
                "modules": ["tools/prod", "imkl/2022.1.0"],
                "nag_license_file": "$HOME/.nag/license.dat",
                "queue": None,
            }
            for refcode in ("FORM01", "FORM02"):
                case = root / "cases" / "Synthetic" / refcode
                case.mkdir(parents=True)
                (case / "cp2_Minimise").write_bytes(executable_bytes)
                (case / "cp2_Minimise").chmod(0o750)
                (case / "run_cp2_local_min.pbs").write_text(
                    "#!/bin/bash\nexit 0\n", encoding="utf-8"
                )
                (case / "run_cp2_local_min.pbs").chmod(0o750)
                (case / "cp2_local_min_manifest.json").write_text(
                    json.dumps(
                        {
                            "system_name": "Synthetic",
                            "refcode": refcode,
                            "supported_scope": "supported_single_component_zprime1",
                            "structurally_supported": True,
                            "mapping_validated": True,
                            "execution": {
                                "staged_cp2_executable": {"sha256": executable_hash},
                                "pbs_settings": pbs_settings,
                            },
                        }
                    ),
                    encoding="utf-8",
                )
                cases.append(case)

            artifacts = build_cp2_local_min_batch(
                root,
                cases,
                per_case_timeout="15m",
            )

            with artifacts.cases_manifest_path.open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(2, len(rows))
            self.assertEqual(
                ["Synthetic__FORM01", "Synthetic__FORM02"],
                [row["case_id"] for row in rows],
            )
            script = artifacts.pbs_script_path.read_text(encoding="utf-8")
            self.assertIn("select=1:ncpus=1:mem=8gb", script)
            self.assertIn("timeout --signal=TERM --kill-after=30s 15m", script)
            self.assertIn("already_converged", script)
            self.assertIn("timed_out", script)
            self.assertNotIn("\nqsub ", script)
            subprocess.run(["bash", "-n", str(artifacts.pbs_script_path)], check=True)

    def test_inventory_prepares_every_case_and_one_batch_job(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            system = root / "SyntheticSystem"
            global_search = system / "5_GlobSrch"
            global_search.mkdir(parents=True)
            (global_search / "input.in").write_text(
                SINGLE_TYPE_INPUT, encoding="utf-8"
            )
            (global_search / "flexible_lam_intra").write_text(
                SYNTHETIC_LAM, encoding="utf-8"
            )
            (global_search / "potential.in").write_text(
                "authoritative potential\n", encoding="utf-8"
            )
            (system / "Zmatrix").write_text(
                CANONICAL_ZMATRIX, encoding="utf-8"
            )
            reference = root / "reference.res"
            experimental = root / "experimental.res"
            reference.write_text(REFERENCE_RES, encoding="utf-8")
            experimental.write_text(EXPERIMENTAL_RES, encoding="utf-8")
            executable = root / "Minimise"
            executable.write_text("#!/bin/bash\nexit 0\n", encoding="utf-8")
            executable.chmod(0o750)
            inventory = root / "inventory.tsv"
            with inventory.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=(
                        "system_name",
                        "refcode",
                        "system_dir",
                        "experimental_structure",
                        "reference_1",
                    ),
                    delimiter="\t",
                )
                writer.writeheader()
                for refcode in ("FORM01", "FORM02"):
                    writer.writerow(
                        {
                            "system_name": "SyntheticSystem",
                            "refcode": refcode,
                            "system_dir": system,
                            "experimental_structure": experimental,
                            "reference_1": reference,
                        }
                    )
            output = root / "batch"

            return_code = prepare_batch_main(
                [
                    str(inventory),
                    str(output),
                    "--cp2-executable",
                    str(executable),
                ]
            )

            self.assertEqual(0, return_code)
            with (output / "cp2_local_min_preparation_status.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                statuses = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                ["prepared_runnable", "prepared_runnable"],
                [row["preparation_status"] for row in statuses],
            )
            self.assertTrue((output / "run_all_cp2_local_min.pbs").is_file())
            self.assertTrue((output / "cp2_local_min_batch_cases.tsv").is_file())

    def test_nonzero_optimizer_info_is_failed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            job = Path(temporary_directory)
            (job / "Minimisation_log.out").write_text(
                "CP2_LOCAL_MIN_RESULT_V1 status=OPTIMIZER_FAILED info=6\n",
                encoding="utf-8",
            )

            status = collect_cp2_local_min_status(job, write_json=False)

            self.assertEqual("failed", status.status)
            self.assertTrue(status.finished)
            self.assertFalse(status.local_minimum_confirmed)
            self.assertEqual(6, status.optimizer_info)


if __name__ == "__main__":
    unittest.main()
