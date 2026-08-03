"""Focused tests for the CrystalPredictor2 local-minimisation adapter."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.cp2_local_min import (  # noqa: E402
    CP2_SUPPORTED_SPACE_GROUPS,
    MappingMetrics,
    collect_cp2_local_min_status,
    discover_global_search_dir,
    parse_cp2_input,
    parse_cp2_lam_topology,
    prepare_cp2_local_min_inputs,
    read_structure,
    resolve_cp2_space_group,
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

    def test_unknown_shelx_space_group_fails_without_override(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "ambiguous.res"
            path.write_text(
                """TITL ordinary SHELX file
CELL 1.0 10.0 11.0 12.0 90.0 100.0 90.0
LATT 1
SYMM -X,0.5+Y,0.5-Z
SFAC C
C1 1 0.1 0.2 0.3 11.0 0.03
END
""",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "pass space_group explicitly"):
                read_structure(path)
            overridden = read_structure(path, space_group="P 21/c")
            self.assertEqual("P 21/c", overridden.space_group)

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

            self.assertEqual("P -1", read_structure(path).space_group)

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
            self.assertIn("NAG_KUSARI_FILE", script)
            self.assertIn("OMP_NUM_THREADS=1", script)
            self.assertIn("MKL_NUM_THREADS=1", script)
            self.assertIn("MKL_DYNAMIC=FALSE", script)
            self.assertIn("exit_status.txt", script)
            self.assertIn("cp2_Minimise", script)
            self.assertNotIn(str(executable.resolve()), script)
            subprocess.run(["bash", "-n", str(artifacts.pbs_script_path)], check=True)

    def test_structured_result_marker_is_collected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            job = Path(temporary_directory)
            (job / "Minimisation_log.out").write_text(
                """IFAIL 1 0
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
            self.assertAlmostEqual(-101.25, status.energies_kj_mol["Utot"]["final"])
            self.assertAlmostEqual(-120.0, status.energies_kj_mol["Utot"]["initial"])
            self.assertAlmostEqual(1200.0, status.result_values["density_kg_m3"])
            self.assertTrue((job / "cp2_local_min_status.json").is_file())

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
