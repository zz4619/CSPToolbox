"""Tests for VASP OUTCAR force parsing."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source.vasp_results import parse_outcar_forces  # noqa: E402


SYNTHETIC_OUTCAR = """
   VRHFIN =C: s2p2
   VRHFIN =H: ultrasoft test
   number of dos      NEDOS =    301   number of ions     NIONS =      3
   ions per type =               1   2

 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      0.00000      0.00000      0.00000         1.000000      0.000000      0.000000
      1.00000      0.00000      0.00000         0.000000      2.000000      0.000000
      0.00000      1.00000      0.00000         0.000000      0.000000      3.000000
 -----------------------------------------------------------------------------------
    total drift:                                1.000000      2.000000      3.000000

 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      0.00000      0.00000      0.00000         0.300000      0.400000      0.000000
      1.00000      0.00000      0.00000         0.000000     -0.120000      0.000000
      0.00000      1.00000      0.00000        -0.010000      0.000000      0.000000
 -----------------------------------------------------------------------------------
    total drift:                                0.290000      0.280000      0.000000
"""


class OutcarForceParsingTests(unittest.TestCase):
    def test_parse_final_force_block_by_default(self) -> None:
        outcar_path = _write_synthetic_outcar()

        parsed = parse_outcar_forces(outcar_path)

        self.assertTrue(parsed.exists)
        self.assertEqual(3, parsed.nions)
        self.assertEqual(("C", "H"), parsed.ion_types)
        self.assertEqual((1, 2), parsed.ion_counts)
        self.assertEqual(2, parsed.force_block_count)
        self.assertEqual(1, parsed.retained_force_block_count)

        block = parsed.final_block
        self.assertIsNotNone(block)
        assert block is not None
        self.assertEqual(2, block.ionic_step)
        self.assertEqual(3, block.atom_count)
        self.assertEqual(("C", "H", "H"), tuple(atom.element for atom in block.atoms))
        self.assertEqual((0.29, 0.28, 0.0), block.total_drift)
        self.assertAlmostEqual(0.5, block.max_force_norm_ev_ang)
        self.assertAlmostEqual(math.sqrt((0.5**2 + 0.12**2 + 0.01**2) / 3), block.rms_force_norm_ev_ang)

    def test_parse_all_force_blocks_when_requested(self) -> None:
        outcar_path = _write_synthetic_outcar()

        parsed = parse_outcar_forces(outcar_path, all_ionic_steps=True)

        self.assertEqual(2, parsed.force_block_count)
        self.assertEqual(2, parsed.retained_force_block_count)
        self.assertEqual((1, 2), tuple(block.ionic_step for block in parsed.blocks))
        self.assertAlmostEqual(3.0, parsed.blocks[0].max_force_norm_ev_ang)
        self.assertAlmostEqual(0.5, parsed.blocks[1].max_force_norm_ev_ang)


def _write_synthetic_outcar() -> Path:
    tmpdir = tempfile.TemporaryDirectory()
    path = Path(tmpdir.name) / "OUTCAR"
    path.write_text(SYNTHETIC_OUTCAR, encoding="utf-8")
    _TEMP_DIRS.append(tmpdir)
    return path


_TEMP_DIRS: list[tempfile.TemporaryDirectory[str]] = []


if __name__ == "__main__":
    unittest.main()
