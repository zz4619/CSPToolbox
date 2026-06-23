from __future__ import annotations

import csv
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Source import CrystalStructure, pdd_distance


DEFAULT_INPUT_ROOT = PROJECT_ROOT.parent / "Calculations" / "01_HA_CSORM"
DEFAULT_OUTPUT_CSV = DEFAULT_INPUT_ROOT / "csorm_res_parsing_emd_k100.csv"


def main() -> int:
    input_root = DEFAULT_INPUT_ROOT
    output_csv = DEFAULT_OUTPUT_CSV

    rows: list[dict[str, object]] = []
    for system_dir in sorted(path for path in input_root.iterdir() if path.is_dir()):
        system_name = system_dir.name
        res_path = system_dir / f"{system_name}.res"
        if not res_path.is_file():
            continue

        row: dict[str, object] = {
            "system_name": system_name,
            "res_path": str(res_path),
            "status": "",
            "normal_atom_count": "",
            "csorm_atom_count": "",
            "emd_k100": "",
            "error": "",
        }

        try:
            structure = CrystalStructure.from_file(res_path, fmt="res")
            normal_expanded = structure.expand_to_explicit_unit_cell()
            csorm_expanded = structure.expand_to_csorm_explicit_unit_cell()
            emd_value = pdd_distance(normal_expanded, csorm_expanded, k=100)

            row["status"] = "ok"
            row["normal_atom_count"] = len(normal_expanded.atoms)
            row["csorm_atom_count"] = len(csorm_expanded.atoms)
            row["emd_k100"] = f"{emd_value:.10f}"
        except Exception as exc:
            row["status"] = "error"
            row["error"] = str(exc)

        rows.append(row)

    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    ok_rows = [row for row in rows if row["status"] == "ok"]
    error_rows = [row for row in rows if row["status"] == "error"]

    print(f"systems={len(rows)}")
    print(f"ok={len(ok_rows)}")
    print(f"errors={len(error_rows)}")
    print(f"output_csv={output_csv}")
    if ok_rows:
        nonzero = sum(float(row["emd_k100"]) > 0.0 for row in ok_rows)
        print(f"nonzero_emd={nonzero}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
