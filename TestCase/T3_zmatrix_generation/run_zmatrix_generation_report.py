"""Write T3 Z-matrix generation example files and reports."""

from __future__ import annotations

from pathlib import Path

from zmatrix_generation_helpers import (
    CASE_ROOT,
    markdown_report,
    write_case_artifacts,
    write_csv_report,
)


GENERATED_DIR = CASE_ROOT / "Generated"
SUMMARY_CSV = CASE_ROOT / "t3_zmatrix_generation_report.csv"
SUMMARY_MD = CASE_ROOT / "t3_zmatrix_generation_report.md"


def main() -> int:
    reports = write_case_artifacts(GENERATED_DIR)
    write_csv_report(SUMMARY_CSV, reports)
    SUMMARY_MD.write_text(markdown_report(reports), encoding="utf-8")

    print(f"generated_dir={GENERATED_DIR}")
    print(f"summary_csv={SUMMARY_CSV}")
    print(f"summary_md={SUMMARY_MD}")

    failed = [report for report in reports if report.validation_error_count]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())

