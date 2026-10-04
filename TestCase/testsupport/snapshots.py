"""Approval ("snapshot") testing for generated text files.

A snapshot is the reviewed, expected text of a generated file, stored next to
the test. When an intended change alters the output, regenerate and review:

    CSPTOOLBOX_UPDATE_SNAPSHOTS=1 python -m pytest TestCase/T11_input_snapshots
    git diff TestCase/T11_input_snapshots/expected

Commit the regenerated files together with the code change that caused them.
"""

from __future__ import annotations

import difflib
import os
from pathlib import Path

UPDATE_ENV = "CSPTOOLBOX_UPDATE_SNAPSHOTS"


def assert_text_snapshot(actual: str, expected_path: Path) -> None:
    """Fail with a unified diff unless ``actual`` equals the stored snapshot."""

    if os.environ.get(UPDATE_ENV) == "1":
        expected_path.parent.mkdir(parents=True, exist_ok=True)
        expected_path.write_text(actual, encoding="utf-8")
        return
    if not expected_path.is_file():
        raise AssertionError(
            f"Missing snapshot {expected_path}. Create it with {UPDATE_ENV}=1 and review it."
        )
    expected = expected_path.read_text(encoding="utf-8")
    if actual != expected:
        diff = "".join(
            difflib.unified_diff(
                expected.splitlines(keepends=True),
                actual.splitlines(keepends=True),
                fromfile=f"expected/{expected_path.name}",
                tofile="actual",
                n=2,
            )
        )
        raise AssertionError(
            f"Output differs from snapshot {expected_path} "
            f"(rerun with {UPDATE_ENV}=1 if the change is intended):\n{diff[:6000]}"
        )
