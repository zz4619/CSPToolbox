"""Shared pytest setup for every TestCase suite.

Plotting code must never open a window during tests, and matplotlib needs a
writable cache directory on clusters with read-only home directories.
"""

from __future__ import annotations

import os
import tempfile

os.environ.setdefault("MPLBACKEND", "Agg")
if not os.access(os.environ.get("MPLCONFIGDIR", os.path.expanduser("~/.matplotlib")), os.W_OK):
    os.environ["MPLCONFIGDIR"] = tempfile.mkdtemp(prefix="csptoolbox-mpl-")
