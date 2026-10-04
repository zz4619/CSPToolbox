"""Compatibility import path for the crystal-structure core.

The implementation lives in the :mod:`Source.crystal` package, split by topic
(file formats, symmetry, connectivity, Z-matrices, plotting). Existing code that
imports from ``Source.crystal_structure`` keeps working through this module;
new code may import from either place.
"""

from .crystal import *  # noqa: F401,F403
from .crystal import __all__  # noqa: F401
