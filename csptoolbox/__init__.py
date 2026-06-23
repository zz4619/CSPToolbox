"""Public package namespace for CSPToolbox.

The historical code lives in the ``Source`` package. This namespace gives new
code a stable import path while legacy scripts continue to work.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

from Source import __all__ as __all__


def __getattr__(name: str) -> Any:
    if name not in __all__:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    return getattr(import_module("Source"), name)
