"""T12: the import paths, exports and packaging that other repositories rely on.

CSPImperial imports ``Source.*``; new code may import ``csptoolbox.*``. Both paths
must expose the same objects, imports must stay free of side effects, and every
console script and packaged asset must exist.
"""

from __future__ import annotations

import importlib
import pkgutil
import subprocess
import sys

import pytest

import Source
from testsupport.paths import REPO_ROOT

SOURCE_MODULES = sorted(
    info.name
    for info in pkgutil.iter_modules(Source.__path__)
    if not info.name.startswith("_") and info.name != "CLI_scripts"
)


def _public_names(module) -> list[str]:
    if hasattr(module, "__all__"):
        return list(module.__all__)
    return [
        name
        for name, value in vars(module).items()
        if not name.startswith("_") and getattr(value, "__module__", "").startswith(module.__name__)
    ]


def test_every_lazy_export_resolves():
    missing = [name for name in Source.__all__ if not hasattr(Source, name)]
    assert not missing


@pytest.mark.parametrize("name", SOURCE_MODULES)
def test_csptoolbox_mirrors_source_module(name):
    source_module = importlib.import_module(f"Source.{name}")
    try:
        mirror = importlib.import_module(f"csptoolbox.{name}")
    except ModuleNotFoundError:
        pytest.fail(f"Add csptoolbox/{name}.py forwarding to Source.{name}")
    differing = [n for n in _public_names(source_module) if getattr(mirror, n, None) is not getattr(source_module, n)]
    assert not differing, f"csptoolbox.{name} does not re-export {differing}"


def test_crystal_structure_facade_matches_package():
    import Source.crystal as package
    import Source.crystal_structure as facade

    assert sorted(facade.__all__) == sorted(package.__all__)
    assert all(getattr(facade, name) is getattr(package, name) for name in package.__all__)


def _run_python(code: str) -> None:
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=REPO_ROOT, capture_output=True, text=True,
        env={"PATH": "", "MPLBACKEND": "svg", "PYTHONPATH": str(REPO_ROOT)},
    )
    assert result.returncode == 0, result.stderr


def test_importing_the_package_root_is_lightweight():
    _run_python(
        "import sys, Source\n"
        "heavy = [m for m in ('numpy', 'rdkit', 'pymatgen', 'ase', 'matplotlib') if m in sys.modules]\n"
        "assert not heavy, heavy\n"
    )


def test_importing_the_core_does_not_touch_matplotlib_state():
    _run_python(
        "import sys, matplotlib\n"
        "before = matplotlib.get_backend()\n"
        "import Source.crystal_structure\n"
        "assert matplotlib.get_backend() == before, matplotlib.get_backend()\n"
        "assert 'matplotlib.pyplot' not in sys.modules\n"
    )


def _pyproject() -> dict:
    tomllib = pytest.importorskip("tomllib")
    return tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))


def test_console_scripts_point_at_callables():
    broken = []
    for command, target in _pyproject()["project"]["scripts"].items():
        module_name, _, attribute = target.partition(":")
        if not callable(getattr(importlib.import_module(module_name), attribute, None)):
            broken.append(command)
    assert not broken


def test_viewer_assets_are_packaged():
    patterns = _pyproject()["tool"]["setuptools"]["package-data"]["Source.zmatrix_viewer"]
    asset_dir = REPO_ROOT / "Source" / "zmatrix_viewer"
    packaged = {path for pattern in patterns for path in asset_dir.glob(pattern)}
    present = {path for path in (asset_dir / "assets").iterdir() if path.is_file()}
    assert present and present <= packaged, sorted(str(p.name) for p in present - packaged)
