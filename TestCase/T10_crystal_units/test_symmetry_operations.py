"""Unit tests for symmetry-operation parsing in Source.crystal.symmetry."""

from __future__ import annotations

import pytest

from Source.crystal.symmetry import apply_symmetry_operation

POINT = (0.1, 0.2, 0.3)


@pytest.mark.parametrize(
    ("operation", "expected"),
    [
        ("x,y,z", (0.1, 0.2, 0.3)),
        ("-x,-y,-z", (-0.1, -0.2, -0.3)),
        ("-x+1/2,y+1/2,-z+1/2", (-0.1 + 1 / 2, 0.2 + 1 / 2, -0.3 + 1 / 2)),
        ("1/2+X, 1/2-Y, Z", (1 / 2 + 0.1, 1 / 2 - 0.2, 0.3)),
        ("x-y,x,z+1/6", (0.1 - 0.2, 0.1, 0.3 + 1 / 6)),
        ("x+0.5,-y+.25,z", (0.1 + 0.5, -0.2 + 0.25, 0.3)),
        ("x+-1/2,y,z", (0.1 - 0.5, 0.2, 0.3)),
    ],
)
def test_operations_give_python_arithmetic_results(operation, expected):
    assert apply_symmetry_operation(operation, POINT) == pytest.approx(expected, abs=0, rel=0)


@pytest.mark.parametrize(
    "operation",
    [
        "__import__('os').system('true'),y,z",
        "().__class__,y,z",
        "2*x,y,z",
        "x+,y,z",
        "(x),y,z",
    ],
)
def test_anything_but_symmetry_arithmetic_is_rejected(operation):
    with pytest.raises(ValueError):
        apply_symmetry_operation(operation, POINT)
