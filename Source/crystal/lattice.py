"""Cell-parameter, fractional and Cartesian coordinate helpers."""

from __future__ import annotations

import math

import numpy as np


def lattice_matrix(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> list[tuple[float, float, float]]:
    """Cell vectors as rows: a along x, b in the xy plane (ASE convention)."""

    alpha_r = math.radians(alpha)
    beta_r = math.radians(beta)
    gamma_r = math.radians(gamma)

    lattice = [
        (a, 0.0, 0.0),
        (b * math.cos(gamma_r), b * math.sin(gamma_r), 0.0),
        (0.0, 0.0, 0.0),
    ]
    c_x = c * math.cos(beta_r)
    c_y = c * (math.cos(alpha_r) - math.cos(beta_r) * math.cos(gamma_r)) / math.sin(gamma_r)
    c_z_sq = c**2 - c_x**2 - c_y**2
    lattice[2] = (c_x, c_y, math.sqrt(max(c_z_sq, 0.0)))
    return lattice


def canonicalize_fractional(
    frac_coords: tuple[float, float, float],
    decimals: int = 8,
) -> tuple[float, float, float]:
    """Wrap fractional coordinates into [0, 1) and round to ``decimals``."""

    normalized: list[float] = []
    tolerance = 10 ** (-decimals)
    for value in frac_coords:
        wrapped = value % 1.0
        if math.isclose(wrapped, 1.0, abs_tol=tolerance):
            wrapped = 0.0
        normalized.append(round(wrapped, decimals))
    return tuple(normalized)  # type: ignore[return-value]


def fractional_minimum_image_distance(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    lattice: list[tuple[float, float, float]],
) -> float:
    """Distance in Å between two fractional positions, nearest image by rounding."""

    delta = np.array(right, dtype=float) - np.array(left, dtype=float)
    delta -= np.rint(delta)
    cartesian = np.dot(delta, np.array(lattice, dtype=float))
    return float(np.linalg.norm(cartesian))


def frac_to_cart(
    frac_coords: tuple[float, float, float],
    lattice: list[tuple[float, float, float]],
) -> tuple[float, float, float]:
    """Convert fractional coordinates to Cartesian Å with row lattice vectors."""

    fx, fy, fz = frac_coords
    a_vec, b_vec, c_vec = lattice
    return (
        fx * a_vec[0] + fy * b_vec[0] + fz * c_vec[0],
        fx * a_vec[1] + fy * b_vec[1] + fz * c_vec[1],
        fx * a_vec[2] + fy * b_vec[2] + fz * c_vec[2],
    )
