"""PNG rendering of the molecules in a unit cell."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.data import atomic_numbers
from ase.data.colors import jmol_colors

from .connectivity import build_connectivity, connected_components, unwrap_component
from .records import DEFAULT_COVALENT_SCALE

if TYPE_CHECKING:
    from .structure import CrystalStructure


def write_unit_cell_molecule_image(
    structure: CrystalStructure,
    destination: str | Path,
    covalent_scale: float = DEFAULT_COVALENT_SCALE,
    *,
    title: str | None = None,
    draw_box: bool = False,
) -> None:
    """Write one PNG image containing all detected molecules in the unit cell."""

    destination_path = Path(destination)
    destination_path.parent.mkdir(parents=True, exist_ok=True)

    adjacency = build_connectivity(structure, covalent_scale)
    components = connected_components(adjacency)
    if not components:
        raise ValueError("No molecular components detected.")

    component_positions: dict[int, np.ndarray] = {}
    for component in components:
        unwrapped = unwrap_component(structure, component, adjacency)
        wrapped_coords = np.array([structure.atoms[index].coordinates for index in component], dtype=float)
        wrapped_centroid = wrapped_coords.mean(axis=0)
        ordered_indices = sorted(component)
        unwrapped_coords = np.vstack([unwrapped[index] for index in ordered_indices])
        translated_coords = unwrapped_coords - unwrapped_coords.mean(axis=0) + wrapped_centroid
        for index, coord in zip(ordered_indices, translated_coords):
            component_positions[index] = coord

    ordered_indices = list(range(len(structure.atoms)))
    coords = np.vstack([component_positions[index] for index in ordered_indices])
    xy, projection_center, projection_basis = _project_for_plot(coords)

    # The Figure API renders off-screen without pyplot, so no global backend is set.
    from matplotlib.figure import Figure

    fig = Figure(figsize=(6, 6), dpi=200)
    ax = fig.subplots()

    if draw_box:
        box_corners = _unit_cell_corners(structure)
        box_xy = _apply_projection(box_corners, projection_center, projection_basis)
        for left, right in _unit_cell_edges():
            ax.plot(
                [box_xy[left, 0], box_xy[right, 0]],
                [box_xy[left, 1], box_xy[right, 1]],
                color="dimgray",
                linewidth=0.8,
                linestyle="--",
                zorder=0,
            )

    for atom_i, neighbors in adjacency.items():
        for edge in neighbors:
            atom_j = edge.neighbor
            if atom_i < atom_j:
                ax.plot(
                    [xy[atom_i, 0], xy[atom_j, 0]],
                    [xy[atom_i, 1], xy[atom_j, 1]],
                    color="black",
                    linewidth=1.2,
                    zorder=1,
                )

    numbers = [int(atomic_numbers[atom.element]) for atom in structure.atoms]
    facecolors = [jmol_colors[number] for number in numbers]
    ax.scatter(
        xy[:, 0],
        xy[:, 1],
        s=220,
        facecolors=facecolors,
        edgecolors="black",
        linewidths=1.0,
        zorder=2,
    )

    span_x = float(np.ptp(xy[:, 0]))
    span_y = float(np.ptp(xy[:, 1]))
    span = max(span_x, span_y, 1.0)
    label_dx = 0.02 * span
    label_dy = 0.02 * span
    font_size = 6 if len(structure.atoms) > 20 else 7

    for atom_index, atom in enumerate(structure.atoms):
        ax.text(
            xy[atom_index, 0] + label_dx,
            xy[atom_index, 1] + label_dy,
            atom.label,
            fontsize=font_size,
            ha="left",
            va="bottom",
            color="black",
            bbox={
                "boxstyle": "round,pad=0.12",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.85,
            },
            zorder=3,
        )

    margin = 0.18 * span
    ax.set_xlim(float(np.min(xy[:, 0]) - margin), float(np.max(xy[:, 0]) + margin))
    ax.set_ylim(float(np.min(xy[:, 1]) - margin), float(np.max(xy[:, 1]) + margin))
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    ax.set_title(title or f"{structure.name} unit-cell molecules", fontsize=9, pad=6)
    fig.tight_layout(pad=0.1)
    fig.savefig(destination_path, bbox_inches="tight", pad_inches=0.05, transparent=False)


def _project_for_plot(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if len(coords) == 0:
        return np.zeros((0, 2)), np.zeros(3), np.eye(3)
    if len(coords) <= 2:
        centered = coords - coords.mean(axis=0)
        basis = np.eye(3)
        return centered[:, :2], coords.mean(axis=0), basis

    centered = coords - coords.mean(axis=0)
    _, _, right_vectors = np.linalg.svd(centered, full_matrices=False)
    rotated = centered @ right_vectors.T
    return rotated[:, :2], coords.mean(axis=0), right_vectors


def _apply_projection(
    coords: np.ndarray,
    center: np.ndarray,
    basis: np.ndarray,
) -> np.ndarray:
    centered = coords - center
    rotated = centered @ basis.T
    return rotated[:, :2]


def _unit_cell_corners(structure: CrystalStructure) -> np.ndarray:
    a_vec, b_vec, c_vec = np.asarray(structure.cell.array)
    origin = np.zeros(3, dtype=float)
    return np.array(
        [
            origin,
            a_vec,
            b_vec,
            c_vec,
            a_vec + b_vec,
            a_vec + c_vec,
            b_vec + c_vec,
            a_vec + b_vec + c_vec,
        ],
        dtype=float,
    )


def _unit_cell_edges() -> list[tuple[int, int]]:
    return [
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 4),
        (1, 5),
        (2, 4),
        (2, 6),
        (3, 5),
        (3, 6),
        (4, 7),
        (5, 7),
        (6, 7),
    ]
