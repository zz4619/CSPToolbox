"""Plot U_tot_DB14_pcm11 vs E_latt from the bookkeeping metadata."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_METADATA_PATH = Path("/Users/zianzhan/Desktop/CSP_sandbox/Calculations/01_results.csv")
DEFAULT_OUTPUT_PNG = Path("/Users/zianzhan/Desktop/CSP_sandbox/Calculations/U_tot_DB14_pcm11_vs_U_latt.png")
DEFAULT_OUTPUT_CSV = Path(
    "/Users/zianzhan/Desktop/CSP_sandbox/Calculations/U_tot_DB14_pcm11_vs_U_latt_points.csv"
)


def _load_points(metadata_path: Path) -> list[tuple[str, float, float, bool]]:
    rows = list(csv.DictReader(metadata_path.open(newline="")))
    points: list[tuple[str, float, float, bool]] = []

    for row in rows:
        identifier = row.get("identifier") or row.get("Identifier") or row.get("name") or ""
        u_latt = (row.get("E_latt") or "").strip()
        u_tot = (row.get("U_tot_DB14_pcm11") or "").strip()
        is_hydrate = (row.get("is_hydrate") or "").strip() == "TRUE"
        if not u_latt or not u_tot:
            continue
        try:
            x_value = float(u_latt)
            y_value = float(u_tot)
        except ValueError:
            continue
        points.append((identifier, x_value, y_value, is_hydrate))

    return points


def _write_points_csv(points: list[tuple[str, float, float, bool]], output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["identifier", "U_latt", "U_tot_DB14_pcm11", "is_hydrate"])
        writer.writerows(points)


def _axis_limits(values: list[float], padding_fraction: float = 0.05) -> tuple[float, float]:
    lower = min(values)
    upper = max(values)
    span = upper - lower
    if span == 0:
        span = max(abs(lower), 1.0) * 0.1
    padding = span * padding_fraction
    return lower - padding, upper + padding


def plot_u_tot_db14_pcm11_vs_u_latt(
    metadata_path: Path = DEFAULT_METADATA_PATH,
    output_png: Path = DEFAULT_OUTPUT_PNG,
    output_csv: Path = DEFAULT_OUTPUT_CSV,
) -> tuple[int, int]:
    points = _load_points(metadata_path)
    _write_points_csv(points, output_csv)

    hydrate_x = [point[1] for point in points if point[3]]
    hydrate_y = [point[2] for point in points if point[3]]
    anhydrate_x = [point[1] for point in points if not point[3]]
    anhydrate_y = [point[2] for point in points if not point[3]]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "mathtext.fontset": "dejavuserif",
        }
    )

    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=220)
    ax.scatter(
        hydrate_x,
        hydrate_y,
        s=22,
        alpha=0.85,
        edgecolors="none",
        color="blue",
        label="Hydrate",
    )
    ax.scatter(
        anhydrate_x,
        anhydrate_y,
        s=22,
        alpha=0.85,
        edgecolors="none",
        color="orange",
        label="Anhydrate",
    )

    all_x = [point[1] for point in points]
    all_y = [point[2] for point in points]
    x_limits = _axis_limits(all_x)
    y_limits = _axis_limits(all_y)
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)

    ref_min = max(x_limits[0], y_limits[0])
    ref_max = min(x_limits[1], y_limits[1])
    if ref_min < ref_max:
        ax.plot([ref_min, ref_max], [ref_min, ref_max], linestyle="--", linewidth=1, color="gray", alpha=0.6)

    ax.set_xlabel(r"$U^{\mathrm{latt}}$", fontsize=16, fontweight="semibold")
    ax.set_ylabel(r"$U^{\mathrm{tot}}_{\mathrm{DB14\_pcm11}}$", fontsize=16, fontweight="semibold")
    ax.tick_params(axis="both", labelsize=11)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=12, markerscale=1.2, handletextpad=0.6)
    plt.tight_layout()

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, bbox_inches="tight")
    plt.close(fig)

    return len(points), len(hydrate_x)


def main() -> None:
    points, hydrate_count = plot_u_tot_db14_pcm11_vs_u_latt()
    print(f"points={points}")
    print(f"hydrate={hydrate_count}")
    print(f"anhydrate={points - hydrate_count}")
    print(DEFAULT_OUTPUT_PNG)
    print(DEFAULT_OUTPUT_CSV)


if __name__ == "__main__":
    main()
