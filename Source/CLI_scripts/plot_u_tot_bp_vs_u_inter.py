"""Plot U_tot_BP vs U_inter from the bookkeeping metadata."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_METADATA_PATH = Path("/Users/zianzhan/Desktop/CSP_sandbox/Calculations/01_results.csv")
DEFAULT_OUTPUT_PNG = Path("/Users/zianzhan/Desktop/CSP_sandbox/Calculations/U_tot_BP_vs_U_inter.png")
DEFAULT_OUTPUT_CSV = Path(
    "/Users/zianzhan/Desktop/CSP_sandbox/Calculations/U_tot_BP_vs_U_inter_points_filtered.csv"
)
OUTLIER_IDENTIFIER = "RAPJIE"
def _load_points(metadata_path: Path) -> tuple[list[tuple[str, float, float, bool]], list[tuple[str, float, float, bool]]]:
    rows = list(csv.DictReader(metadata_path.open(newline="")))
    points: list[tuple[str, float, float, bool]] = []
    removed: list[tuple[str, float, float, bool]] = []

    for row in rows:
        identifier = row.get("identifier") or row.get("Identifier") or row.get("name") or ""
        u_inter = (row.get("U_inter (kJ/mol)") or "").strip()
        u_tot_bp = (row.get("U_tot_BP") or "").strip()
        is_hydrate = (row.get("is_hydrate") or "").strip() == "TRUE"
        if not u_inter or not u_tot_bp:
            continue
        try:
            x_value = float(u_inter)
            y_value = float(u_tot_bp)
        except ValueError:
            continue

        record = (identifier, x_value, y_value, is_hydrate)
        if identifier == OUTLIER_IDENTIFIER:
            removed.append(record)
        else:
            points.append(record)

    return points, removed


def _write_points_csv(points: list[tuple[str, float, float, bool]], output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["identifier", "U_inter_kJ_per_mol", "U_tot_BP", "is_hydrate"])
        writer.writerows(points)


def plot_u_tot_bp_vs_u_inter(
    metadata_path: Path = DEFAULT_METADATA_PATH,
    output_png: Path = DEFAULT_OUTPUT_PNG,
    output_csv: Path = DEFAULT_OUTPUT_CSV,
) -> tuple[int, int, int]:
    points, removed = _load_points(metadata_path)
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
    ax.set_xlabel(r"$U^{\mathrm{inter}}_{\mathrm{VASP}}$", fontsize=16, fontweight="semibold")
    ax.set_ylabel(r"$U^{\mathrm{tot}}_{\mathrm{BP}}$", fontsize=16, fontweight="semibold")
    ax.tick_params(axis="both", labelsize=11)
    ax.set_xlim(-350, 0)
    ax.set_ylim(-350, 0)
    ax.plot([-350, 0], [-350, 0], linestyle="--", linewidth=1, color="gray", alpha=0.6)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=12, markerscale=1.2, handletextpad=0.6)
    plt.tight_layout()

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, bbox_inches="tight")
    plt.close(fig)

    return len(points), len(hydrate_x), len(removed)


def main() -> None:
    points, hydrate_count, removed_count = plot_u_tot_bp_vs_u_inter()
    print(f"points={points}")
    print(f"hydrate={hydrate_count}")
    print(f"anhydrate={points - hydrate_count}")
    print(f"removed={removed_count}")
    print(DEFAULT_OUTPUT_PNG)
    print(DEFAULT_OUTPUT_CSV)


if __name__ == "__main__":
    main()
