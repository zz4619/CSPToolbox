"""Plot U_tot_BP and U_tot_DB14_pcm11 against U_inter."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_METADATA_PATH = Path("/Users/zianzhan/Desktop/CSP_sandbox/Calculations/01_results.csv")
DEFAULT_OUTPUT_PNG = Path(
    "/Users/zianzhan/Desktop/CSP_sandbox/Calculations/U_tot_BP_and_DB14_pcm11_vs_U_inter.png"
)
DEFAULT_OUTPUT_CSV = Path(
    "/Users/zianzhan/Desktop/CSP_sandbox/Calculations/U_tot_BP_and_DB14_pcm11_vs_U_inter_points.csv"
)


def _load_points(
    metadata_path: Path,
) -> tuple[list[tuple[str, float, float]], list[tuple[str, float, float, bool]]]:
    rows = list(csv.DictReader(metadata_path.open(newline="")))
    bp_points: list[tuple[str, float, float]] = []
    pcm11_points: list[tuple[str, float, float, bool]] = []

    for row in rows:
        identifier = row.get("identifier") or row.get("Identifier") or row.get("name") or ""
        u_inter = (row.get("U_inter (kJ/mol)") or "").strip()
        u_tot_bp = (row.get("U_tot_BP") or "").strip()
        u_tot_pcm11 = (row.get("U_tot_DB14_pcm11") or "").strip()
        is_hydrate = (row.get("is_hydrate") or "").strip() == "TRUE"

        if u_inter:
            try:
                x_value = float(u_inter)
            except ValueError:
                x_value = None
        else:
            x_value = None

        if x_value is None:
            continue

        if u_tot_bp:
            try:
                bp_points.append((identifier, x_value, float(u_tot_bp)))
            except ValueError:
                pass

        if u_tot_pcm11:
            try:
                pcm11_points.append((identifier, x_value, float(u_tot_pcm11), is_hydrate))
            except ValueError:
                pass

    return bp_points, pcm11_points


def _write_points_csv(
    bp_points: list[tuple[str, float, float]],
    pcm11_points: list[tuple[str, float, float, bool]],
    output_csv: Path,
) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    pcm11_lookup = {identifier: (x_value, y_value, is_hydrate) for identifier, x_value, y_value, is_hydrate in pcm11_points}

    identifiers = sorted({identifier for identifier, _, _ in bp_points} | set(pcm11_lookup))
    bp_lookup = {identifier: (x_value, y_value) for identifier, x_value, y_value in bp_points}

    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["identifier", "U_inter_kJ_per_mol", "U_tot_BP", "U_tot_DB14_pcm11", "is_hydrate"])
        for identifier in identifiers:
            bp_value = bp_lookup.get(identifier)
            pcm11_value = pcm11_lookup.get(identifier)
            x_value = ""
            bp_y = ""
            pcm11_y = ""
            is_hydrate = ""
            if bp_value is not None:
                x_value = bp_value[0]
                bp_y = bp_value[1]
            if pcm11_value is not None:
                x_value = pcm11_value[0]
                pcm11_y = pcm11_value[1]
                is_hydrate = pcm11_value[2]
            writer.writerow([identifier, x_value, bp_y, pcm11_y, is_hydrate])


def plot_u_tot_bp_and_db14_pcm11_vs_u_inter(
    metadata_path: Path = DEFAULT_METADATA_PATH,
    output_png: Path = DEFAULT_OUTPUT_PNG,
    output_csv: Path = DEFAULT_OUTPUT_CSV,
) -> tuple[int, int, int]:
    bp_points, pcm11_points = _load_points(metadata_path)
    _write_points_csv(bp_points, pcm11_points, output_csv)

    bp_x = [point[1] for point in bp_points]
    bp_y = [point[2] for point in bp_points]
    hydrate_x = [point[1] for point in pcm11_points if point[3]]
    hydrate_y = [point[2] for point in pcm11_points if point[3]]
    anhydrate_x = [point[1] for point in pcm11_points if not point[3]]
    anhydrate_y = [point[2] for point in pcm11_points if not point[3]]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "mathtext.fontset": "dejavuserif",
        }
    )

    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=220)
    ax.scatter(
        bp_x,
        bp_y,
        s=18,
        alpha=0.45,
        edgecolors="none",
        color="gray",
        label="pcm7",
        zorder=1,
    )
    ax.scatter(
        hydrate_x,
        hydrate_y,
        s=24,
        alpha=0.9,
        edgecolors="white",
        linewidths=0.25,
        color="blue",
        label="pcm11 - Hydrate",
        zorder=3,
    )
    ax.scatter(
        anhydrate_x,
        anhydrate_y,
        s=24,
        alpha=0.9,
        edgecolors="white",
        linewidths=0.25,
        color="gold",
        label="pcm11 - Anhydrate",
        zorder=3,
    )
    ax.set_xlabel(r"$U^{\mathrm{inter}}_{\mathrm{VASP}}$", fontsize=16, fontweight="semibold")
    ax.set_ylabel(r"$U^{\mathrm{tot}}$", fontsize=16, fontweight="semibold")
    ax.tick_params(axis="both", labelsize=11)
    ax.set_xlim(-350, 0)
    ax.set_ylim(-350, 0)
    ax.plot([-350, 0], [-350, 0], linestyle="--", linewidth=1, color="gray", alpha=0.6, zorder=0)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=11, markerscale=1.2, handletextpad=0.6)
    plt.tight_layout()

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, bbox_inches="tight")
    plt.close(fig)

    return len(bp_points), len(pcm11_points), len(hydrate_x)


def main() -> None:
    bp_points, pcm11_points, hydrate_count = plot_u_tot_bp_and_db14_pcm11_vs_u_inter()
    print(f"bp_points={bp_points}")
    print(f"pcm11_points={pcm11_points}")
    print(f"pcm11_hydrate={hydrate_count}")
    print(f"pcm11_anhydrate={pcm11_points - hydrate_count}")
    print(DEFAULT_OUTPUT_PNG)
    print(DEFAULT_OUTPUT_CSV)


if __name__ == "__main__":
    main()
