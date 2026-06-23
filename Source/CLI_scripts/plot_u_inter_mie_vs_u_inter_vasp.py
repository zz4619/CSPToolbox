"""Plot U_inter_Mie_pcm11 and U_inter_DB14_pcm11 against U_inter_VASP."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_METADATA_PATH = Path("/Users/zianzhan/Desktop/CSP_sandbox/Calculations/01_results.csv")
DEFAULT_OUTPUT_PNG = Path("/Users/zianzhan/Desktop/CSP_sandbox/Calculations/U_inter_Mie_vs_U_inter_VASP.png")
DEFAULT_OUTPUT_CSV = Path(
    "/Users/zianzhan/Desktop/CSP_sandbox/Calculations/U_inter_Mie_vs_U_inter_VASP_points.csv"
)
OUTLIER_IDENTIFIER = "RAPJIE"


def _load_points(
    metadata_path: Path,
) -> tuple[list[tuple[str, float, float, bool]], list[tuple[str, float, float, bool]]]:
    rows = list(csv.DictReader(metadata_path.open(newline="")))
    db_points: list[tuple[str, float, float, bool]] = []
    mie_points: list[tuple[str, float, float, bool]] = []

    for row in rows:
        identifier = row.get("identifier") or row.get("Identifier") or row.get("name") or ""
        if identifier == OUTLIER_IDENTIFIER:
            continue
        u_inter_vasp = (row.get("U_inter_VASP") or "").strip()
        u_inter_mie = (row.get("U_inter_Mie_pcm11") or "").strip()
        u_inter_db = (row.get("U_inter_DB14_pcm11") or "").strip()
        is_hydrate = (row.get("is_hydrate") or "").strip() == "TRUE"
        if not u_inter_vasp:
            continue
        try:
            x_value = float(u_inter_vasp)
        except ValueError:
            continue

        if u_inter_db:
            try:
                db_points.append((identifier, x_value, float(u_inter_db), is_hydrate))
            except ValueError:
                pass

        if u_inter_mie:
            try:
                mie_points.append((identifier, x_value, float(u_inter_mie), is_hydrate))
            except ValueError:
                pass

    return db_points, mie_points


def _write_points_csv(
    db_points: list[tuple[str, float, float, bool]],
    mie_points: list[tuple[str, float, float, bool]],
    output_csv: Path,
) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    db_lookup = {identifier: (x_value, y_value, is_hydrate) for identifier, x_value, y_value, is_hydrate in db_points}
    mie_lookup = {identifier: (x_value, y_value, is_hydrate) for identifier, x_value, y_value, is_hydrate in mie_points}
    identifiers = sorted(set(db_lookup) | set(mie_lookup))

    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["identifier", "U_inter_VASP", "U_inter_DB14_pcm11", "U_inter_Mie_pcm11", "is_hydrate"])
        for identifier in identifiers:
            x_value = ""
            db_value = ""
            mie_value = ""
            is_hydrate = ""
            if identifier in db_lookup:
                x_value = db_lookup[identifier][0]
                db_value = db_lookup[identifier][1]
                is_hydrate = db_lookup[identifier][2]
            if identifier in mie_lookup:
                x_value = mie_lookup[identifier][0]
                mie_value = mie_lookup[identifier][1]
                is_hydrate = mie_lookup[identifier][2]
            writer.writerow([identifier, x_value, db_value, mie_value, is_hydrate])


def plot_u_inter_mie_vs_u_inter_vasp(
    metadata_path: Path = DEFAULT_METADATA_PATH,
    output_png: Path = DEFAULT_OUTPUT_PNG,
    output_csv: Path = DEFAULT_OUTPUT_CSV,
) -> tuple[int, int]:
    db_points, mie_points = _load_points(metadata_path)
    _write_points_csv(db_points, mie_points, output_csv)

    db_hydrate_x = [point[1] for point in db_points if point[3]]
    db_hydrate_y = [point[2] for point in db_points if point[3]]
    db_anhydrate_x = [point[1] for point in db_points if not point[3]]
    db_anhydrate_y = [point[2] for point in db_points if not point[3]]
    hydrate_x = [point[1] for point in mie_points if point[3]]
    hydrate_y = [point[2] for point in mie_points if point[3]]
    anhydrate_x = [point[1] for point in mie_points if not point[3]]
    anhydrate_y = [point[2] for point in mie_points if not point[3]]

    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "mathtext.fontset": "dejavuserif",
        }
    )

    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=220)
    ax.scatter(
        db_hydrate_x,
        db_hydrate_y,
        s=18,
        alpha=0.65,
        edgecolors="none",
        color="dimgray",
        label="DB-455",
        marker="o",
        zorder=1,
    )
    ax.scatter(
        db_anhydrate_x,
        db_anhydrate_y,
        s=20,
        alpha=0.65,
        edgecolors="none",
        color="dimgray",
        marker="s",
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
        label="CE-755 - Hydrate",
        marker="o",
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
        label="CE-755 - Anhydrate",
        marker="s",
        zorder=3,
    )
    ax.set_xlabel(r"$U^{\mathrm{inter}}_{\mathrm{DFT}}$", fontsize=16, fontweight="semibold")
    ax.set_ylabel(r"$U^{\mathrm{inter}}_{\mathrm{HAIEFF}}$", fontsize=16, fontweight="semibold")
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

    return len(mie_points), len(hydrate_x)


def main() -> None:
    points, hydrate_count = plot_u_inter_mie_vs_u_inter_vasp()
    print(f"points={points}")
    print(f"hydrate={hydrate_count}")
    print(f"anhydrate={points - hydrate_count}")
    print(DEFAULT_OUTPUT_PNG)
    print(DEFAULT_OUTPUT_CSV)


if __name__ == "__main__":
    main()
