"""Summarize UV–vis λ_max from so rich.csv / so rich2.csv and render schematic tables."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "plots"
OUT.mkdir(exist_ok=True)

# (display name, wavelength col, absorbance col, csv path)
WATER_GROUP = [
    (
        r"CuSO$_4\cdot$5H$_2$O (original)",
        "Original CopperSulfatePentaHydrate:Wavelength(nm)",
        "Original CopperSulfatePentaHydrate:Absorbance",
        ROOT / "so rich.csv",
    ),
    (
        "Chloride (water ligand)",
        "Chloride(water ligand):Wavelength(nm)",
        "Chloride(water ligand):Absorbance",
        ROOT / "so rich.csv",
    ),
    (
        "Nitrate (water ligand)",
        "Nitrate (Water ligand):Wavelength(nm)",
        "Nitrate (Water ligand):Absorbance",
        ROOT / "so rich.csv",
    ),
    (
        "Hydroxide (water ligand)",
        "Hydroxide (water ligand):Wavelength(nm)",
        "Hydroxide (water ligand):Absorbance",
        ROOT / "so rich.csv",
    ),
    (
        "Thiocyanate (water ligand)",
        "Thiocyanate (Water ligands):Wavelength(nm)",
        "Thiocyanate (Water ligands):Absorbance",
        ROOT / "so rich.csv",
    ),
]

EDTA_GROUP = [
    (
        r"Cu(EDTA)SO$_4$ (baseline)",
        "Data Set 4:Wavelength(nm)",
        "Data Set 4:Absorbance",
        ROOT / "so rich2.csv",
    ),
    (
        "Chloride (EDTA ligand)",
        "Chloride (EDTA):Wavelength(nm)",
        "Chloride (EDTA):Absorbance",
        ROOT / "so rich.csv",
    ),
    (
        "Nitrate (EDTA ligand)",
        "Nitrate (EDTA):Wavelength(nm)",
        "Nitrate (EDTA):Absorbance",
        ROOT / "so rich.csv",
    ),
    (
        "Hydroxide (EDTA ligand)",
        "Hydroxide (EDTA ligand):Wavelength(nm)",
        "Hydroxide (EDTA ligand):Absorbance",
        ROOT / "so rich.csv",
    ),
    (
        "Thiocyanate (EDTA ligand)",
        "Thiocyanate (EDTA):Wavelength(nm)",
        "Thiocyanate (EDTA):Absorbance",
        ROOT / "so rich.csv",
    ),
]


def peak(wl: pd.Series, ab: pd.Series) -> tuple[float, float]:
    idx = ab.idxmax()
    return float(wl.loc[idx]), float(ab.loc[idx])


def build_table(rows: list[tuple]) -> pd.DataFrame:
    records = []
    for name, w_col, a_col, csv_path in rows:
        df = pd.read_csv(csv_path)
        lam, amp = peak(df[w_col], df[a_col])
        records.append(
            {
                "Complex": name,
                "λ_max (nm)": round(lam, 1),
                "A_max": round(amp, 4),
            }
        )
    return pd.DataFrame(records)


def render_table_image(df: pd.DataFrame, title: str, path: Path) -> None:
    # Matplotlib table needs plain strings (no mathtext in cells)
    display_df = df.copy()
    display_df["Complex"] = display_df["Complex"].str.replace(r"\$", "", regex=True)
    display_df["Complex"] = display_df["Complex"].str.replace("_", "")

    fig, ax = plt.subplots(figsize=(8, 0.45 * len(df) + 1.2))
    ax.axis("off")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=16)

    table = ax.table(
        cellText=display_df.values,
        colLabels=display_df.columns,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.15, 1.6)

    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor("#1f4e79")
            cell.set_text_props(color="white", fontweight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#eef3f8")

    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    water_df = build_table(WATER_GROUP)
    edta_df = build_table(EDTA_GROUP)

    water_csv = OUT / "uvvis_peaks_water_ligand.csv"
    edta_csv = OUT / "uvvis_peaks_edta_ligand.csv"
    water_df.to_csv(water_csv, index=False)
    edta_df.to_csv(edta_csv, index=False)

    render_table_image(
        water_df,
        "UV–vis max absorption — water-ligand Cu(II) complexes",
        OUT / "uvvis_peaks_table_water_ligand.png",
    )
    render_table_image(
        edta_df,
        "UV–vis max absorption — EDTA-ligand Cu(II) complexes",
        OUT / "uvvis_peaks_table_edta_ligand.png",
    )

    print("=== Water ligand group ===")
    print(water_df.to_string(index=False))
    print()
    print("=== EDTA ligand group ===")
    print(edta_df.to_string(index=False))
    print()
    print(f"wrote {water_csv}")
    print(f"wrote {edta_csv}")
    print(f"wrote {OUT / 'uvvis_peaks_table_water_ligand.png'}")
    print(f"wrote {OUT / 'uvvis_peaks_table_edta_ligand.png'}")


if __name__ == "__main__":
    main()
