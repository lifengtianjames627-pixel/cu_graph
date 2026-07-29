"""Plot UV–vis spectra from so rich.csv."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "so rich.csv"
OUT = ROOT / "plots"
OUT.mkdir(exist_ok=True)

# (figure stem, title, water wl, water abs, edta wl, edta abs)
COMPARISONS = [
    (
        "chloride",
        "Chloride",
        "Chloride(water ligand):Wavelength(nm)",
        "Chloride(water ligand):Absorbance",
        "Chloride (EDTA):Wavelength(nm)",
        "Chloride (EDTA):Absorbance",
    ),
    (
        "nitrate",
        "Nitrate",
        "Nitrate (Water ligand):Wavelength(nm)",
        "Nitrate (Water ligand):Absorbance",
        "Nitrate (EDTA):Wavelength(nm)",
        "Nitrate (EDTA):Absorbance",
    ),
    (
        "hydroxide",
        "Hydroxide",
        "Hydroxide (water ligand):Wavelength(nm)",
        "Hydroxide (water ligand):Absorbance",
        "Hydroxide (EDTA ligand):Wavelength(nm)",
        "Hydroxide (EDTA ligand):Absorbance",
    ),
    (
        "thiocyanate",
        "Thiocyanate",
        "Thiocyanate (Water ligands):Wavelength(nm)",
        "Thiocyanate (Water ligands):Absorbance",
        "Thiocyanate (EDTA):Wavelength(nm)",
        "Thiocyanate (EDTA):Absorbance",
    ),
]


def style_ax(ax: plt.Axes) -> None:
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Absorbance")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.25)


def main() -> None:
    df = pd.read_csv(CSV)

    # 1) Original CuSO4·5H2O alone
    wl = df["Original CopperSulfatePentaHydrate:Wavelength(nm)"]
    ab = df["Original CopperSulfatePentaHydrate:Absorbance"]
    fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
    ax.plot(wl, ab, color="#1f4e79", lw=1.8, label=r"CuSO$_4\cdot$5H$_2$O (original)")
    ax.set_title(r"UV–vis: CuSO$_4\cdot$5H$_2$O (original)")
    style_ax(ax)
    path = OUT / "uvvis_original_cuso4.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"wrote {path}")

    # 1b) Cu(EDTA)SO4 from so rich2.csv (same single-trace format)
    df2 = pd.read_csv(ROOT / "so rich2.csv")
    wl2 = df2["Data Set 4:Wavelength(nm)"]
    ab2 = df2["Data Set 4:Absorbance"]
    fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
    ax.plot(wl2, ab2, color="#1f4e79", lw=1.8, label=r"Cu(EDTA)SO$_4$")
    ax.set_title(r"UV–vis: Cu(EDTA)SO$_4$")
    style_ax(ax)
    path = OUT / "uvvis_cu_edta_so4.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"wrote {path}")

    # 2) Water vs EDTA for each anion
    for stem, title, w_wl, w_ab, e_wl, e_ab in COMPARISONS:
        fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
        ax.plot(
            df[w_wl],
            df[w_ab],
            color="#c45c26",
            lw=1.8,
            label=f"{title} (water ligand)",
        )
        ax.plot(
            df[e_wl],
            df[e_ab],
            color="#2a6f6f",
            lw=1.8,
            label=f"{title} (EDTA ligand)",
        )
        ax.set_title(f"UV–vis: {title} — water vs EDTA ligand")
        style_ax(ax)
        path = OUT / f"uvvis_{stem}_water_vs_edta.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
