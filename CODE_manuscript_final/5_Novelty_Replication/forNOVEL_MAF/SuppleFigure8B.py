#!/usr/bin/env python3
"""Cross-population MAF figure for binary ancestry-panel–novel associations.

Same display as Figure4B.py: one point per phenotype–rsID association, colored
by replication class and outlined when the MAF differs by at least two-fold
and |ΔMAF| ≥ 0.01. No Higher-main-candidate or Lower effect-retention class,
and no count line under the panel titles.

Inclusion matches the binary association-level figure: KoGES MAF ≥ 0.01 and a
non-missing allele-matched gnomAD MAF. Manuscript denominators are 55 EAS and
17 EUR associations.

Input (default):
  RESULT_manuscript_final/ExtendedData/ExtendedData2_categorical.xlsx
    sheet Total_summary_Novelty_Repl
  KoGES_gnomAD_EAS_EUR_allele_matching_QC_cate.tsv

Output (default, PNG only):
  RESULT_manuscript_final/SuppleFigure8/SuppleFigure8B.png

Example
-------
python3 SuppleFigure8B.py
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
DEFAULT_XLSX = REPO_ROOT / "RESULT_manuscript_final" / "ExtendedData" / "ExtendedData2_categorical.xlsx"
DEFAULT_CACHE = Path(
    "/data5/joh/openkoges/MAFplot_gnomAD_EAS_EUR/KoGES_gnomAD_EAS_EUR_allele_matching_QC_cate.tsv"
)
DEFAULT_PNG = REPO_ROOT / "RESULT_manuscript_final" / "SuppleFigure8" / "SuppleFigure8B.png"
SHEET = "Total_summary_Novelty_Repl"
FREQ_FLOOR = 1e-3
ABS_DIFF = 0.01
FOLD = 2.0

PANELS = {
    "EAS": {"reference": "gnomAD EAS", "maf_col": "MAF_EAS", "expected_n": 55},
    "EUR": {"reference": "gnomAD EUR", "maf_col": "MAF_EUR", "expected_n": 17},
}


def collapse_replication(nominal: pd.Series, bonferroni: pd.Series) -> pd.Series:
    """Do not separate nominal-only replication from non-replication."""
    nom = nominal.astype(str).str.strip()
    bon = bonferroni.astype(str).str.strip()
    out = pd.Series("Not Bonferroni replicated", index=nom.index, dtype="object")
    out[bon.eq("Replicated")] = "Bonferroni replicated"
    out[nom.eq("Not analyzed") | bon.eq("Not analyzed")] = "Not analyzed"
    out[bon.eq("Replicated")] = "Bonferroni replicated"
    return out


def two_fold_candidate(koges: np.ndarray, ref: np.ndarray) -> np.ndarray:
    fold = koges / np.maximum(ref, FREQ_FLOOR)
    outlier = (fold >= FOLD) | (fold <= 1.0 / FOLD)
    return outlier & (np.abs(koges - ref) >= ABS_DIFF) & np.isfinite(koges) & np.isfinite(ref)


def prep_panel(panel: str, merged: pd.DataFrame) -> pd.DataFrame:
    spec = PANELS[panel]
    nov = merged[merged[f"{panel} locus novelty"].eq("Novel")].copy()
    maf_ok = nov["MAF_KoGES"].ge(0.01)
    ref_ok = nov[spec["maf_col"]].notna()
    keep = nov[maf_ok & ref_ok].copy()
    koges = keep["MAF_KoGES"].to_numpy(dtype=float)
    reference = keep[spec["maf_col"]].to_numpy(dtype=float)
    keep["2-fold candidate"] = two_fold_candidate(koges, reference)
    keep["replication_class_plot"] = collapse_replication(
        keep[f"{panel} meta nominal replication"],
        keep[f"{panel} meta Bonferroni replication"],
    )
    keep["KoGES MAF for plotting"] = np.maximum(koges, FREQ_FLOOR)
    keep["Reference MAF for plotting"] = np.maximum(reference, FREQ_FLOOR)
    return keep


def plot_panel(ax, d: pd.DataFrame, panel: str):
    colors = {
        "Bonferroni replicated": "#E67E22",
        "Not Bonferroni replicated": "#C9C9C9",
    }
    for cls, z in [("Not Bonferroni replicated", 2), ("Bonferroni replicated", 4)]:
        sub = d[d["replication_class_plot"].eq(cls)]
        ax.scatter(
            sub["KoGES MAF for plotting"], sub["Reference MAF for plotting"],
            s=14, c=colors[cls], alpha=0.70, edgecolors="none", zorder=z, clip_on=False,
        )
    na = d[d["replication_class_plot"].eq("Not analyzed")]
    ax.scatter(
        na["KoGES MAF for plotting"], na["Reference MAF for plotting"],
        s=15, facecolors="white", edgecolors="#9A9A9A", linewidths=0.55, alpha=0.95, zorder=3, clip_on=False,
    )
    cand = d[d["2-fold candidate"]]
    ax.scatter(
        cand["KoGES MAF for plotting"], cand["Reference MAF for plotting"],
        s=14, facecolors="none", edgecolors="#6CC5DF", linewidths=0.70, zorder=5, clip_on=False,
    )

    axis_hi = 0.6
    axis_lo = 5e-4
    z = np.logspace(np.log10(axis_lo), np.log10(axis_hi), 300)
    ax.plot(z, z, color="#4D4D4D", linewidth=0.8, zorder=1)
    ax.plot(z, 2 * z, color="#6CC5DF", linewidth=0.85, linestyle=(0, (4, 3)), zorder=1)
    ax.plot(z, z / 2, color="#6CC5DF", linewidth=0.85, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(axis_hi, axis_lo)
    ax.set_ylim(axis_hi, axis_lo)
    ticks = [1e-3, 1e-2, 1e-1, 0.5]
    labels = ["0.001", "0.01", "0.1", "0.5"]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter(labels))
    ax.yaxis.set_major_formatter(FixedFormatter(labels))
    ax.minorticks_off()
    ax.tick_params(labelsize=7, length=2.2, width=0.6)
    ax.grid(False)
    ax.set_aspect("equal", adjustable="box")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.75)
    ax.spines["left"].set_linewidth(0.75)
    ax.set_title(f"Binary {panel}-panel novel associations", fontsize=7.0, fontweight="bold", pad=6)
    ax.text(0.0065, 0.00070, "Higher MAF\nin KoGES", fontsize=4.8, color="#467C8C", ha="center", va="center", clip_on=False, bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=0.8))
    ax.text(0.00105, 0.00105, "Comparable\nMAF", fontsize=4.8, color="#606060", ha="center", va="center", clip_on=False, bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=0.8))
    ax.text(0.00072, 0.0038, "Lower MAF\nin KoGES", fontsize=4.8, color="#467C8C", ha="center", va="center", clip_on=False, bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=0.8))


def make_figure(frames: Dict[str, pd.DataFrame], png_path: Path) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.65), dpi=600)
    for ax, panel, letter in zip(axes, ["EAS", "EUR"], ["A", "B"]):
        plot_panel(ax, frames[panel], panel)
        ax.text(-0.15, 1.08, letter, transform=ax.transAxes, fontsize=8.0, fontweight="bold", va="top")
        ax.set_xlabel("KoGES MAF", fontsize=6.3, labelpad=2)
        ax.set_ylabel(f"{PANELS[panel]['reference']} MAF", fontsize=6.3, labelpad=2)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#E67E22", markeredgecolor="none", markersize=5, label="Bonferroni replicated"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#C9C9C9", markeredgecolor="none", markersize=5, label="Not Bonferroni replicated"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#9A9A9A", markeredgewidth=0.7, markersize=5, label="Not analyzed"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#6CC5DF", markeredgewidth=0.8, markersize=5, label="≥2-fold and |ΔMAF|≥0.01"),
        Line2D([0], [0], color="#4D4D4D", linewidth=0.8, label="Identity"),
        Line2D([0], [0], color="#6CC5DF", linestyle=(0, (4, 3)), linewidth=0.85, label="2-fold boundary"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=5.2, columnspacing=1.05, handletextpad=0.32, bbox_to_anchor=(0.5, 0.005))
    fig.subplots_adjust(left=0.090, right=0.985, bottom=0.205, top=0.88, wspace=0.38)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=600, facecolor="white")
    plt.close(fig)


def load_frames(xlsx: Path, cache_path: Path) -> Dict[str, pd.DataFrame]:
    ann = pd.read_excel(xlsx, sheet_name=SHEET)
    cache = pd.read_csv(cache_path, sep="\t").drop_duplicates("rsID")
    merged = ann.merge(
        cache[["rsID", "MAF_KoGES", "MAF_EAS", "MAF_EUR"]],
        on="rsID",
        how="left",
    )
    frames = {}
    for panel, spec in PANELS.items():
        d = prep_panel(panel, merged)
        if len(d) != spec["expected_n"]:
            raise ValueError(
                f"{panel} association count is {len(d)}; expected {spec['expected_n']}"
            )
        frames[panel] = d
    return frames


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX, help="ExtendedData2_categorical.xlsx")
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE, help="Binary allele-matching cache TSV")
    parser.add_argument("--output", type=Path, default=DEFAULT_PNG, help="PNG path")
    args = parser.parse_args()
    frames = load_frames(args.xlsx, args.cache)
    make_figure(frames, args.output)
    print(args.output)


if __name__ == "__main__":
    main()
