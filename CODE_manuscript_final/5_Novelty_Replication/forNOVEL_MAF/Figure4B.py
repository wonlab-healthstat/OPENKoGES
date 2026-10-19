#!/usr/bin/env python3
"""Cross-population MAF figure for ancestry-panel–novel associations.

Each point is one phenotype–rsID association. rsIDs are not de-duplicated
across phenotypes. The manuscript denominators are 1,185 EAS and 354 EUR rows.
Points are colored by replication class and outlined when the MAF differs by
at least two-fold and |ΔMAF| ≥ 0.01. No separate Higher-main-candidate or
Lower effect-retention class is drawn.

Input (default):
  RESULT_cumulative_inci/Novel_MAF/MAF_outlier_inform_EAS.tsv
  RESULT_cumulative_inci/Novel_MAF/MAF_outlier_inform_EUR.tsv

Output (default, PNG only):
  RESULT_manuscript_final/Figure4/Figure4B.png

Example
-------
python3 Figure4B.py
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
DEFAULT_DATA_DIR = REPO_ROOT / "RESULT_cumulative_inci" / "Novel_MAF"
DEFAULT_PNG = REPO_ROOT / "RESULT_manuscript_final" / "Figure4" / "Figure4B.png"

PANELS = {
    "EAS": {
        "input_name": "MAF_outlier_inform_EAS.tsv",
        "reference": "gnomAD EAS",
        "expected_n": 1185,
    },
    "EUR": {
        "input_name": "MAF_outlier_inform_EUR.tsv",
        "reference": "gnomAD EUR",
        "expected_n": 354,
    },
}


def collapse_replication_class(x: pd.Series) -> pd.Series:
    """Do not separate nominal-only replication from non-replication."""
    s = x.astype(str).str.strip()
    out = pd.Series("Not Bonferroni replicated", index=s.index, dtype="object")
    out[s.eq("Bonferroni replicated")] = "Bonferroni replicated"
    out[s.eq("Not analyzed")] = "Not analyzed"
    return out


def bool_series(x: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(x):
        return x.fillna(False)
    return x.astype(str).str.strip().str.lower().isin(["true", "t", "1", "yes"])


def numeric(x: pd.Series) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def prep_panel(panel: str, data_dir: Path) -> pd.DataFrame:
    spec = PANELS[panel]
    d = pd.read_csv(data_dir / spec["input_name"], sep="\t", low_memory=False)
    for col in ["MAF_KoGES", "reference_MAF"]:
        if col in d:
            d[col] = numeric(d[col])
    d["2-fold candidate"] = bool_series(d["candidate_2fold_absdiff0_01"])
    d["replication_class_plot"] = collapse_replication_class(d["replication_class"])
    d["KoGES MAF for plotting"] = numeric(d["MAF_KoGES_plot"])
    d["Reference MAF for plotting"] = np.maximum(numeric(d["reference_MAF"]), 1e-3)
    return d


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

    ax.set_title(f"{panel}-panel novel associations", fontsize=7.0, fontweight="bold", pad=6)
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


def load_frames(data_dir: Path) -> Dict[str, pd.DataFrame]:
    frames = {}
    for panel, spec in PANELS.items():
        d = prep_panel(panel, data_dir)
        if len(d) != spec["expected_n"]:
            raise ValueError(
                f"{panel} association count is {len(d)}; expected {spec['expected_n']}"
            )
        frames[panel] = d
    return frames


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR, help="Directory with MAF_outlier_inform_{EAS,EUR}.tsv")
    parser.add_argument("--output", type=Path, default=DEFAULT_PNG, help="PNG path")
    args = parser.parse_args()
    frames = load_frames(args.data_dir)
    make_figure(frames, args.output)
    print(args.output)


if __name__ == "__main__":
    main()
