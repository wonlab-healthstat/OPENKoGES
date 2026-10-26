#!/usr/bin/env python3
"""Cross-population MAF figure for EAS-cohort–novel associations (BBJ / TWB / TPMI).

Same display logic as SuppleFigure8B.py: one point per phenotype–lead association,
colored by cohort-level replication class and outlined when the MAF differs by at
least two-fold and |ΔMAF| ≥ 0.01.

Differences from SuppleFigure8B:
  - Continuous and binary leads are both included (2 × 3 layout).
  - X = KoGES MAF; Y = cohort MAF (BBJ, TWB, or TPMI), not gnomAD.
  - Novelty and replication classes come from each cohort's novelty/replication TSV
    (novel_relative_to_{cohort}, replicated_primary, replication_testable).

Binary TWB is unavailable in this manuscript pipeline (no TWB categorical GWAS),
so that panel is left empty with a note.

Inputs (default under DATA_manuscript_final/Extended_Data/):
  KoGES_BOLT_leads_{BBJ,TWB,TPMI}_novelty_replication_cont.tsv
  KoGES_SAIGE_leads_{BBJ,TPMI}_novelty_replication_cate.tsv

Output (default):
  RESULT_manuscript_final/SuppleFigure10/SuppleFigure10.png

Example
-------
python3 SuppleFigure10.py
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
DEFAULT_TSV_DIR = REPO_ROOT / "DATA_manuscript_final" / "Extended_Data"
DEFAULT_PNG = REPO_ROOT / "RESULT_manuscript_final" / "SuppleFigure10" / "SuppleFigure10.png"

FREQ_FLOOR = 1e-3
ABS_DIFF = 0.01
FOLD = 2.0
MAF_MIN = 0.01

# (trait_label, cohort) → TSV file name (None = no data)
PANEL_FILES: Dict[Tuple[str, str], Optional[str]] = {
    ("Continuous", "BBJ"): "KoGES_BOLT_leads_BBJ_novelty_replication_cont.tsv",
    ("Continuous", "TWB"): "KoGES_BOLT_leads_TWB_novelty_replication_cont.tsv",
    ("Continuous", "TPMI"): "KoGES_BOLT_leads_TPMI_novelty_replication_cont.tsv",
    ("Binary", "BBJ"): "KoGES_SAIGE_leads_BBJ_novelty_replication_cate.tsv",
    ("Binary", "TWB"): None,
    ("Binary", "TPMI"): "KoGES_SAIGE_leads_TPMI_novelty_replication_cate.tsv",
}

TRAIT_ORDER = ["Continuous", "Binary"]
COHORT_ORDER = ["BBJ", "TWB", "TPMI"]
COHORT_PREFIX = {"BBJ": "bbj", "TWB": "twb", "TPMI": "tpmi"}


def to_maf(x: pd.Series) -> pd.Series:
    e = pd.to_numeric(x, errors="coerce")
    return np.minimum(e, 1.0 - e)


def as_bool(x: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(x):
        return x.fillna(False)
    s = x.astype(str).str.strip().str.lower()
    return s.isin(["true", "t", "1", "yes"])


def collapse_replication(testable: pd.Series, primary: pd.Series) -> pd.Series:
    """Map cohort replication flags to the same three plot classes as SuppleFigure8B."""
    ok = as_bool(testable)
    bonf = as_bool(primary)
    out = pd.Series("Not Bonferroni replicated", index=testable.index, dtype="object")
    out[~ok] = "Not analyzed"
    out[ok & bonf] = "Bonferroni replicated"
    return out


def two_fold_candidate(koges: np.ndarray, ref: np.ndarray) -> np.ndarray:
    fold = koges / np.maximum(ref, FREQ_FLOOR)
    outlier = (fold >= FOLD) | (fold <= 1.0 / FOLD)
    return outlier & (np.abs(koges - ref) >= ABS_DIFF) & np.isfinite(koges) & np.isfinite(ref)


def prep_cohort_panel(tsv_path: Path, cohort: str) -> pd.DataFrame:
    prefix = COHORT_PREFIX[cohort]
    d = pd.read_csv(tsv_path, sep="\t", low_memory=False)

    novel_col = f"novel_relative_to_{prefix}"
    eaf_col = f"{prefix}_exact_eaf_aligned_a1"
    for col in (novel_col, "koges_eaf", eaf_col, "replication_testable", "replicated_primary"):
        if col not in d.columns:
            raise KeyError(f"{tsv_path.name}: missing required column {col}")

    nov = d[as_bool(d[novel_col])].copy()
    nov["MAF_KoGES"] = to_maf(nov["koges_eaf"])
    nov["MAF_cohort"] = to_maf(nov[eaf_col])

    keep = nov[nov["MAF_KoGES"].ge(MAF_MIN) & nov["MAF_cohort"].notna()].copy()
    koges = keep["MAF_KoGES"].to_numpy(dtype=float)
    reference = keep["MAF_cohort"].to_numpy(dtype=float)
    keep["2-fold candidate"] = two_fold_candidate(koges, reference)
    keep["replication_class_plot"] = collapse_replication(
        keep["replication_testable"],
        keep["replicated_primary"],
    )
    keep["KoGES MAF for plotting"] = np.maximum(koges, FREQ_FLOOR)
    keep["Reference MAF for plotting"] = np.maximum(reference, FREQ_FLOOR)
    keep["cohort"] = cohort
    return keep


def plot_panel(ax, d: Optional[pd.DataFrame], trait: str, cohort: str) -> None:
    if d is None:
        # Leave the panel completely blank (no title / note).
        ax.set_axis_off()
        ax.set_facecolor("white")
        for spine in ax.spines.values():
            spine.set_visible(False)
        return

    colors = {
        "Bonferroni replicated": "#E67E22",
        "Not Bonferroni replicated": "#C9C9C9",
    }
    for cls, z in [("Not Bonferroni replicated", 2), ("Bonferroni replicated", 4)]:
        sub = d[d["replication_class_plot"].eq(cls)]
        ax.scatter(
            sub["KoGES MAF for plotting"],
            sub["Reference MAF for plotting"],
            s=14,
            c=colors[cls],
            alpha=0.70,
            edgecolors="none",
            zorder=z,
            clip_on=False,
        )
    na = d[d["replication_class_plot"].eq("Not analyzed")]
    ax.scatter(
        na["KoGES MAF for plotting"],
        na["Reference MAF for plotting"],
        s=15,
        facecolors="white",
        edgecolors="#9A9A9A",
        linewidths=0.55,
        alpha=0.95,
        zorder=3,
        clip_on=False,
    )
    cand = d[d["2-fold candidate"]]
    ax.scatter(
        cand["KoGES MAF for plotting"],
        cand["Reference MAF for plotting"],
        s=14,
        facecolors="none",
        edgecolors="#6CC5DF",
        linewidths=0.70,
        zorder=5,
        clip_on=False,
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
    ax.tick_params(labelsize=6.2, length=2.2, width=0.6)
    ax.grid(False)
    ax.set_aspect("equal", adjustable="box")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.75)
    ax.spines["left"].set_linewidth(0.75)
    ax.set_title(
        f"{trait} {cohort}-novel associations (n={len(d)})",
        fontsize=6.6,
        fontweight="bold",
        pad=5,
    )
    ax.text(
        0.0065,
        0.00070,
        "Higher MAF\nin KoGES",
        fontsize=4.4,
        color="#467C8C",
        ha="center",
        va="center",
        clip_on=False,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=0.8),
    )
    ax.text(
        0.00105,
        0.00105,
        "Comparable\nMAF",
        fontsize=4.4,
        color="#606060",
        ha="center",
        va="center",
        clip_on=False,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=0.8),
    )
    ax.text(
        0.00072,
        0.0038,
        "Lower MAF\nin KoGES",
        fontsize=4.4,
        color="#467C8C",
        ha="center",
        va="center",
        clip_on=False,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=0.8),
    )


def make_figure(frames: Dict[Tuple[str, str], Optional[pd.DataFrame]], png_path: Path) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(2, 3, figsize=(10.6, 7.0), dpi=600)
    letter_i = 0
    for r, trait in enumerate(TRAIT_ORDER):
        for c, cohort in enumerate(COHORT_ORDER):
            ax = axes[r, c]
            d = frames[(trait, cohort)]
            plot_panel(ax, d, trait, cohort)
            # Skip panel letter on the blank Binary×TWB cell; next panel is "E".
            if d is not None:
                ax.text(
                    -0.14,
                    1.08,
                    chr(ord("A") + letter_i),
                    transform=ax.transAxes,
                    fontsize=8.0,
                    fontweight="bold",
                    va="top",
                )
                letter_i += 1
                ax.set_xlabel("KoGES MAF", fontsize=6.0, labelpad=2)
                ax.set_ylabel(f"{cohort} MAF", fontsize=6.0, labelpad=2)

    handles = [
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor="#E67E22",
            markeredgecolor="none", markersize=5, label="Bonferroni replicated",
        ),
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor="#C9C9C9",
            markeredgecolor="none", markersize=5, label="Not Bonferroni replicated",
        ),
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor="white",
            markeredgecolor="#9A9A9A", markeredgewidth=0.7, markersize=5,
            label="Not analyzed",
        ),
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor="none",
            markeredgecolor="#6CC5DF", markeredgewidth=0.8, markersize=5,
            label="≥2-fold and |ΔMAF|≥0.01",
        ),
        Line2D([0], [0], color="#4D4D4D", linewidth=0.8, label="Identity"),
        Line2D(
            [0], [0], color="#6CC5DF", linestyle=(0, (4, 3)), linewidth=0.85,
            label="2-fold boundary",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=5.4,
        columnspacing=1.05,
        handletextpad=0.32,
        bbox_to_anchor=(0.5, 0.008),
    )
    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.10, top=0.94, wspace=0.32, hspace=0.32)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=600, facecolor="white")
    plt.close(fig)


def load_frames(tsv_dir: Path) -> Dict[Tuple[str, str], Optional[pd.DataFrame]]:
    frames: Dict[Tuple[str, str], Optional[pd.DataFrame]] = {}
    for trait, cohort in PANEL_FILES:
        fname = PANEL_FILES[(trait, cohort)]
        if fname is None:
            frames[(trait, cohort)] = None
            print(f"[skip] {trait} × {cohort}: no input TSV")
            continue
        path = tsv_dir / fname
        if not path.exists():
            raise FileNotFoundError(path)
        d = prep_cohort_panel(path, cohort)
        frames[(trait, cohort)] = d
        n_bonf = int((d["replication_class_plot"] == "Bonferroni replicated").sum())
        n_na = int((d["replication_class_plot"] == "Not analyzed").sum())
        n_fold = int(d["2-fold candidate"].sum())
        print(
            f"[ok] {trait} × {cohort}: n={len(d)} "
            f"(Bonferroni={n_bonf}, not analyzed={n_na}, ≥2-fold={n_fold})"
        )
    return frames


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tsv-dir",
        type=Path,
        default=DEFAULT_TSV_DIR,
        help="Directory with cohort novelty/replication TSVs",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_PNG, help="PNG path")
    args = parser.parse_args()
    frames = load_frames(args.tsv_dir)
    make_figure(frames, args.output)
    print(args.output)


if __name__ == "__main__":
    main()
