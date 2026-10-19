#!/usr/bin/env python3
"""Shared helpers for Cont/Cate external-novelty figures.

Used by Figure4A_SuppleTable6.py (Cont) and SuppleFigure7_SuppleTable7.py (Cate).
"""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path
from typing import Any

if sys.version_info < (3, 8):
    sys.stderr.write(
        "novelty_figure requires Python 3.8+; got %s.%s. Use: python3 ...\n"
        % (sys.version_info[0], sys.version_info[1])
    )
    sys.exit(1)

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

DEFAULT_INPUT_DIR = Path(
    os.environ.get(
        "MERGED_NOVELTY_DIR",
        "/home2/sjim/openKoGES/RESULT_manuscript_final/ExtendedData",
    )
)
DEFAULT_OUTPUT_DIR = Path(
    os.environ.get("OUTPUT_DIR", "/home2/sjim/openKoGES/RESULT_cumulative_inci")
)

SHEET_NAME = "Total_summary_Novelty_Repl"
PHENOTYPE_COL = "Phenotype"
STRICT_COL = "Strict novelty across EAS, EUR, and GWAS Catalog"
INTEGRATED_COL = "Integrated novelty status"
EAS_NOVELTY_COL = "EAS locus novelty"
EUR_NOVELTY_COL = "EUR locus novelty"
CATALOG_NOVELTY_COL = "GWAS Catalog novelty"

STRICT_NOVEL = "Strict level"
INCOMPLETE = "Not strict level"
KNOWN = "Known"
INTEGRATED_NOVEL = "Novel"

OUTPUT_BASENAME = "OpenKoGES_external_novelty_by_phenotype"

SHARED_CATEGORY_ORDER_PREFIX = [
    "EAS + EUR known",
    "EAS known (EUR not known or not analyzed)",
    "EUR known (EAS not known or not analyzed)",
    "GWAS Catalog-only known",
    "EAS + EUR + GWAS catalog novel",
    "EAS + GWAS catalog novel (EUR not analyzed)",
    "EUR + GWAS catalog novel (EAS not analyzed)",
]

SHARED_COLORS = {
    "EAS + EUR known": "#C9D1D9",
    "EAS known (EUR not known or not analyzed)": "#4A5A68",
    "EUR known (EAS not known or not analyzed)": "#9EABB8",
    "GWAS Catalog-only known": "#6D7885",
    "EAS + EUR + GWAS catalog novel": "#174A7E",
    "EAS + GWAS catalog novel (EUR not analyzed)": "#77A7D3",
    "EUR + GWAS catalog novel (EAS not analyzed)": "#77A7D3",
}
SHARED_HATCHES = {
    "EAS + EUR known": None,
    "EAS known (EUR not known or not analyzed)": "///",
    "EUR known (EAS not known or not analyzed)": None,
    "GWAS Catalog-only known": None,
    "EAS + EUR + GWAS catalog novel": None,
    "EAS + GWAS catalog novel (EUR not analyzed)": "///",
    "EUR + GWAS catalog novel (EAS not analyzed)": None,
}
SHARED_HATCH_EDGE_COLORS = {
    "EAS known (EUR not known or not analyzed)": "white",
    "EAS + GWAS catalog novel (EUR not analyzed)": "#4F7FA9",
}

CATALOG_ONLY_NOVEL = "GWAS catalog novel (EAS and EUR not analyzed)"


def default_input_path(cfg: dict[str, Any]) -> Path:
    env_path = os.environ.get("MERGED_NOVELTY_PATH")
    if env_path:
        return Path(env_path)
    return DEFAULT_INPUT_DIR / cfg["default_xlsx"]


def output_stem(cfg: dict[str, Any]) -> str:
    if cfg.get("output_basename"):
        return str(cfg["output_basename"])
    return OUTPUT_BASENAME + cfg["output_suffix"]


def _n_assessed(eas: pd.Series, eur: pd.Series, catalog: pd.Series) -> pd.Series:
    return (
        eas.ne("Not assessed").astype(int)
        + eur.ne("Not assessed").astype(int)
        + catalog.ne("Not assessed").astype(int)
    )


def classify_external_evidence(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.Series:
    """Assign mutually exclusive categories for plotting.

    For known signals, "EAS known (EUR not known or not analyzed)" /
    "EUR known (EAS not known or not analyzed)" mean the named panel is known
    while the counterpart panel is *not known* (Novel) or *not analyzed*
    (Not assessed). GWAS Catalog may still be known in these groups.
    This absorbs EAS/EUR-incomplete known records into the applicable panel
    group, while a signal known only in the GWAS Catalog remains catalog-only.

    Incomplete novel signals are split by which ancestry panel was not analyzed:
    EAS+Catalog novel (EUR not analyzed) vs EUR+Catalog novel (EAS not analyzed),
    and (Cate only) Catalog-only novel (both ancestry panels not analyzed).
    """
    category_order = cfg["category_order"]
    category = pd.Series(pd.NA, index=df.index, dtype="string")
    is_known = df[STRICT_COL].eq(KNOWN)
    eas_status = df[EAS_NOVELTY_COL]
    eur_status = df[EUR_NOVELTY_COL]
    catalog_status = df[CATALOG_NOVELTY_COL]
    n_assessed = _n_assessed(eas_status, eur_status, catalog_status)
    eas_known = eas_status.eq("Known")
    eur_known = eur_status.eq("Known")
    eas_not_known = eas_status.isin(["Novel", "Not assessed"])
    eur_not_known = eur_status.isin(["Novel", "Not assessed"])

    invalid_eas_status = eas_status.loc[~(eas_known | eas_not_known)]
    invalid_eur_status = eur_status.loc[~(eur_known | eur_not_known)]
    if not invalid_eas_status.empty or not invalid_eur_status.empty:
        invalid = pd.concat(
            [
                invalid_eas_status.rename(EAS_NOVELTY_COL),
                invalid_eur_status.rename(EUR_NOVELTY_COL),
            ]
        )
        raise ValueError(
            "Unexpected EAS/EUR locus novelty status(es):\n"
            + invalid.value_counts(dropna=False).to_string()
        )

    category.loc[is_known & eas_known & eur_known] = category_order[0]
    category.loc[is_known & eas_known & eur_not_known] = category_order[1]
    category.loc[is_known & eas_not_known & eur_known] = category_order[2]
    category.loc[is_known & eas_not_known & eur_not_known] = category_order[3]

    incomplete = df[STRICT_COL].eq(INCOMPLETE)
    category.loc[df[STRICT_COL].eq(STRICT_NOVEL)] = category_order[4]
    category.loc[
        incomplete
        & eas_status.eq("Novel")
        & eur_status.eq("Not assessed")
        & catalog_status.eq("Novel")
    ] = category_order[5]
    category.loc[
        incomplete
        & eas_status.eq("Not assessed")
        & eur_status.eq("Novel")
        & catalog_status.eq("Novel")
    ] = category_order[6]
    if cfg["allow_catalog_only_novel"]:
        category.loc[
            incomplete
            & eas_status.eq("Not assessed")
            & eur_status.eq("Not assessed")
            & catalog_status.eq("Novel")
        ] = category_order[7]

    if category.isna().any():
        unexpected = df.loc[
            category.isna(),
            [STRICT_COL, EAS_NOVELTY_COL, EUR_NOVELTY_COL, CATALOG_NOVELTY_COL],
        ]
        raise ValueError(
            "Unexpected novelty classification pattern(s):\n"
            + unexpected.value_counts(dropna=False).to_string()
        )

    catalog_only_not_known = df.loc[
        category.eq(category_order[3]) & ~catalog_status.eq("Known")
    ]
    if not catalog_only_not_known.empty:
        raise ValueError(
            "A GWAS-Catalog-only known record was not known in the trait-matched catalog."
        )

    eas_only = category.eq(category_order[1])
    if not (eas_known.loc[eas_only] & eur_not_known.loc[eas_only]).all():
        raise ValueError(
            "An EAS known (EUR not known or not analyzed) record did not meet "
            "the panel-status definition."
        )
    eur_only = category.eq(category_order[2])
    if not (eas_not_known.loc[eur_only] & eur_known.loc[eur_only]).all():
        raise ValueError(
            "A EUR known (EAS not known or not analyzed) record did not meet "
            "the panel-status definition."
        )

    incomplete_not_novel = df.loc[
        df[STRICT_COL].eq(INCOMPLETE) & ~df[INTEGRATED_COL].eq(INTEGRATED_NOVEL)
    ]
    if not incomplete_not_novel.empty:
        raise ValueError(
            "A Not-strict-level record was not Novel in Integrated novelty status. "
            + cfg["build_final_check_msg"]
        )

    incomplete_n = n_assessed.loc[df[STRICT_COL].eq(INCOMPLETE)]
    allowed = set(cfg["incomplete_n_assessed_allowed"])
    if not incomplete_n.isin(allowed).all():
        raise ValueError(cfg["incomplete_n_assessed_error"])
    if not n_assessed.loc[df[STRICT_COL].eq(STRICT_NOVEL)].eq(3).all():
        raise ValueError(
            "EAS + EUR + GWAS catalog novel records should have all three domains assessed."
        )
    strict_novel = df[STRICT_COL].eq(STRICT_NOVEL)
    if not (
        eas_status.loc[strict_novel].eq("Novel")
        & eur_status.loc[strict_novel].eq("Novel")
        & catalog_status.loc[strict_novel].eq("Novel")
    ).all():
        raise ValueError(
            "EAS + EUR + GWAS catalog novel records should be novel in all three sources."
        )

    if cfg["allow_catalog_only_novel"]:
        catalog_only_novel = category.eq(category_order[7])
        if not (
            eas_status.loc[catalog_only_novel].eq("Not assessed")
            & eur_status.loc[catalog_only_novel].eq("Not assessed")
            & catalog_status.loc[catalog_only_novel].eq("Novel")
            & n_assessed.loc[catalog_only_novel].eq(1)
        ).all():
            raise ValueError(
                "GWAS catalog novel (EAS and EUR not analyzed) records must be "
                "Not assessed in both ancestry panels and Novel in GWAS Catalog."
            )
    return category


def apply_clinical_order(summary: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """Apply the requested biomarker-category and phenotype order."""
    observed = set(summary.index)
    ordered_traits: list[str] = []
    category_by_trait: dict[str, str] = {}
    display_by_trait: dict[str, str] = {}
    for category, traits in cfg["trait_groups"]:
        for source_name, display_name in traits:
            if source_name in observed:
                ordered_traits.append(source_name)
                category_by_trait[source_name] = category
                display_by_trait[source_name] = display_name

    unassigned = observed.difference(ordered_traits)
    if unassigned:
        raise ValueError(
            "Phenotype(s) missing from TRAIT_GROUPS: " + ", ".join(sorted(unassigned))
        )
    if len(ordered_traits) != len(observed):
        raise ValueError("A phenotype was included more than once in TRAIT_GROUPS.")

    summary = summary.loc[ordered_traits].copy()
    summary["Biomarker category"] = [category_by_trait[x] for x in summary.index]
    summary["Phenotype display"] = [display_by_trait[x] for x in summary.index]
    return summary


def make_summary(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    category_order = cfg["category_order"]
    df = df.copy()
    df["External evidence category"] = classify_external_evidence(df, cfg)

    summary = (
        df.groupby([PHENOTYPE_COL, "External evidence category"], observed=False)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=category_order, fill_value=0)
    )
    summary["Known total"] = summary[category_order[:4]].sum(axis=1)
    summary["Novel total"] = summary[category_order[4:]].sum(axis=1)
    summary["Total lead associations"] = summary[category_order].sum(axis=1)

    return apply_clinical_order(summary, cfg)


def draw_figure(summary: pd.DataFrame, output_dir: Path, cfg: dict[str, Any]) -> None:
    # Left panel: stacked absolute counts (known + novel categories).
    # Right panel: novel share of all lead associations for that phenotype (%).
    category_order = cfg["category_order"]
    colors = cfg["colors"]
    hatches = cfg["hatches"]
    hatch_edge_colors = cfg["hatch_edge_colors"]
    category_labels = cfg["category_labels"]
    layout = cfg["layout"]

    n_traits = len(summary)
    fig_height = max(
        layout["fig_height_base"],
        layout["fig_height_per_trait"] * n_traits + layout["fig_height_pad"],
    )
    fig, (ax_group, ax_known, ax_phenotype, ax_novel) = plt.subplots(
        1,
        4,
        figsize=(layout["figsize_width"], fig_height),
        sharey=True,
        gridspec_kw={
            "width_ratios": layout["width_ratios"],
            "wspace": layout["wspace"],
        },
    )

    y = range(n_traits)
    phenotype_labels = summary["Phenotype display"].tolist()

    group_starts: list[tuple[str, int, int]] = []
    start = 0
    for category, rows in summary.groupby("Biomarker category", sort=False):
        end = start + len(rows) - 1
        group_starts.append((category, start, end))
        start = end + 1
    for group_index, (category, start, end) in enumerate(group_starts):
        band_color = "#F3F6F9" if group_index % 2 == 0 else "#FFFFFF"
        for ax in (ax_group, ax_known, ax_phenotype, ax_novel):
            ax.axhspan(start - 0.5, end + 0.5, color=band_color, zorder=0)
        ax_group.text(
            0.5,
            (start + end) / 2,
            category_labels[category],
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            color="#4A5968",
        )
        if start > 0:
            for ax in (ax_group, ax_known, ax_phenotype, ax_novel):
                ax.axhline(start - 0.5, color="#B9C4CF", linewidth=0.8, zorder=3)

    # Optional display remapping (e.g. Cont Figure 4A: all novel -> "Novel 1 level").
    # Categories that share a display label are stacked/legended as one series.
    legend_label_map = cfg.get("legend_label_map") or {}
    plot_groups: list[tuple[str, str, list[str]]] = []
    display_to_idx: dict[str, int] = {}
    for category in category_order:
        display = str(legend_label_map.get(category, category))
        if display in display_to_idx:
            plot_groups[display_to_idx[display]][2].append(category)
        else:
            display_to_idx[display] = len(plot_groups)
            plot_groups.append((display, category, [category]))

    stack_left = pd.Series(0, index=summary.index, dtype=float)
    for display, style_category, members in plot_groups:
        values = summary[members].sum(axis=1)
        ax_known.barh(
            y,
            values,
            left=stack_left,
            color=colors[style_category],
            edgecolor=hatch_edge_colors.get(style_category, "white"),
            linewidth=0.35,
            hatch=hatches[style_category],
            label=display,
            zorder=2,
        )
        stack_left += values

    novel_categories = category_order[4:]
    novel_plot_groups = [
        (display, style_category, members)
        for display, style_category, members in plot_groups
        if any(m in novel_categories for m in members)
    ]
    totals = summary["Total lead associations"].replace(0, pd.NA)
    novel_pct = (100.0 * summary["Novel total"] / totals).fillna(0.0)

    pct_left = pd.Series(0.0, index=summary.index, dtype=float)
    for _display, style_category, members in novel_plot_groups:
        values = (100.0 * summary[members].sum(axis=1) / totals).fillna(0.0)
        ax_novel.barh(
            y,
            values,
            left=pct_left,
            color=colors[style_category],
            edgecolor=hatch_edge_colors.get(style_category, "white"),
            linewidth=0.35,
            hatch=hatches[style_category],
            zorder=2,
        )
        pct_left += values

    max_total = int(summary["Total lead associations"].max())
    if layout["count_xlim_mode"] == "cont":
        count_xlim = max(900, int(math.ceil(max_total / 100.0) * 100) + 100)
        count_label_offset = layout["count_label_offset"]
        known_xticks = [0, count_xlim // 3, 2 * count_xlim // 3, count_xlim]
    else:
        count_xlim = int(math.ceil((max_total * 1.35) / 10.0) * 10)
        count_xlim = max(count_xlim, max_total + 12)
        count_label_offset = max(2.0, count_xlim * 0.035)
        tick_step = max(10, int(math.ceil(count_xlim / 4.0 / 10.0) * 10))
        known_xticks = list(range(0, count_xlim + 1, tick_step))
        if known_xticks[-1] != count_xlim:
            known_xticks.append(count_xlim)

    for i, (_, row) in enumerate(summary.iterrows()):
        total = int(row["Total lead associations"])
        novel_total = int(row["Novel total"])
        pct = float(novel_pct.iloc[i])
        ax_known.text(
            total + count_label_offset,
            i,
            f"{total:,}",
            ha="right",
            va="center",
            fontsize=layout["count_label_fontsize"],
            color="#3F4852",
            clip_on=layout["count_label_clip"],
            zorder=4,
        )
        if novel_total > 0:
            ax_novel.text(
                pct + layout["pct_label_offset"],
                i,
                f"{pct:.1f}% ({novel_total:,}/{total:,})",
                ha="left",
                va="center",
                fontsize=layout["pct_label_fontsize"],
                color="#174A7E",
                clip_on=layout["pct_label_clip"],
                zorder=4,
            )

    for i, phenotype_label in enumerate(phenotype_labels):
        ax_phenotype.text(
            0.5,
            i,
            phenotype_label,
            ha="center",
            va="center",
            fontsize=11,
            color="#1F2328",
            zorder=4,
        )

    ax_known.set_yticks(list(y))
    ax_known.set_yticklabels([])
    ax_known.invert_yaxis()
    ax_known.set_xlim(0, count_xlim)
    ax_known.invert_xaxis()
    ax_known.set_xticks(known_xticks)
    ax_known.set_xlabel(
        "Number of lead variants (known + novel)",
        fontsize=layout["xlabel_fontsize"],
    )
    ax_known.set_title(
        "Phenotype–lead association counts",
        fontsize=layout["title_fontsize"],
        fontweight="bold",
        pad=layout["title_pad"],
    )
    ax_known.xaxis.grid(True, color="#E6E6E6", linewidth=0.7)
    ax_known.set_axisbelow(True)
    ax_known.tick_params(axis="y", length=0, labelleft=False)
    ax_known.tick_params(axis="x", labelsize=layout["tick_labelsize"])

    ax_novel.set_xlim(0, layout["novel_xlim"])
    ax_novel.set_xticks([0, 25, 50, 75, 100])
    ax_novel.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}%"))
    novel_panel_label = "Proportion (%) (Novel/Total)"
    ax_novel.set_xlabel(novel_panel_label, fontsize=layout["xlabel_fontsize"])
    ax_novel.set_title(
        novel_panel_label,
        fontsize=layout["title_fontsize"],
        fontweight="bold",
        pad=layout["title_pad"],
        color="#174A7E",
    )
    ax_novel.xaxis.grid(True, color="#E6E6E6", linewidth=0.7)
    ax_novel.set_axisbelow(True)
    ax_novel.tick_params(axis="y", length=0, labelleft=False)
    ax_novel.tick_params(axis="x", labelsize=layout["tick_labelsize"])

    for ax in (ax_known, ax_novel):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
    ax_group.set_xlim(0, 1)
    ax_group.set_xticks([])
    ax_group.tick_params(axis="y", left=False, labelleft=False)
    for spine in ax_group.spines.values():
        spine.set_visible(False)
    ax_phenotype.set_xlim(0, 1)
    ax_phenotype.set_xticks([])
    ax_phenotype.tick_params(axis="y", left=False, labelleft=False)
    for spine in ax_phenotype.spines.values():
        spine.set_visible(False)

    handles = []
    for display, style_category, members in plot_groups:
        n = int(summary[members].sum().sum())
        handles.append(
            Patch(
                facecolor=colors[style_category],
                edgecolor=hatch_edge_colors.get(style_category, "white"),
                hatch=hatches[style_category],
                label=f"{display} (n={n:,})",
            )
        )

    if layout["legend_style"] == "cont":
        # When novel subcategories are not collapsed, keep the first novel
        # entry on its own legend row above the incomplete-novel labels.
        if len(handles) > 5:
            legend_handles = handles[:5] + [
                Patch(facecolor="none", edgecolor="none", label=" "),
                *handles[5:],
            ]
        else:
            legend_handles = handles
        fig.legend(
            handles=legend_handles,
            loc="upper center",
            ncol=2,
            frameon=False,
            fontsize=10.5,
            bbox_to_anchor=(0.60, 0.985),
        )
        fig.suptitle(
            "External novelty landscape of OpenKoGES lead associations",
            x=layout["suptitle_x"],
            y=layout["suptitle_y"],
            ha=layout["suptitle_ha"],
            fontsize=layout["suptitle_fontsize"],
            fontweight="bold",
        )
        fig.subplots_adjust(**layout["subplots_adjust"])
    else:
        # Two rows of 4: known block, then novel block (avoids title/legend collision).
        fig.subplots_adjust(**layout["subplots_adjust"])
        fig.suptitle(
            "External novelty landscape of OpenKoGES lead associations",
            x=layout["suptitle_x"],
            y=layout["suptitle_y"],
            ha=layout["suptitle_ha"],
            fontsize=layout["suptitle_fontsize"],
            fontweight="bold",
        )
        if layout["subtitle"]:
            fig.text(
                0.5,
                0.955,
                layout["subtitle"],
                ha="center",
                va="top",
                fontsize=11,
                color="#4A5968",
            )
        fig.legend(
            handles=handles,
            loc="upper center",
            ncol=4,
            frameon=False,
            fontsize=8.5,
            columnspacing=1.0,
            handletextpad=0.45,
            labelspacing=0.55,
            bbox_to_anchor=(0.5, 0.93),
        )

    stem = output_stem(cfg)
    png_path = output_dir / f"{stem}.png"
    save_kw: dict[str, Any] = dict(bbox_inches="tight", facecolor="white")
    if layout["savefig_pad_inches"] is not None:
        save_kw["pad_inches"] = layout["savefig_pad_inches"]
    fig.savefig(png_path, dpi=450, **save_kw)
    plt.close(fig)


def format_summary_table(summary: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """Build the export table (e.g. SuppleTable6) from the figure summary.

    If cfg["table_export_columns"] is set, it is a list of
    (source_column_or_None_for_index, export_header) pairs.
    Phenotype comes from the summary index when source is None / "Phenotype".
    """
    export_cols = cfg.get("table_export_columns")
    if not export_cols:
        out = summary.reset_index()
        drop = [c for c in ("Biomarker category", "Phenotype display") if c in out.columns]
        if drop:
            out = out.drop(columns=drop)
        return out

    out = pd.DataFrame(index=summary.index)
    for source, header in export_cols:
        if source is None or source == "Phenotype":
            out[header] = summary.index.astype(str)
        else:
            out[header] = summary[source].values
    return out.reset_index(drop=True)


def run(cfg: dict[str, Any], input_path: Path | None, output_dir: Path | None) -> None:
    """Load Extended Data 2 Excel, write summary table + PNG figure."""
    resolved_input = input_path if input_path is not None else default_input_path(cfg)
    if output_dir is not None:
        resolved_output = Path(output_dir)
    elif cfg.get("default_output_dir") is not None:
        resolved_output = Path(cfg["default_output_dir"])
    else:
        resolved_output = DEFAULT_OUTPUT_DIR

    if not resolved_input.is_file():
        raise FileNotFoundError(f"Input file not found: {resolved_input}")
    resolved_output.mkdir(parents=True, exist_ok=True)

    print(f"Trait group: {cfg['trait_group']}")
    print(f"Input: {resolved_input.resolve()}")
    print(f"Output directory: {resolved_output.resolve()}")

    df = pd.read_excel(resolved_input, sheet_name=SHEET_NAME, engine="openpyxl")
    required = {
        PHENOTYPE_COL,
        STRICT_COL,
        INTEGRATED_COL,
        EAS_NOVELTY_COL,
        EUR_NOVELTY_COL,
        CATALOG_NOVELTY_COL,
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError("Missing required column(s): " + ", ".join(sorted(missing)))

    summary = make_summary(df, cfg)
    if int(summary["Total lead associations"].sum()) != len(df):
        raise ValueError("Trait-level counts do not sum to the input row count.")

    stem = output_stem(cfg)
    if cfg.get("table_output_path") is not None:
        output_table = Path(cfg["table_output_path"])
    else:
        output_table = resolved_output / f"{stem}.tsv"
    output_table.parent.mkdir(parents=True, exist_ok=True)
    sep = "," if output_table.suffix.lower() == ".csv" else "\t"
    table = format_summary_table(summary, cfg)
    table.to_csv(output_table, sep=sep, index=False)
    draw_figure(summary, resolved_output, cfg)

    totals = summary[cfg["category_order"]].sum()
    print(f"Input associations: {len(df):,}")
    print(f"Phenotypes: {summary.shape[0]:,}")
    for category, count in totals.items():
        print(f"{category}: {int(count):,}")
    print(f"Wrote: {output_table}")
    print(f"Wrote: {resolved_output / f'{stem}.png'}")
