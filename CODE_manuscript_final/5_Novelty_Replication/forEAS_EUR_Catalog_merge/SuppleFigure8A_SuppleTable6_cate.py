#!/usr/bin/env python3
"""Create Cate (categorical / binary SAIGE) trait-level external-novelty figure.

Input (default):
  RESULT_manuscript_final/ExtendedData/ExtendedData2_categorical.xlsx
  sheet Total_summary_Novelty_Repl

Output (default):
  RESULT_manuscript_final/SuppleFigure7/SuppleFigure7A.png
  RESULT_manuscript_final/SuppleTable6/SuppleTable6_cate.csv

Novel subcategories (including Cate-only catalog novel) are shown as
one "Novel 1 level" series, matching Figure 4A.

Example
-------
python3 SuppleFigure7A_SuppleTable6_cate.py
python3 SuppleFigure7A_SuppleTable6_cate.py --input /path/ExtendedData2_categorical.xlsx --output-dir /out

Note: requires Python 3 (``python`` on this cluster is often 2.7).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from novelty_figure_common import (
    CATALOG_ONLY_NOVEL,
    SHARED_CATEGORY_ORDER_PREFIX,
    SHARED_COLORS,
    SHARED_HATCH_EDGE_COLORS,
    SHARED_HATCHES,
    run,
)

DEFAULT_FIGURE_DIR = Path(
    "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleFigure7"
)
DEFAULT_TABLE_PATH = Path(
    "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleTable6/SuppleTable6_cate.csv"
)

CFG: dict[str, Any] = {
    "trait_group": "Cate",
    "default_xlsx": "ExtendedData2_categorical.xlsx",
    "default_output_dir": DEFAULT_FIGURE_DIR,
    "output_basename": "SuppleFigure7A",
    "table_output_path": DEFAULT_TABLE_PATH,
    # Same layout as Figure4A_SuppleTable6 / SuppleTable6
    "table_export_columns": [
        (None, "Phenotype"),
        ("EAS + EUR known", "EAS, EUR both known"),
        ("EAS known (EUR not known or not analyzed)", "EAS cohort known"),
        ("EUR known (EAS not known or not analyzed)", "EUR cohort known"),
        ("GWAS Catalog-only known", "GWAS Catalog-only known"),
        ("Novel total", "Novel"),
    ],
    "output_suffix": "_cate",
    "build_final_check_msg": "Check ExtendedData2_cate.R definitions.",
    "incomplete_n_assessed_allowed": (1, 2),
    "incomplete_n_assessed_error": (
        "Not-strict-level records should have one or two assessed domains "
        "(Catalog-only novel, or one ancestry panel + Catalog)."
    ),
    "allow_catalog_only_novel": True,
    "category_order": list(SHARED_CATEGORY_ORDER_PREFIX) + [CATALOG_ONLY_NOVEL],
    "colors": {
        **SHARED_COLORS,
        CATALOG_ONLY_NOVEL: "#AFC9E4",
    },
    "hatches": {
        **SHARED_HATCHES,
        CATALOG_ONLY_NOVEL: "...",
    },
    "hatch_edge_colors": {
        **SHARED_HATCH_EDGE_COLORS,
        CATALOG_ONLY_NOVEL: "#6F93BC",
    },
    # Match Figure 4A: collapse all novel subcategories to one legend/plot series.
    "legend_label_map": {
        "EAS + EUR + GWAS catalog novel": "Novel 1 level",
        "EAS + GWAS catalog novel (EUR not analyzed)": "Novel 1 level",
        "EUR + GWAS catalog novel (EAS not analyzed)": "Novel 1 level",
        CATALOG_ONLY_NOVEL: "Novel 1 level",
    },
    "trait_groups": [
        ("Cardiometabolic", [
            ("Coronary artery disease", "Coronary artery disease"),
            ("Cerebrovascular disease", "Cerebrovascular disease"),
            ("Hypertension", "Hypertension"),
            ("Hyperlipidemia", "Hyperlipidemia"),
            ("Diabetes mellitus", "Diabetes mellitus"),
            ("Fatty liver", "Fatty liver"),
            ("Gout", "Gout"),
            ("Chronic kidney disease (CKD)", "Chronic kidney disease (CKD)"),
        ]),
        ("Cancer", [
            ("Breast cancer", "Breast cancer"),
            ("Gastric cancer", "Gastric cancer"),
            ("Thyroid cancer", "Thyroid cancer"),
            ("Uterine cancer", "Uterine cancer"),
        ]),
        ("Respiratory & Allergy", [
            ("Allergic disease", "Allergic disease"),
            ("Chronic lung disease", "Chronic lung disease"),
        ]),
        ("Gastrointestinal", [
            ("Chronic gastritis", "Chronic gastritis"),
            ("Peptic ulcer", "Peptic ulcer"),
            ("Chronic gastritis / Peptic ulcer", "Chronic gastritis / Peptic ulcer"),
            ("Gallbladder disease", "Gallbladder disease"),
        ]),
        ("Musculoskeletal & Other", [
            ("Osteoporosis", "Osteoporosis"),
            ("Fracture", "Fracture"),
            ("Cataract", "Cataract"),
            ("Periodontal Disease", "Periodontal Disease"),
            ("Thyroid disease", "Thyroid disease"),
            ("Neurosis", "Neurosis"),
        ]),
    ],
    "category_labels": {
        "Cardiometabolic": "Cardiometabolic",
        "Cancer": "Cancer",
        "Respiratory & Allergy": "Respiratory\n& Allergy",
        "Gastrointestinal": "Gastrointestinal",
        "Musculoskeletal & Other": "Musculoskeletal\n& Other",
    },
    "layout": {
        "fig_height_base": 13.5,
        "fig_height_per_trait": 0.34,
        "fig_height_pad": 3.8,
        "figsize_width": 18.5,
        "width_ratios": [1.15, 4.4, 3.5, 3.6],
        "wspace": 0.03,
        "count_xlim_mode": "cate",
        "count_label_offset": None,  # computed from count_xlim
        "pct_label_offset": 1.8,
        "count_label_fontsize": 11,
        "pct_label_fontsize": 10,
        "count_label_clip": False,
        "pct_label_clip": False,
        "xlabel_fontsize": 12,
        "title_fontsize": 13,
        "title_pad": 10,
        "tick_labelsize": 11,
        "novel_xlim": 148,
        "legend_style": "cont",
        "subplots_adjust": dict(left=0.04, right=0.98, top=0.88, bottom=0.05),
        "suptitle_x": 0.14,
        "suptitle_y": 0.995,
        "suptitle_ha": "left",
        "suptitle_fontsize": 18,
        "subtitle": None,
        "savefig_pad_inches": 0.25,
    },
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Merged Cate novelty .xlsx "
        "(default: MERGED_NOVELTY_PATH or ExtendedData2_categorical.xlsx).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=f"Output directory (default: {DEFAULT_FIGURE_DIR})",
    )
    args = parser.parse_args()
    run(CFG, args.input, args.output_dir)


if __name__ == "__main__":
    main()
