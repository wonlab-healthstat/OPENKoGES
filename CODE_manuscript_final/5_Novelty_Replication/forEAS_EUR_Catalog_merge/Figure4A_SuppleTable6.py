#!/usr/bin/env python3
"""Create Cont (continuous) trait-level external-novelty figure (Figure 4A).

Input (default):
  RESULT_manuscript_final/ExtendedData/ExtendedData2_continuous.xlsx
  sheet Total_summary_Novelty_Repl

Output (default):
  RESULT_manuscript_final/Figure4/Figure4A.png
  RESULT_manuscript_final/SuppleTable6/SuppleTable6_cont.csv

Example
-------
python3 Figure4A_SuppleTable6.py
python3 Figure4A_SuppleTable6.py --input /path/ExtendedData2_continuous.xlsx --output-dir /out

Note: requires Python 3 (``python`` on this cluster is often 2.7).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from novelty_figure_common import (
    SHARED_CATEGORY_ORDER_PREFIX,
    SHARED_COLORS,
    SHARED_HATCH_EDGE_COLORS,
    SHARED_HATCHES,
    run,
)

DEFAULT_FIGURE4_DIR = Path(
    "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure4"
)
DEFAULT_TABLE_PATH = Path(
    "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleTable6/SuppleTable6_cont.csv"
)

CFG: dict[str, Any] = {
    "trait_group": "Cont",
    "default_xlsx": "ExtendedData2_continuous.xlsx",
    "default_output_dir": DEFAULT_FIGURE4_DIR,
    "output_basename": "Figure4A",
    "table_output_path": DEFAULT_TABLE_PATH,
    # SuppleTable6 column order / headers (matches manuscript table layout)
    "table_export_columns": [
        (None, "Phenotype"),
        ("EAS + EUR known", "EAS, EUR both known"),
        ("EAS known (EUR not known or not analyzed)", "EAS cohort known"),
        ("EUR known (EAS not known or not analyzed)", "EUR cohort known"),
        ("GWAS Catalog-only known", "GWAS Catalog-only known"),
        ("Novel total", "Novel"),
    ],
    "output_suffix": "",
    "build_final_check_msg": "Check ExtendedData2_cont.R definitions.",
    "incomplete_n_assessed_allowed": (2,),
    "incomplete_n_assessed_error": (
        "Not-strict-level records should have exactly two assessed domains."
    ),
    "allow_catalog_only_novel": False,
    "category_order": list(SHARED_CATEGORY_ORDER_PREFIX),
    "colors": dict(SHARED_COLORS),
    "hatches": dict(SHARED_HATCHES),
    "hatch_edge_colors": dict(SHARED_HATCH_EDGE_COLORS),
    # Figure 4A legend: show every novel subcategory as one "Novel 1 level" entry.
    "legend_label_map": {
        "EAS + EUR + GWAS catalog novel": "Novel 1 level",
        "EAS + GWAS catalog novel (EUR not analyzed)": "Novel 1 level",
        "EUR + GWAS catalog novel (EAS not analyzed)": "Novel 1 level",
    },
    "trait_groups": [
        ("Hematology", [
            ("Red blood cell count (RBC)", "Red blood cell count (RBC)"),
            ("White blood cell count (WBC)", "White blood cell count (WBC)"),
            ("Mean corpuscular volume (MCV)", "Mean corpuscular volume (MCV)"),
            ("Mean corpuscular hemoglobin (MCH)", "Mean corpuscular hemoglobin (MCH)"),
            ("Mean corpuscular hemoglobin concentration (MCHC)", "Mean corpuscular hemoglobin concentration (MCHC)"),
            ("Platelet count (Plt)", "Platelet count (Plt)"),
            ("Hematocrit (HCT)", "Hematocrit (HCT)"),
            ("Hemoglobin (Hb)", "Hemoglobin (Hb)"),
        ]),
        ("Pulmonary Function", [
            ("FEV1", "FEV1"),
            ("FEV1/FVC", "FEV1/FVC"),
            ("FVC", "FVC"),
        ]),
        ("Lipid Profile & Nutrition", [
            ("HDL-cholesterol (HDL-C)", "HDL-cholesterol (HDL-C)"),
            ("LDL-cholesterol (LDL-C)", "LDL-cholesterol (LDL-C)"),
            ("Triglycerides (TG)", "Triglycerides (TG)"),
            ("Total cholesterol (TC)", "Total cholesterol (TC)"),
            ("Apolipoprotein A", "Apolipoprotein A"),
            ("Apolipoprotein B", "Apolipoprotein B"),
            ("Vit.D", "Vit.D"),
            ("Iron", "Iron"),
            ("Folate", "Folate"),
            ("Albumin (Alb)", "Albumin (Alb)"),
            ("Uric acids (UA)", "Uric acids (UA)"),
            ("Homocysteine", "Homocysteine"),
        ]),
        ("Anthropometrics", [
            ("Hip", "Hip"),
            ("Body mass index (Body mass index (BMI))", "Body mass index (BMI)"),
            ("Waist-to-hip ratio (WHR)", "Waist-to-hip ratio (WHR)"),
            ("Waist", "Waist"),
            ("Body weight (BW)", "Body weight (BW)"),
            ("Body height", "Body height"),
        ]),
        ("Liver Function", [
            ("Aspartate transaminase (AST)", "Aspartate transaminase (AST)"),
            ("Alanine aminotransferase (ALT)", "Alanine aminotransferase (ALT)"),
            ("Alkaline phosphatase (ALP)", "Alkaline phosphatase (ALP)"),
            ("r_Gtp (Gamma glutamyl transferase)", "r_Gtp (Gamma glutamyl transferase)"),
            ("Total Bilirubin (TBil)", "Total Bilirubin (TBil)"),
        ]),
        ("Renal Function", [
            ("Blood urea nitrogen (BUN)", "Blood urea nitrogen (BUN)"),
            ("eGFR (CKD-EPI)", "eGFR (CKD-EPI)"),
            ("Cystatin C", "Cystatin C"),
            ("Phosphate", "Phosphate"),
            ("Serum creatinine (sCr)", "Serum creatinine (sCr)"),
        ]),
        ("Glucose Metabolism", [
            ("HbA1C", "HbA1c"),
            ("Insulin", "Insulin"),
            ("Glucose", "Glucose"),
            ("C-peptide", "C-peptide"),
        ]),
        ("Blood Pressure & Pulse", [
            ("Systolic blood pressure (SBP)", "Systolic blood pressure (SBP)"),
            ("Diastolic blood pressure (DBP)", "Diastolic blood pressure (DBP)"),
            ("Pulse", "Pulse"),
        ]),
    ],
    "category_labels": {
        "Hematology": "Hematology",
        "Pulmonary Function": "Pulmonary\nFunction",
        "Lipid Profile & Nutrition": "Lipid Profile\n& Nutrition",
        "Anthropometrics": "Anthropometrics",
        "Liver Function": "Liver\nFunction",
        "Renal Function": "Renal\nFunction",
        "Glucose Metabolism": "Glucose\nMetabolism",
        "Blood Pressure & Pulse": "Blood Pressure\n& Pulse",
    },
    "layout": {
        "fig_height_base": 12.5,
        "fig_height_per_trait": 0.31,
        "fig_height_pad": 2.4,
        "figsize_width": 18.0,
        "width_ratios": [1.2, 4.6, 3.6, 3.2],
        "wspace": 0.02,
        "count_xlim_mode": "cont",
        "count_label_offset": 12,
        "pct_label_offset": 1.2,
        "count_label_fontsize": 12,
        "pct_label_fontsize": 12,
        "count_label_clip": True,
        "pct_label_clip": True,
        "xlabel_fontsize": 13,
        "title_fontsize": 14,
        "title_pad": 16,
        "tick_labelsize": 12,
        "novel_xlim": 100,
        "legend_style": "cont",
        "subplots_adjust": dict(left=0.04, right=0.98, top=0.90, bottom=0.04),
        "suptitle_x": 0.14,
        "suptitle_y": 0.995,
        "suptitle_ha": "left",
        "suptitle_fontsize": 18,
        "subtitle": None,
        "savefig_pad_inches": None,
    },
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Merged Cont novelty .xlsx "
        "(default: MERGED_NOVELTY_PATH or ExtendedData2_continuous.xlsx).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=f"Output directory (default: {DEFAULT_FIGURE4_DIR})",
    )
    args = parser.parse_args()
    run(CFG, args.input, args.output_dir)


if __name__ == "__main__":
    main()
