#!/usr/bin/env python3
"""Generate Supplementary Table 9: Effect-size heterogeneity test (2c).

Direct effect-size heterogeneity between KoGES and EAS panel meta-analysis
for all WC-robust + EAS non-significant (classification 2c) associations.

  Z_het = (β_KoGES − β_EAS_meta) / √(SE²_KoGES + SE²_EAS_meta)
  p_het : two-sided normal

KoGES β is the raw (uncorrected) effect size (EAS_panel_koges_beta_aligned_a1).
Comparison cohort = EAS_contributing_cohorts for that variant.

Outputs:
  RESULT_manuscript_final/SuppleTable9/SuppleTable9.csv
  RESULT_manuscript_final/SuppleTable9/SuppleTable9.xlsx
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path("/home2/sjim/openKoGES")
DEFAULT_OUT_DIR = REPO_ROOT / "RESULT_manuscript_final" / "SuppleTable9"
ED7_SCRIPT = (
    REPO_ROOT
    / "CODE_manuscript_final"
    / "5_Novelty_Replication"
    / "forNOVEL_MAF"
    / "SuppleTable8_ExtendedData7.py"
)


def _load_ed7():
    spec = importlib.util.spec_from_file_location("ed7_suppletable8", ED7_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def fmt_beta(x: float) -> str:
    if pd.isna(x):
        return ""
    # avoid "-0.000" for tiny negatives that round to zero
    v = round(float(x), 3)
    if abs(v) < 5e-4:
        v = 0.0
    return f"{v:.3f}"


def fmt_p_het(p: float) -> str:
    """Scientific notation matching manuscript style, e.g. 3.93E-08."""
    if pd.isna(p):
        return ""
    return f"{p:.2E}"  # e.g. 3.93E-08


def het_table(df: pd.DataFrame) -> pd.DataFrame:
    """Build SuppleTable9 rows for 2c associations (preserve input row order)."""
    s = df[df["classification_code"] == "2c"].copy()
    b_k = pd.to_numeric(s["koges_beta_for_wc"], errors="coerce")
    se_k = pd.to_numeric(s["koges_se"], errors="coerce")
    b_e = pd.to_numeric(s["EAS_panel_beta"], errors="coerce")
    se_e = pd.to_numeric(s["EAS_panel_SE"], errors="coerce")
    z = (b_k - b_e) / np.sqrt(se_k**2 + se_e**2)
    p = 2 * stats.norm.sf(np.abs(z))

    out = pd.DataFrame({
        "Phenotype": s["koges_phenotype"].astype(str).str.strip().to_numpy(),
        "Lead variant": s["koges_lead_variant"].astype(str).str.strip().to_numpy(),
        "KoGES β": [fmt_beta(x) for x in b_k],
        "KoGES SE": [fmt_beta(x) for x in se_k],
        "Comparison cohort": s["EAS_contributing_cohorts"].astype(str).str.strip().to_numpy(),
        "Cohort β": [fmt_beta(x) for x in b_e],
        "Cohort SE": [fmt_beta(x) for x in se_e],
        "p_het": [fmt_p_het(x) for x in p],
        "_koges_beta": np.asarray(b_k, dtype=float),
        "_koges_se": np.asarray(se_k, dtype=float),
        "_cohort_beta": np.asarray(b_e, dtype=float),
        "_cohort_se": np.asarray(se_e, dtype=float),
        "_z_het": np.asarray(z, dtype=float),
        "_p_het": np.asarray(p, dtype=float),
    })
    return out.reset_index(drop=True)


def write_xlsx(path: Path, display: pd.DataFrame, numeric: pd.DataFrame) -> None:
    import openpyxl
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "SuppleTable9"

    headers = [
        "Phenotype",
        "Lead variant",
        "KoGES β",
        "KoGES SE",
        "Comparison cohort",
        "Cohort β",
        "Cohort SE",
        "p_het",
    ]
    thin = Side(style="thin", color="D9D9D9")
    bdr = Border(left=thin, right=thin, top=thin, bottom=thin)
    hdr_fill = PatternFill("solid", fgColor="1F4E79")
    hdr_font = Font(name="Calibri", size=10, bold=True, color="FFFFFF")
    body_font = Font(name="Calibri", size=9)

    for ci, h in enumerate(headers, 1):
        cell = ws.cell(row=1, column=ci, value=h)
        cell.fill = hdr_fill
        cell.font = hdr_font
        cell.border = bdr
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    ws.row_dimensions[1].height = 28

    for ri, (_, r) in enumerate(display.iterrows(), 2):
        for ci, h in enumerate(headers, 1):
            cell = ws.cell(row=ri, column=ci, value=r[h])
            cell.font = body_font
            cell.border = bdr
            cell.alignment = Alignment(
                horizontal="left" if h in ("Phenotype", "Comparison cohort") else "center",
                vertical="center",
            )

    widths = {
        "A": 32, "B": 14, "C": 10, "D": 10,
        "E": 18, "F": 10, "G": 10, "H": 12,
    }
    for col, w in widths.items():
        ws.column_dimensions[col].width = w
    ws.freeze_panes = "A2"

    # numeric sheet for reproducibility
    ws2 = wb.create_sheet("numeric")
    num_cols = [
        "Phenotype", "Lead variant", "Comparison cohort",
        "_koges_beta", "_koges_se", "_cohort_beta", "_cohort_se", "_z_het", "_p_het",
    ]
    rename = {
        "_koges_beta": "KoGES_beta",
        "_koges_se": "KoGES_SE",
        "_cohort_beta": "Cohort_beta",
        "_cohort_se": "Cohort_SE",
        "_z_het": "Z_het",
        "_p_het": "p_het",
    }
    num = numeric[num_cols].rename(columns=rename)
    for ci, h in enumerate(num.columns, 1):
        cell = ws2.cell(row=1, column=ci, value=h)
        cell.font = hdr_font
        cell.fill = hdr_fill
        cell.border = bdr
    for ri, (_, r) in enumerate(num.iterrows(), 2):
        for ci, h in enumerate(num.columns, 1):
            cell = ws2.cell(row=ri, column=ci, value=r[h])
            cell.font = body_font
            cell.border = bdr
            if h in ("KoGES_beta", "KoGES_SE", "Cohort_beta", "Cohort_SE", "Z_het"):
                if isinstance(r[h], float) and not np.isnan(r[h]):
                    cell.number_format = "0.000000"
            if h == "p_het" and isinstance(r[h], float) and not np.isnan(r[h]):
                cell.number_format = "0.00E+00"

    # footnote
    ws3 = wb.create_sheet("footnote")
    footnote = (
        "Direct effect-size heterogeneity test between KoGES and EAS panel "
        "meta-analysis for all WC-robust EAS non-significant (2c) associations. "
        "Z_het = (β_KoGES − β_EAS_meta) / √(SE²_KoGES + SE²_EAS_meta); "
        "p_het: two-sided. Comparison cohort column lists cohorts contributing "
        "to the EAS panel meta-analysis for each variant. KoGES β is the raw "
        "(uncorrected) effect size; all variants are WC-robust (β_WC retains "
        "genome-wide significance after CLA correction). For binary Hypertension, "
        "KoGES N refers to Neff; TPMI N refers to total N (see eas_external_phenotype_N_cate.tsv)."
    )
    ws3["A1"] = footnote
    ws3["A1"].alignment = Alignment(wrap_text=True, vertical="top")
    ws3.column_dimensions["A"].width = 100
    ws3.row_dimensions[1].height = 90

    path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ed-xlsx", type=Path, default=None)
    parser.add_argument("--tsv-dir", type=Path, default=None)
    parser.add_argument("--neff-cate", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    ed7 = _load_ed7()
    ed_xlsx = args.ed_xlsx or ed7.DEFAULT_ED_XLSX
    tsv_dir = args.tsv_dir or ed7.DEFAULT_TSV_DIR
    neff = args.neff_cate or ed7.DEFAULT_NEFF_CATE

    print("Building continuous (for 2c)...")
    cont, *_ = ed7.build(ed_xlsx, tsv_dir, "continuous")
    print("Building binary (for 2c)...")
    cate, *_ = ed7.build(ed_xlsx, tsv_dir, "binary", neff_cate_path=neff)

    tab_cont = het_table(cont)
    tab_cate = het_table(cate)
    full = pd.concat([tab_cont, tab_cate], ignore_index=True)
    print(f"2c rows: continuous={len(tab_cont)}, binary={len(tab_cate)}, total={len(full)}")

    display_cols = [
        "Phenotype", "Lead variant", "KoGES β", "KoGES SE",
        "Comparison cohort", "Cohort β", "Cohort SE", "p_het",
    ]
    display = full[display_cols].copy()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "SuppleTable9.csv"
    xlsx_path = args.out_dir / "SuppleTable9.xlsx"
    display.to_csv(csv_path, index=False)
    write_xlsx(xlsx_path, display, full)
    print(f"Saved: {csv_path}")
    print(f"Saved: {xlsx_path}")
    print(display.to_string(index=False))


if __name__ == "__main__":
    main()
