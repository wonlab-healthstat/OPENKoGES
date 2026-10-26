#!/usr/bin/env python3
"""Generate Extended Data 7A/7B and Supplementary Table 8 (EAS-novel classification).

Classification scheme (revised — WC first, no power step):
  1a : WC-sensitive + EAS GWS     — WC-sensitive but EAS independently GWS
  1b : WC-sensitive + EAS nominal — WC-sensitive, direction-concordant EAS signal
  1c : WC-sensitive + no EAS     — WC-sensitive, no EAS replication
  2a : WC-robust + EAS GWS       — shared signal (EAS panel GWS)
  2b : WC-robust + EAS nominal   — direction-concordant, EAS nominally significant
  2c : WC-robust + EAS non-sig   — unexplained candidates for follow-up
  NA : Not analyzed in EAS panel

Summary table also reports, per group (only cohorts in EAS_contributing_cohorts):
  KoGES>each : n/total where KoGES n > every contributing individual cohort n
  KoGES>meta : n/total where KoGES n > sum of contributing cohort n
  (Binary: KoGES neff and BBJ/TPMI neff used)

Inputs:
  RESULT_manuscript_final/ExtendedData/ExtendedData.xlsx
  DATA_manuscript_final/Extended_Data/
    KoGES_BOLT_leads_EAS_panel_meta.tsv
    KoGES_BOLT_leads_BBJ_novelty_replication_cont.tsv
    KoGES_SAIGE_leads_EAS_panel_meta.tsv
    KoGES_SAIGE_leads_BBJ_novelty_replication_cate.tsv  (or _cate_1.tsv)
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT        = Path("/home2/sjim/openKoGES")
DEFAULT_ED_XLSX  = REPO_ROOT / "RESULT_manuscript_final" / "ExtendedData" / "ExtendedData.xlsx"
DEFAULT_TSV_DIR  = REPO_ROOT / "DATA_manuscript_final" / "Extended_Data"
DEFAULT_ED_OUT    = REPO_ROOT / "RESULT_manuscript_final" / "ExtendedData"
DEFAULT_SUPP_OUT  = REPO_ROOT / "RESULT_manuscript_final" / "SuppleTable8"
DEFAULT_NEFF_CATE = REPO_ROOT / "DATA_manuscript_final" / "Extended_Data" / \
                    "eas_external_phenotype_N_cate.tsv"

# ── Classification labels ─────────────────────────────────────────────────────
CLF = {
    'NA' : 'Not analyzed in EAS panel',
    '1a' : 'WC-sensitive + EAS GWS',
    '1b' : 'WC-sensitive + EAS nominal',
    '1c' : 'WC-sensitive + no EAS signal',
    '2a' : 'WC-robust + EAS GWS',
    '2b' : 'WC-robust + EAS nominal',
    '2c' : 'WC-robust + EAS non-sig',
}

CLF_COLORS = {
    CLF['NA'] : 'F2F2F2',
    CLF['1a'] : 'BDD7EE',
    CLF['1b'] : 'DDEEFF',
    CLF['1c'] : 'C8DCF0',
    CLF['2a'] : 'E2EFDA',
    CLF['2b'] : 'FFF2CC',
    CLF['2c'] : 'FFD966',
}

# ── helpers ──────────────────────────────────────────────────────────────────
def to_maf(s):
    e = pd.to_numeric(s, errors='coerce')
    return np.minimum(e, 1.0 - e)

def maf_cat(k, e, fold=2, thr=0.01):
    if pd.isna(k) or pd.isna(e): return 'Not available'
    if abs(k - e) < thr: return 'Comparable'
    r = k / (e + 1e-8)
    return 'Higher in KoGES' if r >= fold else 'Lower in KoGES' if r <= 1/fold else 'Comparable'

def is_true(v): return v is True or str(v).strip().lower() == 'true'

def classify_row(row):
    """WC first, then EAS result — no power step in classification."""
    if not is_true(row.get('EAS_panel_analyzed')):
        return 'NA'
    # Step 1: WC check
    bwc = row.get('beta_wc')
    se  = row.get('koges_se')
    kb  = row.get('koges_beta_for_wc')
    if pd.notna(bwc) and pd.notna(se) and se > 0 and pd.notna(kb):
        p_wc   = 2 * stats.norm.sf(abs(bwc / se))
        robust = (p_wc < 5e-8) and (np.sign(bwc) == np.sign(kb))
    else:
        robust = False
    # EAS panel result
    ep      = row.get('EAS_panel_P')
    eas_dir = is_true(row.get('EAS_panel_direction_concordant'))
    eas_nom = pd.notna(ep) and ep < 0.05  and eas_dir
    eas_gws = pd.notna(ep) and ep < 5e-8  and eas_dir
    if not robust:
        if eas_gws: return '1a'
        if eas_nom: return '1b'
        return             '1c'
    if eas_gws: return '2a'
    if eas_nom: return '2b'
    return             '2c'

def load_ext_neff_cate(path: Path) -> pd.DataFrame:
    """Load BBJ/TPMI neff table for binary traits.
    Returns df with columns: stem, bbj_neff, tpmi_neff.
    bbj_n / tpmi_n in the file are already 4*case*ctrl/(case+ctrl).
    """
    df = pd.read_csv(path, sep='\t', na_values='NA')
    df = df.rename(columns={'output_stem': 'stem',
                             'bbj_n':       'bbj_neff',
                             'tpmi_n':      'tpmi_neff'})
    df['bbj_neff']  = pd.to_numeric(df['bbj_neff'],  errors='coerce')
    df['tpmi_neff'] = pd.to_numeric(df['tpmi_neff'], errors='coerce')
    return df[['stem', 'bbj_neff', 'tpmi_neff']]

def _contributing_ext_n_cols(row, ext_n_cols):
    """Keep only N columns for cohorts listed in EAS_contributing_cohorts (e.g. BBJ+TWB)."""
    raw = row.get('EAS_contributing_cohorts')
    if pd.isna(raw) or str(raw).strip() == '':
        return []
    cohorts = {c.strip().upper() for c in str(raw).replace(';', '+').split('+') if c.strip()}
    return [c for c in ext_n_cols if c.split('_')[0].upper() in cohorts]

def _koges_gt_indiv(row, ext_n_cols, koges_n_col='koges_n_cmp'):
    kn = row.get(koges_n_col)
    if pd.isna(kn): return False
    cols = _contributing_ext_n_cols(row, ext_n_cols)
    ns = [row.get(c) for c in cols if pd.notna(row.get(c))]
    return all(n < kn for n in ns) if ns else False

def _eas_total_n(row, ext_n_cols):
    cols = _contributing_ext_n_cols(row, ext_n_cols)
    ns = [row.get(c) for c in cols if pd.notna(row.get(c))]
    return sum(ns) if ns else np.nan

# ── data loading ──────────────────────────────────────────────────────────────
def load_ed2_novel(xlsx, sheet):
    df = pd.read_excel(xlsx, sheet_name=sheet, header=1)
    df.columns = [str(c).strip() for c in df.columns]
    df['rsID'] = df['rsID'].astype(str).str.strip()
    df['Phen'] = df['Phenotype'].astype(str).str.strip()
    novel = df[df['EAS locus novelty'] == 'Novel'].copy()
    keys  = set(zip(novel['Phen'], novel['rsID']))
    info  = novel[['Phen','rsID','Nearest gene',
                   'Effect allele','Other allele','KoGES N']].copy()
    return keys, info

def load_ed6(xlsx, sheet, is_binary=False):
    df = pd.read_excel(xlsx, sheet_name=sheet, header=1)
    df.columns = [str(c).strip() for c in df.columns]
    df['rsID'] = df['rsID'].astype(str).str.strip()
    df['Phen'] = df['Phenotype'].astype(str).str.strip()
    df['beta_wc']    = pd.to_numeric(df["Winner's curse corrected beta (β_WC)"], errors='coerce')
    df['shrinkage']  = pd.to_numeric(df['Shrinkage ratio (|β_WC/β_raw|)'],       errors='coerce')
    df['koges_n_ed6']  = pd.to_numeric(df.get('KoGES N',                pd.Series(dtype=float)), errors='coerce')
    df['koges_se_ed6'] = pd.to_numeric(df.get('KoGES standard error',   pd.Series(dtype=float)), errors='coerce')
    df['koges_p_ed6']  = pd.to_numeric(df.get('KoGES P value',          pd.Series(dtype=float)), errors='coerce')
    cols = ['Phen','rsID','beta_wc','shrinkage','koges_n_ed6','koges_se_ed6','koges_p_ed6']
    if is_binary:
        df['koges_case']    = pd.to_numeric(df.get('KoGES case',    pd.Series(dtype=float)), errors='coerce')
        df['koges_control'] = pd.to_numeric(df.get('KoGES control', pd.Series(dtype=float)), errors='coerce')
        df['koges_neff']    = 4 * df['koges_case'] * df['koges_control'] / \
                              (df['koges_case'] + df['koges_control'])
        cols += ['koges_neff']
    return df[cols]

def build(ed_xlsx: Path, tsv_dir: Path, trait_type: str,
          neff_cate_path: Path | None = None):
    is_cont = (trait_type == 'continuous')
    keys, ed2_info = load_ed2_novel(ed_xlsx,
        'Extended Data2A' if is_cont else 'Extended Data2B')
    ed6 = load_ed6(ed_xlsx,
        'Extended Data6A' if is_cont else 'Extended Data6B',
        is_binary=not is_cont)

    if is_cont:
        meta_path  = tsv_dir / 'KoGES_BOLT_leads_EAS_panel_meta.tsv'
        bbj_path   = tsv_dir / 'KoGES_BOLT_leads_BBJ_novelty_replication_cont.tsv'
        ext_n_cols     = ['bbj_external_n','twb_external_n','tpmi_external_n']
        ext_n_cmp_cols = ext_n_cols   # continuous: total N = neff
        eaf_map    = {'bbj_exact_eaf_aligned_a1':'BBJ',
                      'twb_exact_eaf_aligned_a1':'TWB',
                      'tpmi_exact_eaf_aligned_a1':'TPMI'}
        bbj_want   = ['lead_row_id','koges_pheno_code',
                      'koges_eaf','koges_n','koges_se','koges_p']
    else:
        meta_path  = tsv_dir / 'KoGES_SAIGE_leads_EAS_panel_meta.tsv'
        bbj_path   = tsv_dir / 'KoGES_SAIGE_leads_BBJ_novelty_replication_cate.tsv'
        if not bbj_path.exists():
            bbj_path = tsv_dir / 'KoGES_SAIGE_leads_BBJ_novelty_replication_cate_1.tsv'
        ext_n_cols     = ['bbj_external_n','tpmi_external_n']   # display (total N)
        ext_n_cmp_cols = ['bbj_neff','tpmi_neff']              # comparison (neff)
        eaf_map    = {'bbj_exact_eaf_aligned_a1':'BBJ',
                      'tpmi_exact_eaf_aligned_a1':'TPMI'}
        bbj_want   = ['lead_row_id','koges_pheno_code',
                      'koges_eaf','koges_p','koges_beta','koges_se',
                      'koges_n']

    meta = pd.read_csv(meta_path, sep='\t', low_memory=False)
    meta['rsID_m'] = meta['koges_lead_variant'].astype(str).str.strip()
    meta['Phen_m'] = meta['koges_phenotype'].astype(str).str.strip()

    mask = meta.apply(lambda r: (r['Phen_m'], r['rsID_m']) in keys, axis=1)
    d = meta[mask].copy()

    bbj = pd.read_csv(bbj_path, sep='\t', low_memory=False,
                      usecols=lambda c: c in bbj_want)
    d = d.merge(bbj, on=['lead_row_id','koges_pheno_code'], how='left')
    d = d.merge(ed6,      left_on=['Phen_m','rsID_m'],
                right_on=['Phen','rsID'], how='left')
    d = d.merge(ed2_info, left_on=['Phen_m','rsID_m'],
                right_on=['Phen','rsID'], how='left')

    for src, dst in [('koges_n_ed6','koges_n'),
                     ('koges_se_ed6','koges_se'),
                     ('koges_p_ed6','koges_p')]:
        if dst in d.columns:
            d[dst] = pd.to_numeric(d[dst], errors='coerce').fillna(
                     pd.to_numeric(d[src], errors='coerce'))
        else:
            d[dst] = pd.to_numeric(d[src], errors='coerce')

    d['koges_n']  = pd.to_numeric(d['koges_n'],  errors='coerce')
    d['koges_se'] = pd.to_numeric(d['koges_se'], errors='coerce')
    d['koges_p']  = pd.to_numeric(d['koges_p'],  errors='coerce')
    d['EAS_panel_P'] = pd.to_numeric(d['EAS_panel_P'], errors='coerce')
    d['koges_beta_for_wc'] = pd.to_numeric(
        d['EAS_panel_koges_beta_aligned_a1'], errors='coerce')

    # Binary: merge external neff from eas_external_phenotype_N_cate.tsv
    if not is_cont:
        neff_path = neff_cate_path or DEFAULT_NEFF_CATE
        if neff_path.exists():
            neff_df = load_ext_neff_cate(neff_path)
            d['stem'] = d['koges_pheno_code'].str.split('_SAIGE').str[0]
            d = d.merge(neff_df, on='stem', how='left')
            print(f"  Binary neff joined: bbj_neff non-NA={d['bbj_neff'].notna().sum()}, "
                  f"tpmi_neff non-NA={d['tpmi_neff'].notna().sum()}")
        else:
            print(f"  [WARNING] neff file not found: {neff_path}; using total N as fallback")
            d['bbj_neff']  = d.get('bbj_external_n',  pd.Series(dtype=float))
            d['tpmi_neff'] = d.get('tpmi_external_n', pd.Series(dtype=float))

    # KoGES n for comparison (neff for binary, total n for continuous)
    if not is_cont and 'koges_neff' in d.columns:
        d['koges_n_cmp'] = d['koges_neff']
    else:
        d['koges_n_cmp'] = d['koges_n']

    # MAF columns
    if 'koges_eaf' not in d.columns:
        d['koges_eaf'] = np.nan
    d['koges_maf'] = to_maf(d['koges_eaf'])
    for col, cohort in eaf_map.items():
        if col in d.columns:
            cl = cohort.lower()
            d[f'{cl}_maf']     = to_maf(d[col])
            d[f'{cl}_maf_cat'] = d.apply(
                lambda r, c=cl: maf_cat(r['koges_maf'], r.get(f'{c}_maf')), axis=1)

    # Beta ratio
    eb = pd.to_numeric(d['EAS_panel_beta'], errors='coerce')
    d['beta_ratio'] = (
        d['koges_beta_for_wc'].abs() /
        (eb.abs().replace(0, np.nan) + 1e-10)
    ).round(2)

    # Power info columns (descriptive, not classification)
    # binary uses neff-based cols; continuous uses total N
    d['_gt_indiv'] = d.apply(
        lambda r: _koges_gt_indiv(r, ext_n_cmp_cols), axis=1)
    d['_eas_tot_n'] = d.apply(
        lambda r: _eas_total_n(r, ext_n_cmp_cols), axis=1)
    d['_gt_meta'] = d['koges_n_cmp'] > d['_eas_tot_n']

    # Classification
    d['classification_code'] = d.apply(classify_row, axis=1)
    d['classification']      = d['classification_code'].map(CLF)

    return d, ext_n_cols, eaf_map, ext_n_cmp_cols

# ── summary table ─────────────────────────────────────────────────────────────
def frac_fmt(num, denom):
    if denom == 0: return '—'
    return f"{int(num)}/{int(denom)} ({100*num/denom:.1f}%)"

def median_iqr(s, fmt='.3f'):
    s = pd.to_numeric(s, errors='coerce').dropna()
    if len(s) == 0: return '—'
    if 'x' in fmt:
        return f"{s.median():.2f}x ({s.quantile(0.25):.2f}–{s.quantile(0.75):.2f}x)"
    return f"{s.median():{fmt}} ({s.quantile(0.25):{fmt}}–{s.quantile(0.75):{fmt}})"

def make_summary(d, trait_label):
    order    = ['NA','1a','1b','1c','2a','2b','2c']
    total    = len(d)
    analyzed = (d['classification_code'] != 'NA').sum()
    rows     = []

    for code in order:
        sub = d[d['classification_code'] == code]
        n   = len(sub)
        an  = analyzed if code != 'NA' else None
        rows.append({
            'Trait type'                        : trait_label,
            'Classification'                    : CLF[code],
            'N'                                 : n,
            '% of total EAS-novel'              : f'{100*n/total:.1f}%',
            '% of analyzed'                     : f'{100*n/an:.1f}%' if an else '—',
            'Shrinkage ratio (median, IQR)'     : median_iqr(sub.get('shrinkage', pd.Series())),
            'KoGES/EAS beta ratio (median, IQR)': median_iqr(sub.get('beta_ratio', pd.Series()), 'x'),
            'KoGES n > each cohort n'           : frac_fmt(sub['_gt_indiv'].sum(), n) if n else '—',
            'KoGES n > EAS panel total n'       : frac_fmt(sub['_gt_meta'].sum(), n) if n else '—',
        })

    rows.append({
        'Trait type'                        : trait_label,
        'Classification'                    : 'Total EAS-novel',
        'N'                                 : total,
        '% of total EAS-novel'              : '100%',
        '% of analyzed'                     : f'{analyzed}/{total} analyzed',
        **{c: '—' for c in ['Shrinkage ratio (median, IQR)',
                             'KoGES/EAS beta ratio (median, IQR)',
                             'KoGES n > each cohort n',
                             'KoGES n > EAS panel total n']},
    })
    return pd.DataFrame(rows)

# ── Excel output (Extended Data 7) ────────────────────────────────────────────
COL_MAP_BASE = {
    'koges_phenotype'               : 'Phenotype',
    'koges_lead_variant'            : 'Lead variant (rsID)',
    'koges_chr'                     : 'Chr',
    'koges_bp'                      : 'Position (hg19)',
    'Effect allele'                 : 'Effect allele',
    'Other allele'                  : 'Other allele',
    'Nearest gene'                  : 'Nearest gene',
    'koges_n'                       : 'KoGES N',
    'koges_beta_for_wc'             : 'KoGES beta',
    'koges_se'                      : 'KoGES SE',
    'koges_p'                       : 'KoGES P',
    'koges_maf'                     : 'KoGES MAF',
    'beta_wc'                       : 'WC-corrected beta (β_WC)',
    'shrinkage'                     : 'Shrinkage ratio',
    'EAS_panel_analyzed'            : 'EAS panel analyzed',
    'EAS_contributing_cohorts'      : 'EAS contributing cohorts',
    'EAS_panel_beta'                : 'EAS panel beta',
    'EAS_panel_SE'                  : 'EAS panel SE',
    'EAS_panel_P'                   : 'EAS panel P',
    'EAS_panel_I2'                  : 'EAS panel I² (%)',
    'EAS_panel_Q_P'                 : 'Cochran Q P',
    'EAS_panel_direction_concordant': 'Direction concordant',
    'EAS_Q_top_cohort'              : 'Top heterogeneity source',
    'EAS_LOCO_omit_BBJ_I2'         : 'I² omitting BBJ',
    'EAS_LOCO_omit_TWB_I2'         : 'I² omitting TWB',
    'EAS_LOCO_omit_TPMI_I2'        : 'I² omitting TPMI',
    'beta_ratio'                    : '|KoGES beta / EAS panel beta|',
    'bbj_neff'                      : 'BBJ Neff',
    'tpmi_neff'                     : 'TPMI Neff',
    '_gt_indiv'                     : 'KoGES n > each cohort n',
    '_gt_meta'                      : 'KoGES n > EAS meta n',
    'classification'                : 'Analysis classification',
}

WIDTHS = {
    'Phenotype': 28, 'Lead variant (rsID)': 16, 'Chr': 5,
    'Position (hg19)': 12, 'Effect allele': 8, 'Other allele': 8,
    'Nearest gene': 14, 'KoGES N': 10, 'KoGES beta': 10,
    'KoGES SE': 10, 'KoGES P': 12, 'KoGES MAF': 10,
    'WC-corrected beta (β_WC)': 18, 'Shrinkage ratio': 13,
    'EAS panel analyzed': 14, 'EAS contributing cohorts': 20,
    'EAS panel beta': 13, 'EAS panel SE': 13, 'EAS panel P': 13,
    'EAS panel I² (%)': 13, 'Cochran Q P': 14,
    'Direction concordant': 16, 'Top heterogeneity source': 16,
    'I² omitting BBJ': 14, 'I² omitting TWB': 14, 'I² omitting TPMI': 14,
    '|KoGES beta / EAS panel beta|': 20,
    'BBJ Neff': 12, 'TPMI Neff': 12,
    'KoGES n > each cohort n': 18, 'KoGES n > EAS meta n': 18,
    'Analysis classification': 44,
}

def make_col_map(eaf_map, ext_n_cols, ext_n_cmp_cols=None):
    m = dict(COL_MAP_BASE)
    for col, cohort in eaf_map.items():
        cl = cohort.lower()
        m[f'{cl}_maf']     = f'{cohort} MAF'
        m[f'{cl}_maf_cat'] = f'{cohort} MAF category'
    for nc in ext_n_cols:
        m[nc] = f'{nc.split("_")[0].upper()} N (total)'
    # neff display cols already in COL_MAP_BASE; add width entries
    return m

def write_ed7_workbook(path: Path, df, col_map, sheet_name: str):
    import openpyxl
    from openpyxl.styles import PatternFill, Font, Alignment, Border, Side
    from openpyxl.utils import get_column_letter

    wb  = openpyxl.Workbook()
    ws  = wb.active
    ws.title = sheet_name[:31]
    avail = {k: v for k, v in col_map.items() if k in df.columns}
    out   = df[[k for k in avail]].rename(columns=avail).reset_index(drop=True)
    THIN  = Side(style='thin', color='D9D9D9')
    BDR   = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
    fills = {k: PatternFill('solid', fgColor=v) for k, v in CLF_COLORS.items()}
    CLC   = 'Analysis classification'

    for ci, cn in enumerate(out.columns, 1):
        cell = ws.cell(row=1, column=ci, value=cn)
        cell.font      = Font(name='Calibri', size=9, bold=True, color='FFFFFF')
        cell.fill      = PatternFill('solid', fgColor='1F4E79')
        cell.border    = BDR
        cell.alignment = Alignment(wrap_text=True, horizontal='center', vertical='center')
    ws.row_dimensions[1].height = 40

    for ri, row in enumerate(out.itertuples(index=False), 2):
        rd   = dict(zip(out.columns, row))
        fill = fills.get(str(rd.get(CLC, '')))
        for ci, (cn, val) in enumerate(rd.items(), 1):
            cell = ws.cell(row=ri, column=ci)
            try:
                if isinstance(val, bool):
                    cell.value = str(val)
                elif isinstance(val, float) and not np.isnan(val):
                    if any(x in cn for x in (' P', 'Q P', 'KoGES P')):
                        cell.value = val; cell.number_format = '0.00E+00'
                    elif any(x in cn for x in ('beta','SE','MAF','Shrinkage','ratio')):
                        cell.value = round(val, 6); cell.number_format = '0.000000'
                    elif 'I²' in cn or 'beta|' in cn:
                        cell.value = round(val, 2)
                    else:
                        cell.value = val
                elif pd.isna(val):
                    cell.value = ''
                else:
                    cell.value = val
            except Exception:
                cell.value = str(val) if val is not None else ''
            cell.font      = Font(name='Calibri', size=8)
            cell.border    = BDR
            cell.alignment = Alignment(vertical='top')
            if fill: cell.fill = fill

    for ci, cn in enumerate(out.columns, 1):
        ws.column_dimensions[get_column_letter(ci)].width = WIDTHS.get(cn, 12)
    ws.freeze_panes = 'A2'
    n_2c = (out.get(CLC, pd.Series()) == CLF['2c']).sum()
    print(f"  {sheet_name}: {len(out)} EAS-novel | unexplained candidates (2c): {n_2c}")
    path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(path)
    print(f"  Saved: {path}")

def _copy_sheet_values(ws_src, ws_dst):
    for row in ws_src.iter_rows(values_only=True):
        ws_dst.append(list(row))

def update_combined_extended_data(combined_xlsx, ed7a_path, ed7b_path):
    import openpyxl
    if not combined_xlsx.exists():
        raise FileNotFoundError(f"Combined workbook not found: {combined_xlsx}")
    wb = openpyxl.load_workbook(combined_xlsx)
    for name in ("Extended Data7A", "Extended Data7B"):
        if name in wb.sheetnames:
            del wb[name]
    if "Extended Data6B" in wb.sheetnames:
        insert_at = wb.sheetnames.index("Extended Data6B") + 1
    elif "Extended Data8A" in wb.sheetnames:
        insert_at = wb.sheetnames.index("Extended Data8A")
    else:
        insert_at = len(wb.sheetnames)
    for path, sheet_name, offset in (
        (ed7a_path, "Extended Data7A", 0),
        (ed7b_path, "Extended Data7B", 1),
    ):
        src    = openpyxl.load_workbook(path, data_only=True)
        ws_new = wb.create_sheet(title=sheet_name, index=insert_at + offset)
        _copy_sheet_values(src.active, ws_new)
        src.close()
    wb.save(combined_xlsx)
    print(f"  Updated: {combined_xlsx}")

# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ed-xlsx',           type=Path, default=DEFAULT_ED_XLSX)
    parser.add_argument('--tsv-dir',           type=Path, default=DEFAULT_TSV_DIR)
    parser.add_argument('--output-ed7a',       type=Path, default=DEFAULT_ED_OUT / "Extended Data7A.xlsx")
    parser.add_argument('--output-ed7b',       type=Path, default=DEFAULT_ED_OUT / "Extended Data7B.xlsx")
    parser.add_argument('--output-supp-cont',  type=Path, default=DEFAULT_SUPP_OUT / "SuppleTable8_cont.csv")
    parser.add_argument('--output-supp-cate',  type=Path, default=DEFAULT_SUPP_OUT / "SuppleTable8_cate.csv")
    parser.add_argument('--neff-cate',         type=Path, default=DEFAULT_NEFF_CATE,
                        help='TSV with BBJ/TPMI neff for binary traits')
    parser.add_argument('--combined-xlsx',     type=Path, default=DEFAULT_ED_XLSX)
    parser.add_argument('--skip-combined-update', action='store_true')
    args = parser.parse_args()

    print("Building 7A (continuous)...")
    cont, cont_n, cont_eaf, cont_cmp_n = build(
        args.ed_xlsx, args.tsv_dir, 'continuous')
    print("Building 7B (binary)...")
    cate, cate_n, cate_eaf, cate_cmp_n = build(
        args.ed_xlsx, args.tsv_dir, 'binary',
        neff_cate_path=args.neff_cate)

    print("\nWriting Extended Data 7A/7B...")
    write_ed7_workbook(args.output_ed7a, cont,
                       make_col_map(cont_eaf, cont_n, cont_cmp_n), 'Extended Data7A')
    write_ed7_workbook(args.output_ed7b, cate,
                       make_col_map(cate_eaf, cate_n, cate_cmp_n), 'Extended Data7B')

    if not args.skip_combined_update:
        print("\nUpdating combined ExtendedData.xlsx...")
        update_combined_extended_data(
            args.combined_xlsx, args.output_ed7a, args.output_ed7b)

    print("\nWriting Supplementary Table 8...")
    args.output_supp_cont.parent.mkdir(parents=True, exist_ok=True)
    args.output_supp_cate.parent.mkdir(parents=True, exist_ok=True)
    make_summary(cont, 'Continuous').to_csv(args.output_supp_cont, index=False)
    make_summary(cate, 'Binary').to_csv(args.output_supp_cate, index=False)
    print(f"  Saved: {args.output_supp_cont}")
    print(f"  Saved: {args.output_supp_cate}")

    print("\n── Classification summary ──")
    for label, df in [('7A Continuous', cont), ('7B Binary', cate)]:
        print(f"\n{label} (total {len(df)}):")
        for code in ['NA','1a','1b','1c','2a','2b','2c']:
            n = (df['classification_code']==code).sum()
            print(f"  {CLF[code]:<42}: {n}")

if __name__ == '__main__':
    main()
