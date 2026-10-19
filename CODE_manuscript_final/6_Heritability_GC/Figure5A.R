#!/usr/bin/env Rscript
# Continuous-trait SNP-heritability concordance.
# x = KoGES LDSC h2 (EAS LD scores). y = external-cohort LDSC h2.
# Estimates below 0 are plotted as 0. Estimates above 1 are plotted as 1.
# EAS external h2 uses EAS LD scores. EUR external h2 uses EUR LD scores.
#
# Inputs, all under recompute_260921:
#   EAS_LDSC/EAS_h2_rg_summary.tsv          (h2 log prefixes and match keys)
#   EUR_LDSC_h2_99/EUR_h2_summary_99.tsv    (external EUR h2)
#   manifest/EUR_matching_manifest.tsv      (EUR trait -> KoGES phenotype)
#   EAS_LDSC/KoGES_h2_missing10_manifest.tsv (KoGES h2 for 10 EUR-only traits)
# KoGES h2 is read from the logs named by k_h2_prefix.
# External EAS h2 is read from the logs named by ext_h2_prefix.
# No file under /data4/joh/openkoges/Heritability_GC is read.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

RUN_ROOT <- "/data4/sjim/openKoGES_revision/Heritability_GC/recompute_260921"
EAS_SUM <- file.path(RUN_ROOT, "EAS_LDSC", "EAS_h2_rg_summary.tsv")
EUR_SUM <- file.path(RUN_ROOT, "EUR_LDSC_h2_99", "EUR_h2_summary_99.tsv")
EUR_MAN <- file.path(RUN_ROOT, "manifest", "EUR_matching_manifest.tsv")
OUT_DIR <- file.path(RUN_ROOT, "plots_cross_ancestry_260927")
FIG_DIR <- "/data4/joh/openkoges/figure_new0923/output"
ED_DIR <- "/data4/joh/openkoges/figure_new0923/extended_data/cross_ancestry"
FORBIDDEN <- "/data4/joh/openkoges/Heritability_GC"

EAS_COHORTS <- c("BBJ2021", "TWB", "TPMI")
EUR_COHORTS <- c("UKB", "FinnGen", "MVP")
EAS_N <- c(BBJ2021 = 30L, TWB = 28L, TPMI = 22L)
EUR_N <- c(UKB = 37L, FinnGen = 35L, MVP = 27L)
PAPER_LABEL <- c(
  ALBUMIN = "Albumin", ALP = "Alkaline phosphatase", ALT = "Alanine aminotransferase",
  AST = "Aspartate aminotransferase", APOB = "Apolipoprotein B", BMI = "Body mass index",
  BUN = "Blood urea nitrogen", CPEPTIDE = "C-peptide", CREATINE = "Serum creatinine",
  CYSTATINC = "Cystatin C", DBP = "Diastolic blood pressure", eGFR = "eGFR",
  FFRbest = "FEV1/FVC", FEV1best = "FEV1", FEV1pp = "FEV1 % predicted",
  FOLATE = "Folate", FVCbest = "FVC", GLUO = "Glucose", HB = "Hemoglobin",
  HBA1C = "HbA1c", HCT = "Hematocrit", HDL = "HDL cholesterol", HEIGHT = "Height",
  HIP = "Hip circumference", HOMOCYST = "Homocysteine", INSUL = "Insulin",
  IRON = "Iron", LDL = "LDL cholesterol", MCH = "MCH", MCHC = "MCHC", MCV = "MCV",
  PHOSPHATE = "Phosphate", PLAT = "Platelet count", PULSE = "Pulse rate",
  RBC_B = "Red blood cell count", R_GTP = "Gamma-glutamyl transferase",
  SBP = "Systolic blood pressure", T_BIL = "Total bilirubin", TCHL = "Total cholesterol",
  TG = "Triglycerides", URICACID = "Uric acid", WAIST = "Waist circumference",
  WBC_B = "White blood cell count", WEIGHT = "Body weight", WHR = "Waist-to-hip ratio"
)

COHORT_LABEL <- c(
  BBJ2021 = "BBJ 2021", TWB = "TWB", TPMI = "TPMI",
  UKB = "UKB", FinnGen = "FinnGen", MVP = "MVP"
)
COHORT_COLOR <- c(
  BBJ2021 = "#009E73", TWB = "#56B4E9", TPMI = "#E69F00",
  UKB = "#0072B2", FinnGen = "#CC79A7", MVP = "#D55E00"
)
COHORT_SHAPE <- c(
  BBJ2021 = 16, TWB = 17, TPMI = 15,
  UKB = 16, FinnGen = 17, MVP = 15
)

assert_under_recompute <- function(paths, what) {
  paths <- as.character(paths)
  if (any(!nzchar(paths) | is.na(paths))) stop(what, " has an empty path")
  if (any(!startsWith(normalizePath(paths, mustWork = FALSE), RUN_ROOT))) {
    stop(what, " is outside recompute_260921")
  }
  if (any(grepl(FORBIDDEN, paths, fixed = TRUE))) stop(what, " points at a previous result")
}

paper_label <- function(koges_trait) {
  key <- sub("^B_", "", as.character(koges_trait))
  out <- unname(PAPER_LABEL[key])
  out[is.na(out)] <- key[is.na(out)]
  out
}

parse_h2_log <- function(prefix) {
  log <- paste0(prefix, ".log")
  if (!file.exists(log)) stop("Missing LDSC h2 log: ", log)
  lines <- readLines(log, warn = FALSE)
  hit <- grep("Total Observed scale h2:", lines, value = TRUE)
  if (length(hit) != 1L) stop("Expected one h2 line in ", log, "; found ", length(hit))
  m <- regexec("([-+0-9.eE]+)[[:space:]]+\\(([-+0-9.eE]+)\\)", hit)
  reg <- regmatches(hit, m)[[1]]
  if (length(reg) != 3L) stop("Could not parse h2 and SE from ", log, ": ", hit)
  data.table(
    h2 = as.numeric(reg[2]),
    se = as.numeric(reg[3]),
    out_of_bounds = any(grepl("h2[[:space:]]+out of bounds|h2's was out of bounds", lines))
  )
}

parse_h2_string <- function(x) {
  m <- regexec("^([-+0-9.eE]+)[[:space:]]+\\(([-+0-9.eE]+)\\)$", x)
  reg <- regmatches(x, m)
  bad <- vapply(reg, length, integer(1)) != 3L
  if (any(bad)) stop("Could not parse EUR h2 string: ", paste(x[bad], collapse = ", "))
  data.table(
    h2 = vapply(reg, function(z) as.numeric(z[2]), numeric(1)),
    se = vapply(reg, function(z) as.numeric(z[3]), numeric(1))
  )
}

for (f in c(EAS_SUM, EUR_SUM, EUR_MAN)) {
  if (!file.exists(f)) stop("Missing: ", f)
  assert_under_recompute(f, f)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

eas <- fread(EAS_SUM, showProgress = FALSE)
need <- c("trait", "cohort", "koges_trait", "k_h2_prefix", "k_h2_status",
          "ext_h2_prefix", "ext_h2_status", "rg_status")
if (length(setdiff(need, names(eas)))) stop("EAS summary missing columns")
if (nrow(eas) != 80L || any(eas$rg_status != "success") ||
    any(eas$k_h2_status != "success") || any(eas$ext_h2_status != "success")) {
  stop("EAS check failed: expected 80 successful h2 and rg rows")
}
eas_counts <- setNames(eas[, .N, by = cohort][match(names(EAS_N), cohort), N], names(EAS_N))
if (!identical(eas_counts, EAS_N)) stop("EAS cohort counts are not BBJ2021 30, TWB 28, TPMI 22")
assert_under_recompute(eas$k_h2_prefix, "EAS KoGES h2 prefix")
assert_under_recompute(eas$ext_h2_prefix, "EAS external h2 prefix")

koges_key <- unique(eas[, .(koges_trait, trait, k_h2_prefix)])
if (anyDuplicated(koges_key, by = "koges_trait")) stop("A KoGES phenotype maps to more than one h2 log")
koges_h2 <- rbindlist(lapply(seq_len(nrow(koges_key)), function(i) {
  parse_h2_log(koges_key$k_h2_prefix[i])
}))
koges_key <- cbind(koges_key, koges_h2)
setnames(koges_key, c("h2", "se", "out_of_bounds"), c("h2_koges", "se_koges", "koges_oob"))
extra_man <- file.path(RUN_ROOT, "EAS_LDSC", "KoGES_h2_missing10_manifest.tsv")
if (!file.exists(extra_man)) stop("Missing: ", extra_man)
assert_under_recompute(extra_man, extra_man)
extra <- fread(extra_man, showProgress = FALSE)
if (nrow(extra) != 10L || any(extra$h2_status != "success")) {
  stop("Expected 10 successful KoGES h2 logs in KoGES_h2_missing10_manifest.tsv")
}
if (any(extra$koges_trait %chin% koges_key$koges_trait)) stop("A missing-10 trait is already in the EAS KoGES h2 set")
assert_under_recompute(extra$h2_prefix, "extra KoGES h2 prefix")
extra_h2 <- rbindlist(lapply(extra$h2_prefix, parse_h2_log))
extra_key <- cbind(extra[, .(koges_trait, trait = koges_trait, k_h2_prefix = h2_prefix)], extra_h2)
setnames(extra_key, c("h2", "se", "out_of_bounds"), c("h2_koges", "se_koges", "koges_oob"))
koges_key <- rbind(koges_key, extra_key, use.names = TRUE)
if (any(!is.finite(koges_key$h2_koges) | !is.finite(koges_key$se_koges) | koges_key$se_koges <= 0)) {
  stop("KoGES h2 parse produced a non-finite estimate")
}

ext_h2 <- rbindlist(lapply(eas$ext_h2_prefix, parse_h2_log))
eas_pts <- cbind(
  eas[, .(panel = "EAS", cohort, trait, koges_trait, N_koges, N_external)],
  ext_h2
)
eas_pts <- merge(
  eas_pts,
  koges_key[, .(koges_trait, h2_koges, se_koges, koges_oob)],
  by = "koges_trait",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(eas_pts$h2_koges)) stop("An EAS pair has no KoGES h2")

eur <- fread(EUR_SUM, showProgress = FALSE)
eur_need <- c("cohort", "trait", "h2_status", "h2")
if (length(setdiff(eur_need, names(eur)))) stop("EUR h2 summary missing columns")
if (nrow(eur) != 99L || any(eur$h2_status != "success")) {
  stop("EUR check failed: expected 99 successful h2 rows")
}
eur_counts <- setNames(eur[, .N, by = cohort][match(names(EUR_N), cohort), N], names(EUR_N))
if (!identical(eur_counts, EUR_N)) stop("EUR cohort counts are not UKB 37, FinnGen 35, MVP 27")
eur_man <- fread(EUR_MAN, showProgress = FALSE)
eur_man <- unique(eur_man[, .(cohort, trait, koges_trait, N_koges, N_external)])
if (anyDuplicated(eur_man, by = c("cohort", "trait"))) stop("EUR manifest has duplicate cohort-trait keys")
eur <- merge(eur, eur_man, by = c("cohort", "trait"), all.x = TRUE, sort = FALSE)
if (anyNA(eur$koges_trait)) stop("EUR h2 rows missing from the matching manifest")
eur_parsed <- parse_h2_string(eur$h2)
eur_pts <- cbind(
  eur[, .(panel = "EUR", cohort, trait, koges_trait, N_koges, N_external)],
  eur_parsed
)
setnames(eur_pts, c("h2", "se"), c("h2_ext", "se_ext"))
eur_pts[, out_of_bounds := FALSE]
eur_pts <- merge(
  eur_pts,
  koges_key[, .(koges_trait, h2_koges, se_koges, koges_oob)],
  by = "koges_trait",
  all.x = TRUE,
  sort = FALSE
)

eas_pts[, h2_ext := h2]
eas_pts[, se_ext := se]
eas_pts[, reason := fifelse(is.finite(h2_koges) & is.finite(h2_ext) & is.finite(se_koges) & is.finite(se_ext), "plotted", "non_finite_h2")]
eur_pts[, reason := fifelse(is.finite(h2_koges) & is.finite(h2_ext) & is.finite(se_koges) & is.finite(se_ext), "plotted", "no_KoGES_h2")]

pts <- rbind(
  eas_pts[, .(panel, cohort, trait, koges_trait, N_koges, N_external, h2_koges, se_koges, h2_ext, se_ext, reason)],
  eur_pts[, .(panel, cohort, trait, koges_trait, N_koges, N_external, h2_koges, se_koges, h2_ext, se_ext, reason)],
  use.names = TRUE
)
pts[, trait_label := paper_label(koges_trait)]
pts[, cohort_label := unname(COHORT_LABEL[cohort])]
pts[, h2_koges_ldsc := h2_koges]
pts[, h2_ext_ldsc := h2_ext]
pts[, h2_koges := pmin(1, pmax(0, h2_koges))]
pts[, h2_ext := pmin(1, pmax(0, h2_ext))]
pts[, diff := h2_ext - h2_koges]

plotted <- pts[reason == "plotted"]
excluded <- pts[reason != "plotted"]
if (nrow(plotted[panel == "EAS"]) != 80L) {
  stop("Expected 80 plotted EAS pairs, found ", nrow(plotted[panel == "EAS"]))
}
if (nrow(plotted[panel == "EUR"]) != 99L) {
  stop("Expected 99 plotted EUR pairs, found ", nrow(plotted[panel == "EUR"]))
}
if (nrow(excluded) != 0L) stop("Expected no excluded pairs, found ", nrow(excluded))
lim <- c(0, 0.5)

panel_plot <- function(panel_name, cohorts, subtitle) {
  d <- plotted[panel == panel_name]
  d[, cohort := factor(cohort, levels = cohorts)]
  labs_named <- COHORT_LABEL[cohorts]
  ggplot(d, aes(x = h2_koges, y = h2_ext, color = cohort, shape = cohort)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#8A8A8A", linewidth = 0.4) +
    geom_errorbar(
      aes(ymin = h2_ext - se_ext, ymax = h2_ext + se_ext),
      width = 0, linewidth = 0.3, alpha = 0.85, show.legend = FALSE
    ) +
    geom_errorbar(
      aes(xmin = h2_koges - se_koges, xmax = h2_koges + se_koges),
      width = 0, linewidth = 0.3, alpha = 0.85, orientation = "y", show.legend = FALSE
    ) +
    geom_point(size = 2.4, alpha = 0.95) +
    scale_color_manual(values = COHORT_COLOR[cohorts], labels = labs_named, name = "Cohort") +
    scale_shape_manual(values = COHORT_SHAPE[cohorts], labels = labs_named, name = "Cohort") +
    scale_x_continuous(limits = lim, breaks = seq(0, 0.5, by = 0.1), expand = c(0, 0)) +
    scale_y_continuous(limits = lim, breaks = seq(0, 0.5, by = 0.1), expand = c(0, 0)) +
    coord_fixed(clip = "on") +
    labs(
      title = panel_name,
      subtitle = subtitle,
      x = "SNP-based heritability in KoGES",
      y = "SNP-based heritability in the external cohort"
    ) +
    theme_classic(base_family = "Liberation Sans", base_size = 12) +
    theme(
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 11),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 14, hjust = 0, color = "black"),
      plot.subtitle = element_text(size = 10, hjust = 0, color = "#333333", margin = margin(b = 6)),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 12, 4, 8)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

p <- panel_plot("EAS", EAS_COHORTS, "KoGES vs BBJ 2021, TWB, and TPMI") +
  panel_plot("EUR", EUR_COHORTS, "KoGES vs UKB, FinnGen, and MVP") +
  plot_layout(ncol = 2)

stem <- "KoGES_continuous_LDSC_heritability_concordance"
pdf <- file.path(OUT_DIR, paste0(stem, ".pdf"))
png <- file.path(OUT_DIR, paste0(stem, ".png"))
ggsave(pdf, p, width = 12.4, height = 6.35, device = cairo_pdf, bg = "white")
ggsave(png, p, width = 12.4, height = 6.35, dpi = 320, bg = "white")

points_out <- file.path(OUT_DIR, paste0(stem, "_points.tsv"))
excluded_out <- file.path(OUT_DIR, paste0(stem, "_excluded.tsv"))
fwrite(plotted[, .(panel, cohort, trait, koges_trait, trait_label, N_koges, N_external, h2_koges, se_koges, h2_ext, se_ext, h2_koges_ldsc, h2_ext_ldsc, diff)],
       points_out, sep = "\t")
fwrite(excluded[, .(panel, cohort, trait, koges_trait, trait_label, N_koges, N_external, h2_koges_ldsc, se_koges, h2_ext_ldsc, se_ext, reason)],
       excluded_out, sep = "\t")

# Rows drawn in the figure. One row is one cohort-phenotype point.
extended <- plotted[, .(
  panel,
  cohort,
  phenotype = trait_label,
  koges_trait,
  analysis_trait = trait,
  N_koges,
  N_external,
  h2_koges,
  se_koges,
  h2_external = h2_ext,
  se_external = se_ext,
  h2_koges_ldsc,
  h2_external_ldsc = h2_ext_ldsc,
  h2_external_minus_koges = diff,
  koges_ld = "EAS",
  external_ld = fifelse(panel == "EAS", "EAS", "EUR"),
  scale = "LDSC observed scale"
)]
if (anyNA(extended$N_koges) || anyNA(extended$N_external)) stop("Extended data is missing a sample size")
setorder(extended, panel, cohort, phenotype)
if (nrow(extended) != 179L) stop("Extended data should contain the 179 plotted points, found ", nrow(extended))
extended_out <- file.path(OUT_DIR, paste0(stem, "_extended_data.tsv"))
fwrite(extended, extended_out, sep = "\t")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(Sys.getenv("FIGURE5_OUT", unset = "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure5"), recursive = TRUE, showWarnings = FALSE)
ed_copy <- file.path(ED_DIR, basename(extended_out))
copied <- c(
  file.copy(c(pdf, png), FIG_DIR, overwrite = TRUE),
  file.copy(png, file.path(Sys.getenv("FIGURE5_OUT", unset = "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure5"), "Figure5A.png"), overwrite = TRUE),
  file.copy(c(extended_out, points_out, excluded_out), ED_DIR, overwrite = TRUE)
)
if (!all(copied)) stop("Failed to copy the figure or extended data")

cat("Plotted EAS", nrow(plotted[panel == "EAS"]), " EUR", nrow(plotted[panel == "EUR"]), "\n")
print(plotted[, .(n = .N, median_diff = median(diff), mean_diff = mean(diff)), by = .(panel, cohort)])
cat("Excluded\n")
print(excluded[, .(panel, cohort, trait_label, h2_ext, reason)])
cat("Wrote\n", pdf, "\n", png, "\n", extended_out, "\n", ed_copy, "\n", sep = "")
