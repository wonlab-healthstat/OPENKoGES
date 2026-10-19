#!/usr/bin/env Rscript
# Publication heatmaps for the verified recompute_260921 cross-ancestry results.
# EAS LDSC rg and EUR POPCORN pge are different quantities. The separate
# heatmaps keep their own files. The combined figure places them side by side
# on a shared trait axis, with a separate color scale in each panel.
#
# Inputs are only:
#   recompute_260921/EAS_LDSC/EAS_h2_rg_summary.tsv
#   recompute_260921/EUR_POPCORN_pge_99/KoGES_EUR_POPCORN_pge_summary_99.tsv
#   recompute_260921/manifest/EAS_matching_manifest.tsv   (phenotype labels)
#   recompute_260921/manifest/EUR_matching_manifest.tsv   (phenotype labels)
# LDSC point estimates are read from the logs named by rg_prefix in the EAS
# summary. The summary rg / rg_se columns are the section header produced by
# stat() in 01_EAS_all_h2_rg_revision.R, not the numeric estimate.
# No file under /data4/joh/openkoges/Heritability_GC is read.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

RUN_ROOT <- "/data4/sjim/openKoGES_revision/Heritability_GC/recompute_260921"
EAS_SUM <- file.path(RUN_ROOT, "EAS_LDSC", "EAS_h2_rg_summary.tsv")
EUR_SUM <- file.path(RUN_ROOT, "EUR_POPCORN_pge_99", "KoGES_EUR_POPCORN_pge_summary_99.tsv")
EAS_MAN <- file.path(RUN_ROOT, "manifest", "EAS_matching_manifest.tsv")
EUR_MAN <- file.path(RUN_ROOT, "manifest", "EUR_matching_manifest.tsv")
OUT_DIR <- file.path(RUN_ROOT, "plots_cross_ancestry_260927")
FORBIDDEN <- "/data4/joh/openkoges/Heritability_GC"

EAS_COHORTS <- c("BBJ2021", "TWB", "TPMI")
EUR_COHORTS <- c("UKB", "FinnGen", "MVP")
EAS_N <- c(BBJ2021 = 30L, TWB = 28L, TPMI = 22L)
EUR_N <- c(UKB = 37L, FinnGen = 35L, MVP = 27L)

# Manuscript phenotype names, keyed by the KoGES phenotype id without the B_ prefix.
PAPER_LABEL <- c(
  ALBUMIN = "Albumin",
  ALP = "Alkaline phosphatase",
  ALT = "Alanine aminotransferase",
  AST = "Aspartate aminotransferase",
  APOB = "Apolipoprotein B",
  BMI = "Body mass index",
  BUN = "Blood urea nitrogen",
  CPEPTIDE = "C-peptide",
  CREATINE = "Serum creatinine",
  CYSTATINC = "Cystatin C",
  DBP = "Diastolic blood pressure",
  eGFR = "eGFR",
  FFRbest = "FEV1/FVC",
  FEV1best = "FEV1",
  FEV1pp = "FEV1 % predicted",
  FOLATE = "Folate",
  FVCbest = "FVC",
  GLUO = "Glucose",
  HB = "Hemoglobin",
  HBA1C = "HbA1c",
  HCT = "Hematocrit",
  HDL = "HDL cholesterol",
  HEIGHT = "Height",
  HIP = "Hip circumference",
  HOMOCYST = "Homocysteine",
  INSUL = "Insulin",
  IRON = "Iron",
  LDL = "LDL cholesterol",
  MCH = "MCH",
  MCHC = "MCHC",
  MCV = "MCV",
  PHOSPHATE = "Phosphate",
  PLAT = "Platelet count",
  PULSE = "Pulse rate",
  RBC_B = "Red blood cell count",
  R_GTP = "Gamma-glutamyl transferase",
  SBP = "Systolic blood pressure",
  T_BIL = "Total bilirubin",
  TCHL = "Total cholesterol",
  TG = "Triglycerides",
  URICACID = "Uric acid",
  WAIST = "Waist circumference",
  WBC_B = "White blood cell count",
  WEIGHT = "Body weight",
  WHR = "Waist-to-hip ratio"
)

GROUP <- c(
  BMI = "Anthropometric", WEIGHT = "Anthropometric", HEIGHT = "Anthropometric",
  HIP = "Anthropometric", WAIST = "Anthropometric", WHR = "Anthropometric",
  SBP = "Cardiovascular", DBP = "Cardiovascular", PULSE = "Cardiovascular",
  HB = "Hematologic", HCT = "Hematologic", RBC_B = "Hematologic",
  WBC_B = "Hematologic", PLAT = "Hematologic", MCV = "Hematologic",
  MCH = "Hematologic", MCHC = "Hematologic", IRON = "Hematologic",
  FOLATE = "Hematologic",
  BUN = "Kidney", CREATINE = "Kidney", eGFR = "Kidney",
  CYSTATINC = "Kidney", URICACID = "Kidney",
  AST = "Liver", ALT = "Liver", R_GTP = "Liver", ALP = "Liver",
  T_BIL = "Liver", ALBUMIN = "Liver",
  FEV1best = "Lung", FEV1pp = "Lung", FVCbest = "Lung", FFRbest = "Lung",
  GLUO = "Metabolic", HBA1C = "Metabolic", INSUL = "Metabolic",
  CPEPTIDE = "Metabolic", TCHL = "Metabolic", TG = "Metabolic",
  HDL = "Metabolic", LDL = "Metabolic", APOB = "Metabolic",
  HOMOCYST = "Metabolic",
  PHOSPHATE = "Bone"
)
GROUP_ORDER <- c(
  "Anthropometric", "Cardiovascular", "Hematologic", "Kidney",
  "Liver", "Lung", "Metabolic", "Bone", "Other"
)

assert_under_recompute <- function(paths, what) {
  paths <- as.character(paths)
  if (any(!nzchar(paths) | is.na(paths))) stop(what, " has an empty path")
  if (any(!startsWith(normalizePath(paths, mustWork = FALSE), RUN_ROOT))) {
    stop(what, " is outside recompute_260921")
  }
  if (any(grepl(FORBIDDEN, paths, fixed = TRUE))) stop(what, " points at a previous result")
}

require_cols <- function(dt, cols, what) {
  missing <- setdiff(cols, names(dt))
  if (length(missing)) stop(what, " is missing columns: ", paste(missing, collapse = ", "))
}

paper_label <- function(koges_trait, trait) {
  key <- sub("^B_", "", as.character(koges_trait))
  out <- unname(PAPER_LABEL[key])
  unmatched <- is.na(out)
  out[unmatched] <- as.character(trait)[unmatched]
  out
}

trait_group <- function(koges_trait) {
  key <- sub("^B_", "", as.character(koges_trait))
  out <- unname(GROUP[key])
  out[is.na(out)] <- "Other"
  out
}

# LDSC writes the estimate on the results-table row. The prose line
# "Genetic Correlation" is only a header, which is what the summary stored.
parse_ldsc_rg <- function(prefix, koges_munge, ext_munge) {
  log <- paste0(prefix, ".log")
  if (!file.exists(log)) stop("Missing LDSC log: ", log)
  lines <- readLines(log, warn = FALSE)
  tab <- lines[
    grepl(koges_munge, lines, fixed = TRUE) &
      grepl(ext_munge, lines, fixed = TRUE) &
      grepl("sumstats\\.gz[[:space:]]+/", lines)
  ]
  if (length(tab) != 1L) {
    stop("Expected one LDSC summary row in ", log, "; found ", length(tab))
  }
  nums <- strsplit(trimws(tab), "[[:space:]]+")[[1]]
  nums <- nums[!grepl("sumstats\\.gz$", nums)]
  if (length(nums) < 4L) stop("Could not parse rg, se, z, p from ", log)
  to_num <- function(x) if (grepl("^(NA|NaN|nan)$", x)) NA_real_ else as.numeric(x)
  est <- to_num(nums[1])
  rg_oob <- any(grepl("rg out of bounds|rg was out of bounds", lines))
  h2_oob <- any(grepl("h2[[:space:]]+out of bounds|h2's was out of bounds", lines))
  failed <- any(grepl("ERROR computing rg", lines))
  note <- if (rg_oob) "rg_out_of_bounds" else if (h2_oob) "h2_out_of_bounds" else if (failed) "rg_failed" else "ok"
  data.table(
    estimate = est,
    se = to_num(nums[2]),
    z = to_num(nums[3]),
    p_log = to_num(nums[4]),
    out_of_bounds = rg_oob && is.finite(est),
    estimate_note = note
  )
}

attach_labels <- function(dt, manifest) {
  man <- unique(manifest[, .(cohort, trait, koges_manifest = koges_trait)])
  if (anyDuplicated(man, by = c("cohort", "trait"))) {
    stop("Matching manifest has duplicate cohort-trait keys")
  }
  out <- merge(dt, man, by = c("cohort", "trait"), all.x = TRUE, sort = FALSE)
  if (anyNA(out$koges_manifest)) {
    bad <- out[is.na(koges_manifest), paste(cohort, trait, sep = ":")]
    stop("Summary pairs missing from the matching manifest: ", paste(bad, collapse = ", "))
  }
  if ("koges_trait" %in% names(out)) {
    if (any(out$koges_trait != out$koges_manifest)) {
      stop("Summary koges_trait does not match the matching manifest")
    }
    out[, koges_manifest := NULL]
  } else {
    setnames(out, "koges_manifest", "koges_trait")
  }
  out[, trait_label := paper_label(koges_trait, trait)]
  out[, trait_group := trait_group(koges_trait)]
  labels <- unique(out[, .(trait, trait_label)])
  if (anyDuplicated(labels, by = "trait_label")) {
    stop("Two analysis traits collapsed to the same display label")
  }
  out
}

order_traits <- function(dt) {
  info <- unique(dt[, .(trait_label, trait_group)])
  info[, trait_group := factor(trait_group, levels = GROUP_ORDER)]
  setorder(info, trait_group, trait_label)
  info$trait_label
}

complete_grid <- function(dt, cohorts) {
  trait_levels <- order_traits(dt)
  meta <- unique(dt[, .(trait_label, trait, koges_trait, trait_group)])
  if (anyDuplicated(meta, by = "trait_label")) stop("Display label is not unique")
  grid <- CJ(trait_label = trait_levels, cohort = cohorts, unique = TRUE)
  out <- merge(grid, meta, by = "trait_label", all.x = TRUE, sort = FALSE)
  keep <- dt[, .(
    trait_label, cohort, estimate, se, p, out_of_bounds, differs_from_1, estimate_note
  )]
  out <- merge(out, keep, by = c("trait_label", "cohort"), all.x = TRUE, sort = FALSE)
  paired <- unique(dt[, .(trait_label, cohort, paired = TRUE)])
  out <- merge(out, paired, by = c("trait_label", "cohort"), all.x = TRUE, sort = FALSE)
  out[, in_result := fifelse(is.na(paired), FALSE, TRUE)]
  out[in_result == FALSE, estimate_note := "not_matched"]
  out[, paired := NULL]
  out[, trait_label := factor(trait_label, levels = rev(trait_levels))]
  out[, cohort := factor(cohort, levels = cohorts)]
  out[, trait_group := factor(trait_group, levels = GROUP_ORDER)]
  out[]
}

tile_label <- function(estimate, out_of_bounds, in_result) {
  oob <- rep_len(as.logical(out_of_bounds), length(estimate))
  oob[is.na(oob)] <- FALSE
  fifelse(
    !in_result, "",
    fifelse(
      !is.finite(estimate), "NA",
      fifelse(oob, sprintf("%.2f\u2020", estimate), sprintf("%.2f", estimate))
    )
  )
}

base_theme <- function() {
  theme_classic(base_family = "Liberation Sans", base_size = 11) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(face = "bold", size = 12, color = "black"),
      axis.text.y = element_text(size = 8.5, color = "black"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0, color = "#333333", margin = margin(b = 8)),
      plot.caption = element_text(
        size = 8, hjust = 0, color = "#444444", lineheight = 1.2,
        margin = margin(t = 8)
      ),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 14, 8, 8)
    )
}

fill_scale <- function(name) {
  scale_fill_gradient2(
    low = "#C75D5D",
    mid = "white",
    high = "#5A8FD8",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5),
    labels = c("-1", "-0.5", "0", "0.5", "1"),
    oob = squish,
    na.value = "#E6E6E6",
    name = name
  )
}

save_both <- function(plot, stem, width, height) {
  pdf <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  png <- file.path(OUT_DIR, paste0(stem, ".png"))
  ggsave(pdf, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(png, plot, width = width, height = height, dpi = 320, bg = "white")
  c(pdf, png)
}

# ------------------------------------------------------------------------------
# EAS: 80/80 LDSC rg successes
# ------------------------------------------------------------------------------

for (f in c(EAS_SUM, EUR_SUM, EAS_MAN, EUR_MAN)) {
  if (!file.exists(f)) stop("Missing: ", f)
  assert_under_recompute(f, f)
}

eas_raw <- fread(EAS_SUM, showProgress = FALSE)
require_cols(
  eas_raw,
  c("trait", "cohort", "koges_trait", "rg_status", "rg", "rg_se", "rg_p",
    "rg_prefix", "k_munge_prefix", "ext_munge_prefix"),
  "EAS summary"
)
if (nrow(eas_raw) != 80L || any(eas_raw$rg_status != "success")) {
  stop("EAS check failed: expected 80 rows with rg_status success, found ",
       nrow(eas_raw), " rows and ", sum(eas_raw$rg_status == "success"), " successes")
}
eas_counts <- eas_raw[, .N, by = cohort]
if (!identical(setNames(eas_counts[match(names(EAS_N), cohort), N], names(EAS_N)), EAS_N)) {
  stop("EAS cohort counts are not BBJ2021 30, TWB 28, TPMI 22")
}
assert_under_recompute(eas_raw$rg_prefix, "EAS rg_prefix")
assert_under_recompute(eas_raw$k_munge_prefix, "EAS KoGES munge prefix")
assert_under_recompute(eas_raw$ext_munge_prefix, "EAS external munge prefix")

if (any(suppressWarnings(!is.na(as.numeric(eas_raw$rg))))) {
  stop("EAS rg column contains numeric values; this script expects the known non-numeric header and reads estimates from the paired logs")
}

parsed <- rbindlist(lapply(seq_len(nrow(eas_raw)), function(i) {
  parse_ldsc_rg(
    eas_raw$rg_prefix[i],
    eas_raw$k_munge_prefix[i],
    eas_raw$ext_munge_prefix[i]
  )
}))
summary_p <- suppressWarnings(as.numeric(eas_raw$rg_p))
comparable <- is.finite(summary_p)
if (any(comparable & !is.finite(parsed$p_log))) {
  stop("EAS summary rg_p is numeric but the paired log has no p-value")
}
p_gap <- abs(summary_p - parsed$p_log)
if (any(comparable & p_gap > 1e-8 & p_gap / pmax(abs(summary_p), 1e-300) > 1e-4)) {
  stop("Parsed LDSC p-values disagree with rg_p in the EAS summary")
}
if (any(!comparable & is.finite(parsed$estimate) & !parsed$out_of_bounds)) {
  stop("Non-numeric EAS rg_p is not an out-of-bounds estimate")
}

eas <- cbind(
  eas_raw[, .(trait, cohort, koges_trait)],
  parsed[, .(estimate, se, p = p_log, out_of_bounds, estimate_note)]
)
eas[, differs_from_1 := NA]
eas <- attach_labels(eas, fread(EAS_MAN, select = c("cohort", "trait", "koges_trait"), showProgress = FALSE))
eas_plot <- complete_grid(eas, EAS_COHORTS)
eas_plot[, tile_label := tile_label(estimate, out_of_bounds, in_result)]

# ------------------------------------------------------------------------------
# EUR: 99/99 POPCORN fit successes
# ------------------------------------------------------------------------------

eur_raw <- fread(EUR_SUM, showProgress = FALSE)
require_cols(
  eur_raw,
  c("trait", "cohort", "fit_status", "pge", "pge_se", "p_pge_eq_1"),
  "EUR summary"
)
if (nrow(eur_raw) != 99L || any(eur_raw$fit_status != "success")) {
  stop("EUR check failed: expected 99 rows with fit_status success, found ",
       nrow(eur_raw), " rows and ", sum(eur_raw$fit_status == "success"), " successes")
}
eur_counts <- eur_raw[, .N, by = cohort]
if (!identical(setNames(eur_counts[match(names(EUR_N), cohort), N], names(EUR_N)), EUR_N)) {
  stop("EUR cohort counts are not UKB 37, FinnGen 35, MVP 27")
}

eur <- eur_raw[, .(
  trait, cohort,
  estimate = as.numeric(pge),
  se = as.numeric(pge_se),
  p = as.numeric(p_pge_eq_1)
)]
eur[, out_of_bounds := FALSE]
eur[, differs_from_1 := is.finite(p) & p < 0.05]
eur[, estimate_note := fifelse(is.finite(estimate), "ok", "no_finite_pge")]
eur <- attach_labels(eur, fread(EUR_MAN, select = c("cohort", "trait", "koges_trait"), showProgress = FALSE))
eur_plot <- complete_grid(eur, EUR_COHORTS)
eur_plot[, tile_label := tile_label(estimate, FALSE, in_result)]
eur_sig <- eur_plot[differs_from_1 %in% TRUE]

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!capabilities("cairo")) stop("cairo is required to embed Liberation Sans")

eas_caption <- paste(
  "Each tile is the LDSC genetic correlation, printed to two decimal places.",
  "Color runs from \u22121 (red) to 1 (blue). Labels keep the untruncated estimate.",
  "A dagger marks an rg that LDSC reported as out of bounds.",
  "NA is a matched pair with no finite rg. Empty grey tiles were not among the 80 pairs.",
  sep = "\n"
)
eur_caption <- paste(
  "Each tile is the POPCORN genetic-effect correlation, printed to two decimal places.",
  "Open red circles mark estimates that differ from 1 at P < 0.05.",
  "Color runs from \u22121 (red) to 1 (blue). Labels keep the untruncated estimate.",
  "NA is a matched pair with no finite pge. Empty grey tiles were not among the 99 pairs.",
  sep = "\n"
)

p_eas <- ggplot(eas_plot, aes(x = cohort, y = trait_label, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = tile_label), size = 3.1, color = "#1A1A1A", family = "Liberation Sans") +
  fill_scale("rg") +
  scale_x_discrete(position = "top", expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0)) +
  labs(
    title = "KoGES vs East Asian cohorts",
    subtitle = "LDSC genetic correlation (rg)",
    caption = eas_caption
  ) +
  base_theme() +
  guides(fill = guide_colorbar(
    barwidth = unit(6.2, "cm"),
    barheight = unit(0.35, "cm"),
    title.position = "top",
    title.hjust = 0.5
  ))

p_eur <- ggplot(eur_plot, aes(x = cohort, y = trait_label, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = tile_label), size = 3.1, color = "#1A1A1A", family = "Liberation Sans") +
  geom_point(
    data = eur_sig,
    shape = 1,
    color = "white",
    stroke = 1.5,
    size = 2.5,
    position = position_nudge(x = 0.34, y = 0.30)
  ) +
  geom_point(
    data = eur_sig,
    shape = 1,
    color = "#D62728",
    stroke = 0.7,
    size = 2.5,
    position = position_nudge(x = 0.34, y = 0.30)
  ) +
  fill_scale("pge") +
  scale_x_discrete(position = "top", expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0)) +
  labs(
    title = "KoGES vs European cohorts",
    subtitle = "POPCORN genetic-effect correlation (pge)",
    caption = eur_caption
  ) +
  base_theme() +
  guides(fill = guide_colorbar(
    barwidth = unit(6.2, "cm"),
    barheight = unit(0.35, "cm"),
    title.position = "top",
    title.hjust = 0.5
  ))

eas_files <- save_both(p_eas, "KoGES_EAS_LDSC_rg_heatmap", width = 9.0, height = 12.2)
eur_files <- save_both(p_eur, "KoGES_EUR_POPCORN_pge_heatmap", width = 9.0, height = 14.3)

eas_out <- eas_plot[order(trait_group, trait_label, cohort), .(
  trait, koges_trait, trait_label, trait_group, cohort,
  estimate, se, p, out_of_bounds, estimate_note, in_result, tile_label
)]
eur_out <- eur_plot[order(trait_group, trait_label, cohort), .(
  trait, koges_trait, trait_label, trait_group, cohort,
  estimate, se, p, differs_from_1, estimate_note, in_result, tile_label
)]
eas_tsv <- file.path(OUT_DIR, "EAS_LDSC_rg_plot_data.tsv")
eur_tsv <- file.path(OUT_DIR, "EUR_POPCORN_pge_plot_data.tsv")
fwrite(eas_out, eas_tsv, sep = "\t", quote = FALSE, na = "NA")
fwrite(eur_out, eur_tsv, sep = "\t", quote = FALSE, na = "NA")

# ------------------------------------------------------------------------------
# Combined figure: shared trait rows, separate panels and color scales
# ------------------------------------------------------------------------------

COMBINED_KEY_ORDER <- c(
  "BMI", "HEIGHT", "HIP", "WAIST", "WEIGHT", "WHR",
  "PHOSPHATE",
  "DBP", "PULSE", "SBP",
  "FOLATE", "IRON", "HB", "HCT", "MCH", "MCHC", "MCV", "PLAT", "RBC_B", "WBC_B",
  "BUN", "CREATINE", "CYSTATINC", "URICACID", "eGFR",
  "ALBUMIN", "ALP", "ALT", "AST", "R_GTP", "T_BIL",
  "FEV1best", "FEV1pp", "FFRbest", "FVCbest",
  "APOB", "CPEPTIDE", "INSUL", "HOMOCYST", "GLUO", "HBA1C", "HDL", "LDL", "TCHL", "TG"
)
COMBINED_GROUP_ORDER <- c(
  "Anthropometric", "Bone", "Cardiovascular", "Hematologic",
  "Kidney", "Liver", "Lung", "Metabolic", "Other"
)
COMBINED_GROUP_COLOR <- c(
  Anthropometric = "#009E73", Bone = "#F5C710", Cardiovascular = "#D55E00",
  Hematologic = "#CC79A7", Kidney = "#56B4E9", Liver = "#F0E442",
  Lung = "#E69F00", Metabolic = "#0072B2", Other = "#BDBDBD"
)
COMBINED_GROUP_TEXT <- c(
  Anthropometric = "white", Bone = "black", Cardiovascular = "white",
  Hematologic = "white", Kidney = "black", Liver = "black",
  Lung = "black", Metabolic = "white", Other = "black"
)
LABEL_OVERRIDE <- c(
  eGFR = "Estimated glomerular filtration rate",
  ALT = "Alanine transaminase",
  AST = "Aspartate transaminase",
  PULSE = "Pulse"
)

combined_rows <- rbind(
  eas[, .(panel = "EAS", trait, koges_trait, cohort, estimate, estimate_note)],
  eur[, .(panel = "EUR", trait, koges_trait, cohort, estimate, estimate_note)]
)
combined_rows[, usable := estimate_note == "ok" & is.finite(estimate)]
combined_rows[, key := sub("^B_", "", koges_trait)]
combined_rows[, trait_group := trait_group(koges_trait)]
combined_rows[, trait_label := paper_label(koges_trait, trait)]
combined_rows[, override_label := unname(LABEL_OVERRIDE[key])]
combined_rows[!is.na(override_label), trait_label := override_label]
combined_rows[, key_rank := match(key, COMBINED_KEY_ORDER)]
extra_keys <- sort(unique(combined_rows[is.na(key_rank), key]))
if (length(extra_keys)) {
  combined_rows[is.na(key_rank), key_rank := 1000L + match(key, extra_keys)]
}
trait_meta <- unique(combined_rows[, .(key, trait_label, trait_group, key_rank)])
if (anyDuplicated(trait_meta, by = "trait_label")) stop("Combined figure has duplicate trait labels")
trait_meta[, trait_group := factor(trait_group, levels = COMBINED_GROUP_ORDER)]
setorder(trait_meta, trait_group, key_rank, trait_label)
trait_top <- trait_meta$trait_label
n_trait <- length(trait_top)

cohort_label <- c(BBJ2021 = "BBJ 2021", TWB = "TWB", TPMI = "TPMI", UKB = "UKB", FinnGen = "FinnGen", MVP = "MVP")
build_panel <- function(panel_name, cohorts) {
  grid <- CJ(trait_label = trait_top, cohort = cohorts, unique = TRUE)
  keep <- combined_rows[panel == panel_name, .(trait_label, cohort, estimate, usable, estimate_note)]
  out <- merge(grid, keep, by = c("trait_label", "cohort"), all.x = TRUE, sort = FALSE)
  out[, trait_label := factor(trait_label, levels = rev(trait_top))]
  out[, cohort := factor(cohort, levels = cohorts)]
  out[, plot_estimate := fifelse(usable %in% TRUE, estimate, NA_real_)]
  out[, failed := usable %in% FALSE]
  out
}
FAIL_LABEL <- c(
  rg_out_of_bounds = "rg out\nof bounds",
  h2_out_of_bounds = "h2 out\nof bounds",
  rg_failed = "rg failed",
  no_finite_pge = "h2 out\nof bounds"
)
eas_panel <- build_panel("EAS", EAS_COHORTS)
eur_panel <- build_panel("EUR", EUR_COHORTS)
eas_panel[, cell_label := unname(FAIL_LABEL[as.character(estimate_note)])]
eur_panel[, cell_label := unname(FAIL_LABEL[as.character(estimate_note)])]
eas_text <- eas_panel[usable %in% TRUE]
eur_text <- eur_panel[usable %in% TRUE]
eas_text[, label_expr := sprintf("italic(r)[g] == %.2f", estimate)]
eur_text[, label_expr := sprintf("pge == %.2f", estimate)]
eas_text[, text_color := fifelse(estimate > -0.75 & estimate < 0.89, "#1A1A1A", "white")]
eur_text[, text_color := fifelse(estimate < 0.72, "#1A1A1A", "white")]

group_anno <- data.table(
  trait_label = factor(trait_top, levels = rev(trait_top)),
  trait_group = trait_meta$trait_group
)
group_anno[, idx := n_trait - seq_len(n_trait) + 1L]
group_pos <- group_anno[, .(y = mean(idx), n = .N, ymin = min(idx) - 0.5, ymax = max(idx) + 0.5), by = trait_group]
group_pos[, label_color := unname(COMBINED_GROUP_TEXT[as.character(trait_group)])]
group_pos[, label_size := fifelse(n <= 1L, 2.05, 3.1)]

panel_theme <- function(show_y) {
  theme_classic(base_family = "Liberation Sans", base_size = 11) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(face = "bold", size = 12, color = "black"),
      axis.text.y = if (show_y) element_text(size = 9, color = "black") else element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0, color = "black"),
      plot.subtitle = element_text(size = 10, hjust = 0, color = "#4D4D4D", margin = margin(b = 6)),
      plot.title.position = "panel",
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(6, 8, 4, if (show_y) 2 else 8)
    )
}
legend_bar <- guide_colorbar(
  barwidth = unit(4.6, "cm"), barheight = unit(0.32, "cm"),
  title.position = "top", title.hjust = 0.5
)

p_group <- ggplot() +
  geom_rect(
    data = group_pos,
    aes(xmin = 0.5, xmax = 1.5, ymin = ymin, ymax = ymax, fill = trait_group),
    color = NA
  ) +
  geom_text(
    data = group_pos,
    aes(x = 1, y = y, label = trait_group, color = label_color, size = label_size),
    inherit.aes = FALSE,
    angle = 90,
    fontface = "bold",
    family = "Liberation Sans"
  ) +
  scale_size_identity() +
  scale_fill_manual(values = COMBINED_GROUP_COLOR) +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, n_trait + 0.5), expand = c(0, 0)) +
  theme_void(base_family = "Liberation Sans") +
  theme(legend.position = "none", plot.margin = margin(28, 0, 36, 4))

p_eas_c <- ggplot(eas_panel, aes(x = cohort, y = trait_label, fill = plot_estimate)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_tile(
    data = eas_panel[failed %in% TRUE],
    aes(x = cohort, y = trait_label),
    inherit.aes = FALSE,
    fill = "#D5D5D5",
    color = "white",
    linewidth = 0.7
  ) +
  geom_text(
    data = eas_panel[!is.na(cell_label)],
    aes(x = cohort, y = trait_label, label = cell_label),
    inherit.aes = FALSE,
    color = "#3A3A3A",
    size = 2.15,
    lineheight = 0.92,
    family = "Liberation Sans"
  ) +
  geom_text(
    data = eas_text,
    aes(label = label_expr, color = text_color),
    parse = TRUE,
    size = 2.8,
    fontface = "bold",
    family = "Liberation Sans"
  ) +
  scale_color_identity() +
  scale_fill_gradientn(
    colours = c("#B2182B", "#F4A582", "#F7F7F7", "#DEEBF7", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C"),
    values = rescale(c(-1, -0.5, 0, 0.60, 0.75, 0.85, 0.90, 0.95, 1)),
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    oob = squish,
    na.value = "white",
    name = expression(italic(r)[g])
  ) +
  scale_x_discrete(position = "top", expand = expansion(add = 0), labels = cohort_label) +
  scale_y_discrete(expand = expansion(add = 0), drop = FALSE) +
  labs(
    title = "EAS external cohorts",
    subtitle = "KoGES vs BBJ 2021, TWB, and TPMI"
  ) +
  panel_theme(TRUE) +
  guides(fill = legend_bar)

p_eur_c <- ggplot(eur_panel, aes(x = cohort, y = trait_label, fill = plot_estimate)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_tile(
    data = eur_panel[failed %in% TRUE],
    aes(x = cohort, y = trait_label),
    inherit.aes = FALSE,
    fill = "#D5D5D5",
    color = "white",
    linewidth = 0.7
  ) +
  geom_text(
    data = eur_panel[!is.na(cell_label)],
    aes(x = cohort, y = trait_label, label = cell_label),
    inherit.aes = FALSE,
    color = "#3A3A3A",
    size = 2.15,
    lineheight = 0.92,
    family = "Liberation Sans"
  ) +
  geom_text(
    data = eur_text,
    aes(label = label_expr, color = text_color),
    parse = TRUE,
    size = 2.8,
    fontface = "bold",
    family = "Liberation Sans"
  ) +
  scale_color_identity() +
  scale_fill_gradientn(
    colours = c("#FCFBFD", "#F4EAF8", "#E1C4EE", "#C49AD6", "#A45BBF", "#8E44AD", "#6C3483", "#4A235A"),
    values = rescale(c(0, 0.25, 0.45, 0.60, 0.70, 0.80, 0.90, 1)),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    oob = squish,
    na.value = "white",
    name = "pge"
  ) +
  scale_x_discrete(position = "top", expand = expansion(add = 0), labels = cohort_label) +
  scale_y_discrete(expand = expansion(add = 0), drop = FALSE) +
  labs(
    title = "EUR external cohorts",
    subtitle = "KoGES vs UKB, FinnGen, and MVP"
  ) +
  panel_theme(FALSE) +
  guides(fill = legend_bar)

p_combined <- p_group + p_eas_c + p_eur_c +
  plot_layout(widths = c(0.42, 4.6, 4.6))

combined_height <- 1.35 + 0.31 * n_trait
combined_files <- save_both(
  p_combined,
  "KoGES_EAS_EUR_cross_ancestry_heatmap",
  width = 16.4,
  height = combined_height
)
fig_root <- "/data4/joh/openkoges/figure_new0923"
fig_out <- file.path(fig_root, "output")
fig_panel <- file.path(fig_root, "panel")
fig_ed <- file.path(fig_root, "extended_data", "cross_ancestry")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)
dir.create(Sys.getenv("FIGURE5_OUT", unset = "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure5"), recursive = TRUE, showWarnings = FALSE)
dir.create(fig_panel, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_ed, recursive = TRUE, showWarnings = FALSE)
pdf_only <- function(paths) paths[grepl("\\.pdf$", paths)]
copied <- c(
  file.copy(combined_files, fig_out, overwrite = TRUE),
  file.copy(
    combined_files[grepl("\\.png$", combined_files)],
    file.path(Sys.getenv("FIGURE5_OUT", unset = "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure5"), "Figure5B.png"),
    overwrite = TRUE
  ),
  file.copy(pdf_only(c(eas_files, eur_files)), fig_panel, overwrite = TRUE),
  file.copy(c(eas_tsv, eur_tsv), fig_ed, overwrite = TRUE)
)
if (!all(copied)) stop("Could not copy the figure into output, panel, and extended_data")

cat("EAS pairs: ", nrow(eas), " rg_status success\n", sep = "")
print(eas[, .N, by = estimate_note])
cat("EUR pairs: ", nrow(eur), " fit_status success; finite pge: ", sum(is.finite(eur$estimate)),
    "; P<0.05 vs 1: ", sum(eur$differs_from_1, na.rm = TRUE), "\n", sep = "")
cat("Unmapped EAS labels: ", sum(!eas$trait_label %in% PAPER_LABEL), "\n", sep = "")
cat("Unmapped EUR labels: ", sum(!eur$trait_label %in% PAPER_LABEL), "\n", sep = "")
cat(paste(c(eas_files, eur_files, combined_files, fig_out, fig_panel, fig_ed, eas_tsv, eur_tsv), collapse = "\n"), "\n", sep = "")
