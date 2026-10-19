#!/usr/bin/env Rscript
# ============================================================
# Figure3C.r
#
# Genetic-correlation heatmap for continuous traits (Figure 3C).
# Parses within-KoGES LDSC rg logs and plots the pairwise rg
# matrix ordered by clinical group.
#
# Adapted from:
#   /data4/joh/openkoges/figure_new0923/code/figure3/05_extended_gc.R
#
# Inputs:
#   /data4/joh/openkoges/figure_new0923/intermediate/ldsc_rg_data5/*.log
#   (from ldsc_rg.R)
#
# Outputs:
#   DATA_manuscript_final/Figure3/ldsc_rg_results.tsv
#   DATA_manuscript_final/Figure3/genetic_correlation_continuous.tsv
#   RESULT_manuscript_final/Figure3/Figure3C.png
#
# Usage:
#   Rscript Figure3C.r
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

###############################################################################
# Paths
###############################################################################

LDSC_LOG_DIR <- "/data4/joh/openkoges/figure_new0923/intermediate/ldsc_rg_data5"
DATA_DIR <- "/home2/sjim/openKoGES/DATA_manuscript_final/Figure3"
OUT_PNG  <- "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure3/Figure3C.png"
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Font
###############################################################################

setup_figure_font <- function() {
  if (requireNamespace("showtext", quietly = TRUE) &&
      requireNamespace("sysfonts", quietly = TRUE)) {
    arial_reg  <- "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"
    arial_bold <- "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"
    if (file.exists(arial_reg)) {
      sysfonts::font_add(
        "Arial",
        regular = arial_reg,
        bold = if (file.exists(arial_bold)) arial_bold else arial_reg
      )
      showtext::showtext_auto()
      showtext::showtext_opts(dpi = 300)
    }
  }
  invisible(TRUE)
}
setup_figure_font()

###############################################################################
# Phenotype lists / labels / groups / colors
###############################################################################

CONT_PHENOS <- c(
  "ALBUMIN","ALP","ALT","APOA","APOB","AST","BMI","BUN",
  "CPEPTIDE","CREATINE","CYSTATINC","DBP","eGFR",
  "FEV1best","FEV1pp","FFRbest","FOLATE","FVCbest","FVCpp",
  "GLU0","HB","HBA1C","HCT","HDL","HEIGHT","HIP","HOMOCYST",
  "INSUL","IRON","LDL","MCH","MCHC","MCV","PHOSPHATE","PLAT",
  "PULSE","R_GTP","RBC_B","SBP","T_BIL","TCHL","TG","URICACID",
  "VITD","WAIST","WBC_B","WEIGHT","WHR"
)

BIN_PHENOS <- c(
  "HTN2","DM2","LIP","THY","ULCER","LIV","FLIV","STOMUL",
  "THYCA","NOI","MI","UTCA","CLD","GCA","COLCA","LCA","GT",
  "ANG","BRCA","CAD","OSTE","TB","COPD2","FRAC","GLAU","CVA",
  "BRON","STR","ASTH","CKD","KD","GB","UB","UT","ALLER","CATA",
  "GASTRO","ARTH","PER","POL","BPH"
)

pheno_label_cont <- c(
  RBC_B = "Red blood cell count", WBC_B = "White blood cell count",
  MCV = "MCV", MCH = "MCH", MCHC = "MCHC",
  PLAT = "Platelet count", HCT = "Hematocrit", HB = "Hemoglobin",
  FEV1best = "FEV1", FEV1pp = "FEV1 %pred", FFRbest = "FEV1/FVC",
  FVCbest = "FVC", FVCpp = "FVC %pred",
  HDL = "HDL cholesterol", LDL = "LDL cholesterol",
  TG = "Triglycerides", TCHL = "Total cholesterol",
  APOA = "Apolipoprotein A", APOB = "Apolipoprotein B",
  VITD = "Vitamin D", IRON = "Iron", FOLATE = "Folate",
  ALBUMIN = "Albumin", URICACID = "Uric acid", HOMOCYST = "Homocysteine",
  HIP = "Hip circumference", BMI = "Body mass index",
  WHR = "Waist to hip ratio", WAIST = "Waist circumference",
  WEIGHT = "Body weight", HEIGHT = "Body height",
  AST = "Aspartate transaminase", ALT = "Alanine aminotransferase",
  ALP = "Alkaline phosphatase", R_GTP = "Gamma glutamyl transferase",
  T_BIL = "Total bilirubin",
  BUN = "Blood urea nitrogen",
  eGFR = "Estimated glomerular filtration rate",
  CYSTATINC = "Cystatin C", PHOSPHATE = "Phosphate",
  CREATINE = "Serum creatinine",
  HBA1C = "HbA1c", INSUL = "Insulin", GLU0 = "Glucose", CPEPTIDE = "C peptide",
  SBP = "Systolic blood pressure", DBP = "Diastolic blood pressure",
  PULSE = "Pulse"
)

get_group_cont <- function(ph) {
  if (ph %in% c("BMI","HEIGHT","WEIGHT","WAIST","HIP","WHR")) return("Anthropometric")
  if (ph %in% c("VITD","PHOSPHATE")) return("Bone")
  if (ph %in% c("SBP","DBP","PULSE")) return("Cardiovascular")
  if (ph %in% c("HB","HCT","RBC_B","WBC_B","PLAT","MCV","MCH","MCHC","IRON","FOLATE"))
    return("Hematologic")
  if (ph %in% c("BUN","CREATINE","eGFR","CYSTATINC","URICACID")) return("Kidney")
  if (ph %in% c("AST","ALT","R_GTP","ALP","T_BIL","ALBUMIN")) return("Liver")
  if (ph %in% c("FEV1best","FEV1pp","FVCbest","FFRbest","FVCpp")) return("Lung")
  if (ph %in% c("GLU0","HBA1C","INSUL","CPEPTIDE","TCHL","TG","HDL","LDL","APOA","APOB","HOMOCYST"))
    return("Metabolic")
  "Other"
}

group_colors_cont <- c(
  Anthropometric = "#009E73", Bone = "#F5C710",
  Cardiovascular = "#D55E00", Hematologic = "#CC79A7",
  Kidney = "#56B4E9", Liver = "#F0E442",
  Lung = "#E69F00", Metabolic = "#0072B2"
)

###############################################################################
# Parse LDSC rg logs
###############################################################################

parse_trait <- function(p) {
  p <- basename(gsub("\\\\", "/", p))
  p <- sub("\\.sumstats\\.gz$", "", p)
  p <- sub("\\.sumstats$", "", p)
  p <- sub("_binary$", "", p)
  p <- sub("_noNC$", "", p)
  p <- sub("^RIN_", "", p)
  p <- sub("^B_", "", p)
  p <- sub("clean$", "", p)
  p
}

parse_one_log <- function(f) {
  txt <- readLines(f, warn = FALSE)
  hit <- grep("^Summary of Genetic Correlation Results", txt)
  if (length(hit) == 0) return(NULL)
  rest <- txt[(hit[1] + 1):length(txt)]
  rest <- rest[nzchar(trimws(rest))]
  if (length(rest) < 2) return(NULL)
  hdr <- rest[1]
  dat <- rest[2]
  if (!grepl("\\srg\\s", hdr) && !grepl("^\\s*p1\\s", hdr)) return(NULL)

  parts <- strsplit(trimws(dat), "\\s+")[[1]]
  if (length(parts) < 6) return(NULL)
  data.table(
    pheno1 = parse_trait(parts[1]),
    pheno2 = parse_trait(parts[2]),
    rg = as.numeric(parts[3]),
    se = as.numeric(parts[4]),
    z  = as.numeric(parts[5]),
    p  = as.numeric(parts[6])
  )
}

RG_TSV <- file.path(DATA_DIR, "ldsc_rg_results.tsv")
RG_TSV_FALLBACK <- "/data4/joh/openkoges/figure_new0923/extended_data/figure3/ldsc_rg_results.tsv"

if (file.exists(RG_TSV)) {
  message("Loading pre-parsed rg table: ", RG_TSV)
  res <- fread(RG_TSV)
  n_parsed <- nrow(res)
} else if (file.exists(RG_TSV_FALLBACK)) {
  message("Loading pre-parsed rg table: ", RG_TSV_FALLBACK)
  res <- fread(RG_TSV_FALLBACK)
  n_parsed <- nrow(res)
  fwrite(res, RG_TSV, sep = "\t")
} else {
  logs <- list.files(LDSC_LOG_DIR, pattern = "\\.log$", full.names = TRUE)
  if (length(logs) == 0) stop("No LDSC logs found in: ", LDSC_LOG_DIR, call. = FALSE)
  message("Parsing ", length(logs), " LDSC logs from ", LDSC_LOG_DIR, " ...")
  res_raw <- rbindlist(lapply(logs, parse_one_log), fill = TRUE)
  n_parsed <- nrow(res_raw)
  res <- res_raw[!is.na(rg) & !is.na(se) & is.finite(rg) & is.finite(se)]
  res[, pheno1 := sub("_binary$", "", pheno1)]
  res[, pheno2 := sub("_binary$", "", pheno2)]
  res[, pair := paste(pmin(pheno1, pheno2), pmax(pheno1, pheno2), sep = "|")]
  res <- res[!duplicated(pair)]
  res[, pair := NULL]
  fwrite(res, RG_TSV, sep = "\t")
}

res <- res[!is.na(rg) & !is.na(se) & is.finite(rg) & is.finite(se)]
n_cont <- sum(res$pheno1 %in% CONT_PHENOS & res$pheno2 %in% CONT_PHENOS)
message("rg pairs finite: ", nrow(res), " / parsed ", n_parsed,
        "  (cont-cont ", n_cont, " / ", choose(length(CONT_PHENOS), 2), ")")

###############################################################################
# Heatmap
###############################################################################

build_gc_heatmap <- function(res, traits_keep, get_group_fun, group_order,
                             label_map, group_colors, x_angle = 40) {
  res_use <- copy(res)
  res_use <- res_use[pheno1 %in% traits_keep & pheno2 %in% traits_keep]
  res_use[, z := rg / se]
  res_use[, p := 2 * pnorm(-abs(z))]

  res_sym <- rbind(
    res_use[, .(trait_x = pheno1, trait_y = pheno2, rg, p)],
    res_use[, .(trait_x = pheno2, trait_y = pheno1, rg, p)]
  )
  traits <- traits_keep
  n_pair <- nrow(unique(res_use[, .(pmin(pheno1, pheno2), pmax(pheno1, pheno2))]))
  message("continuous GC: ", length(traits), " traits, ", n_pair, " off-diagonal pairs with rg")

  diag_dt <- data.table(trait_x = traits, trait_y = traits, rg = 1, p = 0)
  final_df <- rbindlist(list(res_sym, diag_dt), fill = TRUE)
  final_df <- final_df[trait_x %in% traits & trait_y %in% traits]
  final_df <- final_df[!duplicated(final_df, by = c("trait_x", "trait_y"))]
  final_df[, label_text :=
             ifelse(is.na(p), "",
                    ifelse(p < 1e-4, "***",
                           ifelse(p < 1e-3, "**",
                                  ifelse(p < 0.05, "*", ""))))]

  trait_info <- data.table(
    trait = traits,
    group = vapply(traits, get_group_fun, character(1))
  )
  trait_info <- trait_info[!is.na(group)]
  trait_info[, group := factor(group, levels = group_order)]
  trait_info <- trait_info[order(group, trait)]
  trait_order <- trait_info$trait

  grid <- CJ(trait_x = trait_order, trait_y = trait_order, unique = TRUE)
  final_df <- merge(grid, final_df, by = c("trait_x", "trait_y"), all.x = TRUE)
  final_df[is.na(label_text), label_text := ""]
  final_df[trait_x == trait_y & is.na(rg), `:=`(rg = 1, p = 0, label_text = "")]
  final_df[, trait_x := factor(trait_x, levels = trait_order)]
  final_df[, trait_y := factor(trait_y, levels = rev(trait_order))]

  trait_info[, idx := seq_len(.N)]
  n_trait <- nrow(trait_info)
  group_label_pos <- trait_info[, .(ymin = min(idx), ymax = max(idx)), by = group]
  group_label_pos[, y_center := n_trait - (ymin + ymax) / 2 + 1]
  anno_left <- trait_info[, .(
    trait = factor(trait, levels = rev(trait_order)),
    group = group
  )]
  group_colors_alpha <- scales::alpha(group_colors, 0.5)

  p_left <- ggplot(anno_left, aes(x = 1, y = trait, fill = group)) +
    geom_tile() +
    geom_text(
      data = group_label_pos, inherit.aes = FALSE,
      aes(x = 1, y = y_center, label = group),
      angle = 90, fontface = "bold", family = "Arial", size = 6
    ) +
    scale_fill_manual(values = group_colors_alpha) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, 1.5)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = ggplot2::margin(0, -50, 0, 0, "pt"))

  p_heat <- ggplot(final_df, aes(x = trait_x, y = trait_y, fill = rg)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label_text), size = 8, fontface = "bold", family = "Arial") +
    scale_fill_gradient2(
      low = "#C75D5D", mid = "white", high = "#5A8FD8",
      midpoint = 0, limits = c(-1, 1), oob = squish,
      na.value = "grey90",
      name = expression(italic(r[g]))
    ) +
    scale_x_discrete(labels = label_map, position = "top", expand = c(0, 0)) +
    scale_y_discrete(labels = label_map, expand = c(0, 0)) +
    coord_fixed() +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 26, margin = ggplot2::margin(r = 2)),
      axis.text.x = element_text(size = 26, angle = x_angle, hjust = 0, vjust = 0),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(-20, 0, 0, -1, "pt"),
      legend.title = element_text(size = 24, face = "bold"),
      legend.text  = element_text(size = 22),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width  = unit(3, "cm")
    )

  final_plot <- p_left + p_heat + plot_layout(widths = c(0.5, 12))
  list(plot = final_plot, data = final_df)
}

gc_cont <- build_gc_heatmap(
  res = res,
  traits_keep = CONT_PHENOS,
  get_group_fun = get_group_cont,
  group_order = c("Anthropometric","Bone","Cardiovascular",
                  "Hematologic","Kidney","Liver","Lung","Metabolic"),
  label_map = pheno_label_cont,
  group_colors = group_colors_cont,
  x_angle = 40
)
fwrite(gc_cont$data, file.path(DATA_DIR, "genetic_correlation_continuous.tsv"), sep = "\t")

ggsave(
  filename = OUT_PNG,
  plot = gc_cont$plot,
  width = 36,
  height = 30,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote: ", normalizePath(OUT_PNG, winslash = "/"))
message("Wrote: ", file.path(DATA_DIR, "ldsc_rg_results.tsv"))
message("Wrote: ", file.path(DATA_DIR, "genetic_correlation_continuous.tsv"))
