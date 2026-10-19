#!/usr/bin/env Rscript
# ============================================================
# SuppleFigure5C.r
#
# Genetic-correlation heatmap for binary traits
# (Supplementary Figure 5C). Parses within-KoGES LDSC rg logs
# and plots the pairwise rg matrix ordered by disease group.
# Traits with no finite rg (MI, KD, UB, UT) are dropped.
#
# Adapted from:
#   /data4/joh/openkoges/figure_new0923/code/figure3/05_extended_gc.R
#
# Inputs:
#   /data4/joh/openkoges/figure_new0923/intermediate/ldsc_rg_data5/*.log
#   (from ldsc_rg.R)
#
# Outputs:
#   DATA_manuscript_final/SuppleFigure5/genetic_correlation_binary.tsv
#   RESULT_manuscript_final/SuppleFigure5/SuppleFigure5C.png
#
# Usage:
#   Rscript SuppleFigure5C.r
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
DATA_DIR <- "/home2/sjim/openKoGES/DATA_manuscript_final/SuppleFigure5"
OUT_PNG  <- "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleFigure5/SuppleFigure5C.png"
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

BIN_PHENOS <- c(
  "HTN2","DM2","LIP","THY","ULCER","LIV","FLIV","STOMUL",
  "THYCA","NOI","MI","UTCA","CLD","GCA","COLCA","LCA","GT",
  "ANG","BRCA","CAD","OSTE","TB","COPD2","FRAC","GLAU","CVA",
  "BRON","STR","ASTH","CKD","KD","GB","UB","UT","ALLER","CATA",
  "GASTRO","ARTH","PER","POL","BPH"
)

# drop traits with no finite rg (grey rows): MI, KD, UB, UT
BIN_PHENOS_GC <- setdiff(BIN_PHENOS, c("MI", "KD", "UB", "UT"))

pheno_label_bin <- c(
  HTN2 = "Hypertension", DM2 = "Diabetes mellitus",
  LIP = "Hyperlipidemia", THY = "Thyroid disease", GT = "Gout",
  OSTE = "Osteoporosis", BPH = "Benign prostatic hyperplasia",
  MI = "Myocardial infarction", ANG = "Angina pectoris",
  CAD = "Coronary artery disease", CVA = "Cerebrovascular disease",
  STR = "Stroke",
  CLD = "Chronic lung disease", TB = "Pulmonary tuberculosis",
  COPD2 = "Chronic obstructive pulmonary disease",
  BRON = "Chronic bronchitis", ASTH = "Asthma",
  ULCER = "Peptic ulcer", LIV = "Acute liver disease",
  FLIV = "Fatty liver", STOMUL = "Chronic gastritis and peptic ulcer",
  GASTRO = "Chronic gastritis", POL = "Intestinal polyp",
  GB = "Gallbladder disease",
  CKD = "Chronic kidney disease", KD = "Kidney disease",
  UB = "Cystitis", UT = "Urinary tract infection",
  THYCA = "Thyroid cancer", UTCA = "Uterine cancer",
  GCA = "Gastric cancer", COLCA = "Colorectal cancer",
  LCA = "Lung cancer", BRCA = "Breast cancer",
  NOI = "Neurosis", FRAC = "Fracture", GLAU = "Glaucoma",
  ALLER = "Allergic disease", CATA = "Cataract",
  ARTH = "Arthritis", PER = "Periodontal disease"
)

get_group_bin <- function(ph) {
  ph0 <- sub("2$", "", ph)
  if (ph %in% c("POL","FLIV","GASTRO","ULCER","LIV","STOMUL","GB") ||
      ph0 %in% c("POL","FLIV","GASTRO","ULCER","LIV","STOMUL","GB"))
    return("Digestive")
  if (ph %in% c("HTN2","HTN","ANG","MI","CAD","CVA","STR") ||
      ph0 %in% c("HTN","ANG","MI","CAD","CVA","STR"))
    return("Cardiovascular")
  if (ph %in% c("GCA","LCA","BRCA","COLCA","UTCA","THYCA"))
    return("Cancer")
  if (ph %in% c("ASTH","TB","COPD2","COPD","BRON","CLD") ||
      ph0 %in% c("ASTH","TB","COPD","BRON","CLD"))
    return("Respiratory")
  if (ph %in% c("UT","UB","KD","CKD","BPH"))
    return("Renal_Urologic")
  if (ph %in% c("DM2","DM","LIP","THY","GT") ||
      ph0 %in% c("DM","LIP","THY","GT"))
    return("Endocrine_Metabolic")
  if (ph %in% c("ARTH","OSTE","FRAC"))
    return("Musculoskeletal")
  if (ph %in% c("CATA","GLAU"))
    return("Ophthalmologic")
  if (ph %in% c("ALLER","NOI","PER"))
    return("Other")
  "Other"
}

get_group_bin_gc <- function(ph) {
  gsub("_", " ", get_group_bin(ph))
}

group_colors_bin <- c(
  Digestive = "#009E73", Cardiovascular = "#F5C710",
  Cancer = "#D55E00", Respiratory = "#CC79A7",
  Renal_Urologic = "#56B4E9", Endocrine_Metabolic = "#F0E442",
  Musculoskeletal = "#E69F00", Ophthalmologic = "#0072B2",
  Other = "#999999"
)
group_colors_bin_gc <- group_colors_bin
names(group_colors_bin_gc) <- gsub("_", " ", names(group_colors_bin_gc))

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

RG_TSV <- "/home2/sjim/openKoGES/DATA_manuscript_final/Figure3/ldsc_rg_results.tsv"
RG_TSV_FALLBACK <- "/data4/joh/openkoges/figure_new0923/extended_data/figure3/ldsc_rg_results.tsv"

if (file.exists(RG_TSV)) {
  message("Loading pre-parsed rg table: ", RG_TSV)
  res <- fread(RG_TSV)
  n_parsed <- nrow(res)
} else if (file.exists(RG_TSV_FALLBACK)) {
  message("Loading pre-parsed rg table: ", RG_TSV_FALLBACK)
  res <- fread(RG_TSV_FALLBACK)
  n_parsed <- nrow(res)
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
}

res <- res[!is.na(rg) & !is.na(se) & is.finite(rg) & is.finite(se)]
n_bin <- sum(res$pheno1 %in% BIN_PHENOS & res$pheno2 %in% BIN_PHENOS)
message("rg pairs finite: ", nrow(res), " / parsed ", n_parsed,
        "  (bin-bin ", n_bin, " / ", choose(length(BIN_PHENOS), 2), ")")
message("binary GC traits: ", length(BIN_PHENOS_GC),
        " (dropped ", paste(setdiff(BIN_PHENOS, BIN_PHENOS_GC), collapse = ", "), ")")

###############################################################################
# Heatmap
###############################################################################

build_gc_heatmap <- function(res, traits_keep, get_group_fun, group_order,
                             label_map, group_colors, x_angle = 45) {
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
  message("binary GC: ", length(traits), " traits, ", n_pair, " off-diagonal pairs with rg")

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

gc_bin <- build_gc_heatmap(
  res = res,
  traits_keep = BIN_PHENOS_GC,
  get_group_fun = get_group_bin_gc,
  group_order = c(
    "Digestive","Cardiovascular","Cancer","Respiratory",
    "Renal Urologic","Endocrine Metabolic","Musculoskeletal",
    "Ophthalmologic","Other"
  ),
  label_map = pheno_label_bin,
  group_colors = group_colors_bin_gc,
  x_angle = 45
)
fwrite(gc_bin$data, file.path(DATA_DIR, "genetic_correlation_binary.tsv"), sep = "\t")

ggsave(
  filename = OUT_PNG,
  plot = gc_bin$plot,
  width = 36,
  height = 30,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote: ", normalizePath(OUT_PNG, winslash = "/"))
message("Wrote: ", file.path(DATA_DIR, "genetic_correlation_binary.tsv"))
