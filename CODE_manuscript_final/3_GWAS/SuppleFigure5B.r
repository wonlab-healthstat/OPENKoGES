#!/usr/bin/env Rscript
# ============================================================
# SuppleFigure5B.r
#
# SNP-based heritability forest plot for binary traits
# (Supplementary Figure 5B). Negative lower CI bounds are clipped
# at zero for display.
#
# Inputs:
#   RESULT_manuscript_final/SuppleTable3/SuppleTable3B.csv
#   (build with SuppleTable3.r after heritability.r)
#
# Output:
#   RESULT_manuscript_final/SuppleFigure5/SuppleFigure5B.png
#
# Usage:
#   Rscript SuppleFigure5B.r
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})



# Resolve openKoGES repo root without hard-coding host absolute paths.
.openkoges_repo_root <- function() {
  env <- Sys.getenv("OPENKOGES_ROOT", unset = "")
  if (nzchar(env)) return(normalizePath(env, winslash = "/", mustWork = FALSE))
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[startsWith(args, "--file=")])
  start <- if (length(fa)) dirname(normalizePath(fa[[1]], winslash = "/", mustWork = FALSE)) else getwd()
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(12L)) {
    if (dir.exists(file.path(d, "CODE_manuscript_final")) &&
        dir.exists(file.path(d, "DATA_manuscript_final"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate openKoGES root; set OPENKOGES_ROOT.", call. = FALSE)
}
REPO_ROOT <- .openkoges_repo_root()

.require_env_dir <- function(name) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) {
    stop(sprintf("Set environment variable %s (directory path).", name), call. = FALSE)
  }
  normalizePath(val, winslash = "/", mustWork = FALSE)
}



OPENKOGES_AUX_DIR <- .require_env_dir("OPENKOGES_AUX_DIR")

###############################################################################
# Paths
###############################################################################

IN_CSV  <- file.path(REPO_ROOT, "RESULT_manuscript_final", "SuppleTable3", "SuppleTable3B.csv")
MANIFEST <- file.path(OPENKOGES_AUX_DIR, "table_new0923", "intermediate", "suppletable3", "ldsc_binary_260924", "manifest.tsv")
OUT_PNG <- file.path(REPO_ROOT, "RESULT_manuscript_final", "SuppleFigure5", "SuppleFigure5B.png")

if (!file.exists(IN_CSV)) {
  stop(
    "Heritability table not found: ", IN_CSV,
    "\nRun SuppleTable3.r first.",
    call. = FALSE
  )
}
if (!file.exists(MANIFEST)) {
  stop("Binary LDSC manifest not found: ", MANIFEST, call. = FALSE)
}

###############################################################################
# Manuscript phenotype set and display labels (same as SuppleFigure 5A)
###############################################################################

BIN_PHENOS <- c(
  "HTN2","DM2","LIP","THY","ULCER","LIV","FLIV","STOMUL",
  "THYCA","NOI","MI","UTCA","CLD","GCA","COLCA","LCA","GT",
  "ANG","BRCA","CAD","OSTE","TB","COPD2","FRAC","GLAU","CVA",
  "BRON","STR","ASTH","CKD","KD","GB","UB","UT","ALLER","CATA",
  "GASTRO","ARTH","PER","POL","BPH"
)

pheno_label <- c(
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
# Plot
###############################################################################

h2b <- fread(IN_CSV)
man <- fread(MANIFEST)
h2b <- merge(h2b, man[, .(Phenotype = phenotype, variable)], by = "Phenotype", all.x = TRUE)
h2_col <- "Observed-scale h2 (SE)"
h2b[, h2 := as.numeric(sub(" \\(.*", "", get(h2_col)))]
h2b[, se := as.numeric(sub(".*\\((.*)\\)", "\\1", get(h2_col)))]
h2b[, trait_short := sub("^B_", "", variable)]

koges <- h2b[trait_short %in% BIN_PHENOS]
koges[, trait_full := pheno_label[trait_short]]
koges <- koges[!is.na(h2) & !is.na(se) & !is.na(trait_full)]
setorder(koges, h2)
koges[, trait_full := factor(trait_full, levels = trait_full)]
koges[, row_id := as.numeric(trait_full)]
koges[, ci_lo := h2 - se]
koges[, ci_hi := h2 + se]
# Clip negative lower bounds for binary observed-scale display.
koges[, h2 := pmax(h2, 0)]
koges[, ci_lo := pmax(ci_lo, 0)]
koges[, ci_hi := pmax(ci_hi, 0)]
bg_data <- koges[row_id %% 2 == 1]

x_max <- max(c(koges$h2, koges$ci_hi), na.rm = TRUE)

p <- ggplot(koges, aes(x = h2, y = trait_full)) +
  geom_rect(
    data = bg_data,
    aes(ymin = row_id - 0.5, ymax = row_id + 0.5, xmin = -Inf, xmax = Inf),
    inherit.aes = FALSE, fill = "grey97"
  ) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi),
    height = 0.25, color = "#0072B2", linewidth = 0.7
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_x_continuous(
    limits = c(0, x_max * 1.08),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = expression("SNP-based heritability"), y = NULL) +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = OUT_PNG,
  plot = p,
  width = 10,
  height = 11,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote: ", normalizePath(OUT_PNG, winslash = "/"),
        "  n=", nrow(koges))
