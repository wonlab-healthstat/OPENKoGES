#!/usr/bin/env Rscript
# ============================================================
# Figure3B.r
#
# SNP-based heritability forest plot for continuous traits
# (Figure 3B).
#
# Inputs:
#   RESULT_manuscript_final/SuppleTable3/SuppleTable3A.csv
#   (build with SuppleTable3.r; continuous LDSC manifest for variable IDs)
#
# Output:
#   RESULT_manuscript_final/Figure3/Figure3B.png
#
# Usage:
#   Rscript Figure3B.r
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

###############################################################################
# Paths
###############################################################################

IN_CSV  <- "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleTable3/SuppleTable3A.csv"
MANIFEST <- "/data4/joh/openkoges/table_new0923/intermediate/suppletable3/ldsc_260923/manifest.tsv"
OUT_PNG <- "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure3/Figure3B.png"

if (!file.exists(IN_CSV)) {
  stop(
    "Heritability table not found: ", IN_CSV,
    "\nRun SuppleTable3.r first.",
    call. = FALSE
  )
}
if (!file.exists(MANIFEST)) {
  stop("Continuous LDSC manifest not found: ", MANIFEST, call. = FALSE)
}

###############################################################################
# Manuscript phenotype set and display labels (same as Figure 3A)
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

pheno_label <- c(
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

h2c <- fread(IN_CSV)
man <- fread(MANIFEST)
h2c <- merge(h2c, man[, .(Phenotype = phenotype, variable)], by = "Phenotype", all.x = TRUE)
h2_col <- "Observed-scale h2 (SE)"
h2c[, h2 := as.numeric(sub(" \\(.*", "", get(h2_col)))]
h2c[, se := as.numeric(sub(".*\\((.*)\\)", "\\1", get(h2_col)))]
h2c[, trait_short := sub("^B_", "", variable)]
h2c[trait_short == "LDLclean", trait_short := "LDL"]

koges <- h2c[trait_short %in% CONT_PHENOS]
koges[, trait_full := pheno_label[trait_short]]
koges <- koges[!is.na(h2) & !is.na(se) & !is.na(trait_full)]
setorder(koges, h2)
koges[, trait_full := factor(trait_full, levels = trait_full)]
koges[, row_id := as.numeric(trait_full)]
koges[, ci_lo := h2 - se]
koges[, ci_hi := h2 + se]
bg_data <- koges[row_id %% 2 == 1]

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
