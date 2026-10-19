#!/usr/bin/env Rscript
# ============================================================
# Figure3A.r
#
# Continuous-trait multi-phenotype locus plot (Figure 3A).
# X = chromosome; Y = phenotype (clinical group); points = lead
# SNPs shared across traits at P < 5e-6. Gene labels are from
# ANNOVAR annotation of group-shared coding/near-gene variants.
#
# Inputs:
#   DATA_manuscript_final/Figure3/lead_snps_continuous.tsv
#   DATA_manuscript_final/Figure3/sumstat_lead_continuous_all.tsv
#   ANNOVAR continuous annotation table (ANNOT_FILE_CONT)
#
# Output:
#   RESULT_manuscript_final/Figure3/Figure3A.png
#
# Usage:
#   Rscript Figure3A.r
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

###############################################################################
# Paths
###############################################################################

DATA_DIR <- "/home2/sjim/openKoGES/DATA_manuscript_final/Figure3"
OUT_PNG <- "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure3/Figure3A.png"
ANNOT_FILE_CONT <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/ANNOVAR/final_annotation_bolt.tsv"

lead_file <- file.path(DATA_DIR, "lead_snps_continuous.tsv")
sum_file  <- file.path(DATA_DIR, "sumstat_lead_continuous_all.tsv")
if (!file.exists(lead_file)) stop("Lead file not found: ", lead_file, call. = FALSE)
if (!file.exists(sum_file)) stop("Sumstat file not found: ", sum_file, call. = FALSE)
message("Lead file: ", lead_file)
message("Sumstat file: ", sum_file)

###############################################################################
# Phenotype labels / groups / colors
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

clean_cont_pheno <- function(x) {
  x <- as.character(x)
  x <- sub("_BOLT_result$", "", x)
  x <- sub("^RIN_", "", x)
  x <- sub("^B_", "", x)
  x <- sub("_noNC$", "", x)
  x <- sub("clean$", "", x)
  x
}

clean_gene <- function(gene) {
  if (is.na(gene) || gene == "." || gene == "") return(NA_character_)
  genes <- unlist(strsplit(gene, "[;,]"))
  genes <- trimws(genes)
  if (length(genes) == 0) return(NA_character_)
  genes[1]
}

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
# 1. Lead SNPs × P<5e-6 sumstats
###############################################################################

lead_cont <- fread(lead_file)
sum_all   <- fread(sum_file)

plot_dt <- merge(lead_cont, sum_all, by = "SNP", allow.cartesian = TRUE)
plot_dt <- plot_dt[!is.na(CHR) & !is.na(POS) & !is.na(pheno)]
plot_dt[, POS := as.numeric(POS)]
plot_dt[, CHR := as.integer(CHR)]
if (nrow(plot_dt) == 0) stop("plot_dt is empty for continuous traits.", call. = FALSE)

###############################################################################
# 2. Y axis (gap row between groups)
###############################################################################

pheno_list <- data.table(pheno = CONT_PHENOS)
pheno_list[, pheno_group := factor(
  vapply(pheno, get_group_cont, character(1)),
  levels = c("Anthropometric","Bone","Cardiovascular",
             "Hematologic","Kidney","Liver","Lung","Metabolic")
)]
setorder(pheno_list, pheno_group, pheno)

pheno_list <- rbindlist(
  lapply(split(pheno_list, pheno_list$pheno_group), function(dt)
    rbind(dt, data.table(pheno_group = unique(dt$pheno_group), pheno = NA_character_))
  ),
  fill = TRUE
)
pheno_list[, y_idx := .I]
max_y <- nrow(pheno_list)
pheno_list[, y_coord := max_y - y_idx + 1]
pheno_list[, display_label := fifelse(is.na(pheno), "", pheno_label_cont[pheno])]

group_labels <- pheno_list[!is.na(pheno), .(g_mid = mean(y_coord)), by = pheno_group]
plot_dt[pheno_list, y_coord := i.y_coord, on = "pheno"]
plot_dt <- plot_dt[!is.na(y_coord)]

###############################################################################
# 3. X axis (chromosome)
###############################################################################

chr_info <- plot_dt[, .(chr_len = as.numeric(max(POS))), by = CHR][order(as.integer(CHR))]
chr_info[, chr_offset := cumsum(as.numeric(shift(chr_len, fill = 0)))]
chr_info[, chr_mid := chr_offset + chr_len / 2]
plot_dt[chr_info, chr_offset := i.chr_offset, on = "CHR"]
plot_dt[, x_global := POS + chr_offset]
max_x <- max(plot_dt$x_global)

chr_blocks <- chr_info[as.integer(CHR) %% 2 == 1]
zebra_dt <- CJ(
  chr_idx = seq_len(nrow(chr_blocks)),
  pheno_idx = seq_len(nrow(pheno_list[!is.na(pheno)]))
)
zebra_dt <- cbind(
  zebra_dt,
  chr_blocks[zebra_dt$chr_idx],
  pheno_list[!is.na(pheno)][zebra_dt$pheno_idx, .(ymin = y_coord - 0.5, ymax = y_coord + 0.5)]
)

###############################################################################
# 4. Gene annotation (ANNOVAR)
###############################################################################

if (!file.exists(ANNOT_FILE_CONT)) {
  stop("ANNOVAR file not found: ", ANNOT_FILE_CONT, call. = FALSE)
}
anno <- fread(ANNOT_FILE_CONT)
anno[, pheno := clean_cont_pheno(Phenotype)]
anno <- anno[pheno %in% CONT_PHENOS]
anno <- anno[Func.refGene %in% c("exonic","UTR3","UTR5","upstream","downstream","upstream;downstream")]
if (nrow(anno) > 0) {
  if ("Gene_clean" %in% names(anno)) {
    anno[, Gene := fifelse(
      is.na(Gene_clean) | Gene_clean == "" | Gene_clean == ".",
      vapply(Gene.refGene, clean_gene, character(1)),
      as.character(Gene_clean)
    )]
  } else {
    anno[, Gene := vapply(Gene.refGene, clean_gene, character(1))]
  }
  ann_all <- anno[, .(
    CHR = as.integer(Chr),
    POS = as.numeric(Start),
    Gene,
    pheno,
    pheno_group = vapply(pheno, get_group_cont, character(1))
  )]
} else {
  ann_all <- data.table(
    CHR = integer(), POS = numeric(), Gene = character(),
    pheno = character(), pheno_group = character()
  )
}

gene_counts <- ann_all[!is.na(Gene), .(N_pheno = uniqueN(pheno)), by = .(pheno_group, Gene)]
pheno_group_n <- pheno_list[!is.na(pheno), .N, by = pheno_group]
gene_counts2 <- merge(gene_counts, pheno_group_n, by = "pheno_group")
gene_counts2[, prop := N_pheno / N]
gene_counts_filt <- gene_counts2[prop >= 0.6]
top_genes <- gene_counts_filt[, .SD[N_pheno == max(N_pheno)], by = pheno_group]

gene_coords <- merge(
  ann_all[, .(CHR, POS, Gene, pheno_group)],
  unique(plot_dt[, .(CHR, POS, x_global)]),
  by = c("CHR", "POS"),
  allow.cartesian = TRUE
)
gene_annots <- merge(top_genes, gene_coords, by = c("pheno_group", "Gene"))[
  , .(CHR = unique(CHR)[1], x_pos = median(x_global)),
  by = .(pheno_group, Gene)
]
gene_annots <- merge(gene_annots, group_labels, by = "pheno_group", all.x = TRUE)
gene_annots[, label_text := Gene]
setorder(gene_annots, CHR, x_pos)

if (nrow(gene_annots) > 0) {
  zigzag_levels <- c(+1, -1, +2, -2)
  gene_annots[, zigzag_idx := seq_len(.N), by = CHR]
  gene_annots[, zigzag_id := (zigzag_idx - 1) %% length(zigzag_levels) + 1]
  label_y_base <- -3
  label_y_step <- 1.5
  gene_annots[, label_y := label_y_base + zigzag_levels[zigzag_id] * label_y_step]
  y_tol <- label_y_step * 0.5
  x_tol <- max_x * 0.00015
  gene_annots[, label_y_adj := label_y]
  if (nrow(gene_annots) >= 2) {
    for (i in 2:nrow(gene_annots)) {
      same_chr <- gene_annots$CHR[i] == gene_annots$CHR[i - 1]
      close_x  <- abs(gene_annots$x_pos[i] - gene_annots$x_pos[i - 1]) < x_tol
      same_y   <- abs(gene_annots$label_y_adj[i] - gene_annots$label_y_adj[i - 1]) < y_tol
      if (same_chr && close_x && same_y) {
        gene_annots$label_y_adj[i] <-
          gene_annots$label_y_adj[i] + sign(gene_annots$label_y_adj[i]) * label_y_step * 0.8
      }
    }
  }
  final_gene_key <- unique(gene_annots[, .(pheno_group, Gene)])
  gene_locus_pts <- merge(ann_all, final_gene_key, by = c("pheno_group", "Gene"))
  gene_locus_pts <- merge(
    gene_locus_pts,
    plot_dt[, .(CHR, POS, pheno, x_global, y_coord)],
    by = c("CHR", "POS", "pheno"),
    allow.cartesian = TRUE
  )
} else {
  gene_annots <- data.table()
  gene_locus_pts <- data.table()
}

###############################################################################
# 5. Plot → PNG
###############################################################################

p <- ggplot() +
  geom_rect(
    data = zebra_dt,
    aes(xmin = chr_offset, xmax = chr_offset + chr_len, ymin = ymin, ymax = ymax),
    fill = "grey94"
  ) +
  geom_rect(
    data = pheno_list[!is.na(pheno)],
    aes(xmin = 0, xmax = max_x, ymin = y_coord - 0.5, ymax = y_coord + 0.5, fill = pheno_group),
    alpha = 0.2
  ) +
  geom_hline(yintercept = seq(0.5, max_y + 0.5, 1), color = "white", linewidth = 0.5) +
  geom_rect(
    data = chr_info,
    aes(xmin = chr_offset, xmax = chr_offset + chr_len, ymin = max_y + 1.2, ymax = max_y + 2.5),
    fill = "white", color = "black", linewidth = 0.4
  ) +
  geom_text(
    data = chr_info, aes(x = chr_mid, y = max_y + 1.85, label = CHR),
    size = 7, fontface = "bold", family = "Arial"
  ) +
  geom_point(
    data = plot_dt, aes(x = x_global, y = y_coord),
    shape = 21, size = 1.2, alpha = 0.7
  ) +
  geom_text(
    data = group_labels, aes(x = -max_x * 0.13, y = g_mid, label = pheno_group),
    hjust = 0, fontface = "bold", size = 7.5, family = "Arial"
  ) +
  geom_text(
    data = pheno_list, aes(x = -max_x * 0.01, y = y_coord, label = display_label),
    hjust = 1, size = 6, family = "Arial"
  )

if (nrow(gene_locus_pts) > 0) {
  p <- p +
    geom_point(
      data = gene_locus_pts,
      aes(x = x_global, y = y_coord, color = pheno_group),
      shape = 21, size = 1.2, fill = "white", stroke = 1.1
    )
}
if (nrow(gene_annots) > 0) {
  p <- p +
    geom_label(
      data = gene_annots,
      aes(x = x_pos, y = label_y_adj, label = label_text, fill = pheno_group),
      size = 6, fontface = "bold", family = "Arial", label.size = 0.35, alpha = 0.45
    ) +
    geom_segment(
      data = gene_annots,
      aes(x = x_pos, xend = x_pos, y = g_mid, yend = label_y_adj, color = pheno_group),
      linetype = "dotted", linewidth = 0.9, alpha = 0.8
    )
}

p <- p +
  scale_fill_manual(values = group_colors_cont) +
  scale_color_manual(values = group_colors_cont) +
  scale_x_continuous(limits = c(-max_x * 0.15, max_x), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-10, max_y + 4), expand = c(0, 0)) +
  theme_void() +
  theme(
    text = element_text(family = "Arial"),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = OUT_PNG,
  plot = p,
  width = 48,
  height = 15,
  units = "in",
  dpi = 300,
  bg = "white"
)

message(
  "Wrote: ", normalizePath(OUT_PNG, winslash = "/"),
  "  points=", nrow(plot_dt),
  "  genes=", if (nrow(gene_annots)) nrow(gene_annots) else 0
)
