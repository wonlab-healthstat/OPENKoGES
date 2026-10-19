#!/usr/bin/env Rscript
# ============================================================
# SuppleFigure5A.r
#
# Binary-trait multi-phenotype locus plot (Supplementary Figure 5A).
# Same layout as Figure 3A for disease phenotypes: X = chromosome;
# Y = phenotype (disease group); points = lead SNPs shared across
# traits at P < 5e-6. Gene labels are from ANNOVAR annotation of
# group-shared coding/near-gene variants.
#
# Inputs:
#   DATA_manuscript_final/SuppleFigure5/lead_snps_binary.tsv
#   DATA_manuscript_final/SuppleFigure5/sumstat_lead_binary_all.tsv
#   ANNOVAR binary annotation table (ANNOT_FILE_BIN)
#
# Output:
#   RESULT_manuscript_final/SuppleFigure5/SuppleFigure5A.png
#
# Usage:
#   Rscript SuppleFigure5A.r
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

###############################################################################
# Paths
###############################################################################

DATA_DIR <- "/home2/sjim/openKoGES/DATA_manuscript_final/SuppleFigure5"
OUT_PNG <- "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleFigure5/SuppleFigure5A.png"
ANNOT_FILE_BIN <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/ANNOVAR/final_annotation_binary.tsv"

lead_file <- file.path(DATA_DIR, "lead_snps_binary.tsv")
sum_file  <- file.path(DATA_DIR, "sumstat_lead_binary_all.tsv")
if (!file.exists(lead_file)) stop("Lead file not found: ", lead_file, call. = FALSE)
if (!file.exists(sum_file)) stop("Sumstat file not found: ", sum_file, call. = FALSE)
message("Lead file: ", lead_file)
message("Sumstat file: ", sum_file)

###############################################################################
# Phenotype labels / groups / colors
###############################################################################

BIN_PHENOS <- c(
  "HTN2","DM2","LIP","THY","ULCER","LIV","FLIV","STOMUL",
  "THYCA","NOI","MI","UTCA","CLD","GCA","COLCA","LCA","GT",
  "ANG","BRCA","CAD","OSTE","TB","COPD2","FRAC","GLAU","CVA",
  "BRON","STR","ASTH","CKD","KD","GB","UB","UT","ALLER","CATA",
  "GASTRO","ARTH","PER","POL","BPH"
)

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

group_colors_bin <- c(
  Digestive = "#009E73", Cardiovascular = "#F5C710",
  Cancer = "#D55E00", Respiratory = "#CC79A7",
  Renal_Urologic = "#56B4E9", Endocrine_Metabolic = "#F0E442",
  Musculoskeletal = "#E69F00", Ophthalmologic = "#0072B2",
  Other = "#999999"
)

clean_bin_pheno <- function(x) {
  x <- as.character(x)
  x <- sub("_SAIGE_cumulative_inci\\.SAIGE\\.txt$", "", x)
  x <- sub("_SAIGE_cumulative_inci.*$", "", x)
  x <- sub("^B_", "", x)
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

lead_bin <- fread(lead_file)
sum_all  <- fread(sum_file)

plot_dt <- merge(lead_bin, sum_all, by = "SNP", allow.cartesian = TRUE)
plot_dt <- plot_dt[!is.na(CHR) & !is.na(POS) & !is.na(pheno)]
plot_dt[, POS := as.numeric(POS)]
plot_dt[, CHR := as.integer(CHR)]
if (nrow(plot_dt) == 0) stop("plot_dt is empty for binary traits.", call. = FALSE)

###############################################################################
# 2. Y axis
###############################################################################

group_levels <- c(
  "Digestive","Cardiovascular","Cancer","Respiratory",
  "Renal_Urologic","Endocrine_Metabolic","Musculoskeletal",
  "Ophthalmologic","Other"
)
pheno_list <- data.table(pheno = BIN_PHENOS)
pheno_list[, pheno_group := factor(
  vapply(pheno, get_group_bin, character(1)),
  levels = group_levels
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
pheno_list[, display_label := fifelse(is.na(pheno), "", pheno_label_bin[pheno])]

group_labels <- pheno_list[!is.na(pheno), .(g_mid = mean(y_coord)), by = pheno_group]
plot_dt[pheno_list, y_coord := i.y_coord, on = "pheno"]
plot_dt <- plot_dt[!is.na(y_coord)]

###############################################################################
# 3. X axis
###############################################################################

chr_info <- plot_dt[, .(chr_len = as.numeric(max(POS))), by = CHR][order(as.integer(CHR))]
chr_info[, chr_offset := cumsum(as.numeric(shift(chr_len, fill = 0)))]
chr_info[, chr_mid := chr_offset + chr_len / 2]
plot_dt[chr_info, chr_offset := i.chr_offset, on = "CHR"]
plot_dt[, x_global := POS + chr_offset]
max_x <- max(plot_dt$x_global)

chr_blocks <- chr_info[as.integer(CHR) %% 2 == 1]
if (nrow(chr_blocks) == 0) chr_blocks <- chr_info
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

if (!file.exists(ANNOT_FILE_BIN)) {
  stop("ANNOVAR file not found: ", ANNOT_FILE_BIN, call. = FALSE)
}
anno <- fread(ANNOT_FILE_BIN)
anno[, pheno := clean_bin_pheno(Phenotype)]
anno <- anno[pheno %in% BIN_PHENOS]
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
    pheno_group = vapply(pheno, get_group_bin, character(1))
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
if (nrow(gene_counts_filt) == 0) {
  gene_counts_filt <- gene_counts2[N_pheno >= 2]
}
top_genes <- gene_counts_filt[, .SD[N_pheno == max(N_pheno)], by = pheno_group]

gene_coords <- merge(
  ann_all[, .(CHR, POS, Gene, pheno_group)],
  unique(plot_dt[, .(CHR, POS, x_global)]),
  by = c("CHR", "POS"),
  allow.cartesian = TRUE
)
gene_annots <- merge(top_genes, gene_coords, by = c("pheno_group", "Gene"))
if (nrow(gene_annots) > 0) {
  gene_annots <- gene_annots[
    , .(CHR = unique(CHR)[1], x_pos = median(x_global)),
    by = .(pheno_group, Gene)
  ]
  gene_annots <- merge(gene_annots, group_labels, by = "pheno_group", all.x = TRUE)
  gene_annots[, label_text := Gene]
  setorder(gene_annots, CHR, x_pos)
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
  scale_fill_manual(values = group_colors_bin) +
  scale_color_manual(values = group_colors_bin) +
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
