#!/usr/bin/env Rscript
# SuppleFigure4.R
#
# Purpose:
#   Plot KoGES samples projected onto 1000G Phase 3 ALL-population
#   reference PCs (Supplementary Figure 4: PC1–PC2 and PC1–PC3).
#   Points are coloured by continental super-population (AFR/AMR/EAS/
#   EUR/SAS) plus KoGES.
#
# Upstream analysis:
#   bash 1000G_projection_PCA.sh all
#   -> .../DATA_projectionPCA/1kg_reference_projection_all/
#        pca_1kg_reference.coordinates.tsv
#
# Inputs:
#   pca_1kg_reference.coordinates.tsv
#     (FID IID SOURCE POP SUPERPOP PC1..PC20;
#      SOURCE=1KG or KoGES)
# Outputs:
#   RESULT_manuscript_final/SuppleFigure4/SuppleFigure4_pc1_2.png
#   RESULT_manuscript_final/SuppleFigure4/SuppleFigure4_pc1_3.png
#
# Usage:
#   Rscript SuppleFigure4.R
#   Rscript SuppleFigure4.R all <coordinates.tsv>

suppressMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

mode <- "all"
if (length(args) >= 1 && args[[1]] %in% c("kchip", "all")) {
  mode <- args[[1]]
  args <- args[-1]
}

if (mode == "kchip") {
  default_coord <- "/data4/sjim/openKoGES_revision/DATA_projectionPCA/1kg_reference_projection_kchip/pca_1kg_reference.coordinates.tsv"
} else {
  default_coord <- "/data4/sjim/openKoGES_revision/DATA_projectionPCA/1kg_reference_projection_all/pca_1kg_reference.coordinates.tsv"
}

coord_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  default_coord
}

out_dir <- "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleFigure4"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_pc12 <- file.path(out_dir, "SuppleFigure4_pc1_2.png")
out_pc13 <- file.path(out_dir, "SuppleFigure4_pc1_3.png")

expected_pngs <- c(out_pc12, out_pc13)
if (all(file.exists(expected_pngs))) {
  message("[INFO] All figure PNGs already exist under ", out_dir, "; skip.")
  quit(save = "no", status = 0)
}

# ============================================================
# Style (Nature Comms-safe): colors & shapes
# ============================================================

KOGES_COL   <- "#6F359B"
KOGES_ALPHA <- 0.5
KOGES_SIZE  <- 1.6

# Continental super-population colours (Okabe–Ito) + KoGES
pal_super <- c(
  "AFR" = "#000000",
  "AMR" = "#E69F00",
  "EAS" = "#56B4E9",
  "EUR" = "#009E73",
  "SAS" = "#CC79A7",
  "KoGES" = "#6F359B"
)

shape_super <- c(
  "AFR" = 22,
  "AMR" = 23,
  "EAS" = 21,
  "EUR" = 24,
  "SAS" = 25,
  "KoGES" = 21
)

# ============================================================
# (1) Load projection coordinates
# ============================================================
if (!file.exists(coord_file)) {
  stop("Coordinates not found: ", coord_file)
}
message("[INFO] coordinates: ", coord_file)
message("[INFO] out_dir:     ", out_dir)

E <- fread(coord_file, header = TRUE)

need_cols <- c("IID", "SOURCE", "POP", "SUPERPOP", paste0("PC", 1:3))
missing_cols <- need_cols[!(need_cols %in% names(E))]
if (length(missing_cols) > 0L) {
  stop(
    "Missing columns in coordinates: ",
    paste(missing_cols, collapse = ", "),
    "\nAll columns: ", paste(names(E), collapse = ", ")
  )
}

if (!all(E$SOURCE %chin% c("1KG", "KoGES"))) {
  stop(
    "Unexpected SOURCE values: ",
    paste(unique(E$SOURCE), collapse = ", ")
  )
}

if (anyDuplicated(E$IID)) stop("Duplicate IID in coordinates.")

n_1kg   <- sum(E$SOURCE == "1KG")
n_koges <- sum(E$SOURCE == "KoGES")
if (n_1kg != 2504L) {
  stop("Expected 2,504 1KG samples, found: ", n_1kg)
}

# Remap labels for plotting (KoGES TSV rows keep SUPERPOP="EAS")
E[, super_pop := fifelse(SOURCE == "KoGES", "KoGES", SUPERPOP)]
E[, pop       := fifelse(SOURCE == "KoGES", "KoGES", POP)]

message("[INFO] samples: total=", nrow(E),
        " KoGES=", n_koges,
        " 1KG=", n_1kg)

# ============================================================
# (2) Plot function (ALL reference, continental colours)
# ============================================================
plot_pair <- function(d, x, y, file) {
  dd <- copy(d)
  dd[, grp := super_pop]
  dd[, is_koges := (grp == "KoGES")]

  x_range <- range(dd[[x]], na.rm = TRUE)
  y_range <- range(dd[[y]], na.rm = TRUE)
  x_limits <- c(x_range[1] - 0.003, x_range[2] + 0.003)
  y_limits <- c(y_range[1] - 0.003, y_range[2] + 0.003)

  p <- ggplot() +
    geom_point(
      data = dd[is_koges == TRUE],
      aes_string(x, y),
      shape = 16,
      color = KOGES_COL,
      size  = KOGES_SIZE,
      alpha = KOGES_ALPHA,
      show.legend = FALSE
    ) +
    geom_point(
      data = dd[is_koges == TRUE][1],
      aes_string(x, y, fill = "grp", shape = "grp"),
      size = 0, alpha = 0, stroke = 0
    ) +
    geom_point(
      data = dd[is_koges == FALSE],
      aes_string(x, y, fill = "grp", shape = "grp"),
      size  = 1.6,
      alpha = 0.90,
      color = "white",
      stroke = 0.2
    ) +
    scale_fill_manual(
      values = pal_super, drop = FALSE,
      breaks = names(pal_super),
      guide = guide_legend(
        override.aes = list(
          alpha = 1, size = 3, colour = "white", stroke = 0.35
        )
      )
    ) +
    scale_shape_manual(
      values = shape_super, drop = FALSE,
      breaks = names(shape_super)
    ) +
    labs(fill = NULL, shape = NULL) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    theme_classic(base_size = 16) +
    labs(x = x, y = y) +
    theme(
      aspect.ratio = 1,
      legend.position = "top",
      legend.direction = "horizontal",
      legend.key = element_blank(),
      legend.text = element_text(size = 11, family = "Arial"),
      axis.title = element_text(size = 18, family = "Arial"),
      axis.text  = element_text(size = 14, family = "Arial"),
      text = element_text(family = "Arial")
    )

  if (file.exists(file)) {
    message("[INFO] Skip (exists): ", file)
    return(invisible(NULL))
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, p, width = 7.2, height = 5.6, dpi = 300)
  message("[WROTE] ", file)
}

# ============================================================
# (3) Produce plots (PC1–PC2 / PC1–PC3)
# ============================================================
plot_pair(E, "PC1", "PC2", out_pc12)
plot_pair(E, "PC1", "PC3", out_pc13)

message("[DONE] Outputs: ", out_pc12, ", ", out_pc13)
