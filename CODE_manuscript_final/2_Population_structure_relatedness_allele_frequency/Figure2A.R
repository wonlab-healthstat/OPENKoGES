#!/usr/bin/env Rscript
# Figure2A.R
#
# Purpose:
#   Plot KoGES samples projected onto 1000G Phase 3 EAS-only reference PCs
#   (Figure 2A: PC1–PC2 and PC1–PC3).
#
# Upstream analysis:
#   bash 1000G_projection_PCA.sh eas
#   -> .../DATA_projectionPCA/1kg_eas_reference_projection_all/
#        pca_1kg_eas_reference.coordinates.tsv
#
# Inputs:
#   pca_1kg_eas_reference.coordinates.tsv
#     (FID IID SOURCE POP SUPERPOP PC1..PC20;
#      SOURCE=1KG_EAS or KoGES)
# Outputs:
#   RESULT_manuscript_final/Figure2/Figure2A_pc1_2.png
#   RESULT_manuscript_final/Figure2/Figure2A_pc1_3.png
#
# Usage:
#   Rscript Figure2A.R
#   Rscript Figure2A.R all <coordinates.tsv>

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
  default_coord <- "/data4/sjim/openKoGES_revision/DATA_projectionPCA/1kg_eas_reference_projection_kchip/pca_1kg_eas_reference.coordinates.tsv"
} else {
  default_coord <- "/data4/sjim/openKoGES_revision/DATA_projectionPCA/1kg_eas_reference_projection_all/pca_1kg_eas_reference.coordinates.tsv"
}

coord_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  default_coord
}

out_dir <- "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_pc12 <- file.path(out_dir, "Figure2A_pc1_2.png")
out_pc13 <- file.path(out_dir, "Figure2A_pc1_3.png")

expected_pngs <- c(out_pc12, out_pc13)
if (all(file.exists(expected_pngs))) {
  message("[INFO] All figure PNGs already exist under ", out_dir, "; skip.")
  quit(save = "no", status = 0)
}

EXPECTED_1KG_EAS_N <- 504L

# ============================================================
# Style (Nature Comms-safe): colors & shapes
# ============================================================

KOGES_COL   <- "#6F359B"
KOGES_ALPHA <- 0.5
KOGES_SIZE  <- 2.1

pal_eas <- c(
  "CDX"   = "#D55E00",
  "CHB"   = "#009E73",
  "CHS"   = "#E69F00",
  "JPT"   = "#56B4E9",
  "KHV"   = "#0072B2",
  "KoGES" = "#6F359B"
)

shape_eas <- c(
  "CDX"   = 22,
  "CHB"   = 23,
  "CHS"   = 24,
  "JPT"   = 25,
  "KHV"   = 22,
  "KoGES" = 21
)

eas_pops <- c("CHB", "CHS", "JPT", "CDX", "KHV")

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

if (!all(E$SOURCE %chin% c("1KG_EAS", "KoGES"))) {
  stop(
    "Unexpected SOURCE values: ",
    paste(unique(E$SOURCE), collapse = ", ")
  )
}

if (anyDuplicated(E$IID)) stop("Duplicate IID in coordinates.")

n_1kg   <- sum(E$SOURCE == "1KG_EAS")
n_koges <- sum(E$SOURCE == "KoGES")
if (n_1kg != EXPECTED_1KG_EAS_N) {
  stop("Expected ", EXPECTED_1KG_EAS_N, " 1KG EAS samples, found: ", n_1kg)
}

# Remap labels for plotting (KoGES TSV rows keep SUPERPOP="EAS")
E[, is_ref := (SOURCE == "1KG_EAS")]
E[, pop := fifelse(SOURCE == "KoGES", "KoGES", POP)]

if (any(E[is_ref == TRUE, SUPERPOP] != "EAS")) {
  stop("Non-EAS SUPERPOP found among 1KG_EAS reference rows.")
}
if (any(!E[is_ref == TRUE, pop] %chin% eas_pops)) {
  stop(
    "Unexpected EAS population label among reference rows: ",
    paste(setdiff(unique(E[is_ref == TRUE, pop]), eas_pops), collapse = ", ")
  )
}

E[, grp := fifelse(SOURCE == "KoGES", "KoGES", pop)]

message("[INFO] samples: total=", nrow(E),
        " KoGES=", n_koges,
        " 1KG_EAS=", n_1kg)

# ============================================================
# (2) Plot function (EAS population colours)
# ============================================================
plot_pair <- function(d, x, y, file) {
  dd <- copy(d)
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
      size   = 2.1,
      alpha  = 0.95,
      color  = "white",
      stroke = 0.35
    ) +
    scale_fill_manual(
      values = pal_eas, drop = FALSE,
      breaks = names(pal_eas),
      guide = guide_legend(
        override.aes = list(
          alpha = 1, size = 3, colour = "white", stroke = 0.35
        )
      )
    ) +
    scale_shape_manual(
      values = shape_eas, drop = FALSE,
      breaks = names(shape_eas)
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
# (3) Produce plots
# ============================================================
plot_pair(E, "PC1", "PC2", out_pc12)
plot_pair(E, "PC1", "PC3", out_pc13)

message("[DONE] Outputs: ", out_pc12, ", ", out_pc13)
