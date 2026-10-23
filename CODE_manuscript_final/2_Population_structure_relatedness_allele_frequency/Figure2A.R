#!/usr/bin/env Rscript
# Figure2A.R
#
# Purpose:
#   Plot KoGES samples projected onto 1000 Genomes Phase 3 EAS-only
#   reference PCs with OADP correction (Figure 2A).
#   Panels: PC1-PC2 and PC1-PC3, combined into one figure with a shared
#   top-centre legend.
#
# KoGES point sizes (KoGES samples only; 1KG EAS sizes unchanged):
#   - KoGES (Non-outlier): 50% of the baseline size
#   - KoGES (Outlier)    : baseline size
#
# KoGES outlier rule (plot-specific; KoGES rows only):
#   One-pass Tukey fence with k = 5:
#     For each plotted axis among KoGES,
#       L = Q1 - 5*IQR ,  U = Q3 + 5*IQR
#     Flag as outlier if outside [L, U] on either axis of the panel.
#   Outlier status is recomputed separately for each panel.
#
# Upstream analysis:
#   bash 1000G_projection_OADP.sh eas
#   -> .../DATA_projectionPCA/1kg_eas_reference_projection_all/
#        pca_1kg_eas_reference.oadp.coordinates.tsv
#
# Inputs:
#   pca_1kg_eas_reference.oadp.coordinates.tsv
#     (FID IID SOURCE POP SUPERPOP PC1..PC20;
#      SOURCE = 1KG_EAS or KoGES)
#
# Output:
#   RESULT_manuscript_final/Figure2/Figure2A.png
#
# Usage:
#   Rscript Figure2A.R
#   Rscript Figure2A.R all <coordinates.tsv>

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)

mode <- "all"
if (length(args) >= 1 && args[[1]] %in% c("kchip", "all")) {
  mode <- args[[1]]
  args <- args[-1]
}

if (mode == "kchip") {
  default_coord <- "/data4/sjim/openKoGES_revision/DATA_projectionPCA/1kg_eas_reference_projection_kchip/pca_1kg_eas_reference.oadp.coordinates.tsv"
} else {
  default_coord <- "/data4/sjim/openKoGES_revision/DATA_projectionPCA/1kg_eas_reference_projection_all/pca_1kg_eas_reference.oadp.coordinates.tsv"
}

coord_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  default_coord
}

out_dir <- "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "Figure2A.png")

EXPECTED_1KG_EAS_N <- 504L
IQR_K <- 5

# ============================================================
# Style (Nature Comms-safe): colors & shapes
# ============================================================

KOGES_COL       <- "#6F359B"
KOGES_ALPHA     <- 0.5
KOGES_SIZE      <- 2.1                 # Outlier / baseline size
KOGES_SIZE_CORE <- KOGES_SIZE * 0.5    # Non-outlier = 50%

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

tukey_outlier <- function(x, k = IQR_K) {
  qs <- as.numeric(quantile(x, c(0.25, 0.75), names = FALSE, na.rm = TRUE))
  iqr <- qs[2] - qs[1]
  (x < qs[1] - k * iqr) | (x > qs[2] + k * iqr)
}

# ============================================================
# (1) Load projection coordinates
# ============================================================
if (!file.exists(coord_file)) {
  stop("Coordinates not found: ", coord_file)
}
message("[INFO] coordinates: ", coord_file)
message("[INFO] out_file:    ", out_file)
message("[INFO] KoGES outlier rule: ", IQR_K, "xIQR Tukey fence on plotted axes")
message("[INFO] KoGES sizes: non-outlier=", KOGES_SIZE_CORE,
        " outlier=", KOGES_SIZE)

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
plot_pair <- function(d, x, y) {
  dd <- copy(d)
  dd[, is_koges := (grp == "KoGES")]

  koges_idx <- which(dd$is_koges)
  out_flag <- rep(FALSE, nrow(dd))
  if (length(koges_idx) > 0L) {
    out_x <- tukey_outlier(dd[[x]][koges_idx])
    out_y <- tukey_outlier(dd[[y]][koges_idx])
    out_flag[koges_idx] <- out_x | out_y
  }
  dd[, is_outlier := out_flag]

  n_out <- sum(dd$is_koges & dd$is_outlier)
  n_core <- sum(dd$is_koges & !dd$is_outlier)
  message("[INFO] ", x, "-", y, ": KoGES non-outlier=", n_core,
          " outlier=", n_out, " (", IQR_K, "xIQR)")

  x_range <- range(dd[[x]], na.rm = TRUE)
  y_range <- range(dd[[y]], na.rm = TRUE)
  x_limits <- c(x_range[1] - 0.003, x_range[2] + 0.003)
  y_limits <- c(y_range[1] - 0.003, y_range[2] + 0.003)

  koges_core <- dd[is_koges == TRUE & is_outlier == FALSE]
  koges_out  <- dd[is_koges == TRUE & is_outlier == TRUE]
  size_levels <- c("KoGES (Non-outlier)", "KoGES (Outlier)")
  koges_plot <- rbindlist(list(koges_core, koges_out), use.names = TRUE, fill = TRUE)
  koges_plot[, koges_lab := factor(
    fifelse(is_outlier, "KoGES (Outlier)", "KoGES (Non-outlier)"),
    levels = size_levels
  )]
  pop_breaks <- setdiff(names(pal_eas), "KoGES")
  size_vals <- setNames(c(KOGES_SIZE_CORE, KOGES_SIZE), size_levels)
  if (!all(size_levels %in% as.character(unique(koges_plot$koges_lab)))) {
    miss <- setdiff(size_levels, as.character(unique(koges_plot$koges_lab)))
    pad <- koges_plot[rep(1L, length(miss))]
    pad[, koges_lab := factor(miss, levels = size_levels)]
    koges_plot <- rbindlist(list(koges_plot, pad), use.names = TRUE)
  }

  ggplot() +
    geom_point(
      data = koges_plot,
      aes_string(x, y, colour = "koges_lab", size = "koges_lab"),
      shape = 16,
      alpha = KOGES_ALPHA
    ) +
    geom_point(
      data = dd[is_koges == FALSE],
      aes_string(x, y, fill = "grp", shape = "grp"),
      size   = 2.1,
      alpha  = 0.95,
      color  = "white",
      stroke = 0.35
    ) +
    # Nearly-identical colours keep points looking the same but prevent
    # ggplot2 from collapsing discrete colour keys that share one hex.
    scale_colour_manual(
      name = NULL,
      values = setNames(c("#6F359B", "#6F359C"), size_levels),
      breaks = size_levels,
      limits = size_levels,
      drop = FALSE,
      guide = guide_legend(
        order = 1,
        override.aes = list(
          shape = 16,
          size = c(2.0, 4.0),
          alpha = 1
        )
      )
    ) +
    scale_size_manual(
      values = size_vals,
      breaks = size_levels,
      limits = size_levels,
      drop = FALSE,
      guide = "none"
    ) +
    scale_fill_manual(
      name = NULL,
      values = pal_eas, drop = FALSE,
      breaks = pop_breaks,
      guide = guide_legend(
        order = 2,
        override.aes = list(
          alpha = 1,
          size = 3,
          colour = "white",
          stroke = 0.35,
          shape = unname(shape_eas[pop_breaks])
        )
      )
    ) +
    scale_shape_manual(
      name = NULL,
      values = shape_eas, drop = FALSE,
      breaks = pop_breaks,
      guide = "none"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    theme_classic(base_size = 16) +
    labs(x = x, y = y) +
    theme(
      aspect.ratio = 1,
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.key = element_blank(),
      legend.text = element_text(size = 11, family = "Arial"),
      axis.title = element_text(size = 18, family = "Arial"),
      axis.text  = element_text(size = 14, family = "Arial"),
      text = element_text(family = "Arial")
    )
}

# ============================================================
# (3) Combined figure (shared legend once at top-centre)
# ============================================================
p12 <- plot_pair(E, "PC1", "PC2")
p13 <- plot_pair(E, "PC1", "PC3")

# cowplot 1.0 + ggplot2 3.5: get_legend() may return zeroGrob; pick the
# non-empty guide-box component instead.
guide_boxes <- get_plot_component(
  p12 + theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center"
  ),
  "guide-box",
  return_all = TRUE
)
shared_legend <- guide_boxes[[
  which(!vapply(guide_boxes, inherits, logical(1), "zeroGrob"))[1L]
]]
panels <- plot_grid(
  p12 + theme(legend.position = "none"),
  p13 + theme(legend.position = "none"),
  ncol = 2, align = "hv", axis = "tblr"
)
combined <- plot_grid(
  shared_legend,
  panels,
  ncol = 1,
  rel_heights = c(0.12, 1)
)
ggsave(out_file, combined, width = 12.5, height = 6.6, dpi = 300)
message("[WROTE] ", out_file)
message("[DONE] ", out_file)
