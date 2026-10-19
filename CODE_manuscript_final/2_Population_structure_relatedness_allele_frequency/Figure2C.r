#!/usr/bin/env Rscript
# Figure2C.r
#
# Purpose:
#   Genome-wide MAF spectrum of KoGES autosomal variants (Figure 2C):
#   - Main panel: MAF histogram for QC-retained autosomal markers
#   - Inset: MAC distribution among low-MAF markers (MAF < 0.05)
#     with dashed connectors from the low-MAF zoom box to the inset
#
# Upstream analysis:
#   bash Allele_frequency.sh
#
# Inputs:
#   DATA_manuscript_final/Figure2/KoGES_autosome_MAF.afreq
# Outputs:
#   RESULT_manuscript_final/Figure2/Figure2C.png

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(grid)
  library(cowplot)
})

pref <- "/home2/sjim/openKoGES/DATA_manuscript_final/Figure2/KoGES_autosome_MAF"
outp <- "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure2/Figure2C.png"

# -----------------------------------------------------------------------------
# 0) Load PLINK2 allele-frequency table
# -----------------------------------------------------------------------------
f_in <- paste0(pref, ".afreq")
if (!file.exists(f_in)) stop("Frequency file not found: ", f_in)

dt <- fread(f_in)
names_l <- tolower(names(dt))

get_col <- function(possible) {
  idx <- match(tolower(possible), names_l)
  idx <- idx[!is.na(idx)][1]
  if (!is.na(idx)) names(dt)[idx] else NA_character_
}

col_maf <- get_col(c("MAF", "maf"))
col_af  <- get_col(c("ALT_FREQS", "alt_freqs", "ALT_FREQ", "alt_freq", "ALT_FREQS1", "alt_freqs1"))
col_obs <- get_col(c("OBS_CT", "obs_ct", "NOBS", "nobs", "N", "n"))

# -----------------------------------------------------------------------------
# 1) Derive MAF and OBS_CT
# -----------------------------------------------------------------------------
if (!is.na(col_maf)) {
  dt[, maf := as.numeric(get(col_maf))]
} else if (!is.na(col_af)) {
  # PLINK2 may write comma-separated ALT frequencies; use the first ALT
  af_raw <- dt[[col_af]]
  af_num <- suppressWarnings(as.numeric(tstrsplit(as.character(af_raw), ",", fixed = TRUE)[[1]]))
  dt[, maf := pmin(af_num, 1 - af_num)]
} else {
  stop("No MAF/ALT_FREQS column found in ", f_in)
}

if (is.na(col_obs)) {
  stop(
    "No OBS_CT (allele observation count) column found in ", f_in,
    "\nOBS_CT is required to compute MAC for the inset."
  )
}
dt[, obs_ct := as.numeric(get(col_obs))]

dt <- dt[
  is.finite(maf) & !is.na(maf) & maf >= 0 & maf <= 0.5 &
    is.finite(obs_ct) & obs_ct > 0
]

# -----------------------------------------------------------------------------
# 2) Main MAF histogram (all QC-retained variants)
# -----------------------------------------------------------------------------
binw <- 0.01

main <- ggplot(dt, aes(x = maf)) +
  geom_histogram(
    binwidth = binw, boundary = 0, closed = "left",
    color = "white", linewidth = 0.30, fill = "grey30"
  ) +
  coord_cartesian(xlim = c(0, 0.5), expand = FALSE) +
  scale_x_continuous(
    breaks = seq(0, 0.5, by = 0.1),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_y_continuous(
    labels = function(x) label_comma(accuracy = 1)(x / 1000),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = "MAF",
    y = expression(paste("Number of markers (", "\u00D7 ", 10^3, ")"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    axis.text.x = element_text(margin = margin(t = 2, unit = "pt")),
    axis.text.y = element_text(margin = margin(r = 5, unit = "pt")),
    axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
    plot.margin = margin(6, 10, 12, 6, unit = "pt")
  )

# Zoom box over MAF 0–0.05 (height matches the tallest bin in that range)
x0 <- 0.00
x1 <- 0.05
y0 <- 0

h <- dt[, .N, by = .(bin = cut(maf, breaks = seq(0, 0.5, by = binw), right = FALSE))]
h[, bin_left := as.numeric(sub("^\\[([^,]+),.*$", "\\1", as.character(bin)))]
h_zoom <- h[bin_left < x1]
y1 <- max(h_zoom$N) * 1.0001

main_zoom <- main +
  annotate(
    "rect",
    xmin = x0, xmax = x1, ymin = y0, ymax = y1,
    fill = NA, color = "grey30", linewidth = 0.12,
    linetype = "dashed"
  )

# -----------------------------------------------------------------------------
# 3) Inset: MAC among low-MAF variants (MAF < 0.05)
# -----------------------------------------------------------------------------
dt_inset <- dt[maf < 0.05]
dt_inset[, mac := round(2 * obs_ct * maf)]  # diploid autosome approximation

mac_breaks <- c(-Inf, 0, 10, 100, 1000, 10000, Inf)
mac_labels <- c("0", "(1–10)", "(10–100)", "(100–1,000)", "(1,000–10,000)", "(>10,000)")

mac_dt <- dt_inset[, .N, by = .(
  mac_bin = cut(mac, breaks = mac_breaks, labels = mac_labels, right = TRUE)
)]
mac_dt[, mac_bin := factor(mac_bin, levels = mac_labels)]
mac_dt <- mac_dt[order(mac_bin)]

inset <- ggplot(mac_dt, aes(x = mac_bin, y = N)) +
  geom_col(width = 0.70, fill = "#6F359B", color = "#6F359B") +
  scale_y_continuous(
    labels = function(x) label_comma(accuracy = 1)(x / 1000),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Minor allele count",
    y = expression(paste("Number of markers (", "\u00D7 ", 10^3, ")"))
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title.x = element_text(size = 11, margin = margin(t = 7, unit = "pt")),
    axis.title.y = element_text(size = 11, margin = margin(r = 7, unit = "pt")),
    axis.text.x = element_text(
      size = 10, angle = 70, hjust = 1, vjust = 1,
      margin = margin(t = 4, unit = "pt")
    ),
    axis.text.y = element_text(size = 10, margin = margin(r = 3, unit = "pt")),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(4, 6, 10, 4, unit = "pt")
  )

# -----------------------------------------------------------------------------
# 4) Compose main + inset with dashed connectors
# -----------------------------------------------------------------------------
ix <- 0.38
iy <- 0.35
iw <- 0.60
ih <- 0.63

# Map data coordinates of the zoom box into the composed figure space
gb <- ggplot_build(main_zoom)
panel_info <- gb$layout$panel_params[[1]]
xr <- panel_info$x.range
yr <- panel_info$y.range

width_inch <- 7.2
height_inch <- 4.6
margin_npc_x_left <- 6 / (width_inch * 72)
margin_npc_x_right <- 6 / (width_inch * 72)
margin_npc_y_bottom <- 6 / (height_inch * 72)
margin_npc_y_top <- 12 / (height_inch * 72)

main_height <- 0.80
main_y <- 0.0

panel_x0 <- margin_npc_x_left
panel_x1 <- 1 - margin_npc_x_right
panel_y0 <- margin_npc_y_bottom
panel_y1 <- 1 - margin_npc_y_top

to_panel_x <- function(x) panel_x0 + (x - xr[1]) / (xr[2] - xr[1]) * (panel_x1 - panel_x0)
to_panel_y <- function(y) panel_y0 + (y - yr[1]) / (yr[2] - yr[1]) * (panel_y1 - panel_y0)

zx_right <- to_panel_x(x1)
zy_top <- to_panel_y(y1) * main_height + main_y

fig <- ggdraw() +
  draw_plot(main_zoom, x = 0, y = main_y, width = 1, height = main_height) +
  draw_plot(inset, x = ix, y = iy, width = iw, height = ih) +
  draw_grob(
    grid::rectGrob(
      x = ix + iw / 2,
      y = iy + ih / 2,
      width = iw,
      height = ih,
      gp = grid::gpar(col = "grey20", fill = NA, lwd = 0.5, lty = "longdash")
    )
  ) +
  # Connectors from zoom-box top-right toward inset left corners
  draw_line(
    x = c(zx_right + 0.13, ix),
    y = c(zy_top + 0.013, iy + ih),
    color = "grey30",
    size = 0.13,
    linetype = "longdash"
  ) +
  draw_line(
    x = c(zx_right + 0.13, ix),
    y = c(zy_top + 0.013, iy),
    color = "grey30",
    size = 0.13,
    linetype = "longdash"
  )

# -----------------------------------------------------------------------------
# 5) Save
# -----------------------------------------------------------------------------
ggsave(outp, fig, width = 5.5, height = 5.5, units = "in", dpi = 300, bg = "white")
message("[WROTE] ", outp)
