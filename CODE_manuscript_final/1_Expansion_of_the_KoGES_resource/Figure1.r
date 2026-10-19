# Figure 1: Geographic distribution and demographic composition of the KoGES cohort
#
# This script generates the main manuscript figures for the Korean Genome and
# Epidemiology Study (KoGES) overview:
#   - Figure 1A: Choropleth map of participant residential distribution by
#     province (SIDO), with recruitment site markers by cohort type
#     (Community-based, Rural, Urban).
#   - Figure 1B: Stacked bar charts of sex-by-cohort composition for the
#     total KoGES cohort and the genotyped subset.
#
# Data workflow (sample-level data are NOT for public release):
#   1) Read sample-level inputs from DATA_.../Figure1/raw/
#   2) Write publishable summary tables to DATA_.../Figure1/
#   3) Read those tables (plus public shapefiles/icons/recruit sites) and draw figures
#
# Inputs (raw, private):
#   DATA_manuscript_final/Figure1/raw/samples_recruit_region.rda
#   DATA_manuscript_final/Figure1/raw/df_juso_pc.rda
# Intermediate (publishable):
#   DATA_manuscript_final/Figure1/sido_sample_counts.csv
#   DATA_manuscript_final/Figure1/sex_cohort_freq.csv
# Public assets (kept in DATA_.../Figure1/):
#   bnd_sido_*.shp, bnd_sigungu_*.shp, recruit_sites.csv, icon_*.{svg,png}
# Outputs:
#   RESULT_manuscript_final/Figure1/Figure1A.png, Figure1B.png

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(cowplot)
  library(grid)
  library(png)
  library(rsvg)
})

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflict_prefer("filter", "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("select", "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("first", "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("margin", "ggplot2", quiet = TRUE)
}

data_dir   <- "/home2/sjim/openKoGES/DATA_manuscript_final/Figure1"
raw_dir    <- file.path(data_dir, "raw")
output_dir <- "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sido_count_path <- file.path(data_dir, "sido_sample_counts.csv")
sex_cohort_path <- file.path(data_dir, "sex_cohort_freq.csv")

shape_breaks <- c("Community-based", "Rural", "Urban")
cohort_outline <- c(
  "Community-based" = "#26958F",
  "Rural" = "#3B5B92",
  "Urban" = "#60BE70"
)
cohort_levels <- c("AS", "NC", "CT")
cohort_cols <- setNames(unname(cohort_outline[shape_breaks]), cohort_levels)

# -----------------------------------------------------------------------------
# 1) Build publishable summary tables from private sample-level raw data
# -----------------------------------------------------------------------------
if (file.exists(file.path(raw_dir, "samples_recruit_region.rda"))) {
  load(file.path(raw_dir, "samples_recruit_region.rda"))  # df
  df_sido_count <- df %>%
    mutate(SIDO_NM = na_if(trimws(as.character(SIDO_NM)), "")) %>%
    filter(!is.na(SIDO_NM)) %>%
    count(SIDO_NM, name = "n_total")
  write_csv(df_sido_count, sido_count_path)
  message("Wrote: ", sido_count_path)
}

if (file.exists(file.path(raw_dir, "df_juso_pc.rda"))) {
  load(file.path(raw_dir, "df_juso_pc.rda"))  # df_pc
  sex_cohort_freq <- bind_rows(
    df_pc %>% mutate(set = "All"),
    df_pc %>% filter(gene == 1) %>% mutate(set = "gene=1")
  ) %>%
    mutate(
      sex = case_when(
        sex == 1 ~ "Male",
        sex == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      cohort = as.character(cohort)
    ) %>%
    filter(!is.na(sex), !is.na(cohort), cohort %in% cohort_levels) %>%
    count(set, sex, cohort, name = "n")
  write_csv(sex_cohort_freq, sex_cohort_path)
  message("Wrote: ", sex_cohort_path)
}

if (!file.exists(sido_count_path) || !file.exists(sex_cohort_path)) {
  stop(
    "Summary tables not found. Place sample-level data under raw/ and re-run, ",
    "or provide sido_sample_counts.csv and sex_cohort_freq.csv."
  )
}

# -----------------------------------------------------------------------------
# Figure 1A: geographic map (from sido counts + public shapefiles/recruit sites)
# -----------------------------------------------------------------------------
df_sido_count <- read_csv(sido_count_path, show_col_types = FALSE)

sido <- st_read(file.path(data_dir, "bnd_sido_00_2024_2Q.shp"), quiet = TRUE) %>%
  select(SIDO_NM, geometry)
sgg <- st_read(file.path(data_dir, "bnd_sigungu_00_2024_2Q.shp"), quiet = TRUE) %>%
  select(SIGUNGU_NM, geometry)

sido_map <- sido %>%
  left_join(df_sido_count, by = "SIDO_NM") %>%
  mutate(n_total = replace_na(n_total, 0))

sido_map_main <- sido_map %>%
  filter(!SIDO_NM %in% c("제주특별자치도", "울릉군"))

recruit <- read_csv(file.path(data_dir, "recruit_sites.csv"), show_col_types = FALSE) %>%
  mutate(
    row_id   = row_number(),
    cohort   = as.character(cohort),
    followup = as.integer(followup),
    SIDO_NM  = na_if(trimws(as.character(SIDO_NM)), "")
  )

recruit$cohort <- ifelse(
  recruit$cohort == "urban", "Urban",
  ifelse(recruit$cohort == "rural", "Rural", "Community-based")
)

sig_cols <- grep("^SIGUNGU_NM", names(recruit), value = TRUE)
recruit <- recruit %>%
  mutate(across(all_of(sig_cols), ~ na_if(trimws(as.character(.x)), "")))

rec_sig_long <- recruit %>%
  select(row_id, cohort, followup, SIDO_NM, all_of(sig_cols)) %>%
  pivot_longer(cols = all_of(sig_cols), names_to = "sig_col", values_to = "SIGUNGU_NM") %>%
  filter(!is.na(SIGUNGU_NM)) %>%
  distinct(row_id, cohort, followup, SIDO_NM, SIGUNGU_NM)

sgg_joined <- sgg %>%
  right_join(rec_sig_long, by = "SIGUNGU_NM") %>%
  st_as_sf()

pts_sgg_sf <- sgg_joined %>%
  group_by(row_id) %>%
  summarise(
    cohort   = dplyr::first(cohort),
    followup = dplyr::first(followup),
    label    = paste(unique(SIGUNGU_NM), collapse = " + "),
    geometry = st_point_on_surface(st_union(geometry)),
    .groups  = "drop"
  )

rows_with_sig <- unique(rec_sig_long$row_id)

rec_sido_only <- recruit %>%
  filter(!row_id %in% rows_with_sig) %>%
  filter(!is.na(SIDO_NM))

sido_joined <- sido %>%
  right_join(rec_sido_only, by = "SIDO_NM") %>%
  st_as_sf()

pts_sido_sf <- sido_joined %>%
  group_by(row_id) %>%
  summarise(
    cohort   = dplyr::first(cohort),
    followup = dplyr::first(followup),
    label    = dplyr::first(SIDO_NM),
    geometry = st_point_on_surface(st_union(geometry)),
    .groups  = "drop"
  )

pts_recruit_sf <- bind_rows(pts_sgg_sf, pts_sido_sf) %>%
  mutate(
    cohort   = factor(cohort, levels = shape_breaks),
    followup = factor(followup, levels = c(0, 1), labels = c("Baseline only", "Follow-up"))
  )

pts_recruit_sf_main <- pts_recruit_sf %>%
  filter(!(label %in% c("제주특별자치도", "울릉군"))) %>%
  st_transform(st_crs(sido_map_main)) %>%
  filter(!st_is_empty(geometry))

final_plot <- ggplot() +
  geom_sf(data = sido_map_main, aes(fill = n_total), color = "white", linewidth = 0.25) +
  scale_fill_gradient(
    low = "grey95", high = "grey50", trans = "log1p",
    breaks = c(1000, 10000), labels = c("1,000", "10,000"),
    name = "Number of samples"
  ) +
  geom_sf(
    data = pts_recruit_sf_main,
    aes(shape = cohort, colour = cohort),
    size = 2.6, stroke = 1.0, fill = "white"
  ) +
  scale_shape_manual(
    values = c("Community-based" = 21, "Rural" = 24, "Urban" = 22),
    name = " "
  ) +
  scale_colour_manual(values = cohort_outline, name = " ") +
  coord_sf(datum = NA) +
  theme_void(base_size = 8) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "top",
    legend.box = "vertical",
    legend.justification = "center",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.height = grid::unit(16, "pt"),
    legend.key.width = grid::unit(26, "pt"),
    legend.box.spacing = grid::unit(10, "pt"),
    legend.spacing.x = grid::unit(10, "pt"),
    legend.spacing.y = grid::unit(8, "pt"),
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      nrow = 1,
      title.position = "left",
      label.position = "right",
      override.aes = list(
        colour = unname(cohort_outline[shape_breaks]),
        fill = "white",
        size = 3.5,
        stroke = 1.4
      )
    ),
    colour = "none",
    fill = guide_colorbar(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      label.theme = element_text(size = 10),
      title.theme = element_text(size = 14)
    )
  )

ggsave(
  filename = file.path(output_dir, "Figure1A.png"),
  plot = final_plot,
  width = 7.5,
  height = 7.0,
  units = "in",
  dpi = 300
)

# -----------------------------------------------------------------------------
# Figure 1B: sex x cohort composition (from frequency table)
# -----------------------------------------------------------------------------
dat_1b <- read_csv(sex_cohort_path, show_col_types = FALSE) %>%
  mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    cohort = factor(cohort, levels = c("AS", "NC", "CT")),
    set = factor(set, levels = c("All", "gene=1")),
    label = format(n, big.mark = ",")
  )

icon_male   <- rsvg_nativeraster(file.path(data_dir, "icon_male.svg"),   width = 200)
icon_female <- rsvg_nativeraster(file.path(data_dir, "icon_female.svg"), width = 200)
icon_dna    <- readPNG(file.path(data_dir, "icon_dna.png"))

n_all  <- dat_1b %>% filter(set == "All") %>% summarise(n = sum(n)) %>% pull(n)
n_gene <- dat_1b %>% filter(set == "gene=1") %>% summarise(n = sum(n)) %>% pull(n)
y_max  <- max(dat_1b %>% group_by(set, sex) %>% summarise(n = sum(n), .groups = "drop") %>% pull(n))

make_fig1b_panel <- function(dat, title, with_dna = FALSE) {
  title_y   <- y_max * 1.28
  bracket_y <- y_max * 1.18

  title_df <- data.frame(x = 1.5, y = title_y, label = title)
  bracket_df <- data.frame(
    x    = c(1.05, 1.05, 1.95),
    xend = c(1.95, 1.05, 1.95),
    y    = c(bracket_y, bracket_y, bracket_y),
    yend = c(bracket_y, bracket_y - y_max * 0.03, bracket_y - y_max * 0.03)
  )

  tops <- dat %>%
    group_by(sex) %>%
    summarise(y = sum(n), .groups = "drop") %>%
    mutate(x = as.numeric(sex))

  p <- ggplot(dat, aes(x = sex, y = n, fill = cohort)) +
    geom_col(width = 0.62, color = "white", linewidth = 0.35) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 14 / .pt,
      color = "white",
      family = "Arial",
      fontface = "bold"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.45,
      color = "grey35",
      linetype = "dotted"
    ) +
    geom_text(
      data = title_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 14 / .pt,
      fontface = "bold",
      family = "Arial",
      color = "black"
    ) +
    scale_fill_manual(values = cohort_cols, drop = FALSE) +
    scale_x_discrete(expand = expansion(add = 0.55)) +
    scale_y_continuous(limits = c(0, y_max * 1.38), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 12) +
    theme(
      text = element_text(family = "Arial"),
      legend.position = "none",
      plot.margin = margin(8, 12, 4, 12)
    )

  for (i in seq_len(nrow(tops))) {
    sex_i <- as.character(tops$sex[i])
    x_i   <- tops$x[i]
    y_i   <- tops$y[i]
    icon  <- if (sex_i == "Male") icon_male else icon_female
    icon_ymin <- y_i + y_max * 0.012
    icon_ymax <- icon_ymin + y_max * 0.095

    if (with_dna) {
      p <- p +
        annotation_custom(
          rasterGrob(icon, interpolate = TRUE),
          xmin = x_i - 0.26, xmax = x_i - 0.02,
          ymin = icon_ymin, ymax = icon_ymax
        ) +
        annotation_custom(
          rasterGrob(icon_dna, interpolate = TRUE),
          xmin = x_i + 0.00, xmax = x_i + 0.16,
          ymin = icon_ymin, ymax = icon_ymax
        )
    } else {
      p <- p +
        annotation_custom(
          rasterGrob(icon, interpolate = TRUE),
          xmin = x_i - 0.14, xmax = x_i + 0.14,
          ymin = icon_ymin, ymax = icon_ymax
        )
    }
  }
  p
}

p_all <- make_fig1b_panel(
  dat_1b %>% filter(set == "All"),
  title = sprintf("Total KoGES cohort (N=%s)", format(n_all, big.mark = ",")),
  with_dna = FALSE
)
p_gene <- make_fig1b_panel(
  dat_1b %>% filter(set == "gene=1"),
  title = sprintf("Genotyped subset (N=%s)", format(n_gene, big.mark = ",")),
  with_dna = TRUE
)

legend_order <- c("Community-based", "Rural", "Urban")
legend_cols  <- setNames(unname(cohort_outline[legend_order]), legend_order)

drawDetails.cohortLegend <- function(x, recording = TRUE) {
  labels <- x$labels
  colors <- x$colors
  fontsize <- x$fontsize
  swatch_w_pt <- x$swatch_w_pt
  swatch_h_pt <- x$swatch_h_pt
  label_gap_pt <- x$label_gap_pt
  item_gap_pt <- x$item_gap_pt
  fontfamily <- x$fontfamily

  gp_txt <- gpar(fontsize = fontsize, fontfamily = fontfamily, col = "black")
  label_w <- vapply(labels, function(lab) {
    convertWidth(grobWidth(textGrob(lab, gp = gp_txt)), "pt", valueOnly = TRUE)
  }, numeric(1))

  n <- length(labels)
  item_w <- swatch_w_pt + label_gap_pt + label_w + c(rep(item_gap_pt, n - 1), 0)
  x0 <- -sum(item_w) / 2
  xs <- x0 + c(0, cumsum(item_w[-n]))

  for (i in seq_len(n)) {
    grid.rect(
      x = unit(0.5, "npc") + unit(xs[i] + swatch_w_pt / 2, "pt"),
      y = unit(0.5, "npc"),
      width = unit(swatch_w_pt, "pt"),
      height = unit(swatch_h_pt, "pt"),
      gp = gpar(fill = colors[i], col = NA)
    )
    grid.text(
      labels[i],
      x = unit(0.5, "npc") + unit(xs[i] + swatch_w_pt + label_gap_pt, "pt"),
      y = unit(0.5, "npc"),
      just = c("left", "center"),
      gp = gp_txt
    )
  }
}

legend_centered <- gTree(
  labels = legend_order,
  colors = unname(legend_cols[legend_order]),
  fontfamily = "Arial",
  fontsize = 14,
  swatch_w_pt = 26,
  swatch_h_pt = 16,
  label_gap_pt = 6,
  item_gap_pt = 10,
  cl = "cohortLegend"
)

figure1B <- plot_grid(
  legend_centered,
  plot_grid(p_all, p_gene, nrow = 1, align = "h", rel_widths = c(1, 1)),
  ncol = 1,
  rel_heights = c(0.09, 1)
)

ggsave(
  filename = file.path(output_dir, "Figure1B.png"),
  plot = figure1B,
  width = 7.5,
  height = 7.0,
  units = "in",
  dpi = 300,
  bg = "white"
)
