# Supplementary Figure 1 (new layout): Baseline vs cumulative disease prevalence
#
# Adapted from CODE_PCA_kinship_maf/supplefigure1.R into the same data/code
# workflow used by CODE_manuscript_final/SuppleFigure1.r.
#
# This script compares baseline (*_base) vs cumulative follow-up (*_FUnew)
# case counts for major diseases, separately for:
#   - Total KoGES cohort (epidemiology)
#   - Genotyped subset
#
# Data workflow (sample-level data are NOT for public release):
#   1) Read sample-level inputs from DATA_.../SuppleFigure1/raw/
#   2) Write a publishable frequency table to DATA_.../SuppleFigure1/
#   3) Read that table and draw the figures
#
# Inputs (raw, private):
#   DATA_manuscript_final/SuppleFigure1/raw/disease_data.RData
#   DATA_manuscript_final/SuppleFigure1/raw/genome_ids.txt
# Intermediate (publishable):
#   DATA_manuscript_final/SuppleFigure1/disease_freq_baseline_cumulative.csv
# Outputs:
#   RESULT_manuscript_final/SuppleFigure1/SuppleFigure1_new_epidemiology.png
#   RESULT_manuscript_final/SuppleFigure1/SuppleFigure1_new_genome.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(grid)
})

data_dir   <- "/home2/sjim/openKoGES/DATA_manuscript_final/SuppleFigure1"
raw_dir    <- file.path(data_dir, "raw")
output_dir <- "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleFigure1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

disease_name_order <- c(
  "Hypertension", "Diabetes", "Hyperlipidemia", "Cerebrovascular Disease",
  "Thyroid Disease", "Allergic Disease", "Asthma", "Gastritis/Peptic Ulcer"
)
base_vars_order  <- c("HTN_base", "DM_base", "LIP_base", "CEVA_base", "THY_base", "ALLER_base", "ASTH_base", "ULCER_base")
funew_vars_order <- c("HTN_FUnew", "DM_FUnew", "LIP_FUnew", "CEVA_FUnew", "THY_FUnew", "ALLER_FUnew", "ASTH_FUnew", "ULCER_FUnew")

freq_path <- file.path(data_dir, "disease_freq_baseline_cumulative.csv")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
# Count control(0)/case(1) by disease for one timepoint and population
count_disease_freq <- function(dat, var_names, disease_names, timepoint, population) {
  rows <- list()
  for (i in seq_along(var_names)) {
    x <- dat[[var_names[i]]]
    rows[[length(rows) + 1]] <- data.frame(
      disease = disease_names[i],
      timepoint = timepoint,
      population = population,
      outcome = "control",
      n = sum(x == 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      disease = disease_names[i],
      timepoint = timepoint,
      population = population,
      outcome = "case",
      n = sum(x == 1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# Build case-only plot table for one population:
# rows 1:8 = baseline cases; rows 9:16 = cumulative cases
freq_to_case_plot_table <- function(freq, population) {
  disease_name <- disease_name_order[disease_name_order %in% unique(freq$disease)]
  sub <- freq[freq$population == population & freq$outcome == "case", ]
  sub$disease <- factor(sub$disease, levels = disease_name)
  sub <- sub[order(sub$disease), ]

  baseline <- sub[sub$timepoint == "baseline", ]
  cumulative <- sub[sub$timepoint == "cumulative", ]

  dat44 <- rbind(
    data.frame(
      Group = factor(as.character(baseline$disease), levels = disease_name),
      Values = as.numeric(baseline$n),
      Type = "baseline_1",
      stringsAsFactors = FALSE
    ),
    data.frame(
      Group = factor(as.character(cumulative$disease), levels = disease_name),
      Values = as.numeric(cumulative$n),
      Type = "cumulative_1",
      stringsAsFactors = FALSE
    )
  )
  dat44$Type <- factor(dat44$Type, levels = c("baseline_1", "cumulative_1"))
  list(disease_name = disease_name, dat44 = dat44)
}

# Draw one population panel (Baseline vs Cumulative case bars)
# Colors: #99C3CC = Baseline, #B699CC = Cumulative
create_disease_plot <- function(dat44, subtitle, output_filename) {
  options(scipen = 999)

  # Short bars get labels near x=0 instead of mid-bar
  short_bar_indices <- c(4, 7)  # Cerebrovascular Disease, Asthma
  baseline_values <- dat44$Values[1:8]
  cumulative_values <- dat44$Values[9:16]

  baseline_text_df <- data.frame(
    x = 1:8 + 0.15,
    y = ifelse(1:8 %in% short_bar_indices,
               max(dat44$Values) * 0.02,
               baseline_values / 2),
    label = scales::comma(baseline_values),
    hjust = ifelse(1:8 %in% short_bar_indices, 0, 0.5)
  )

  ggp42 <- ggplot() +
    geom_segment(
      aes(x = 0.5, xend = 8.5,
          y = seq(0, max(dat44$Values), 10000),
          yend = seq(0, max(dat44$Values), 10000)),
      color = "grey70"
    ) +
    geom_col(aes(x = 1:8 - 0.15, y = Values), dat44[9:16, ], width = 0.5, fill = "#B699CC") +
    geom_col(aes(x = 1:8 + 0.15, y = Values), dat44[1:8, ], width = 0.5, fill = "#99C3CC") +
    geom_text(aes(x = x, y = y, label = label, hjust = hjust), baseline_text_df, family = "Arial", size = 4) +
    geom_text(
      aes(x = 1:8 - 0.15, y = dat44$Values[9:16], label = scales::comma(dat44$Values[9:16])),
      hjust = -0.1, nudge_y = max(dat44$Values) * 0.01, family = "Arial", size = 4
    ) +
    geom_text(aes(x = 1:8, y = 0, label = dat44$Group[1:8]), hjust = 1.1, family = "Arial", size = 4) +
    coord_flip(ylim = c(0, 1.2) * max(dat44$Values), xlim = c(0.4, 8.6), expand = FALSE) +
    labs(title = "", y = NULL, x = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      text = element_text(size = 13, family = "Arial"),
      plot.title = element_text(size = 13, hjust = 0.5),
      axis.text.x = element_text(size = 10)
    )

  # Center disease-name strip (kept for layout parity with source script)
  ggp43 <- ggplot() +
    geom_text(aes(x = 0, y = 1:8, label = dat44$Group[1:8]), family = "Arial", size = 4) +
    coord_cartesian(ylim = c(0.4, 8.6), expand = FALSE) +
    labs(title = "Case", y = NULL, x = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      text = element_text(size = 13, family = "Arial"),
      plot.title = element_text(hjust = 0.5, color = "transparent"),
      axis.text.x = element_text(size = 10, color = "transparent"),
      plot.margin = unit(c(5.5, -5.5, 5.5, -5.5), "pt")
    )

  legd1 <- get_legend(
    ggplot() +
      geom_col(aes(x = rep(1:8, each = 2), y = Values, fill = Type), dat44, width = 0.5) +
      scale_fill_manual(values = c("#99C3CC", "#B699CC"), labels = c("Baseline", "Cumulative")) +
      guides(fill = guide_legend(NULL, nrow = 1, byrow = TRUE)) +
      theme(text = element_text(size = 13, family = "Arial"))
  )

  ggp4 <- plot_grid(ggp43, ggp42, nrow = 1, rel_widths = c(0.3, 1))
  final_plot <- plot_grid(
    NULL, ggp4, plot_grid(NULL, legd1, nrow = 1, rel_widths = c(1, 9)),
    nrow = 3, rel_heights = c(0.7, 11, 1)
  )

  final_plot <- ggdraw(final_plot) +
    draw_label("   Prevalence of Major Diseases", x = 0, y = 0.98, hjust = 0, vjust = 1,
               fontfamily = "Arial", size = 13, fontface = "bold") +
    draw_label(paste0("   ", subtitle), x = 0, y = 0.95, hjust = 0, vjust = 1,
               fontfamily = "Arial", size = 13, fontface = "plain")

  ggsave(
    filename = file.path(output_dir, output_filename),
    plot = final_plot,
    width = 9, height = 5.6, dpi = 300, bg = "white"
  )
  invisible(final_plot)
}

# -----------------------------------------------------------------------------
# 1) Build publishable frequency table from private sample-level raw data
# -----------------------------------------------------------------------------
if (file.exists(file.path(raw_dir, "disease_data.RData")) &&
    file.exists(file.path(raw_dir, "genome_ids.txt"))) {

  load(file.path(raw_dir, "disease_data.RData"))  # dat_disease (list by cohort)
  dat_genome <- read.table(file.path(raw_dir, "genome_ids.txt"),
                           stringsAsFactors = FALSE, header = TRUE)
  dat_disease <- do.call(rbind, dat_disease)

  base_vars <- base_vars_order[base_vars_order %in% colnames(dat_disease)]
  funew_vars <- funew_vars_order[funew_vars_order %in% colnames(dat_disease)]
  disease_name <- disease_name_order[base_vars_order %in% colnames(dat_disease)]

  # Keep only diseases that have both baseline and cumulative columns
  keep <- base_vars_order %in% colnames(dat_disease) & funew_vars_order %in% colnames(dat_disease)
  base_vars <- base_vars_order[keep]
  funew_vars <- funew_vars_order[keep]
  disease_name <- disease_name_order[keep]

  dat_all <- dat_disease[, c("ID", base_vars, funew_vars)]
  dat_g <- subset(dat_all, ID %in% dat_genome$V1)

  freq <- rbind(
    count_disease_freq(dat_all, base_vars,  disease_name, "baseline",   "epidemiology"),
    count_disease_freq(dat_all, funew_vars, disease_name, "cumulative", "epidemiology"),
    count_disease_freq(dat_g,   base_vars,  disease_name, "baseline",   "genotyped"),
    count_disease_freq(dat_g,   funew_vars, disease_name, "cumulative", "genotyped")
  )
  write.csv(freq, freq_path, row.names = FALSE)
  message("Wrote: ", freq_path)
}

if (!file.exists(freq_path)) {
  stop("Frequency table not found. Place raw sample-level data under raw/ and re-run, ",
       "or provide disease_freq_baseline_cumulative.csv.")
}

# -----------------------------------------------------------------------------
# 2) Draw figures from the frequency table
# -----------------------------------------------------------------------------
freq <- read.csv(freq_path, stringsAsFactors = FALSE)

plot_epi <- freq_to_case_plot_table(freq, "epidemiology")
create_disease_plot(
  plot_epi$dat44,
  "(Total KoGES cohort)",
  "SuppleFigure1A.png"
)

plot_g <- freq_to_case_plot_table(freq, "genotyped")
create_disease_plot(
  plot_g$dat44,
  "(Genotyped subset)",
  "SuppleFigure1B.png"
)

