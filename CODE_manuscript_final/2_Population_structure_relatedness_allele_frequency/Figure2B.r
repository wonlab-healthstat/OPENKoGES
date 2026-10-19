#!/usr/bin/env Rscript
# Figure2B.r
#
# Purpose:
#   Summarise KING relatedness in the KoGES cohort (Figure 2B):
#   donut of relationship categories and kinship-coefficient histogram.
#   Sample-level .kin0 is private; only frequency tables are publishable.
#
# Upstream analysis:
#   bash Kinship.sh
#
# Data workflow:
#   1) If private KING .kin0 is present, build publishable frequency tables
#   2) Read those frequency tables and draw figures / print summaries
#
# Inputs (raw, private):
#   DATA_manuscript_final/Figure2/raw/..._pr0.2_king.kin0
#   DATA_manuscript_final/Figure2/raw/...king.cutoff.{in,out}.id
# Intermediate (publishable):
#   DATA_manuscript_final/Figure2B/kinship_meta.csv
#   DATA_manuscript_final/Figure2B/kinship_donut_freq.csv
#   DATA_manuscript_final/Figure2B/kinship_category_freq.csv
#   DATA_manuscript_final/Figure2B/kinship_hist_freq.csv
# Outputs:
#   RESULT_manuscript_final/Figure2/Figure2B_kinship_donut.png
#   RESULT_manuscript_final/Figure2/Figure2B_kinship_hist.png
#   (frequency summary is printed to the console / log)

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

raw_dir   <- "/home2/sjim/openKoGES/DATA_manuscript_final/Figure2/raw"
kin0_path <- file.path(raw_dir, "KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8_pr0.2_king.kin0")
data_dir  <- "/home2/sjim/openKoGES/DATA_manuscript_final/Figure2B"
out_dir   <- "/home2/sjim/openKoGES/RESULT_manuscript_final/Figure2"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta_path  <- file.path(data_dir, "kinship_meta.csv")
donut_path <- file.path(data_dir, "kinship_donut_freq.csv")
cat_path   <- file.path(data_dir, "kinship_category_freq.csv")
hist_path  <- file.path(data_dir, "kinship_hist_freq.csv")

# Analysis thresholds (method parameters; not cohort counts)
cut_2nd <- 0.0884
cut_1st <- 0.177
cut_dup <- 0.354
hist_bins <- 100L

cutoff_1st_out <- file.path(
  raw_dir,
  "KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8_pr0.2_king_1stdeg_2.king.cutoff.out.id"
)
cutoff_1st_in <- file.path(
  raw_dir,
  "KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8_pr0.2_king_1stdeg_2.king.cutoff.in.id"
)
cutoff_2nd_out <- file.path(
  raw_dir,
  "KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8_pr0.2_king_2nddeg_2.king.cutoff.out.id"
)
cutoff_2nd_in <- file.path(
  raw_dir,
  "KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8_pr0.2_king_2nddeg_2.king.cutoff.in.id"
)

# PLINK king.cutoff.*.id has a header row
count_id_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required ID file not found: ", path)
  }
  as.integer(max(0L, length(readLines(path, warn = FALSE)) - 1L))
}

# -----------------------------------------------------------------------------
# 1) Build publishable frequency tables from private raw outputs (if present)
# -----------------------------------------------------------------------------
if (file.exists(kin0_path)) {
  n_removed_1st <- count_id_file(cutoff_1st_out)
  n_kept_1st <- count_id_file(cutoff_1st_in)
  n_removed_2nd <- count_id_file(cutoff_2nd_out)
  n_kept_2nd <- count_id_file(cutoff_2nd_in)
  n_total <- n_removed_1st + n_kept_1st
  if ((n_removed_2nd + n_kept_2nd) != n_total) {
    stop(
      "Cohort size mismatch between 1st- and 2nd-degree KING cutoff files: ",
      n_total, " vs ", n_removed_2nd + n_kept_2nd
    )
  }

  pairs <- fread(kin0_path)
  id1_col <- intersect(c("ID1", "IID1", "i1", "IID", "#ID1"), names(pairs))[1]
  id2_col <- intersect(c("ID2", "IID2", "i2"), names(pairs))[1]
  kin_col <- intersect(c("KINSHIP", "Kinship", "kinship"), names(pairs))[1]
  if (is.na(id1_col) || is.na(id2_col) || is.na(kin_col)) {
    stop("Could not detect ID1/ID2/KINSHIP columns in kin0 file.")
  }
  setnames(pairs, c(id1_col, id2_col, kin_col), c("IID1", "IID2", "KINSHIP"))
  pairs[, `:=`(
    IID1 = as.character(IID1),
    IID2 = as.character(IID2),
    KINSHIP = as.numeric(KINSHIP)
  )]
  pairs <- pairs[!is.na(IID1) & !is.na(IID2) & !is.na(KINSHIP)]
  pairs <- pairs[KINSHIP >= cut_2nd]
  if (nrow(pairs) == 0L) stop("No pairs with KINSHIP >= cut_2nd.")

  n_unique_related <- uniqueN(c(pairs$IID1, pairs$IID2))
  n_not_in_related <- n_total - n_unique_related
  n_unique_1st_or_closer <- uniqueN(c(
    pairs[KINSHIP >= cut_1st, IID1],
    pairs[KINSHIP >= cut_1st, IID2]
  ))

  pairs[, rel_cat := fifelse(
    KINSHIP > cut_dup, "Duplicate/MZ twins",
    fifelse(KINSHIP > cut_1st, "1st-degree", "2nd-degree")
  )]
  rel_levels <- c("2nd-degree", "1st-degree", "Duplicate/MZ twins")

  cat_summary <- pairs[, .(n_pairs = .N), by = rel_cat]
  cat_unique <- rbindlist(lapply(rel_levels, function(cc) {
    sub <- pairs[rel_cat == cc]
    data.table(
      rel_cat = cc,
      n_unique_individuals = uniqueN(c(sub$IID1, sub$IID2))
    )
  }))
  cat_summary <- merge(cat_summary, cat_unique, by = "rel_cat", all = TRUE)
  cat_summary[is.na(n_pairs), n_pairs := 0L]
  cat_summary[is.na(n_unique_individuals), n_unique_individuals := 0L]
  cat_summary[, rel_cat := factor(rel_cat, levels = rel_levels)]
  setorder(cat_summary, rel_cat)
  cat_summary[, rel_cat := as.character(rel_cat)]

  donut_tbl <- data.table(
    group = c("not_in_related_pairs", "in_related_pairs"),
    label = c(
      "Individuals not in related pairs",
      sprintf("Individuals in related pairs (kinship \u2265 %.4f)", cut_2nd)
    ),
    n = c(n_not_in_related, n_unique_related)
  )

  xmax <- min(max(as.numeric(quantile(pairs$KINSHIP, 0.9995, na.rm = TRUE)), 0.28), 0.50)
  breaks <- seq(cut_2nd, xmax, length.out = hist_bins + 1L)
  pairs_hist <- pairs[KINSHIP >= cut_2nd & KINSHIP <= xmax]
  h <- hist(pairs_hist$KINSHIP, breaks = breaks, plot = FALSE, right = TRUE)
  hist_freq <- data.table(
    bin_left = head(h$breaks, -1L),
    bin_right = tail(h$breaks, -1L),
    bin_mid = h$mids,
    n = as.integer(h$counts)
  )

  meta <- data.table(
    param = c(
      "cut_2nd", "cut_1st", "cut_dup",
      "N_TOTAL",
      "N_REMOVED_1ST", "N_REMOVED_2ND",
      "n_pairs_total", "n_unique_related", "n_not_in_related",
      "n_unique_1st_or_closer",
      "hist_bins", "hist_xmin", "hist_xmax"
    ),
    value = c(
      cut_2nd, cut_1st, cut_dup,
      n_total,
      n_removed_1st, n_removed_2nd,
      nrow(pairs), n_unique_related, n_not_in_related,
      n_unique_1st_or_closer,
      hist_bins, cut_2nd, xmax
    )
  )

  fwrite(meta, meta_path)
  fwrite(donut_tbl, donut_path)
  fwrite(cat_summary, cat_path)
  fwrite(hist_freq, hist_path)
  message("Wrote frequency tables under: ", data_dir)
}

needed <- c(meta_path, donut_path, cat_path, hist_path)
missing <- needed[!file.exists(needed)]
if (length(missing) > 0L) {
  stop(
    "Frequency tables not found:\n  ", paste(missing, collapse = "\n  "),
    "\nPlace private .kin0 at:\n  ", kin0_path,
    "\nand re-run, or provide the four CSV tables under:\n  ", data_dir
  )
}

# -----------------------------------------------------------------------------
# 2) Read frequency tables and plot
# -----------------------------------------------------------------------------
meta <- fread(meta_path)
get_meta <- function(k) as.numeric(meta[param == k, value][1])
cut_2nd <- get_meta("cut_2nd")
cut_1st <- get_meta("cut_1st")
cut_dup <- get_meta("cut_dup")
N_TOTAL <- as.integer(get_meta("N_TOTAL"))
N_REMOVED_1ST <- as.integer(get_meta("N_REMOVED_1ST"))
N_REMOVED_2ND <- as.integer(get_meta("N_REMOVED_2ND"))
n_pairs_total <- as.integer(get_meta("n_pairs_total"))
n_unique_related <- as.integer(get_meta("n_unique_related"))
n_not_in_related <- as.integer(get_meta("n_not_in_related"))
n_unique_1st_or_closer <- as.integer(get_meta("n_unique_1st_or_closer"))
hist_xmin <- get_meta("hist_xmin")
hist_xmax <- get_meta("hist_xmax")

donut_raw <- fread(donut_path)
cat_summary <- fread(cat_path)
hist_freq <- fread(hist_path)

required_meta <- c(
  "N_TOTAL", "N_REMOVED_1ST", "N_REMOVED_2ND",
  "n_pairs_total", "n_unique_related", "n_not_in_related",
  "n_unique_1st_or_closer", "cut_2nd", "cut_1st", "cut_dup"
)
missing_meta <- required_meta[vapply(required_meta, function(k) is.na(get_meta(k)), logical(1))]
if (length(missing_meta) > 0L) {
  stop(
    "Incomplete kinship_meta.csv (missing: ", paste(missing_meta, collapse = ", "),
    "). Re-run with raw .kin0 and KING cutoff ID files present under:\n  ", raw_dir
  )
}

pct_of_cohort <- function(n) 100 * as.numeric(n) / as.numeric(N_TOTAL)

COL_REL  <- "#6F359B"
COL_IND  <- "#E2D7EB"
COL_HIST <- "grey70"

out_png1 <- file.path(out_dir, "Figure2B_kinship_donut.png")
out_png2 <- file.path(out_dir, "Figure2B_kinship_hist.png")

lab_not <- donut_raw[group == "not_in_related_pairs", label][1]
lab_in  <- donut_raw[group == "in_related_pairs", label][1]
donut <- data.table(
  group = factor(c(lab_not, lab_in), levels = c(lab_not, lab_in)),
  n = c(
    donut_raw[group == "not_in_related_pairs", n][1],
    donut_raw[group == "in_related_pairs", n][1]
  )
)
donut[, pct := n / sum(n)]
donut_cols <- setNames(c(COL_IND, COL_REL), c(lab_not, lab_in))

p_donut <- ggplot(donut, aes(x = 2.0, y = pct, fill = group)) +
  geom_col(width = 0.92, color = "white", linewidth = 2.2) +
  coord_polar(theta = "y", clip = "off") +
  xlim(0.55, 3.90) +
  scale_fill_manual(values = donut_cols) +
  theme_void(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "none",
    plot.margin = margin(20, 55, 20, 55)
  ) +
  geom_text(
    data = donut[group == lab_not],
    aes(
      label = paste0(
        "Individuals not in related pairs\n",
        comma(n), " (", percent(pct, accuracy = 0.1), ")"
      )
    ),
    position = position_stack(vjust = 0.5),
    size = 4.8,
    color = "black",
    lineheight = 1.05
  ) +
  annotate(
    "segment",
    x = 2.08, xend = 2.60, y = 0.060, yend = 0.060,
    color = "black", linewidth = 0.7, linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 2.95, y = 0.028,
    label = paste0(
      lab_in, "\n",
      comma(donut[group == lab_in, n]),
      " (", percent(donut[group == lab_in, pct], accuracy = 0.1), ")"
    ),
    size = 4.1, color = "black", lineheight = 1.05, family = "Arial"
  )

ggsave(out_png1, p_donut, width = 8.5, height = 6.2, dpi = 300, bg = "white")
message("[WROTE] ", out_png1)

bin_w <- if (nrow(hist_freq) > 1L) {
  median(diff(hist_freq$bin_left), na.rm = TRUE)
} else {
  hist_freq$bin_right[1] - hist_freq$bin_left[1]
}

p_hist <- ggplot(hist_freq, aes(x = bin_mid, y = n)) +
  geom_col(width = bin_w, fill = COL_HIST, color = "white", linewidth = 0.15) +
  geom_vline(xintercept = cut_2nd, linetype = "dashed", linewidth = 0.8, color = "red") +
  geom_vline(xintercept = cut_1st, linetype = "dashed", linewidth = 0.8, color = "red") +
  geom_vline(xintercept = cut_dup, linetype = "dashed", linewidth = 0.8, color = "red") +
  annotate(
    "text", x = cut_2nd, y = Inf, vjust = 2.0, hjust = -0.05,
    label = sprintf("2nd-degree\n(K = %.4f)", cut_2nd),
    size = 3.8, family = "Arial", color = "red"
  ) +
  annotate(
    "text", x = cut_1st, y = Inf, vjust = 2.0, hjust = -0.05,
    label = sprintf("1st-degree\n(K = %.3f)", cut_1st),
    size = 3.8, family = "Arial", color = "red"
  ) +
  annotate(
    "text", x = cut_dup, y = Inf, vjust = 2.0, hjust = -0.05,
    label = sprintf("Duplicate/MZ\n(K = %.3f)", cut_dup),
    size = 3.8, family = "Arial", color = "red"
  ) +
  coord_cartesian(xlim = c(hist_xmin, hist_xmax), clip = "on") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.00, 0.05))) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(size = 12, hjust = 0, color = "grey30"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.margin = margin(t = 40, r = 10, b = 30, l = 10)
  ) +
  labs(
    title = "Distribution of kinship coefficients",
    subtitle = sprintf("Filtered pairs with kinship coefficient \u2265 %.4f", cut_2nd),
    x = "KING kinship coefficient",
    y = "Number of pairs"
  )

ggsave(out_png2, p_hist, width = 10.0, height = 6.0, dpi = 300, bg = "white")
message("[WROTE] ", out_png2)

summary_lines <- c(
  "Relatedness summary (Figure 2B)",
  "===============================",
  sprintf("Total cohort samples (N): %s", comma(N_TOTAL)),
  sprintf(
    "Related pairs with kinship >= %.4f: %s",
    cut_2nd, comma(n_pairs_total)
  ),
  sprintf(
    "Unique individuals in related pairs: %s (%.1f%% of cohort)",
    comma(n_unique_related), pct_of_cohort(n_unique_related)
  ),
  sprintf(
    "Individuals not in related pairs: %s (%.1f%% of cohort)",
    comma(n_not_in_related), pct_of_cohort(n_not_in_related)
  ),
  sprintf(
    "Unique individuals in >=1 first-degree pair (kinship >= %.3f): %s (%.1f%% of cohort)",
    cut_1st, comma(n_unique_1st_or_closer), pct_of_cohort(n_unique_1st_or_closer)
  ),
  "",
  "Pair counts by relationship category:",
  paste0(
    "  - ", cat_summary$rel_cat, ": ", comma(cat_summary$n_pairs),
    " pairs; ", comma(cat_summary$n_unique_individuals), " unique individuals"
  ),
  "",
  "KING cutoff removals (unrelated-set construction):",
  sprintf(
    "  - Removed by 1st-degree cutoff (kinship >= %.3f): %s",
    cut_1st, comma(N_REMOVED_1ST)
  ),
  sprintf(
    "  - Removed by 2nd-degree cutoff (kinship >= %.4f): %s",
    cut_2nd, comma(N_REMOVED_2ND)
  )
)
cat(paste(summary_lines, collapse = "\n"), "\n", sep = "")
