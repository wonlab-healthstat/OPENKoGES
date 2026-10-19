# ==============================================================================
# Shared helpers for Cont/Cate ExtendedData2_*.R
# Sourced by ExtendedData2_cont.R and ExtendedData2_cate.R
#
# Excel outputs (set in the calling scripts):
#   Cont -> RESULT_manuscript_final/ExtendedData/ExtendedData2_continuous.xlsx
#   Cate -> RESULT_manuscript_final/ExtendedData/ExtendedData2_categorical.xlsx
# ==============================================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("package 'data.table' is required")
}
suppressPackageStartupMessages(library(data.table))

require_cols <- function(dt, cols, label) {
  miss <- setdiff(cols, names(dt))
  if (length(miss) > 0L) {
    stop(label, " missing required column(s): ", paste(miss, collapse = ", "))
  }
}

## Summary_* sheets may store R-style dotted headers even when sep.names=" ".
restore_space_colnames <- function(dt) {
  nms <- names(dt)
  if (any(grepl(".", nms, fixed = TRUE))) {
    setnames(dt, nms, gsub(".", " ", nms, fixed = TRUE))
  }
  invisible(dt)
}

as_flag <- function(x) {
  if (is.logical(x)) return(x)
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
  out[y %in% c("FALSE", "F", "0", "NO", "N")] <- FALSE
  out
}

norm_chr <- function(x) {
  y <- trimws(as.character(x))
  y <- sub("^chr", "", y, ignore.case = TRUE)
  y
}

norm_pos <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  out <- trimws(as.character(x))
  ok <- is.finite(y)
  out[ok] <- format(round(y[ok]), scientific = FALSE, trim = TRUE)
  out
}

clean_status <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(y) | y == "" | y %in% c("—", "–")] <- "Not assessed"
  y
}

normalize_novelty <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep("Not assessed", length(y))
  out[which(!is.na(y) & y == "novel")] <- "Novel"
  out[which(!is.na(y) & y == "known")] <- "Known"
  out
}

replication_status <- function(flag, analyzed) {
  flag <- as_flag(flag)
  analyzed <- as_flag(analyzed)
  out <- rep("Not analyzed", length(flag))
  out[analyzed %in% TRUE & flag %in% TRUE] <- "Replicated"
  out[analyzed %in% TRUE & flag %in% FALSE] <- "Not replicated"
  out[analyzed %in% TRUE & is.na(flag)] <- "Not available"
  out
}

direction_status <- function(flag, analyzed) {
  flag <- as_flag(flag)
  analyzed <- as_flag(analyzed)
  out <- rep("Not analyzed", length(flag))
  out[analyzed %in% TRUE & flag %in% TRUE] <- "Concordant"
  out[analyzed %in% TRUE & flag %in% FALSE] <- "Discordant"
  out[analyzed %in% TRUE & is.na(flag)] <- "Not available"
  out
}

yes_no <- function(x) {
  f <- as_flag(x)
  fifelse(f %in% TRUE, "Yes", fifelse(f %in% FALSE, "No", "Not available"))
}

as_num <- function(x) suppressWarnings(as.numeric(x))

make_identity_key <- function(phenotype, chr, bp, rsid) {
  paste(
    trimws(as.character(phenotype)),
    norm_chr(chr),
    norm_pos(bp),
    trimws(as.character(rsid)),
    sep = "\t"
  )
}

append_dt_sheet <- function(wb, dt, dest_sheet, tab_colour = NULL) {
  if (dest_sheet %in% wb$sheet_names) {
    stop("Destination sheet already exists: ", dest_sheet)
  }
  openxlsx::addWorksheet(wb, dest_sheet, tabColour = tab_colour)
  openxlsx::writeData(wb, dest_sheet, as.data.frame(dt), withFilter = TRUE)
  openxlsx::freezePane(wb, dest_sheet, firstRow = TRUE, firstCol = TRUE)
  invisible(list(sheet = dest_sheet, nrow = nrow(dt), ncol = ncol(dt)))
}

make_meta_block <- function(d, panel) {
  analyzed <- as_flag(d$panel_analyzed)
  nominal <- as_flag(d$replication_nominal)
  bonf <- as_flag(d$replication_Bonferroni)
  signed_bonf <- as_flag(d$signedZ_replication_Bonferroni)
  direction <- as_flag(d$direction_concordant_with_KoGES)

  out <- data.table(lead_row_id = d$lead_row_id)

  out[[paste0(panel, " panel meta analyzed")]] <-
    fifelse(analyzed %in% TRUE, "Analyzed", "Not analyzed")

  out[[paste0(panel, " panel contributing cohorts")]] <-
    as.character(d$contributing_cohorts)

  out[[paste0(panel, " panel meta beta")]] <-
    as_num(d$external_panel_beta)

  out[[paste0(panel, " panel meta standard error")]] <-
    as_num(d$external_panel_SE)

  out[[paste0(panel, " panel meta P value")]] <-
    as_num(d$external_panel_P)

  out[[paste0(panel, " meta direction vs KoGES")]] <-
    direction_status(direction, analyzed)

  out[[paste0(panel, " meta nominal replication")]] <-
    replication_status(nominal, analyzed)

  out[[paste0(panel, " meta Bonferroni replication")]] <-
    replication_status(bonf, analyzed)

  out[[paste0(panel, " signed-Z Bonferroni sensitivity")]] <-
    replication_status(signed_bonf, analyzed)

  out
}
