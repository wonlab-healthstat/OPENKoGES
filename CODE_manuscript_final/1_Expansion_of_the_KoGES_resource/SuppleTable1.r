# Supplementary Table 1: Summary of data availability in the KoGES cohort
#
# This script generates a cohort-by-follow-up summary of:
#   - survey period
#   - number of participants
#   - number of genotyped participants
#   - number of variables
#
# Data workflow (sample-level data are NOT for public release):
#   1) Read sample-level IDs from DATA_.../SuppleTable1/raw/
#   2) Merge with publishable codebook metadata and write a summary table
#   3) Format and write the manuscript Supplementary Table 1 CSV
#
# Inputs (raw, private):
#   DATA_manuscript_final/SuppleTable1/raw/cohort_participant_ids.RData
#   DATA_manuscript_final/SuppleTable1/raw/genome_ids.txt
# Inputs (publishable metadata):
#   DATA_manuscript_final/SuppleTable1/codebook_meta.csv
# Intermediate (publishable):
#   DATA_manuscript_final/SuppleTable1/data_availability_summary.csv
# Outputs:
#   RESULT_manuscript_final/SuppleTable1/SuppleTable1.csv

suppressPackageStartupMessages({
  library(dplyr)
})

data_dir   <- "/home2/sjim/openKoGES/DATA_manuscript_final/SuppleTable1"
raw_dir    <- file.path(data_dir, "raw")
output_dir <- "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleTable1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta_path    <- file.path(data_dir, "codebook_meta.csv")
summary_path <- file.path(data_dir, "data_availability_summary.csv")
out_path     <- file.path(output_dir, "SuppleTable1.csv")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
ordinal_suffix <- function(n) {
  if (n %% 10 == 1 && n %% 100 != 11) return("st")
  if (n %% 10 == 2 && n %% 100 != 12) return("nd")
  if (n %% 10 == 3 && n %% 100 != 13) return("rd")
  "th"
}

follow_up_label <- function(fu) {
  if (fu == 0) "Baseline" else paste0(fu, ordinal_suffix(fu), " Follow-up")
}

# -----------------------------------------------------------------------------
# 1) Build publishable summary from private sample-level IDs + codebook meta
# -----------------------------------------------------------------------------
if (file.exists(file.path(raw_dir, "cohort_participant_ids.RData")) &&
    file.exists(file.path(raw_dir, "genome_ids.txt")) &&
    file.exists(meta_path)) {

  load(file.path(raw_dir, "cohort_participant_ids.RData"))  # dat1: named list of IDs
  genome_ids_df <- read.table(file.path(raw_dir, "genome_ids.txt"),
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  genome_ids <- unique(genome_ids_df[[colnames(genome_ids_df)[1]]])

  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  meta$participants <- NA_integer_
  meta$genome_participants <- NA_integer_

  for (i in seq_len(nrow(meta))) {
    cohort_name <- meta$cohort[i]
    if (!cohort_name %in% names(dat1)) next
    sample_ids <- dat1[[cohort_name]]
    meta$participants[i] <- length(sample_ids)
    meta$genome_participants[i] <- sum(sample_ids %in% genome_ids, na.rm = TRUE)
  }

  write.csv(meta, summary_path, row.names = FALSE)
  message("Wrote: ", summary_path)
}

if (!file.exists(summary_path)) {
  stop("Summary table not found. Place raw sample-level IDs under raw/ ",
       "(with codebook_meta.csv) and re-run, or provide data_availability_summary.csv.")
}

# -----------------------------------------------------------------------------
# 2) Format Supplementary Table 1
# -----------------------------------------------------------------------------
summary_tbl <- read.csv(summary_path, stringsAsFactors = FALSE)

TableS1 <- summary_tbl %>%
  mutate(
    cohort_type = substr(cohort, 1, 2),
    FU = as.numeric(substr(cohort, 3, 4)) - 1,
    follow_up_label = vapply(FU, follow_up_label, character(1)),
    `Survey Period` = survey_period,
    `Number of Participants(n)` = participants,
    `Number of Genotyped Participants(n)` = genome_participants,
    `Number of Variables(n)` = variable_count
  ) %>%
  arrange(cohort_type, FU) %>%
  select(
    cohort_type,
    follow_up_label,
    `Survey Period`,
    `Number of Participants(n)`,
    `Number of Genotyped Participants(n)`,
    `Number of Variables(n)`
  )

write.csv(TableS1, out_path, row.names = FALSE)
message("Wrote: ", out_path)
