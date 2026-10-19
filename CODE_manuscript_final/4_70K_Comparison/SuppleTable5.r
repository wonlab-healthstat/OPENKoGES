#!/usr/bin/env Rscript
# ============================================================
# SuppleTable5.r
#
# Supplementary Table 5. Comparison of genome-wide significant loci
# in KoGES 160K and KoGES 70K across 34 shared traits.
#
# The source workbook stores the 1 Mb comparison in seven columns.
# This table keeps the five columns used in the manuscript:
#   Significant loci in KoGES160K
#   Significant loci in KoGES70K
#   Overlapping loci          = 160K loci with a 70K GWS variant within 1 Mb
#   Novel loci in KoGES160K   = 160K loci with no 70K GWS variant within 1 Mb
#   Loci only in KoGES70K     = 70K loci with no 160K GWS variant within 1 Mb
#
# Output:
#   RESULT_manuscript_final/SuppleTable5/SuppleTable5_cont.csv
#   RESULT_manuscript_final/SuppleTable5/SuppleTable5_cate.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

SRC <- "/data4/joh/openkoges/table_new0923/output/Supplementary_Table_4.xlsx"
OUT_DIR <- Sys.getenv(
  "SUPPLETABLE5_OUT",
  unset = "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleTable5"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

to_five <- function(sheet) {
  d <- as.data.table(read.xlsx(SRC, sheet = sheet, check.names = FALSE))
  if (ncol(d) != 7L) stop(sheet, " expected 7 columns, found ", ncol(d))
  out <- data.table(
    Phenotype = d[[1]],
    `Significant loci in KoGES160K` = d[[2]],
    `Significant loci in KoGES70K` = d[[5]],
    `Overlapping loci` = d[[3]],
    `Novel loci in KoGES160K` = d[[4]],
    `Loci only in KoGES70K` = d[[7]]
  )
  out
}

cont <- to_five("4a_continuous")
cate <- to_five("4b_binary")
if (nrow(cont) != 25L) stop("Continuous traits should be 25, found ", nrow(cont))
if (nrow(cate) != 9L) stop("Binary traits should be 9, found ", nrow(cate))

alb <- cont[Phenotype == "Albumin (Alb)"]
if (nrow(alb) != 1L || alb[[2]] != 123 || alb[[3]] != 63 || alb[[4]] != 86 || alb[[5]] != 37 || alb[[6]] != 18) {
  stop("Albumin row does not match 123, 63, 86, 37, 18")
}

cont_csv <- file.path(OUT_DIR, "SuppleTable5_cont.csv")
cate_csv <- file.path(OUT_DIR, "SuppleTable5_cate.csv")
fwrite(cont, cont_csv)
fwrite(cate, cate_csv)
message("Wrote ", cont_csv, " n=", nrow(cont))
message("Wrote ", cate_csv, " n=", nrow(cate))
