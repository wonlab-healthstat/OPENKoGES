#!/usr/bin/env Rscript
# ==============================================================================
# Purpose
#   Lift GWAS summary / lead-variant TSV coordinates from hg19 (GRCh37) to
#   GRCh38 using UCSC liftOver and the hg19ToHg38 chain file. The same table
#   schema is retained; CHR/BP columns are overwritten with GRCh38 positions
#   (optionally keeping the original position in BP_GRCh37).
#
#   Downstream GWAS Catalog novelty scripts expect GRCh38 lead TSVs
#   (e.g. Cate_lead_variants_hg38.tsv / Cont_lead_variants_hg38.tsv).
#
# Run order (typical manuscript pipeline)
#   1. Prepare lead-variant TSVs on GRCh37 (Extended_Data/{Cate,Cont}_lead_variants.tsv)
#   2. Rscript liftover_sumstat.R <input.tsv> <output_hg38.tsv>
#   3. (Optional) Fill remaining unmapped loci via <stem>_unmapped_add_info.tsv
#   4. Run gwas_catalog_*_novelty scripts on the hg38 output
#
# Inputs
#   - TSV with chromosome and base-pair columns (auto-detected or via --chr/--bp)
#   - UCSC liftOver binary (default: /data/bin/liftOver)
#   - Chain: hg19ToHg38.over.chain.gz
#   - Optional companion: <output_stem>_unmapped_add_info.tsv next to the output
#       (liftover_note column = manual integer GRCh38 BP for failed loci only)
#
# Outputs
#   - Main TSV with CHR/BP updated to GRCh38 (plus BP_GRCh37 unless --no-keep-hg19-bp)
#   - Companion <output>_unmapped.tsv for loci that failed liftOver (unless
#     --no-unmapped-out); includes liftover_note when available from liftOver
#
# Usage
#   Rscript liftover_sumstat.R input.tsv output_hg38.tsv
#   Rscript liftover_sumstat.R input.tsv output_hg38.tsv --chr CHR --bp BP
#   Rscript liftover_sumstat.R input.tsv out.tsv --no-keep-hg19-bp
#   Rscript liftover_sumstat.R input.tsv out.tsv --unmapped-out /path/unmapped.tsv
#
# Examples (manuscript Extended Data leads)
#   Rscript liftover_sumstat.R \
#     .../Extended_Data/Cont_lead_variants.tsv \
#     .../Extended_Data/Cont_lead_variants_hg38.tsv
#   Rscript liftover_sumstat.R \
#     .../Extended_Data/Cate_lead_variants.tsv \
#     .../Extended_Data/Cate_lead_variants_hg38.tsv
#
# Notes
#   - BED conversion uses 0-based half-open intervals (start = BP-1, end = BP).
#   - Exit status is 0 even when some loci remain unmapped (written to companion TSV).
#   - Manual add_info updates only __lift_fail rows; successful automatic mappings
#     are never overwritten.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

DEFAULT_LIFTOVER <- "/data/bin/liftOver"
DEFAULT_CHAIN <- "/home2/sjim/openKoGES/DATA_baseline/forGWAS_annot/hg19ToHg38.over.chain.gz"

parse_args <- function(argv) {
  if (length(argv) < 2L) {
    cat("Usage: Rscript liftover_sumstat.R <input.tsv> <output.tsv> [options]\n",
        "Options:\n",
        "  --chr COL              CHR column name (default: CHR, Chromosome, ...)\n",
        "  --bp COL               BP/POS column name (default: BP, Base pair position, ...)\n",
        "  --liftover PATH        liftOver binary (default: ", DEFAULT_LIFTOVER, ")\n",
        "  --chain PATH           chain file (default: ", DEFAULT_CHAIN, ")\n",
        "  --no-keep-hg19-bp      do not add BP_GRCh37 column\n",
        "  --unmapped-out PATH    write liftOver-failed rows (NA) to this TSV\n",
        "                         (default: <output>_unmapped.tsv next to main output)\n",
        "  --no-unmapped-out      do not write unmapped locus file\n",
        "  --temp-dir DIR         temp directory (default: same dir as output)\n",
        sep = "")
    quit(status = 1L)
  }
  inp <- argv[1L]
  out <- argv[2L]
  rest <- argv[-(1:2)]
  opt <- list(
    input = inp,
    output = out,
    chr = NA_character_,
    bp = NA_character_,
    liftover = DEFAULT_LIFTOVER,
    chain = DEFAULT_CHAIN,
    keep_hg19_bp = TRUE,
    unmapped_out = NA_character_,
    write_unmapped = TRUE,
    temp_dir = NA_character_
  )
  i <- 1L
  while (i <= length(rest)) {
    a <- rest[i]
    if (a == "--chr") {
      opt$chr <- rest[i + 1L]; i <- i + 2L
    } else if (a == "--bp") {
      opt$bp <- rest[i + 1L]; i <- i + 2L
    } else if (a == "--liftover") {
      opt$liftover <- rest[i + 1L]; i <- i + 2L
    } else if (a == "--chain") {
      opt$chain <- rest[i + 1L]; i <- i + 2L
    } else if (a == "--no-keep-hg19-bp") {
      opt$keep_hg19_bp <- FALSE
      i <- i + 1L
    } else if (a == "--unmapped-out") {
      opt$unmapped_out <- rest[i + 1L]; i <- i + 2L
    } else if (a == "--no-unmapped-out") {
      opt$write_unmapped <- FALSE
      i <- i + 1L
    } else if (a == "--temp-dir") {
      opt$temp_dir <- rest[i + 1L]; i <- i + 2L
    } else {
      stop("Unknown argument: ", a)
    }
  }
  opt
}

# Normalize chromosome labels to UCSC BED style (chr1, chrX, chrM, ...).
norm_chr_bed <- function(chr) {
  x <- trimws(as.character(chr))
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x %in% c("23", "X")] <- "X"
  x[x %in% c("24", "Y")] <- "Y"
  x[x %in% c("25", "MT", "M")] <- "M"
  paste0("chr", x)
}

resolve_col <- function(names_vec, primary, fallbacks) {
  if (!is.na(primary) && primary %in% names_vec) return(primary)
  for (f in fallbacks) {
    if (f %in% names_vec) return(f)
  }
  NULL
}

# Parse liftOver unmapped.bed: "#reason" line followed by BED + row id (4th field).
parse_liftover_unmap_notes <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L) {
    return(data.table(`__row` = integer(), liftover_note = character()))
  }
  txt <- readLines(path, warn = FALSE)
  out <- list()
  i <- 1L
  n <- length(txt)
  while (i <= n) {
    line <- txt[i]
    if (nzchar(line) && startsWith(line, "#")) {
      note <- sub("^#\\s*", "", line)
      i <- i + 1L
      if (i <= n) {
        parts <- strsplit(txt[i], "\t", fixed = TRUE)[[1L]]
        if (length(parts) >= 4L) {
          rid <- suppressWarnings(as.integer(parts[4L]))
          if (!is.na(rid)) {
            out[[length(out) + 1L]] <- data.table(`__row` = rid, liftover_note = note)
          }
        }
      }
    }
    i <- i + 1L
  }
  if (length(out)) rbindlist(out) else data.table(`__row` = integer(), liftover_note = character())
}

default_unmapped_path <- function(main_out) {
  if (grepl("\\.[Tt][Ss][Vv]$", main_out)) {
    sub("\\.[Tt][Ss][Vv]$", "_unmapped.tsv", main_out)
  } else {
    paste0(main_out, "_unmapped.tsv")
  }
}

normalize_chr_for_join <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x %in% c("23", "X")] <- "X"
  x[x %in% c("24", "Y")] <- "Y"
  x[x %in% c("25", "MT", "M")] <- "M"
  x
}

# Locate <output_stem>_unmapped_add_info.tsv next to the main output
# (or the sole *_unmapped_add_info.tsv in that directory).
find_add_info_path <- function(out_path) {
  out_dir <- dirname(normalizePath(out_path, mustWork = FALSE))
  if (!nzchar(out_dir) || !dir.exists(out_dir)) return(NA_character_)
  p1 <- if (grepl("\\.[Tt][Ss][Vv]$", out_path)) {
    sub("\\.[Tt][Ss][Vv]$", "_unmapped_add_info.tsv", out_path)
  } else {
    paste0(out_path, "_unmapped_add_info.tsv")
  }
  if (file.exists(p1)) return(p1)
  g <- Sys.glob(file.path(out_dir, "*_unmapped_add_info.tsv"))
  if (length(g) == 1L) return(g[1L])
  if (length(g) > 1L) {
    message("[liftover_sumstat] multiple *_unmapped_add_info.tsv in dir; using: ", g[1L])
    return(g[1L])
  }
  NA_character_
}

# Manual GRCh38 BP rescue from add_info:
#   liftover_note = integer GRCh38 BP; non-numeric notes are skipped.
#   Matching keys are GRCh37 (hg19) CHR/BP; only __lift_fail rows are updated
#   (successful automatic liftOver rows are never overwritten).
apply_manual_hg38_add_info <- function(dt, path, chr_col, bp_col, keep_hg19) {
  ai <- fread(path, sep = "auto", header = TRUE, data.table = TRUE, showProgress = FALSE)
  chr_ai <- resolve_col(names(ai), chr_col,
                        c("CHR", "#CHR", "Chromosome", "chrom", "CHROM"))
  bp_ai <- resolve_col(names(ai), bp_col,
                       c("BP", "POS", "pos", "GENPOS", "Base pair position"))
  if (is.null(chr_ai)) {
    stop("[add_info] could not find CHR/Chromosome column (same names as main table preferred)")
  }
  if (is.null(bp_ai)) {
    stop("[add_info] could not find BP/Base pair position column")
  }
  if (!"liftover_note" %in% names(ai)) {
    stop("[add_info] missing column: liftover_note")
  }
  ai[, `__kchr` := normalize_chr_for_join(get(chr_ai))]
  ai[, `__kbp` := suppressWarnings(as.integer(get(bp_ai)))]
  ai[, `__hg38` := suppressWarnings(as.integer(trimws(as.character(liftover_note))))]
  v <- ai[!is.na(`__hg38`) & !is.na(`__kbp`) & !is.na(`__kchr`)]
  if (nrow(v) == 0L) {
    message("[liftover_sumstat] add_info: no rows with numeric liftover_note (hg38 BP); skipped")
    return(invisible(0L))
  }
  n_ai <- nrow(v)
  v <- unique(v, by = c("__kchr", "__kbp"))
  if (nrow(v) < n_ai) {
    message("[liftover_sumstat] add_info: duplicate CHR/BP keys; using first per key")
  }
  dt[, `__kchr` := normalize_chr_for_join(get(chr_col))]
  dt[, `__kbp` := `__bp_raw`]
  lf <- dt[`__lift_fail` == TRUE]
  if (nrow(lf) == 0L) {
    message("[liftover_sumstat] add_info: no liftover-failed rows to update; skipped")
    return(invisible(0L))
  }
  merged <- lf[v[, .(`__kchr`, `__kbp`, `__hg38`)],
                 on = c("__kchr", "__kbp"), nomatch = 0L]
  if (nrow(merged) == 0L) {
    message("[liftover_sumstat] add_info: no CHR/BP overlap with current unmapped rows; skipped")
    return(invisible(0L))
  }
  n_before <- nrow(dt[`__lift_fail` == TRUE])
  if (isTRUE(keep_hg19)) {
    dt[merged, BP_GRCh37 := i.__bp_raw, on = "__row"]
  }
  dt[merged, (bp_col) := i.__hg38, on = "__row"]
  dt[merged, `__lift_fail` := FALSE, on = "__row"]
  n_after <- nrow(dt[`__lift_fail` == TRUE])
  message(sprintf(
    "[liftover_sumstat] add_info applied from %s | manual fixes: %d | still unmapped: %d",
    path, n_before - n_after, n_after
  ))
  invisible(n_before - n_after)
}

argv <- commandArgs(trailingOnly = TRUE)
opt <- parse_args(argv)

if (!file.exists(opt$input)) stop("Input not found: ", opt$input)
if (!file.exists(opt$liftover)) stop("liftOver not found: ", opt$liftover)
if (!file.exists(opt$chain)) stop("Chain file not found: ", opt$chain)

dt <- fread(opt$input, sep = "auto", header = TRUE, data.table = TRUE, showProgress = FALSE)
n_in <- nrow(dt)
if (n_in == 0L) stop("Input has 0 rows.")

chr_col <- resolve_col(names(dt), opt$chr,
                       c("CHR", "#CHR", "Chromosome", "chrom", "CHROM"))
bp_col <- resolve_col(names(dt), opt$bp,
                      c("BP", "POS", "pos", "GENPOS", "Base pair position"))
if (is.null(chr_col)) stop("Could not find CHR column. Use --chr.")
if (is.null(bp_col)) stop("Could not find BP/POS column. Use --bp.")

chr_type <- typeof(dt[[chr_col]])

dt[, `__row` := .I]
dt[, `__chr_raw` := get(chr_col)]
dt[, `__bp_raw` := suppressWarnings(as.integer(get(bp_col)))]
dt[, `__ok` := !is.na(`__bp_raw`) & !is.na(`__chr_raw`) & nzchar(trimws(as.character(`__chr_raw`)))]

dt_ok <- dt[`__ok` == TRUE]
n_skip <- n_in - nrow(dt_ok)
if (nrow(dt_ok) == 0L) stop("No rows with valid CHR and integer BP.")

out_dir <- dirname(normalizePath(opt$output, mustWork = FALSE))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
tmp_base <- if (!is.na(opt$temp_dir) && nzchar(opt$temp_dir)) opt$temp_dir else out_dir
bed_in <- tempfile("sumstat_", fileext = ".bed", tmpdir = tmp_base)
bed_out <- sub("\\.bed$", "_hg38.bed", bed_in)
bed_unmap <- sub("\\.bed$", "_unmapped.bed", bed_in)

chr_bed <- norm_chr_bed(dt_ok[["__chr_raw"]])
bpv <- dt_ok[["__bp_raw"]]
# BED: 0-based start, end exclusive (same as liftover.R awk pipeline: BP-1, BP)
bed_lines <- paste(chr_bed, pmax(0L, bpv - 1L), bpv, dt_ok[["__row"]], sep = "\t")
writeLines(bed_lines, bed_in)

rc <- system2(opt$liftover,
              args = c(bed_in, opt$chain, bed_out, bed_unmap),
              stdout = TRUE, stderr = TRUE)
if (!is.null(attr(rc, "status")) && attr(rc, "status") != 0L) {
  stop("liftOver failed:\n", paste(rc, collapse = "\n"))
}

if (!file.exists(bed_out) || file.info(bed_out)$size == 0L) {
  stop("liftOver produced no mapped output. Check coordinates and chain file.")
}

mapped <- fread(bed_out, header = FALSE, sep = "\t", data.table = TRUE, showProgress = FALSE)
# V1=chr, V2=start0, V3=end, V4=row id
setnames(mapped, c("V1", "V2", "V3", "V4"), c("bed_chr", "bed_start0", "bed_end", "map_row"))
mapped[, `__row` := as.integer(map_row)]
mapped[, new_bp := ifelse(bed_end - bed_start0 == 1L, as.integer(bed_start0) + 1L,
                          as.integer(round((bed_start0 + bed_end) / 2)))]
mapped[, new_chr := sub("^chr", "", bed_chr, ignore.case = TRUE)]

# Prefer first mapping when a locus splits into multiple intervals
mapped <- unique(mapped, by = "__row")

dt[mapped, c("BP_new", "CHR_new") := list(i.new_bp, i.new_chr), on = "__row"]

dt[, `__lift_fail` := `__ok` & is.na(BP_new)]

n_mapped <- dt[!is.na(BP_new), .N]
n_unmapped <- dt_ok[, .N] - mapped[, uniqueN(`__row`)]

unmap_note <- parse_liftover_unmap_notes(bed_unmap)

if (opt$keep_hg19_bp) {
  dt[`__ok` & !is.na(BP_new), BP_GRCh37 := get(bp_col)]
}
dt[`__ok` & !is.na(BP_new), (bp_col) := BP_new]
if (chr_type == "integer") {
  dt[`__ok` & !is.na(BP_new), (chr_col) := suppressWarnings(as.integer(CHR_new))]
} else if (chr_type == "double") {
  dt[`__ok` & !is.na(BP_new), (chr_col) := suppressWarnings(as.numeric(CHR_new))]
} else {
  dt[`__ok` & !is.na(BP_new), (chr_col) := CHR_new]
}

add_path <- find_add_info_path(opt$output)
if (!is.na(add_path) && nzchar(add_path) && file.exists(add_path)) {
  apply_manual_hg38_add_info(dt, add_path, chr_col, bp_col, opt$keep_hg19_bp)
}

dt_unmapped <- copy(dt[`__lift_fail` == TRUE])
if (nrow(dt_unmapped) > 0L && nrow(unmap_note) > 0L) {
  dt_unmapped[unmap_note, liftover_note := i.liftover_note, on = "__row"]
}
if (nrow(dt_unmapped) > 0L && !"liftover_note" %in% names(dt_unmapped)) {
  dt_unmapped[, liftover_note := NA_character_]
} else if (nrow(dt_unmapped) > 0L) {
  dt_unmapped[is.na(liftover_note), liftover_note := NA_character_]
}

n_unmapped_final <- dt[`__lift_fail` == TRUE, .N]

drop_internal <- intersect(c("__row", "__chr_raw", "__bp_raw", "__ok", "BP_new", "CHR_new",
                             "__lift_fail", "__kchr", "__kbp"), names(dt))
if (length(drop_internal)) dt[, (drop_internal) := NULL]

fwrite(dt, opt$output, sep = "\t", quote = FALSE, na = "NA")

unmapped_path <- NA_character_
if (isTRUE(opt$write_unmapped) && nrow(dt_unmapped) > 0L) {
  unmapped_path <- if (!is.na(opt$unmapped_out) && nzchar(opt$unmapped_out)) {
    opt$unmapped_out
  } else {
    default_unmapped_path(opt$output)
  }
  udir <- dirname(normalizePath(unmapped_path, mustWork = FALSE))
  if (!dir.exists(udir)) dir.create(udir, recursive = TRUE)
  du_drop <- intersect(c("__row", "__chr_raw", "__bp_raw", "__ok", "BP_new", "CHR_new",
                         "__lift_fail", "__kchr", "__kbp"), names(dt_unmapped))
  if (length(du_drop)) dt_unmapped[, (du_drop) := NULL]
  fwrite(dt_unmapped, unmapped_path, sep = "\t", quote = FALSE, na = "NA")
}

msg_end <- normalizePath(opt$output, mustWork = FALSE)
if (!is.na(unmapped_path) && nzchar(unmapped_path)) {
  msg_end <- paste0(msg_end, " | unmapped: ", normalizePath(unmapped_path, mustWork = FALSE))
}
message(sprintf(
  "[liftover_sumstat] rows in: %d | skipped (bad CHR/BP): %d | mapped: %d | unmapped (after add_info): %d | written: %s",
  n_in, n_skip, n_mapped, n_unmapped_final, msg_end
))

unlink(c(bed_in, bed_out, bed_unmap), force = TRUE)

# Successful completion should return exit status 0 even when some loci remain
# unmapped and are written to the companion TSV.
quit(save = "no", status = 0L)
