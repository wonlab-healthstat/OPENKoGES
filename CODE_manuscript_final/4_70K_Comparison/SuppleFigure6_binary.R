library(data.table)
library(parallel)
library(magick)
library(grid)

# ==============================================================================
# KoGES 70K vs 160K disease Manhattan overlay plot
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 경로 설정
# ------------------------------------------------------------------------------

# SE 보정 후 160K disease SAIGE summary statistics
# 예: B_DM2_SAIGE_cumulative_inci.SAIGE.SEfixed.txt
GWAS160K_DIR <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/GWAS_sumstat"

# 70K categorical GWAS summary statistics 경로
GWAS70K_DIR <- "/data6/OpenGWAS/Sumstats/KBA70K/Cate"

# 새 longitudinal clumping 및 ANNOVAR 결과 경로
CLUMP_ANNO_DIR <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/Clumping"

# Figure 출력 경로
OUT_DIR <- "/data4/joh/openkoges/figure/Overlaid_comparision_70Kvs160K/parallel_disease_260928_SEfixed"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. 비교 phenotype
# ------------------------------------------------------------------------------

pheno_map <- data.table(
  pheno_160K = c(
    "B_DM2", "B_HTN2", "B_LIP", "B_THY",
    "B_GCA", "B_BRCA", "B_THYCA"
  ),
  pheno_70K = c(
    "DM", "HTN", "LIP", "THYD",
    "GCA", "BRCA", "THYCA"
  )
)

# ------------------------------------------------------------------------------
# 2. phenotype 표시명
# ------------------------------------------------------------------------------

pheno_map2 <- c(
  DM    = "Diabetes mellitus",
  HTN   = "Hypertension",
  LIP   = "Hyperlipidemia",
  THYD  = "Thyroid disease",
  GCA   = "Gastric cancer",
  BRCA  = "Breast cancer",
  THYCA = "Thyroid cancer"
)

# ------------------------------------------------------------------------------
# 3. gene annotation 보조 함수
# ------------------------------------------------------------------------------

clean_gene <- function(gene) {
  if (is.na(gene) || gene == ".") {
    return(NA_character_)
  }
  
  genes <- unlist(strsplit(gene, "[;,]"))
  genes <- trimws(genes)
  
  genes[1]
}

make_gene_pos <- function(anno_dt) {
  gene_pos <- unique(
    anno_dt[
      Gene.refGene != "." & !is.na(Gene.refGene),
      .(
        Gene  = sapply(Gene.refGene, clean_gene),
        CHR   = as.numeric(Chr),
        START = as.numeric(Start),
        END   = as.numeric(End)
      )
      ]
  )
  
  gene_pos[!is.na(Gene)]
}

nearest_gene_factory <- function(gene_pos) {
  force(gene_pos)
  
  function(chr, pos) {
    sub <- gene_pos[CHR == chr]
    
    if (nrow(sub) == 0) {
      return(NA_character_)
    }
    
    dist <- pmin(
      abs(pos - sub$START),
      abs(pos - sub$END)
    )
    
    sub$Gene[which.min(dist)]
  }
}

# ------------------------------------------------------------------------------
# 4. phenotype별 데이터 불러오기
# ------------------------------------------------------------------------------

load_all_data <- function(pheno160, pheno70) {
  
  cat("Loading:", pheno160, "/", pheno70, "\n")
  
  # ------------------ 160K disease SAIGE summary statistics ------------------
  
  d1_path <- file.path(
    GWAS160K_DIR,
    paste0(pheno160, "_SAIGE_cumulative_inci.SAIGE.SEfixed.txt")
  )
  
  if (!file.exists(d1_path)) {
    message("160K GWAS file not found: ", d1_path)
    return(NULL)
  }
  
  d1 <- fread(d1_path)
  
  if (!"P" %in% names(d1) && "p.value" %in% names(d1)) {
    setnames(d1, "p.value", "P")
  }
  
  d1[, P := as.numeric(P)]
  d1[, CHR := as.numeric(CHR)]
  d1[, POS := as.numeric(POS)]
  
  d1 <- d1[
    CHR %in% 1:22 &
      !is.na(POS) &
      !is.na(P)
    ]
  
  setorder(d1, CHR, POS)
  
  # Manhattan plot x축 누적 순서
  d1[, ss := .I]
  
  chr_n <- table(d1$CHR)
  CM <- cumsum(chr_n)
  n.chr <- length(chr_n)
  
  # ------------------ 70K categorical GWAS ------------------
  
  d2_path <- file.path(
    GWAS70K_DIR,
    paste0("phenocode-KoGES_", pheno70, "_unique.txt")
  )
  
  if (!file.exists(d2_path)) {
    message("70K GWAS file not found: ", d2_path)
    return(NULL)
  }
  
  d2 <- fread(d2_path)
  
  if (!"P" %in% names(d2) && "p" %in% names(d2)) {
    setnames(d2, "p", "P")
  }
  
  d2[, P := as.numeric(P)]
  d2[, CHR := as.numeric(CHR)]
  d2[, POS := as.numeric(POS)]
  
  d2 <- d2[
    CHR %in% 1:22 &
      !is.na(POS) &
      !is.na(P)
    ]
  
  # 160K 좌표에 70K P-value 결합
  combined <- merge(
    d1[, .(CHR, POS, ss)],
    d2[, .(CHR, POS, P)],
    by = c("CHR", "POS"),
    all.x = TRUE
  )
  
  setorder(combined, ss)
  
  # ------------------ 새 longitudinal clumping 결과 ------------------
  #
  # 예:
  # B_DM2_SAIGE_cumulative_inci.SAIGE.txt_saige_clump_5e-8_0.1_1000.clumped
  #
  
  clump_prefix <- paste0(
    pheno160,
    "_SAIGE_cumulative_inci.SAIGE.txt_saige_clump_5e-8_0.1_1000"
  )
  
  clump_path <- file.path(
    CLUMP_ANNO_DIR,
    paste0(clump_prefix, ".clumped")
  )
  
  if (!file.exists(clump_path)) {
    message("Clump file not found: ", clump_path)
    return(NULL)
  }
  
  clump <- fread(clump_path)
  
  # ------------------ 새 longitudinal ANNOVAR 결과 ------------------
  
  anno_path <- file.path(
    CLUMP_ANNO_DIR,
    paste0(clump_prefix, ".hg19_multianno.txt")
  )
  
  if (!file.exists(anno_path)) {
    message("ANNOVAR file not found: ", anno_path)
    return(NULL)
  }
  
  anno <- fread(anno_path)
  
  list(
    pheno    = pheno70,
    pheno160 = pheno160,
    d1       = d1,
    combined = combined,
    clump    = clump,
    anno     = anno,
    CM       = CM,
    n.chr    = n.chr
  )
}

# ------------------------------------------------------------------------------
# 5. 데이터 불러오기
# ------------------------------------------------------------------------------

N_CORES <- 25

pheno_data_list <- mclapply(
  seq_len(nrow(pheno_map)),
  function(i) {
    load_all_data(
      pheno160 = pheno_map$pheno_160K[i],
      pheno70 = pheno_map$pheno_70K[i]
    )
  },
  mc.cores = N_CORES,
  mc.preschedule = FALSE
)

pheno_data_list <- pheno_data_list[
  sapply(
    pheno_data_list,
    function(x) is.list(x) && !inherits(x, "try-error")
  )
  ]

# panel label
letters_vec <- letters[seq_along(pheno_data_list)]

# ------------------------------------------------------------------------------
# 6. Manhattan plot 함수
# ------------------------------------------------------------------------------

plot.manhattan.fast <- function(
  d1,
  combined,
  clump_top20,
  CM,
  n.chr,
  title_text
) {
  p1 <- as.numeric(d1$P)
  p2 <- as.numeric(combined$P)
  
  y1 <- -log10(pmax(p1, 1e-300))
  y2 <- -log10(pmax(p2, 1e-300))
  
  genome_wide_line <- -log10(5e-8)
  
  ymax_signal <- max(c(y1, y2), na.rm = TRUE)
  
  ymax_annot <- if (!is.null(clump_top20) && nrow(clump_top20) > 0) {
    max(clump_top20$y_label, na.rm = TRUE)
  } else {
    -Inf
  }
  
  ymax_val <- max(
    ymax_signal * 1.15,
    ymax_annot + 12,
    genome_wide_line + 5,
    na.rm = TRUE
  )
  
  y_ticks <- pretty(c(0, ymax_val), n = 6)
  ymax_val <- max(y_ticks)
  
  plot(
    d1$ss,
    y1,
    type = "n",
    xlab = "Chromosome",
    ylab = expression(-log[10](p)),
    ylim = c(0, ymax_val),
    xaxt = "n",
    yaxt = "n",
    xaxs = "i",
    yaxs = "i",
    main = "",
    bty = "l"
  )
  
  par(xpd = FALSE)
  
  # KoGES 160K
  points(
    d1$ss,
    y1,
    col = "#0072B2",
    pch = 20,
    cex = 0.35
  )
  
  # KoGES 70K
  points(
    combined$ss,
    y2,
    col = "#999999",
    pch = 20,
    cex = 0.35
  )
  
  if (!is.null(clump_top20) && nrow(clump_top20) > 0) {
    
    # 160K lead SNP
    points(
      clump_top20$ss,
      clump_top20$y,
      pch = 4,
      col = "red",
      cex = 1,
      lwd = 1.5
    )
    
    # 70K에서도 매칭되는 동일 위치 SNP
    m <- match(clump_top20$ss, combined$ss)
    p_rep <- combined$P[m]
    
    if (any(!is.na(p_rep))) {
      points(
        clump_top20$ss[!is.na(p_rep)],
        -log10(pmax(as.numeric(p_rep[!is.na(p_rep)]), 1e-300)),
        pch = 1,
        col = "black",
        cex = 0.9
      )
    }
    
    # SNP와 gene label 연결
    segments(
      x0 = clump_top20$ss,
      y0 = clump_top20$y,
      x1 = clump_top20$x_annot,
      y1 = clump_top20$y_label - 1.5,
      col = "#4682B4",
      lwd = 0.5,
      lty = 3
    )
    
    x_limit <- max(d1$ss, na.rm = TRUE)
    
    for (i in seq_len(nrow(clump_top20))) {
      curr_x <- clump_top20$x_annot[i]
      
      adj_v <- 0.5
      if (curr_x < x_limit * 0.05) adj_v <- 0
      if (curr_x > x_limit * 0.95) adj_v <- 1
      
      text(
        x = curr_x,
        y = clump_top20$y_label[i],
        labels = clump_top20$Gene[i],
        cex = 0.7,
        font = 3,
        col = "#4682B4",
        adj = adj_v
      )
    }
  }
  
  CM1 <- c(0, CM)
  midpoints <- diff(CM1) / 2 + CM1[-(n.chr + 1)]
  
  axis(1, at = midpoints, labels = FALSE)
  mtext(1:n.chr, side = 1, line = 0.5, at = midpoints, cex = 0.8)
  
  axis(2, at = y_ticks, labels = y_ticks, cex.axis = 0.8, las = 1)
  
  segments(
    x0 = 0,
    y0 = genome_wide_line,
    x1 = max(d1$ss),
    y1 = genome_wide_line,
    col = "red",
    lty = 2
  )
  
  mtext(
    title_text,
    side = 3,
    line = 2,
    adj = 0,
    cex = 1.5,
    font = 2
  )
}

# ------------------------------------------------------------------------------
# 7. phenotype별 plot 생성
# ------------------------------------------------------------------------------

plot_one_parallel <- function(obj, idx) {
  
  if (!is.list(obj) || nrow(obj$d1) == 0 || nrow(obj$clump) == 0) {
    return(NULL)
  }
  
  gene_pos <- make_gene_pos(obj$anno)
  nearest_gene <- nearest_gene_factory(gene_pos)
  
  clump <- copy(obj$clump)
  setorder(clump, P)
  
  clump_top20 <- clump[1:min(20, .N), 1:5]
  
  # clumped file의 CHR, BP와 160K summary statistics의 CHR, POS 매칭
  clump_top20 <- merge(
    clump_top20,
    obj$d1[, .(CHR, POS, ss)],
    by.x = c("CHR", "BP"),
    by.y = c("CHR", "POS"),
    all.x = TRUE
  )
  
  clump_top20 <- clump_top20[!is.na(ss)]
  
  if (nrow(clump_top20) == 0) {
    return(NULL)
  }
  
  clump_top20[, Gene := mapply(nearest_gene, CHR, BP)]
  clump_top20[, y := -log10(pmax(as.numeric(P), 1e-300))]
  
  setorder(clump_top20, ss)
  
  # 동일 gene은 가장 유의한 SNP 하나만 label
  clump_top20 <- clump_top20[
    ,
    .SD[which.max(y)],
    by = Gene
    ]
  
  setorder(clump_top20, ss)
  
  # ------------------ label fan-out ------------------
  
  x_limit <- max(obj$d1$ss, na.rm = TRUE)
  x_mid <- x_limit / 2
  min_gap <- x_limit * 0.045
  
  new_x <- clump_top20$ss
  
  if (nrow(clump_top20) > 1) {
    for (iter in 1:200) {
      for (i in seq_len(nrow(clump_top20) - 1)) {
        distance <- new_x[i + 1] - new_x[i]
        
        if (distance < min_gap) {
          bias <- if (new_x[i] < x_mid) 0.6 else 0.4
          push_force <- (min_gap - distance) * 0.3
          
          new_x[i] <- new_x[i] - push_force * bias
          new_x[i + 1] <- new_x[i + 1] +
            push_force * (1 - bias)
        }
      }
      
      new_x <- pmax(
        pmin(new_x, x_limit * 0.97),
        x_limit * 0.03
      )
    }
  }
  
  clump_top20[, x_annot := new_x]
  clump_top20[, y_label := y + 10 + (seq_len(.N) %% 2) * 12]
  
  # ------------------ PNG 저장 ------------------
  
  prefix <- paste0(letters_vec[idx], ". ")
  
  out_file <- file.path(
    OUT_DIR,
    paste0(prefix, obj$pheno, "_supple.png")
  )
  
  png(
    filename = out_file,
    width = 5000,
    height = 2700,
    res = 300
  )
  
  par(
    mar = c(5, 5, 10, 15),
    family = "Arial"
  )
  
  trait_label <- ifelse(
    is.na(pheno_map2[obj$pheno]),
    obj$pheno,
    pheno_map2[obj$pheno]
  )
  
  title_text <- paste0(
    letters_vec[idx],
    ". ",
    trait_label
  )
  
  plot.manhattan.fast(
    d1 = obj$d1,
    combined = obj$combined,
    clump_top20 = clump_top20,
    CM = obj$CM,
    n.chr = obj$n.chr,
    title_text = title_text
  )
  
  par(xpd = TRUE)
  usr <- par("usr")
  
  legend(
    x = usr[2] + (usr[2] - usr[1]) * 0.02,
    y = usr[4] * 0.9,
    legend = c(
      "KoGES160K",
      "KoGES70K",
      "Lead SNP in 160K",
      "Same variant in KoGES 70K"
    ),
    col = c("#0072B2", "#999999", "red", "black"),
    pch = c(20, 20, 4, 1),
    pt.cex = c(1.5, 1.5, 1.2, 1.2),
    bty = "n",
    cex = 1.0
  )
  
  dev.off()
  
  message("Saved: ", out_file)
  
  invisible(out_file)
}

# ------------------------------------------------------------------------------
# 8. PNG 생성
# ------------------------------------------------------------------------------

plot_files <- mclapply(
  seq_along(pheno_data_list),
  function(i) plot_one_parallel(pheno_data_list[[i]], i),
  mc.cores = N_CORES,
  mc.preschedule = FALSE
)

# ------------------------------------------------------------------------------
# 9. PNG → PDF 변환
# ------------------------------------------------------------------------------

png_files <- list.files(
  OUT_DIR,
  pattern = "_supple\\.png$",
  full.names = TRUE,
  ignore.case = TRUE
)

png_files <- sort(png_files)

if (length(png_files) == 0) {
  stop("No PNG files found in: ", OUT_DIR)
}

pheno_names <- basename(png_files)
pheno_names <- sub(
  "_supple\\.png$",
  "",
  pheno_names,
  ignore.case = TRUE
)

start_num <- 5

fig_numbers <- seq(
  from = start_num,
  length.out = length(png_files)
)

pdf_dir <- file.path(OUT_DIR, "pdffile")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(png_files)) {
  
  output_path <- file.path(
    pdf_dir,
    paste0(
      "FigureS",
      fig_numbers[i],
      "_",
      pheno_names[i],
      ".pdf"
    )
  )
  
  img <- image_read(png_files[i])
  info <- image_info(img)
  
  pdf(
    file = output_path,
    width = info$width / 300,
    height = info$height / 300
  )
  
  grid::grid.raster(as.raster(img))
  dev.off()
  
  message("Saved: ", output_path)
}