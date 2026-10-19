library(data.table)
library(parallel)
library(magick)

# ==============================================================================
# KoGES 70K vs 160K Manhattan overlay plot
#
# 변경된 clumping / ANNOVAR 경로:
# /data6/OpenGWAS/KBA160K_longitudinal/Result/Clumping_BOLT/
#
# 파일명 예:
# RIN_B_ALBUMIN_bolt_clump_5e-8_0.1_1000.clumped
# RIN_B_ALBUMIN_bolt_clump_5e-8_0.1_1000.hg19_multianno.txt
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 경로 설정
# ------------------------------------------------------------------------------

GWAS160K_DIR <- "/data6/OpenGWAS/KBA160K/Result/BOLT-LMM"
GWAS70K_DIR  <- "/data6/OpenGWAS/Sumstats/KBA70K/Cont"

# 새로 생성한 longitudinal clumping + ANNOVAR 결과 폴더
CLUMP_ANNO_DIR <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/Clumping_BOLT"

OUT_DIR <- "/data4/joh/openkoges/figure/Overlaid_comparision_70Kvs160K/parallel_260907"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. 비교할 phenotype 목록
# ------------------------------------------------------------------------------

pheno_map <- data.table(
  pheno_160K = c(
    "B_ALBUMIN", "B_ALT", "B_AST", "B_BMI", "B_BUN", "B_CREATINE",
    "B_DBP", "B_GLU0", "B_HB", "B_HBA1C", "B_HDL",
    "B_HEIGHT", "B_HIP", "B_LDLclean", "B_PLAT", "B_PULSE",
    "B_R_GTP", "B_RBC_B", "B_SBP", "B_T_BIL",
    "B_TCHL", "B_TG", "B_WAIST", "B_WBC_B", "B_WEIGHT"
  ),
  pheno_70K = c(
    "ALB", "ALT", "AST", "BMI", "BUN", "CREATINE",
    "DBP", "Glucose", "HB", "HbA1C", "HDL",
    "HEIGHT", "HIP", "LDL", "PLT", "PULSE",
    "GGTP", "RBC", "SBP", "TBIL",
    "TC", "TG", "WAIST", "WBC", "WEIGHT"
  )
)

# 테스트 시 아래 phenotype 목록만 사용하려면 위 pheno_map 대신 사용
#
# pheno_map <- data.table(
#   pheno_160K = c("B_DBP", "B_PULSE", "B_SBP", "B_WAIST"),
#   pheno_70K = c("DBP", "PULSE", "SBP", "WAIST")
# )

# ------------------------------------------------------------------------------
# 2. phenotype 표시명
# ------------------------------------------------------------------------------

pheno_map2 <- c(
  RBC       = "Red blood cell count (RBC)",
  WBC       = "White blood cell count (WBC)",
  MCV       = "Mean corpuscular volume (MCV)",
  MCH       = "Mean corpuscular hemoglobin (MCH)",
  MCHC      = "Mean corpuscular hemoglobin concentration (MCHC)",
  PLT       = "Platelet count (Plt)",
  HCT       = "Hematocrit (HCT)",
  HB        = "Hemoglobin (Hb)",
  FEV1best  = "FEV1",
  FEV1pred  = "FEV1_%pred",
  FFRbest   = "FEV1/FVC",
  FVCbest   = "FVC",
  FVCpred   = "FVC_%pred",
  HDL       = "HDL-cholesterol (HDL-C)",
  LDL       = "LDL-cholesterol (LDL-C)",
  TG        = "Triglycerides (TG)",
  TC        = "Total cholesterol (TC)",
  APOA      = "Apolipoprotein A",
  APOB      = "Apolipoprotein B",
  VITD      = "Vit.D",
  IRON      = "Iron",
  FOLATE    = "Folate",
  ALB       = "Albumin (Alb)",
  URICACID  = "Uric acids (UA)",
  HOMOCYST  = "Homocysteine",
  HIP       = "Hip",
  BMI       = "Body mass index (BMI)",
  WHR       = "Waist-to-hip ratio (WHR)",
  WAIST     = "Waist",
  WEIGHT    = "Body weight (BW)",
  HEIGHT    = "Body height",
  AST       = "Aspartate transaminase (AST)",
  ALT       = "Alanine aminotransferase (ALT)",
  ALP       = "Alkaline phosphatase (ALP)",
  GGTP      = "Gamma glutamyl transferase (r_GTP)",
  TBIL      = "Total bilirubin (TBil)",
  BUN       = "Blood urea nitrogen (BUN)",
  eGFR      = "eGFR (CKD-EPI)",
  CYSTATINC = "Cystatin C",
  PHOSPHATE = "Phosphate",
  CREATINE  = "Serum creatinine (sCr)",
  HBA1C     = "HbA1C",
  INSUL     = "Insulin",
  GLU0      = "Glucose",
  CPEPTIDE  = "C-peptide",
  SBP       = "Systolic blood pressure (SBP)",
  DBP       = "Diastolic blood pressure (DBP)",
  PULSE     = "Pulse"
)

# ------------------------------------------------------------------------------
# 3. 유전자 annotation 보조 함수
# ------------------------------------------------------------------------------

clean_gene <- function(gene) {
  if (is.na(gene) || gene == ".") return(NA_character_)
  
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
    
    if (nrow(sub) == 0) return(NA_character_)
    
    dist <- pmin(abs(pos - sub$START), abs(pos - sub$END))
    
    sub$Gene[which.min(dist)]
  }
}

# ------------------------------------------------------------------------------
# 4. phenotype별 160K / 70K / clump / ANNOVAR 데이터 불러오기
# ------------------------------------------------------------------------------

load_all_data <- function(pheno160, pheno70) {
  
  cat("Loading:", pheno160, "/", pheno70, "\n")
  
  # ------------------ 160K GWAS summary statistics ------------------
  
  d1_path <- file.path(
    GWAS160K_DIR,
    paste0("RIN_", pheno160, "_BOLT_result")
  )
  
  if (!file.exists(d1_path)) {
    message("160K GWAS file not found: ", d1_path)
    return(NULL)
  }
  
  d1 <- fread(d1_path)
  
  if ("BP" %in% names(d1)) {
    setnames(d1, "BP", "POS")
  }
  
  if (!"P" %in% names(d1) && "P_BOLT_LMM_INF" %in% names(d1)) {
    setnames(d1, "P_BOLT_LMM_INF", "P")
  }
  
  d1[, P := as.numeric(P)]
  d1[, CHR := as.numeric(CHR)]
  d1[, POS := as.numeric(POS)]
  
  d1 <- d1[CHR %in% 1:22 & !is.na(POS) & !is.na(P)]
  
  setorder(d1, CHR, POS)
  
  # Manhattan plot을 위한 누적 SNP 순서
  d1[, ss := .I]
  
  chr_n <- table(d1$CHR)
  CM <- cumsum(chr_n)
  n.chr <- length(chr_n)
  
  # ------------------ 70K GWAS summary statistics ------------------
  
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
  
  d2 <- d2[CHR %in% 1:22 & !is.na(POS) & !is.na(P)]
  
  # 160K SNP 좌표 기준으로 70K P-value 매칭
  combined <- merge(
    d1[, .(CHR, POS, ss)],
    d2[, .(CHR, POS, P)],
    by = c("CHR", "POS"),
    all.x = TRUE
  )
  
  setorder(combined, ss)
  
  # ------------------ 신규 longitudinal clump 결과 ------------------
  
  clump_path <- file.path(
    CLUMP_ANNO_DIR,
    paste0(
      "RIN_", pheno160,
      "_bolt_clump_5e-8_0.1_1000.clumped"
    )
  )
  
  if (!file.exists(clump_path)) {
    message("Clump file not found: ", clump_path)
    return(NULL)
  }
  
  clump <- fread(clump_path)
  
  # ------------------ 신규 longitudinal ANNOVAR annotation ------------------
  
  anno_path <- file.path(
    CLUMP_ANNO_DIR,
    paste0(
      "RIN_", pheno160,
      "_bolt_clump_5e-8_0.1_1000.hg19_multianno.txt"
    )
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
# 5. 데이터 병렬 로딩
# ------------------------------------------------------------------------------

N_CORES <- 25

pheno_data_list <- mclapply(
  seq_len(nrow(pheno_map)),
  function(i) {
    load_all_data(
      pheno160 = pheno_map$pheno_160K[i],
      pheno70  = pheno_map$pheno_70K[i]
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
    
    # 70K에서 동일 좌표 SNP 표시
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
    
    # lead SNP와 gene label 연결선
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
# 7. phenotype별 plot 생성 함수
# ------------------------------------------------------------------------------

plot_one_parallel <- function(obj, idx) {
  
  if (!is.list(obj) || nrow(obj$d1) == 0 || nrow(obj$clump) == 0) {
    return(NULL)
  }
  
  # ANNOVAR gene 좌표 정보
  gene_pos <- make_gene_pos(obj$anno)
  nearest_gene <- nearest_gene_factory(gene_pos)
  
  clump <- copy(obj$clump)
  setorder(clump, P)
  
  # clumped 파일의 첫 5개 열이 CHR, SNP, BP, P 등을 포함한다는 기존 구조 유지
  clump_top20 <- clump[1:min(20, .N), 1:5]
  
  # clump SNP를 160K Manhattan 좌표에 연결
  clump_top20 <- merge(
    clump_top20,
    obj$d1[, .(CHR, POS, ss)],
    by.x = c("CHR", "BP"),
    by.y = c("CHR", "POS"),
    all.x = TRUE
  )
  
  # 160K Manhattan 좌표와 연결되지 않은 SNP 제거
  clump_top20 <- clump_top20[!is.na(ss)]
  
  if (nrow(clump_top20) == 0) {
    return(NULL)
  }
  
  # 가장 가까운 gene 이름
  clump_top20[, Gene := mapply(nearest_gene, CHR, BP)]
  
  # lead SNP P-value
  clump_top20[, y := -log10(pmax(as.numeric(P), 1e-300))]
  
  # 같은 gene에 여러 lead SNP가 있으면 가장 유의한 SNP 하나만 사용
  setorder(clump_top20, ss)
  clump_top20 <- clump_top20[
    ,
    .SD[which.max(y)],
    by = Gene
    ]
  setorder(clump_top20, ss)
  
  # ------------------ Gene label fan-out 위치 계산 ------------------
  
  x_limit <- max(obj$d1$ss, na.rm = TRUE)
  x_mid <- x_limit / 2
  min_gap <- x_limit * 0.045
  
  if (nrow(clump_top20) > 0) {
    
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
        
        # 양쪽 plot 경계에서 3% 여백 유지
        new_x <- pmax(
          pmin(new_x, x_limit * 0.97),
          x_limit * 0.03
        )
      }
    }
    
    clump_top20[, x_annot := new_x]
  }
  
  # 위아래 지그재그 gene label
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
# 8. PNG figure 생성
# ------------------------------------------------------------------------------

letters_vec <- letters

plot_files <- mclapply(
  seq_along(pheno_data_list),
  function(i) plot_one_parallel(pheno_data_list[[i]], i),
  mc.cores = N_CORES,
  mc.preschedule = FALSE
)

# ------------------------------------------------------------------------------
# 9. 생성된 PNG를 Supplementary PDF로 변환
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

# FigureS5부터 번호를 시작하려면 5 유지
start_num <- 5

fig_numbers <- seq(
  from = start_num,
  length.out = length(png_files)
)

pdf_dir <- file.path(OUT_DIR, "pdffile")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(png_files)) {
  
  new_filename <- paste0(
    "FigureS",
    fig_numbers[i],
    "_",
    pheno_names[i],
    ".pdf"
  )
  
  output_path <- file.path(pdf_dir, new_filename)
  
  img <- image_read(png_files[i])
  info <- image_info(img)
  
  pdf(
    file = output_path,
    width = info$width / 300,
    height = info$height / 300
  )
  
  grid::grid.raster(img)
  dev.off()
  
  message("Saved: ", output_path)
}