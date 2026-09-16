###############################################################################
# Mfuzz pipeline (3 weeks input + 3 differential input)
# - Mfuzz clustering: FeatureID 기준 (기존 그대로)
# - 초록선(trajectories): Q_CUT & LFC_CUT 통과한 FeatureID만 그림
# - 우측 라벨: Interested 열이 비어있지 않은 것만 "Name"으로 라벨 (Name 비면 스킵)
###############################################################################

# install.packages("BiocManager")
# BiocManager::install(c("Biobase","Mfuzz"))
# install.packages(c("readr","dplyr","stringr"))
# install.packages("readxl")  # 필요

rm(list = ls())

suppressPackageStartupMessages({
  library(Biobase)
  library(Mfuzz)
  library(readr)
  library(dplyr)
  library(stringr)
  library(readxl)
})

setwd('C:/Github_repository/Metabolomics_Tutorial/Data')

#####################################################################################
# -----------------------
# 0) 사용자 설정
# -----------------------

FILES_WEEK <- c(
  Week0 = 'fecal_neg_0w',
  Week1 = 'fecal_neg_1w',
  Week2 = 'fecal_neg_2w'
)
WEEK_FILE <- 'dataset_filtered_scaled_merged.xlsx'

## (선택1) 여러 엑셀 파일에서 데이터를 가져올 때
#FILES_DIFF <- c(
#  Week0 = 'MA_stats+annot_0w.xlsx',
#  Week1 = 'MA_stats+annot_1w.xlsx',
#  Week2 = 'MA_stats+annot_2w.xlsx'
#)

## (선택2) 하나의 엑셀 파일의 여러 시트에서 데이터를 가져올 때
FILES_DIFF <- c(
  Week0 = 'fecal_neg_0w',  #Edit
  Week1 = 'fecal_neg_1w',  #Edit
  Week2 = 'fecal_neg_2w'   #Edit
)
DIFF_FILE <- 'stats_merged.xlsx'


Q_CUT   <- 0.05  #Edit: default 0.05
LFC_CUT <- 0  #Edit: default 1
SIG_MODE <- "union"  # "union" or "intersection"

C <- 3  #Edit: 클러스터 수 조정 (일반적으로 3 or 4)
SEED <- 42

OUTDIR <- "mfuzz_out"
dir.create(OUTDIR, showWarnings = FALSE)

TOP_N_NAME <- 15

#####################################################################################
# -----------------------
# 1) I/O 유틸
# -----------------------

read_table_auto <- function(path, sheet = 1){
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx","xls")) return(readxl::read_excel(path, sheet = sheet))
  if (ext %in% c("tsv","txt"))  return(readr::read_tsv(path, show_col_types = FALSE))
  if (ext %in% c("csv"))       return(readr::read_csv(path, show_col_types = FALSE))
  return(readr::read_csv(path, show_col_types = FALSE))
}

set_rownames_from_first_col <- function(df){
  df <- as.data.frame(df)
  if (ncol(df) < 2) stop("입력 파일 컬럼이 너무 적습니다.")
  rn <- df[[1]] |> as.character() |> stringr::str_trim()
  keep <- !is.na(rn) & rn != ""
  df <- df[keep, , drop = FALSE]
  rn <- rn[keep]
  rn <- make.unique(rn)
  df <- df[, -1, drop = FALSE]
  rownames(df) <- rn
  df
}

find_qcol <- function(df){
  nms <- colnames(df)
  cand <- nms[str_detect(tolower(nms), "bky")]  #Edit: 다른 fdr method를 쓰는 경우 변경
  if (length(cand) == 0) stop("q-val/FDR 컬럼을 찾지 못했습니다.")
  cand[1]
}

find_lfccol <- function(df){
  nms <- colnames(df)
  cand <- nms[str_detect(tolower(nms), "lfc|log2fc")]
  if (length(cand) == 0) stop("lfc/logFC 컬럼을 찾지 못했습니다. (lfc, log2FC 등)")
  cand[1]
}

find_namecol <- function(df){
  nms <- colnames(df)
  cand <- nms[str_detect(tolower(nms), "^name$|metabolite|compound|annotation")]
  if (length(cand) == 0) return(NA_character_)
  cand[1]
}

find_interestedcol <- function(df){
  nms <- colnames(df)
  cand <- nms[str_detect(tolower(nms), "interested")]
  if (length(cand) == 0) return(NA_character_)
  cand[1]
}

#####################################################################################
# -----------------------
# 2) differential 3개에서 유의 FeatureID set + Name/Interested 매핑 생성
# -----------------------

sig_lists <- list()

# Name, Interested 매핑은 "diff 파일 전체에서" 만들고, 나중에 sig_set에 대해 인덱싱
name_map_all <- character(0)        # names = FeatureID, value = Name
interested_all <- character(0)      # names = FeatureID, value = Interested (비어있지 않으면 저장)

for (w in names(FILES_DIFF)){
  ## (선택1) 여러 엑셀 파일에서 데이터를 가져올 때
#  dfa <- read_table_auto(FILES_DIFF[[w]])
  ## (선택2) 하나의 엑셀 파일의 여러 시트에서 데이터를 가져올 때
  dfa <- readxl::read_excel(DIFF_FILE, sheet = FILES_DIFF[[w]], col_types = 'text')
  dfa <- set_rownames_from_first_col(dfa)  # rownames = FeatureID
  
  # 유의 set
  qcol  <- find_qcol(dfa)
  lfcol <- find_lfccol(dfa)
  dfa[[qcol]]  <- suppressWarnings(as.numeric(dfa[[qcol]]))
  dfa[[lfcol]] <- suppressWarnings(as.numeric(dfa[[lfcol]]))
  
  sig <- dfa %>%
    filter(!is.na(.data[[qcol]]), !is.na(.data[[lfcol]])) %>%
    filter(.data[[qcol]] < Q_CUT, abs(.data[[lfcol]]) > LFC_CUT)
  
  sig_lists[[w]] <- rownames(sig)
  cat(sprintf("[DIFF %s] significant = %d\n", w, length(sig_lists[[w]])))
  
  # Name 매핑(첫 등장 우선)
  ncol <- find_namecol(dfa)
  if (!is.na(ncol)){
    nm <- as.character(dfa[[ncol]])
    nm <- stringr::str_trim(nm)
    nm[nm == ""] <- NA_character_
    names(nm) <- rownames(dfa)
    
    new_ids <- setdiff(names(nm)[!is.na(nm)], names(name_map_all))
    if (length(new_ids) > 0){
      name_map_all[new_ids] <- nm[new_ids]
    }
  }
  
  # Interested 매핑(비어있지 않은 것만, 첫 등장 우선)
  icol <- find_interestedcol(dfa)
  if (!is.na(icol)){
    iv <- as.character(dfa[[icol]])
    iv <- stringr::str_trim(iv)
    iv[iv == ""] <- NA_character_
    names(iv) <- rownames(dfa)
    
    new_ids2 <- setdiff(names(iv)[!is.na(iv)], names(interested_all))
    if (length(new_ids2) > 0){
      interested_all[new_ids2] <- iv[new_ids2]
    }
  }
}

if (SIG_MODE == "union"){
  sig_set <- Reduce(union, sig_lists)
} else if (SIG_MODE == "intersection"){
  sig_set <- Reduce(intersect, sig_lists)
} else {
  stop("SIG_MODE는 'union' 또는 'intersection'만 허용")
}
cat(sprintf("Significant FeatureID set (%s) = %d\n", SIG_MODE, length(sig_set)))


# -----------------------
# 수동 제외 FeatureID (나중에 결과 피규어 보고 중복 Name 있으면 여기서 제거)
# -----------------------
REMOVE_FEATURES <- c(
#  'FeatureID_삭제할값1',
  ## Fecal_NEG
  '188.03554-5.585',  ## Kynurenic acid
  '188.03554-5.821',  ## Kynurenic acid (Intensity 큰놈)
  '202.05126-17.793',  ## Indole-3-pyruvic acid
  ## Fecal_POS
  '190.0499-5.9',  ## Kynurenic acid (Intensity 큰놈)
  '190.04982-8.565',  ## Kynurenic acid
  '190.04988-5.63',  ## Kynurenic acid
  '190.0499-6.464',  ## Kynurenic acid
  '190.0498-7.785',  ## Kynurenic acid
  '206.08109-7.276',  ## Indole-3-lactic acid
  '146.11763-1.456',  ## Acetylcholine
  '192.06553-7.642',  ## 5-Hydroxyindole-3-acetic acid
  '192.06546-6.725',  ## 5-Hydroxyindole-3-acetic acid
  '221.09201-1.691'  ## 5-Hydroxy-DL-tryptophan
)

sig_set <- setdiff(sig_set, REMOVE_FEATURES)

cat(sprintf("After manual removal = %d\n", length(sig_set)))


##**##
# 우측 라벨링용: Interested가 비어있지 않은 실제 FeatureID
interested_feature_set <- intersect(sig_set, names(interested_all))
cat(sprintf("Label-eligible FeatureIDs (Interested non-empty) = %d\n", length(interested_feature_set)))
##**##

##**##
# -----------------------
# Name -> FeatureID 역매핑
# 같은 Name이 여러 FeatureID에 걸릴 수 있으므로 "첫 등장 1개"만 사용
# -----------------------
name_to_feature <- names(name_map_all)
names(name_to_feature) <- unname(name_map_all)

# 중복 Name이면 첫 번째만 유지
name_to_feature <- name_to_feature[!duplicated(names(name_to_feature))]
##**##


# -----------------------
# 수동 Name 변경
# -----------------------
RENAME_FEATURES <- c(
#  'FeatureID_1' = '표기할_Name_1',
  ## Fecal_NEG
  '173.00927-1.126' = 'Aconitic acid',
  '173.00929-1.445' = 'Aconitic acid',
  '202.05128-16.408' = 'Indole-3-pyruvic acid',
  '219.07772-2.061' = '5-Hydroxy-tryptophan',
  '524.27889-13.665' = 'LysoPE(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0)',
  '554.34702-15.063' = 'LysoPC',
  ## Fecal_POS
  '205.09707-5.741' = 'Tryptophan',
  '206.08106-8.623' = 'Indole-3-lactic acid',
  '466.29327-14.068' = 'LysoPC(14:1(9Z)/0:0)',
  '221.09201-2.066' = '5-Hydroxy-tryptophan'
)

name_map_all[names(RENAME_FEATURES)] <- RENAME_FEATURES


#####################################################################################

week_mats <- list()

for (w in names(FILES_WEEK)){

  ## (선택1) Control + VNAM 합친 데이터 사용  
#  dfw <- read_table_auto(FILES_WEEK[[w]])
#  dfw <- as.data.frame(dfw)
  
  # 첫 컬럼 = FeatureID
#  feature_ids <- dfw[[1]] |> as.character() |> stringr::str_trim()
#  dfw <- dfw[, -1, drop = FALSE]
#  rownames(dfw) <- make.unique(feature_ids)
  
  ## (선택2) VNAM 데이터만 사용
#  dfw <- read_table_auto(FILES_WEEK[[w]])
  dfw <- read_table_auto(WEEK_FILE, sheet = FILES_WEEK[[w]])
  dfw <- as.data.frame(dfw)
  
  # 첫 컬럼 = FeatureID
  feature_ids <- dfw[[1]] |> as.character() |> stringr::str_trim()
  
  # 그룹 정보 (Label row)
  group_row <- as.character(dfw[1, -1])
  
  # 데이터 부분만 분리
  dfw <- dfw[-1, -1, drop = FALSE]
  
  # rownames 설정
  rownames(dfw) <- make.unique(feature_ids[-1])
  
  # 숫자화
  dfw[] <- lapply(dfw, function(x){
    suppressWarnings(as.numeric(x))
  })
  
  # -----------------------
  # VNAM 샘플만 선택
  # -----------------------
  keep_cols <- group_row == "VNAM"
  dfw <- dfw[, keep_cols, drop = FALSE]
  
  # QC 제거
  dfw <- dfw[, !grepl("QC", colnames(dfw), ignore.case = TRUE), drop = FALSE]
  
  # abundance matrix
  mat_abund <- as.matrix(dfw)
  
  # NA/Inf 제거
  ok_feat <- apply(mat_abund, 1, function(v) all(is.finite(v)))
  mat_abund <- mat_abund[ok_feat, , drop = FALSE]
  
  week_mats[[w]] <- mat_abund
  
  cat(sprintf("[WEEK %s] features=%d, samples=%d\n",
              w, nrow(mat_abund), ncol(mat_abund)))
  
  
  # 숫자화
#  dfw[] <- lapply(dfw, function(x){
#    if (is.character(x)) suppressWarnings(as.numeric(x)) else x
#  })
  
  # -----------------------
  # 1) QC 샘플 제거
  # -----------------------
#  keep_cols <- !grepl("QC", colnames(dfw), ignore.case = TRUE)
#  dfw <- dfw[, keep_cols, drop = FALSE]
  
  # -----------------------
  # 2) autoscaling (feature-wise z-score)
  # -----------------------
#  dfw_scaled <- t(scale(t(as.matrix(dfw)), center = TRUE, scale = TRUE))
  
  ##**##
#  mat_rel <- dfw_scaled
  ##**##
  
  # 혹시 남은 NaN/NA 제거
#  ok_feat <- apply(mat_rel, 1, function(v) all(is.finite(v)))
#  mat_rel <- mat_rel[ok_feat, , drop = FALSE]
  
#  week_mats[[w]] <- mat_rel
#  cat(sprintf("[WEEK %s] features=%d, samples=%d\n",
#              w, nrow(mat_rel), ncol(mat_rel)))
}

# -----------------------
# 4) 공통 feature + 유의 feature만 유지
# -----------------------
common_features <- Reduce(intersect, lapply(week_mats, rownames))
keep <- intersect(common_features, sig_set)

cat(sprintf("Common features = %d\n", length(common_features)))
cat(sprintf("Keep (common ∩ significant) = %d\n", length(keep)))

if (length(keep) < 5)
  stop("남는 metabolite가 너무 적습니다.")

# Week별 평균으로 timepoint matrix 구성
WEEKS_ORDER <- names(FILES_WEEK)

mat <- sapply(WEEKS_ORDER, function(w){
  rowMeans(week_mats[[w]][keep, , drop = FALSE], na.rm = TRUE)
})

mat <- as.matrix(mat)
rownames(mat) <- keep
colnames(mat) <- WEEKS_ORDER

#####################################################################################
# -----------------------
# 4) Mfuzz 실행
# -----------------------

eset <- ExpressionSet(assayData = mat)
eset <- standardise(eset)

# 안전장치: NA/NaN 있는 feature 제거
x <- exprs(eset)
bad <- which(apply(x, 1, function(v) any(!is.finite(v))))
if (length(bad) > 0) eset <- eset[-bad, ]

m <- mestimate(eset)
cat(sprintf("Estimated m = %.3f\n", m))

set.seed(SEED)
cl <- mfuzz(eset, c = C, m = m)

assign_df <- data.frame(
  FeatureID = featureNames(eset),
  Cluster   = cl$cluster,
  cl$membership,
  check.names = FALSE
)
write.csv(assign_df, file.path(OUTDIR, "Mfuzz_cluster_membership.csv"), row.names = FALSE)

##**##
# 클러스터별로 오른쪽 라벨에 실제 표시된 metabolite 이름 저장
cluster_label_names <- vector("list", C)
names(cluster_label_names) <- paste0("Cluster", seq_len(C))
##**##

#####################################################################################
#####################################################################################

# -----------------------
# 5) 플롯 함수만 교체 (그라데이션 + 박스 분리 + 라벨 위정렬 + 간격)
# -----------------------

plot_cluster_panel <- function(k, min_mem = 0.5){
  
  mfuzz.plot2(
    eset,
    cl = cl,
    mfrow = c(1, 1),
    single = k,
#    colo = "fancy",
    x11 = FALSE,
    time.labels = colnames(exprs(eset)),
    min.mem = min_mem,
    xlab = "",
    ylab = "Metabolite changes",
    cex.main = 1.2
  )
  
  # center line 추가
  lines(
    x = seq_len(ncol(exprs(eset))),
    y = cl$centers[k, ],
    lwd = 3,
    col = "black"
  )
  
  # 오른쪽 라벨용 이름 저장
  idx <- which(cl$cluster == k)
  feat_ids <- featureNames(eset)[idx]
  memk <- cl$membership[idx, k]
  ord <- order(memk, decreasing = TRUE)
  
  ##**##
  feat_order <- feat_ids[ord]
  nm_vec <- name_map_all[feat_order]
  
  # FeatureID 기준으로 Interested 있는 것만 유지
  ok <- feat_order %in% interested_feature_set
  ok <- ok & !is.na(nm_vec) & nm_vec != ""
  
  lab_df <- data.frame(
    FeatureID = feat_order[ok],
    Name = unname(nm_vec[ok]),
    stringsAsFactors = FALSE
  )
  
  # 같은 Name이 여러 FeatureID에 있더라도, 실제 Interested 있는 FeatureID만 남은 상태에서
  # cluster 내 membership 높은 것 1개만 유지
  lab_df <- lab_df[!duplicated(lab_df$Name), , drop = FALSE]
  
  if (nrow(lab_df) > TOP_N_NAME) {
    lab_df <- lab_df[1:TOP_N_NAME, , drop = FALSE]
  }
  
  cluster_label_names[[k]] <<- lab_df$Name
  ##**##
}

#####################################################################################

for(k in 1:C){
  png(file.path(OUTDIR, sprintf("Mfuzz_Cluster_%d.png", k)),
      width = 3200, height = 2200, res = 500)
  
  par(mar = c(4, 5, 3, 2), oma = c(0, 0, 0, 0))
  
  plot_cluster_panel(k, min_mem = 0)  #Edit: min_mem 부분
  
  dev.off()
}

#####################################################################################
## 검토 (아래 3개의 결과가 같아야함)

length(featureNames(eset))

sum(sapply(1:C, function(k) sum(cl$cluster == k)))

sum(table(cl$cluster))

#####################################################################################
#####################################################################################

##**##
get_cluster_lab_df <- function(k){
  idx <- which(cl$cluster == k)
  feat_ids <- featureNames(eset)[idx]
  memk <- cl$membership[idx, k]
  ord <- order(memk, decreasing = TRUE)
  
  feat_order <- feat_ids[ord]
  nm_vec <- name_map_all[feat_order]
  
  # FeatureID 기준으로 Interested 있는 것만
  ok <- feat_order %in% interested_feature_set
  ok <- ok & !is.na(nm_vec) & nm_vec != ""
  
  out <- data.frame(
    FeatureID = feat_order[ok],
    Name = unname(nm_vec[ok]),
    stringsAsFactors = FALSE
  )
  
  # 같은 Name 중복 제거: membership 높은 순으로 이미 정렬되어 있으므로 첫 번째만 유지
  out <- out[!duplicated(out$Name), , drop = FALSE]
  
  if (nrow(out) > TOP_N_NAME){
    out <- out[1:TOP_N_NAME, , drop = FALSE]
  }
  
  out
}
##**##

##**##
plot_right_name_panel <- function(k){
  lab_df <- get_cluster_lab_df(k)
  
  if (nrow(lab_df) == 0){
    plot.new()
    title(main = paste0("Cluster", k))
    text(0.5, 0.5, "No labeled metabolites")
    return(invisible(NULL))
  }
  
  feat_ids  <- lab_df$FeatureID
  lab_names <- lab_df$Name
  
  sub <- exprs(eset)[feat_ids, , drop = FALSE]
  x <- seq_len(ncol(sub))
  
  yr <- range(sub, finite = TRUE)
  ypad <- diff(yr) * 0.15
  if (!is.finite(ypad) || ypad == 0) ypad <- 0.3
  yr <- c(yr[1] - ypad, yr[2] + ypad)
  
  cols <- grDevices::hcl.colors(nrow(sub), palette = "Dark 3")
  
  plot(x, sub[1, ],
       type = "n",
       xlim = c(1, ncol(sub)),
       ylim = yr,
       xaxt = "n",
       xlab = "",
       ylab = "Metabolite changes",
       main = paste0("Cluster", k),
       bty = "o")
  
  axis(1, at = x, labels = colnames(sub))
  abline(h = 0, lty = 2, col = "grey70")
  grid(nx = NA, ny = NULL, col = "grey90", lty = 1)
  
  for(i in seq_len(nrow(sub))){
    lines(x, sub[i, ],
          type = "b",
          pch = 19,
          cex = 1.4,
          lwd = 3.1,
          col = cols[i])
  }
  
  x_leg <- par("usr")[2] + 0.05 * diff(par("usr")[1:2])
  
  legend(x = x_leg,
         y = mean(par("usr")[3:4]),
         legend = lab_names,
         col = cols,
         lty = 1,
         lwd = 1.4,
         pch = 16,
         pt.cex = 1.2,
         bty = "n",
         cex = 0.8,
         title = "Metabolites",
         xpd = NA,
         xjust = 0,
         yjust = 0.5)
}
##**##

plot_cluster_with_right_panel <- function(k, min_mem = 0){
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.6, 1))
  
  par(mar = c(4, 5, 3, 2))
  mfuzz.plot2(
    eset,
    cl = cl,
    mfrow = c(1, 1),
    single = k,
    #   colo = "fancy",
    x11 = FALSE,
    time.labels = colnames(exprs(eset)),
    min.mem = min_mem,
    xlab = "",
    ylab = "Metabolite changes",
    cex.main = 1.2
  )
  
  lines(
    x = seq_len(ncol(exprs(eset))),
    y = cl$centers[k, ],
    lwd = 3,
    col = "black"
  )
  
  par(mar = c(3, 4.5, 3, 22.5))  #Edit: 아래, 왼, 위, 오른
  plot_right_name_panel(k)
  
  layout(1)
}

#####################################################################################

for(k in 1:C){
  png(file.path(OUTDIR, sprintf("Mfuzz_Cluster_%d_with_names.png", k)),
      width = 4200, height = 3200, res = 500)
  
  plot_cluster_with_right_panel(k, min_mem = 0)
  
  dev.off()
}

#####################################################################################
