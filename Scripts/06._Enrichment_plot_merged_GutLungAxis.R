rm(list = ls())

library(tidyverse)
library(viridis)
library(mutoss)

# =========================
# 0) 기본 경로
# =========================
BASE_DIR <- "C:/Github_repository/Metabolomics_Tutorial/Data"

# =========================
# 1) 같이 표시할 데이터 조합
# =========================
input_info <- tibble(
  STYPE = c("fecal_pos", "fecal_pos", "serum_pos"),
  WEEK  = c(1, 2, 2),
  Label = c("Fecal 1w", "Fecal 2w", "Serum 2w")
)

TOP_N_PER_GROUP <- 999
Q_CUTOFF <- 0.05


# =========================
# 2) 파일 읽기 함수
# =========================
read_pathway_file <- function(stype, week, label) {
  
  file <- file.path(
    BASE_DIR,
    stype,
    paste0("pathway_results_", week, "w.csv")
  )
  
  pw <- read.csv(file, check.names = FALSE)
  
  names(pw)[names(pw) == ""] <- "Pathway"
  
  pw <- pw %>%
    rename(
      TotalCmpd = any_of("Total Cmpd"),
      RawP = any_of("Raw p"),
      HolmAdjust = any_of("Holm adjust")
    )
  
  if (!"RawP" %in% names(pw)) {
    stop(paste0(file, " 파일에 Raw p 컬럼이 없습니다."))
  }
  
  if (!"Impact" %in% names(pw)) {
    stop(paste0(file, " 파일에 Impact 컬럼이 없습니다."))
  }
  
  pw <- pw %>%
    mutate(
      RawP = as.numeric(RawP),
      Impact = as.numeric(Impact),
      BKY = suppressWarnings(
        tryCatch(
          mutoss::p.adjust(RawP, method = "BKY"),
          error = function(e) stats::p.adjust(RawP, method = "BY")
        )
      ),
      minuslog10BKY = -log10(BKY),
      STYPE = stype,
      Week = paste0(week, "w"),
      Label = label
    ) %>%
    filter(is.finite(minuslog10BKY))
  
  return(pw)
}


# =========================
# 3) 여러 파일 병합
# =========================
pw_all <- pmap_dfr(
  input_info,
  function(STYPE, WEEK, Label) {
    read_pathway_file(STYPE, WEEK, Label)
  }
)

pw_all <- pw_all %>%
  mutate(Label = factor(Label, levels = input_info$Label))

write.csv(
  pw_all,
  "pathway_results_fecal1w_fecal2w_serum2w_merged_BKY.csv",
  row.names = FALSE
)


# =========================
# 4) 표시할 pathway 선정
#    각 Label별 q-value 상위 pathway를 합침
# =========================
selected_pathways <- pw_all %>%
  filter(BKY < Q_CUTOFF) %>%
  group_by(Label) %>%
  arrange(BKY, .by_group = TRUE) %>%
  slice_head(n = TOP_N_PER_GROUP) %>%
  ungroup() %>%
  pull(Pathway) %>%
  unique()


# pathway 순서: 전체에서 가장 작은 q-value 기준
pathway_order <- pw_all %>%
  filter(Pathway %in% selected_pathways) %>%
  group_by(Pathway) %>%
  summarise(best_q = min(BKY, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_q) %>%
  pull(Pathway)


pw_plot <- pw_all %>%
  filter(
    Pathway %in% selected_pathways,
    BKY < Q_CUTOFF
  ) %>%
  mutate(
    Pathway = factor(Pathway, levels = rev(pathway_order))
  )


# =========================
# 5) Figure
# =========================
fig_height <- max(4, 0.28 * length(selected_pathways) + 1.5)

p <- ggplot(
  pw_plot,
  aes(x = Label, y = Pathway)
) +
  geom_point(
    aes(size = Impact, color = minuslog10BKY),
    alpha = 0.95
  ) +
  scale_color_viridis_c(
    option = "plasma",
    direction = -1,
    name = expression(-Log[10](q-value))
  ) +
  scale_size_continuous(
    name = "Pathway Impact",
    range = c(2.5, 7)
  ) +
  labs(
    title = "Pathway enrichment analysis",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0, size = 13),
    
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 8),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  )

p

ggsave(
  "enrichment_fecal1w_fecal2w_serum2w.png",
  p,
  width = 7.5,
  height = fig_height,
  dpi = 500
)
