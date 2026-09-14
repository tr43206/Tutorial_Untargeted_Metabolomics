rm(list = ls())

library(tidyverse)
library(viridis)
library(mutoss)

STYPE = 'fecal_pos'
WEEK = '2'

setwd(paste0('C:/Github_repository/Metabolomics_Tutorial/Data',STYPE))

# 1) CSV 로드 (공백 있는 컬럼명 때문에 check.names=FALSE)
pw <- read.csv(paste0('pathway_results_',WEEK,'w.csv'), check.names = FALSE)

print(names(pw))

names(pw)[names(pw) == ""] <- "Pathway"


# 2) 컬럼 정리 + 지표 생성
pw <- pw %>%
  rename(TotalCmpd = `Total Cmpd`,
         RawP = `Raw p`,
         HolmAdjust = `Holm adjust`) %>%
  mutate(
    RawP = as.numeric(RawP),
    BKY = suppressWarnings(tryCatch(
      mutoss::p.adjust(RawP, method = "BKY"),
      error = function(e) stats::p.adjust(RawP, method = "BY")
    )),
    minuslog10BKY = -log10(BKY),
    
    FoldEnrichment = Impact
  ) %>%
  filter(is.finite(minuslog10BKY))

write.csv(pw, paste0('pathway_results_',WEEK,'w_BKY.csv'), row.names = FALSE)


# p-value < 0.05 전부 표시
pw_plot <- pw %>%
  filter(BKY < 0.05) %>%
  arrange(desc(minuslog10BKY)) %>%
  slice_head(n = 15) %>%
  mutate(Pathway = factor(Pathway, levels = rev(Pathway)))


# 4) Figure
p <- ggplot(pw_plot, aes(x = minuslog10BKY, y = Pathway)) +
  geom_point(aes(size = FoldEnrichment, color = minuslog10BKY)) +
  scale_color_viridis_c(option = "plasma", direction = -1,
                        name = expression(-Log[10](q-value))) +
  scale_size_continuous(name = "Pathway Impact", range = c(2.5, 5.5)) +
#  scale_y_discrete(position = "right") +
  labs(title = paste0('Fecal  (Positive ion mode)  -  Week ',WEEK),  #Edit
       x = expression(-Log[10](q-value)), y = NULL) +
#  guides(
#    size  = guide_legend(order = 2),
#    color = guide_colorbar(order = 1, barheight = unit(35, "mm"))
#  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0),  ## 왼쪽 정렬
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
#    axis.text.y.left  = element_blank(),
#    axis.ticks.y.left = element_blank(),
    
#    legend.position = c(0.77, 0.03),
#    legend.justification = c(0, 0),
#    legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.3),
#    legend.background = element_blank(),
#    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white", color = NA),
  )

p

ggsave(paste0('enrichment_',WEEK,'w.png'), p, width = 7, height = 4, dpi = 500)
