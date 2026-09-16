rm(list = ls())

library(readxl)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(mutoss)

STYPE = 'serum_pos'
WEEK = '2'
ABS_Log2FC = 1
q_val = 0.05

setwd('C:/Github_repository/Metabolomics_Tutorial/Data')

df <- read_excel('stats_merged.xlsx', sheet = paste0(STYPE,'_',WEEK,'w'))
head(df)

#####################################################################################

df$q_val <- as.numeric(df$`q-val_BKY`)  #Edit: p-val or q-val 실제 열이름


df$diffexpressed <- 'NO'
df$diffexpressed[df$lfc > ABS_Log2FC & df$q_val < q_val] <- 'UP'
df$diffexpressed[df$lfc < -ABS_Log2FC & df$q_val < q_val] <- 'DOWN'
head(df[order(df$q_val) & df$diffexpressed == 'DOWN', ])

df_qfilter <- df[df$q_val < q_val,]


targets <- c("Tryptophan metabolism", "Interested")
sig_idx <- df$q_val < q_val & abs(df$lfc) > ABS_Log2FC

## 기본값: 전부 NA
df$delabel <- NA_character_
## (선택1) Interested도 채워져 있고, 유의성 기준도 만족하는 경우만 annotation
idx_label <- sig_idx & (df$Interested %in% targets)
## (선택2) 유의성 기준 만족하는 경우만 annotation
#idx_label <- sig_idx

#df$delabel[idx_label] <- df$Name[idx_label]
df$delabel[idx_label] <- ifelse(
  !is.na(df$ReName[idx_label]) & df$ReName[idx_label] != "",
  df$ReName[idx_label],
  df$Name[idx_label]
)


df$color_group <- "Unannotated metabolites"  # 기본값: 전부 회색

# 1) 비유의: Annotated(o)만 색(짙은 회색)
df$color_group[!sig_idx & df$Annotated == "o"] <- "Annotated metabolites"

# 2) 유의: (a) Annotated도 색, (b) targets는 각 색 (targets가 Annotated보다 우선)
df$color_group[sig_idx & df$Annotated == "o"] <- "Annotated metabolites"

## (선택1) targets 이름으로 그룹명 지정
df$color_group[sig_idx & (df$Interested %in% targets)] <- df$Interested[sig_idx & (df$Interested %in% targets)]
## (선택2) 유의한 차이 (lfc, p-val)가 나는 것만 지정
#df$color_group[sig_idx] <- "Differentially abundant metabolites"

df$color_group <- factor(df$color_group,
                         levels = c("Unannotated metabolites", "Annotated metabolites", "Tryptophan metabolism", "SCFAs", "Interested"#,
#                                    "Differentially abundant metabolites"
                                    ))


ymin <- floor(min(-log10(df$q_val), na.rm = TRUE))
ymax <- ceiling(max(-log10(df$q_val), na.rm = TRUE))
xmin <- floor(min(df$lfc, na.rm = TRUE))
xmax <- ceiling(max(df$lfc, na.rm = TRUE))

######################################################################################

theme_set(theme_classic())

volcano <- 
  ggplot(data = df, aes(x = lfc, y = -log10(q_val), col = color_group, label = delabel)) +  #Edit: -log10() p_val로 할지 q_val로 할지 결정해야됨
#  geom_vline(xintercept = c(-ABS_Log2FC, ABS_Log2FC), col = "red", linetype = 'dashed', linewidth = 1.0) +
#  geom_hline(yintercept = -log10(q_val), col = "red", linetype = 'dashed', linewidth = 1.0) +
  
  # 1) Others(연한 회색)
  geom_point(data = subset(df, color_group == "Unannotated metabolites"), size = 2) +
  # 2) Annotated(짙은 회색)
  geom_point(data = subset(df, color_group == "Annotated metabolites"), size = 2) +
  # 3) 타겟 컬러
  geom_point(data = subset(df, !(color_group %in% c("Unannotated metabolites", "Annotated metabolites"))), size = 3) +
  
  scale_color_manual(values = c("Unannotated metabolites"  = "grey80",  ## 기준 안맞거나, Interested가 아닌 것들
                                "Annotated metabolites" = "grey70",
                                "Tryptophan metabolism" = "royalblue",
#                                "Glycerophospholipid metabolism" = "darkorange",
#                                "TCA cycle"    = "forestgreen",
                                "SCFAs"        = "purple",
                                "Interested"  = "firebrick"#,
#                                "Differentially abundant metabolites" = "#00A19C"
                                ),
                     breaks = c("Tryptophan metabolism", "Glycerophospholipid metabolism", "TCA cycle", "SCFAs", "Interested",
#                                "Differentially abundant metabolites",
                                "Annotated metabolites", "Unannotated metabolites")) +
  labs(color = NULL,
       x = expression("log"[2]*"FC (VNAM/Control)"),  #Edit
       y = expression("-log"[10]*" (q-value)")) +  #Edit: p-val or q-val
  
  coord_cartesian(xlim = c(xmin-(xmax+xmin)-1.0, xmax+1.0),
                  ylim = c(1.0, ymax+0.2)) +
  
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = seq(0, ceiling(ymax), by = 0.5),
                     expand = expansion(mult = c(0, 0.05))) +
  
  ggtitle(paste0('Serum  (Positive ion mode)  -  Week ',WEEK)) +  #Edit: ggtitle 수정
  
  geom_text_repel(max.overlaps = Inf, fontface = "bold",
                  ## Annotation label끼리 겹칠 땐 아래의 세팅값들 조정
                  size = 7,  #Edit: annotation label 글자크기
                  point.size = 10,  #Edit: annotation label끼리의 거리
                  box.padding = 0.5,
                  point.padding = 0.4,
                  segment.size = 0.7,  #Edit: 선 두께
                  force = 5,  #Edit: annotation label을 점에서 밀어내는 힘
                  force_pull = 0.1  #Edit: 점 쪽으로 당기는 힘
                  ) +
  geom_vline(xintercept = c(-ABS_Log2FC, ABS_Log2FC), col = "red", linetype = 'dashed', linewidth = 1.0) +
  geom_hline(yintercept = -log10(q_val), col = "red", linetype = 'dashed', linewidth = 1.0) +
  
  theme(legend.position = "right",
        legend.text = element_text(size = 20),   # legend 글자크기
        legend.key.size = unit(1.5, "cm"),       # legend 점/박스 간격 크기)
        axis.title.x = element_text(size = 22, face = "bold", color = "black"),
        axis.title.y = element_text(size = 22, face = "bold", color = "black"),
        
        axis.text.x = element_text(size = 20, color = "black"),  # xticks 글자크기
        axis.text.y = element_text(size = 20, color = "black"),  # yticks 글자크기
        
        plot.title = element_text(size = 22, face = "bold"))

print(volcano)

ggsave(paste0('volcano_',WEEK,'w_Trp.png'), width = 12, height = 9, dpi = 500)

#####################################################################################
