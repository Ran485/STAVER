rm(list = ls())
# 先加载包
library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)

data <- read.csv("~/DIA-QC/results/covid19-2/qc_orignal.csv")
colnames(data)
# data[data == 0] = NA
df <- data
df$Tech_Repeats_1 <- log2(df$Tech_Repeats_1)
df$Tech_Repeats_2 <- log2(df$Tech_Repeats_2)

ggplot(df, aes(x = Tech_Repeats_1, y = Tech_Repeats_2)) +
  geom_point() +
  geom_smooth(method = lm, color = "black") +
  labs(
    title = "Correlation for STAVER Tech.Repeats",
    x = "Tech. Repeats", y = "Tech. Repeats"
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic()

plot <- ggscatter(df,
  x = "Tech_Repeats_1",
  y = "Tech_Repeats_2",
  color = "#B2DF89",
  add = "reg.line"
) +
  # geom_point()+
  # geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(
    title = paste(colnames(gene_data)[2],
      colnames(gene_data)[3],
      sep = "_"
    ),
    x = "Tech_Repeats_1",
    y = "Tech_Repeats_2"
  ) +
  # scale_color_brewer(palette="paired") +
  theme_minimal() +
  stat_cor(aes(color = "red", alpha = 0.3), method = "spearman", alpha = 1) +
  ggtitle("TMB correlation with signature") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.6, colour = "black", family = "Arial", size = 16),
    axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", family = "Arial", size = 16),
    axis.text.y = element_text(family = "Arial", size = 16, face = "plain"),
    axis.title.y = element_text(family = "Arial", size = 20, face = "plain"),
    panel.border = element_blank(), axis.line = element_line(colour = "black", size = 0.8),
    legend.text = element_text(
      face = "plain", family = "Arial", colour = "black",
      size = 12
    ),
    legend.title = element_text(
      face = "plain", family = "Arial", colour = "black",
      size = 12
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
