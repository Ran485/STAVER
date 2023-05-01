rm(list = ls())
# 先加载包
library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)

data <- read.csv("/Users/ranpeng/Desktop/DIA-QC/results/covid19-2/qc_orignal.csv")
# data$Group <- as.factor(data$Group)
# colnames(data)
# data[data == 0] = NA
df <- data
df$tQC_NA_NA_114 <- log2(df$tQC_NA_NA_114)
df$tQC_NA_NA_216 <- log2(df$tQC_NA_NA_216)

ggplot(df, aes(x = tQC_NA_NA_114, y = tQC_NA_NA_114)) +
  geom_point() +
  geom_smooth(method = lm, color = "black") +
  labs(
    title = "Miles per gallon \n according to the weight",
    x = "Weight (lb/1000)", y = "Miles/(US) gallon"
  ) +
  xlim(0, 15) +
  ylim(0, 15) +
  theme_classic()

p1 <- ggscatter(df,
  x = "tQC_NA_NA_114",
  y = "tQC_NA_NA_216",
  color = "#B2DF89",
  add = "reg.line"
) +
  # geom_point()+
  # geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs( # title= paste(colnames(gene_data)[2],colnames(gene_data)[3],sep = '_'),
    x = "PLT_log", y = "GMP"
  ) +
  # scale_color_brewer(palette="paired") +
  theme_minimal() +
  stat_cor(aes(color = "red", alpha = 0.3), method = "spearman", alpha = 1) + ### 添加相关性 pearson spearman
  theme( # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    # axis.line = element_line(colour = "black")
  ) +
  ggtitle("TMB correlation with signature")

p1 + theme(
  panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black")
)

ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  geom_smooth(method = lm, color = "black") +
  labs(
    title = "Miles per gallon \n according to the weight",
    x = "Weight (lb/1000)", y = "Miles/(US) gallon"
  ) +
  theme_classic()

p1 + theme_define2

theme_define1 <- theme(
  panel.border = element_blank(),
  # panel.grid.major = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black")
)

theme_define2 <- theme(
  axis.text.x = element_text(angle = 0, hjust = 0.6, colour = "black", family = "Arial", size = 16), # 设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
  axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", family = "Arial", size = 16),
  axis.text.y = element_text(family = "Arial", size = 16, face = "plain"), # 设置y轴刻度标签的字体簇，字体大小，字体样式为plain
  axis.title.y = element_text(family = "Arial", size = 20, face = "plain"), # 设置y轴标题的字体属性
  panel.border = element_blank(), axis.line = element_line(colour = "black", size = 0.8), # 去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
  legend.text = element_text(
    face = "plain", family = "Arial", colour = "black", # 设置图例的子标题的字体属性
    size = 12
  ),
  legend.title = element_text(
    face = "plain", family = "Arial", colour = "black", # 设置图例的总标题的字体属性
    size = 12
  ),
  panel.grid.major = element_blank(), # 不显示网格线
  panel.grid.minor = element_blank()
)
