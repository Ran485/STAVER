library(ggplot2)
library(scales)


data <- read.csv("~/results/plasma_data_QC_proteins.csv", row.names = 1)
data$LFQ <- data$LFQ_intensity
protein_hihglight <- c("C1R", "C8G", "F5 ")
data <- data[data$QC_proteins %in% protein_hihglight, ]
# Change color/shape by groups
# Remove confidence bands
p <- ggplot(data, aes(x = index, y = LFQ, color = QC_proteins)) +
  geom_point(shape = 1) + # 使用空心圆
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
  # xlim(10000, 1000000000000)+
  labs(
    title = "Analytical Variability",
    x = "Sample ID",
    y = "LFQ intensity"
  )

# Continuous colors
p + scale_color_brewer(palette = "Paired") +
  theme_classic() +
  scale_y_continuous(
    trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.6, colour = "black", family = "Arial", size = 0),
    axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", family = "Arial", size = 16),
    axis.text.y = element_text(family = "Arial", size = 16, face = "plain"),
    axis.title.y = element_text(family = "Arial", size = 20, face = "plain"),
    panel.border = element_blank(), axis.line = element_line(colour = "black", size = 0.8),
    legend.text = element_text(
      face = "plain", family = "Arial", colour = "black", size = 12
    ),
    legend.title = element_text(
      face = "plain", family = "Arial", colour = "black", size = 14
    ),
    panel.grid.major = element_blank(), # 不显示网格线
    panel.grid.minor = element_blank()
  )
