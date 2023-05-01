library(ggplot2)
library(ggbreak)

df2 <- data.frame(
  supp = rep(c("VC", "OJ"), each = 3),
  dose = rep(c("D0.5", "D1", "D2"), 2),
  len = c(6.8, 15, 33, 4.2, 10, 29.5)
)
head(df2)

df2 <- read.csv("/Users/ranpeng/Library/Mobile Documents/com~apple~CloudDocs/Desktop/CTC项目/processed_data/bar_plot.csv")
# Stacked barplot with multiple groups
ggplot(data = df2, aes(x = dose, y = len, fill = supp)) +
  geom_bar(stat = "identity")
# Use position=position_dodge()
ggplot(data = df2, aes(x = dose, y = len, fill = supp)) +
  geom_bar(stat = "identity", position = position_dodge())

ggplot(data = df2, aes(x = dose, y = len, fill = supp)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = len),
    vjust = 1.6, color = "white",
    position = position_dodge(0.9), size = 3.5
  ) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_break(c(1200, 1400), scales = 1) +
  theme_classic()
