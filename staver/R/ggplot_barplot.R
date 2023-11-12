library("ggpubr")
library(ggplot2)
library(ggbreak)

# -------------  test data  ----------------
# Calculate the z-score of the mpg data
dfm$mpg_z <- (dfm$mpg - mean(dfm$mpg)) / sd(dfm$mpg)
dfm$mpg_grp <- factor(ifelse(dfm$mpg_z < 0, "low", "high"),
    levels = c("low", "high")
)
# Inspect the data
head(dfm[, c("name", "wt", "mpg", "mpg_z", "mpg_grp", "cyl")])

# ------------- load data ----------------
df <- read.csv("~/DIANN_Data/merge_data/merge_barplot.csv")
colnames(df)


ggbarplot(df,
    x = "Gene", y = "logFC_IPX",
    fill = "cluster", # change fill color by mpg_level
    color = "white", # Set bar border colors to white
    palette = "jco", # jco journal color palett. see ?ggpar
    sort.val = "asc", # Sort the value in ascending order
    sort.by.groups = FALSE, # Don't sort inside each group
    x.text.angle = 90, # Rotate vertically x axis texts
    ylab = "MPG z-score",
    xlab = FALSE,
    rotate = TRUE,
    legend.title = "MPG Group"
)


df2 <- read.csv("~/processed_data/bar_plot.csv")
ggplot(data = df2, aes(x = dose, y = len, fill = supp)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = len),
        vjust = 1.6, color = "white",
        position = position_dodge(0.9), size = 3.5
    ) +
    scale_fill_brewer(palette = "Paired") +
    scale_y_break(c(1200, 1400), scales = 1) +
    theme_classic()
