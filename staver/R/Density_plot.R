library("ggplot2")
library(plyr)

df <- read.csv("~/STAVER-nc-revised/Revised-Tables/Table_metadata/HEK-293T_cv_compare.csv")
colnames(df)
df$CVs <- as.numeric(df$CVs)
df$CVs <- df$CVs * 100
mu <- ddply(df, "Type", summarise, grp.mean = median(CVs))
head(mu)


# Basic density
ggplot(df, aes(x = CVs, fill = Type)) +
  geom_density(fill = "gray") +
  geom_vline(aes(xintercept = mean(CVs)),
    color = "blue",
    linetype = "dashed"
  ) +
  labs(title = "CV density curve", x = "P-value", y = "Density") +
  theme_classic()
# Change line colors by groups
p <- ggplot(df, aes(x = CVs, color = Type)) +
  geom_density() +
  geom_vline(
    data = mu, aes(xintercept = grp.mean, color = Type),
    linetype = "dashed"
  ) +
  labs(title = "CVs density curve", x = "Coefficient of Variation [%]", y = "Density")

p + scale_color_manual(values = c("#56B4E9", "#E69F00", "#56B4E9")) +
  theme_classic() +
  xlim(0, 300)
