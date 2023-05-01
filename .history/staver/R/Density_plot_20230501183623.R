library('ggplot2')
library(plyr)

df = read.csv("/Volumes/T7_Shield/staver/results/DIA_repeat20/processed_data/Density_plot_fig2.csv")

df$weight = as.numeric(df$weight)
mu <- ddply(df, "sex", summarise, grp.mean=median(weight))
head(mu)


# Basic density
ggplot(df, aes(x=weight, fill=sex)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=mean(weight)), color="blue",
             linetype="dashed")+
  labs(title="CV density curve",x="P-value", y = "Density")+
  theme_classic()
# Change line colors by groups
p<- ggplot(df, aes(x=weight, color=sex)) +
  geom_density()+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")+
  labs(title="Weight density curve",x="Coefficient of Variation [%]", y = "Density")

p + scale_color_manual(values=c("#56B4E9", "#E69F00", "#56B4E9"))+
  theme_classic()

