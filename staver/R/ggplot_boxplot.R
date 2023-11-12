library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)
library("openxlsx")

## ---------------- Defined functions ------------------------

# Function to handle outliers using IQR method
capping_outliers_IQR <- function(df, category_var, factor = 4) {
  qnt <- quantile(df[, category_var], probs = c(.25, .75), na.rm = T)
  caps <- quantile(df[, category_var], probs = c(.10, .90), na.rm = T)
  H <- factor * IQR(df[, category_var], na.rm = T)
  lower_bound <- qnt[1] - H
  upper_bound <- qnt[2] + H
  df[, category_var] <- ifelse(df[, category_var] < lower_bound, caps[1], df[, category_var])
  df[, category_var] <- ifelse(df[, category_var] > upper_bound, caps[2], df[, category_var])
  return(df)
}

# Function to create a customized boxplot or violin plot
customed_boxplot <- function(data, target, category_var_name, log2_transform = FALSE, capping_outliers = FALSE, violin_plot = FALSE) {
  df <- data[, c(target, category_var_name)]
  if (capping_outliers) {
    df <- capping_outliers_IQR(df, category_var = category_var_name)
  }
  if (log2_transform) {
    df[, category_var_name] <- log2(df[, category_var_name] + 1)
  }
  plot <- ggplot(df, aes_string(x = target, y = category_var_name, color = target))
  if (violin_plot) {
    plot <- plot + geom_violin(trim = FALSE) +
      geom_boxplot(width = 0.6, position = position_dodge(0.9), outlier.colour = "red") +
      theme_bw() +
      stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, position = position_dodge(width = 0.9)) +
      scale_color_manual(values = c("#0E9F87", "#3C5588")) +
      ylab("Relative exprssion [log2]") +
      xlab("") +
      stat_compare_means(aes_string(group = target), label = "p.signif", method = "t.test") +
      ggtitle(category_var_name) +
      # ylim(23,27)+
      theme(
        axis.text.x = element_text(angle = 30, hjust = 0.6, colour = "black", family = "ArialMT", size = 16),
        axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", family = "ArialMT", size = 16),
        axis.text.y = element_text(family = "ArialMT", size = 16, face = "plain"),
        axis.title.y = element_text(family = "ArialMT", size = 16, face = "plain"),
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = 0.8),
        legend.text = element_text(
          face = "plain", family = "ArialMT", colour = "black",
          size = 12
        ),
        legend.title = element_text(
          face = "plain", family = "ArialMT", colour = "black",
          size = 12
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  } else {
    plot <- ggboxplot(df,
      x = target,
      y = category_var_name,
      color = target,
      fill = target,
      alpha = 0.12,
      width = 0.5,
      palette = c("#0E9F87", "#3C5588", "#FF9E29", "#86AA00", "#F94F21", "#916CA0", "#599BAD", "#DBD289"),
      x.text.angle = 45,
    ) +
      ylim(0, 2000) +
      rotate_x_text(angle = 45, hjust = .5, vjust = 1) +
      labs(title = "", x = "", y = "Relative exprssion [log2]") +
      geom_jitter(alpha = 0.7, aes_string(colour = target), width = 0.15, height = 0) +
      theme_test() +
      stat_compare_means(label.y = c(), paired = FALSE, method = "t.test") +
      # coord_flip()
      theme(
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
      ) +
      ggtitle(category_var_name) + # 设置变量名
      theme(
        axis.text.x = element_text(angle = 18, hjust = 0.6, colour = "black", family = "ArialMT", size = 16),
        axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", family = "ArialMT", size = 16),
        axis.text.y = element_text(family = "ArialMT", size = 16, face = "plain"),
        axis.title.y = element_text(family = "ArialMT", size = 16, face = "plain"),
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = 0.8),
        legend.text = element_text(
          face = "plain", family = "ArialMT", colour = "black",
          size = 12
        ),
        legend.title = element_text(
          face = "plain", family = "ArialMT", colour = "black",
          size = 14
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  }
  return(plot)
}

## ------------------------ Example Usage -----------------------------
data <- read.csv("/Volumes/Disk3/STAVER-revised/validation20_293T/results/diann_protein_num_compare.csv")
customed_boxplot(data,
  target = "Type",
  category_var_name = "Protein.Numbers",
  log2_transform = FALSE,
  capping_outliers = TRUE,
  violin_plot = FALSE
)


## --------- Analyzes the specified `category_var_colname` list in bulk. -----------
colnames(data)
category_list <- colnames(data)[-c(1:1)]
# category_list = c('Basophil_count', 'HBP', 'Plasma_fibrinogen', 'Chlorine')
category_var_plots <- list()

for (category_var in category_list) {
  category_var_plots[[category_var]] <- customed_boxplot(data, target = "Feature", category_var_name = category_var, log2_transform = F)
  print(category_var_plots[[category_var]])
  # ggsave(category_var_plots[[category_var]], file=paste0("plot_", category_var,".png"), width = 44.45, height = 27.78, units = "cm", dpi=300)
}

## A critical step in generating multiple diagrams
plot_grid(plotlist = category_var_plots)
