library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(tidyr)
library(openxlsx)

# Capping Outliers using Interquartile Range (IQR)
#
# This function caps outliers in a dataframe based on the IQR method.
#
# Args:
#   df: A dataframe containing the data.
#   category_var: The name of the variable/column to apply the outlier capping.
#   factor: A numeric value representing the factor to calculate the bounds for outlier capping (default is 4).
#
# Returns:
#   A dataframe with outliers capped in the specified column.
capping_outliers_IQR <- function(df, category_var, factor = 4) {
    qnt <- quantile(df[, category_var], probs = c(.25, .75), na.rm = TRUE)
    caps <- quantile(df[, category_var], probs = c(.10, .90), na.rm = TRUE)
    H <- factor * IQR(df[, category_var], na.rm = TRUE)
    lower_bound <- qnt[1] - H
    upper_bound <- qnt[2] + H
    df[, category_var] <- ifelse(df[, category_var] < lower_bound, caps[1], df[, category_var])
    df[, category_var] <- ifelse(df[, category_var] > upper_bound, caps[2], df[, category_var])
    return(df)
}

# Producing Summary Statistics
#
# This function calculates the mean, standard deviation, and confidence intervals for a given numeric vector.
#
# Args:
#   x: A numeric vector.
#
# Returns:
#   A numeric vector containing the mean, lower bound, and upper bound.
data_summary <- function(x) {
    m <- mean(x)
    ymin <- m - sd(x)
    ymax <- m + sd(x)
    return(c(y = m, ymin = ymin, ymax = ymax))
}

# Extract Rows by Specified Groups
#
# Extracts rows from a dataframe based on specified group values.
#
# Args:
#   dataframe: A dataframe from which to extract rows.
#   groups: A vector of group values to filter the dataframe.
#
# Returns:
#   A dataframe containing only the rows matching the specified groups.
extract_rows_by_group <- function(dataframe, groups) {
    return(dataframe[dataframe$group %in% groups, ])
}

# Custom Boxplot Function
#
# Creates a customized boxplot for given data, allowing various options like outlier capping, log2 transformation, and violin plot.
#
# Args:
#   data: Dataframe containing the data.
#   target: The target column for plotting (default is "Feature").
#   category_var_name: The name of the variable/column to plot (default is "SOD1").
#   capping_outliers: Boolean to indicate if outliers should be capped (default is FALSE).
#   log2_transform: Boolean to indicate if log2 transformation is needed (default is FALSE).
#   group_vs_rest: Specific group to compare against the rest (default is NULL).
#   desired_groups: Vector of groups to include in the plot (default is NULL).
#   violin_plot: Boolean to indicate if a violin plot should be drawn (default is FALSE).
#
# Returns:
#   A ggplot object representing the customized boxplot.
customed_boxplot <- function(data,
                             target = "Feature",
                             category_var_name = "SOD1",
                             capping_outliers = FALSE,
                             log2_transform = FALSE,
                             drop_col = NULL,
                             group_vs_rest = NULL,
                             desired_groups = NULL,
                             method = "t.test",
                             handel_na = FALSE,
                             handle_target_na = FALSE,
                             sort_mean = FALSE,
                             violin_plot = FALSE) {
    df <- data[, c(target, category_var_name)]
    if (capping_outliers) {
        df1 <- capping_outliers_IQR(df, category_var = category_var_name)
    } else {
        df1 <- df
    }
    if (log2_transform) {
        df1[, category_var_name] <- log2(df1[, category_var_name] + 1)
    }
    if (exists("drop_col") && !is.null(drop_col)) {
        df1 <- df1[df1[, target] != drop_col, ]
    }
    # Processing for group_vs_rest and desired_groups...
    if (exists("group_vs_rest") && !is.null(group_vs_rest)) {
        print("group_vs_rest is defined and not NULL!")
        print(paste("The group_vs_rest is:", group_vs_rest))
        df1[, target] <- ifelse(df1[, target] == group_vs_rest, group_vs_rest, "rest")
        print(table(df1[, target]))
    } else {
        print("group_vs_rest is either not defined or NULL!")
    }

    if (exists("desired_groups") && !is.null(desired_groups)) {
        df1 <- df1[df1[, target] %in% desired_groups, ]
        df1[, category_var_name][df1[, category_var_name] == 0 & df1[, target] == desired_groups[1]] <- NA
        print(table(df1[, target]))
    }

    if (handel_na) {
        df1[df1[, category_var_name] == 0, category_var_name] <- NA
    }

    if (handle_target_na) {
        df1[, category_var_name][df1[, category_var_name] == 0 & df1[, target] == group_vs_rest] <- NA
    }
    # Calculating the mean for each group ...
    group_means <- aggregate(get(category_var_name) ~ get(target), data = df1, FUN = mean)
    # Sorting the groups by mean ...
    sorted_groups <- group_means[[1]][order(group_means[[2]])]
    if (sort_mean) {
        df1[, target] <- factor(df1[, target], levels = sorted_groups)
    }

    # Plotting logic for boxplot and violin plot...
    # Note: The actual plotting code (using ggplot2 functions) should be here.
    if (violin_plot) {
        p <- ggplot(df1, aes_string(x = target, y = category_var_name, color = target)) +
            geom_violin(trim = FALSE) +
            geom_boxplot(width = 0.3, position = position_dodge(0.8)) +
            scale_color_manual(values = c("#0E9F87", "#3C5588", "#FF9E29", "#86AA00", "#F94F21", "#916CA0", "#599BAD", "#DBD289", "#EFC000", "#00AEBA")) +
            theme_bw() +
            stat_summary(fun.y = mean, geom = "point", shape = 23, size = 4, position = position_dodge(width = 0.9)) +
            theme(
                panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
            ) +
            # ggtitle("Cis efffect on Protein")+
            # ylim(0,2000)+
            theme(
                axis.text.x = element_text(angle = 30, hjust = 0.6, colour = "black", family = "ArialMT", size = 16),
                axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", family = "ArialMT", size = 16),
                axis.text.y = element_text(family = "ArialMT", size = 16, face = "plain"),
                axis.title.y = element_text(family = "ArialMT", size = 16, face = "plain"),
                panel.border = element_blank(), axis.line = element_line(colour = "black", size = 0.8),
                legend.text = element_text(face = "plain", family = "ArialMT", colour = "black", size = 12),
                legend.title = element_text(face = "plain", family = "ArialMT", colour = "black", size = 12),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
            ) +
            ylab("Relative exprssion [log2]") +
            xlab("") +
            stat_compare_means(aes_string(group = target), method = method, label = "p.format", label.y = c())
        if (exists("group_vs_rest") && !is.null(group_vs_rest)) {
            p <- p + ggtitle(paste0(paste0(paste0(group_vs_rest, " vs. rest [ "), category_var_name), " ]"))
        } else {
            p <- p + ggtitle(category_var_name)
        } # Pairwise comparison against all t.test wilcox
        print(p)
    } else {
        p <- ggboxplot(df1,
            x = target, y = category_var_name, combine = F,
            color = target,
            fill = target,
            alpha = 0.12,
            # order =  c("WT","Mutation"),
            order = c(group_vs_rest, "rest"),
            # palette = "Paired",
            palette = c("#0E9F87", "#3C5588", "#FF9E29", "#86AA00", "#F94F21", "#916CA0", "#599BAD", "#DBD289"),
            width = 0.5, x.text.angle = 45
        ) +
            rotate_x_text(angle = 45, hjust = .5, vjust = 1) +
            # Global p-value
            stat_compare_means(label.y = c(), paired = FALSE) +
            labs(title = "", x = "", y = "Relative exprssion [log2]") +
            geom_jitter(alpha = 0.7, aes_string(colour = target), width = 0.15, height = 0) +
            # ylim(0,0.2)+
            theme_test() +
            # coord_flip()
            theme(
                panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
            ) +

            theme(
                axis.text.x = element_text(angle = 18, hjust = 0.6, colour = "black", family = "ArialMT", size = 16),
                axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", family = "ArialMT", size = 16),
                axis.text.y = element_text(family = "ArialMT", size = 16, face = "plain"),
                axis.title.y = element_text(family = "ArialMT", size = 16, face = "plain"),
                panel.border = element_blank(), axis.line = element_line(colour = "black", size = 0.8),
                legend.text = element_text(face = "plain", family = "ArialMT", colour = "black", size = 12),
                legend.title = element_text(face = "plain", family = "ArialMT", colour = "black", size = 14),
                panel.grid.major = element_blank(), # 不显示网格线
                panel.grid.minor = element_blank()
            )
        if (exists("group_vs_rest") && !is.null(group_vs_rest)) {
            p <- p + ggtitle(paste0(paste0(paste0(group_vs_rest, " vs. rest [ "), category_var_name), " ]"))
        } else {
            p <- p + ggtitle(category_var_name)
        }
        print(p)
    }
    return(p)
}

# Load and process data from CSV files
#
# This function reads data from a CSV file and its corresponding annotation file,
# transposes the data, binds annotation to it, and ensures no negative or NA values.
#
# Args:
#   data_path: Path to the data CSV file.
#   annotation_path: Path to the annotation CSV file.
#
# Returns:
#   A transposed and annotated dataframe with no negative or NA values.
load_data <- function(data_path, annotation_path) {
    # Read the data and annotation files
    if (!file.exists(data_path)) {
        stop("Data file does not exist: ", data_path)
    }
    if (!file.exists(annotation_path)) {
        stop("Annotation file does not exist: ", annotation_path)
    }

    data <- read.csv(data_path, row.names = 1)
    annotation <- read.csv(annotation_path, row.names = 1)

    # Transpose the data and bind with annotation
    data_t <- t(data)
    data_t <- cbind(data_t, annotation)

    # Replace NA and negative values with zero
    data_t[is.na(data_t) | data_t < 0] <- 0

    return(data_t)
}


## -------------------------  Load Data  -----------------------------
STAVER_data <- load_data("~/STAVER-revised/PCA/STAVER_data.csv", "~/STAVER-revised/PCA/PCA_anatation_scanpy.csv")
raw_data <- load_data("~/STAVER-revised/PCA/raw_data.csv", "~/STAVER-revised/PCA/PCA_anatation_scanpy.csv")


## ------------------------ Example Usage -----------------------------
# desired_groups <- c("Non-Small Cell Lung Carcinoma", "Small Cell Lung Carcinoma")
customed_boxplot(STAVER_data,
    target = "Cancer_type",
    group_vs_rest = "Hepatocellular Carcinoma",
    # desired_groups = desired_groups,
    category_var_name = "CTSE",
    drop_col = "Other Solid Carcinomas",
    method = "t.test", # "anova", "wilcox.test", "t.test"
    log2_transform = F,
    # capping_outliers = T,
    violin_plot = T,
    sort_mean = F
)
