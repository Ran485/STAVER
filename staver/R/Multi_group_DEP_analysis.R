# =========================================================================
#         Two groups analysis of differentially expressed proteins
# =========================================================================
library(dplyr)
library(greenbrown)
library("openxlsx")
library(xlsx)
### —————————————— wilcox.test and p.val-adjust function  ——————————————————————
# handle outliers
capping_outliers_IQR <- function(df, category_var, factor = 1.5) {
  qnt <- quantile(df[, category_var], probs = c(.25, .75), na.rm = T)
  caps <- quantile(df[, category_var], probs = c(.10, .90), na.rm = T)
  H <- factor * IQR(df[, category_var], na.rm = T)
  lower_bound <- qnt[1] - H
  upper_bound <- qnt[2] + H
  df[, category_var] <- ifelse(df[, category_var] < lower_bound, caps[1], df[, category_var])
  df[, category_var] <- ifelse(df[, category_var] > upper_bound, caps[2], df[, category_var])
  return(df)
}


group_DEP <- function(data, ana_index, capping_outliers = FALSE, group1, group2, method) {
  if (capping_outliers) {
    data <- apply(data, 1, capping_outliers_IQR)
  }
  index <- ana_index[ana_index[, 2] == group1 | ana_index[, 2] == group2, ]
  print(index)
  group_info <- data.frame(table(index$type))
  rownames(index) <- index[, 1]
  print(group_info)
  sep1 <- group_info[, group_info$Var1 == group1][2]
  sep2 <- sep1 + 1
  cat("sep1: ", sep1, "\n")
  cat("sep2: ", sep2, "\n")
  # library(dplyr)
  index_col <- intersect(index$id, colnames(data))
  expset <- data[, index_col]
  index_1 <- group1
  index_2 <- group2
  index_group1 <- index[index[, 2] == group1, ]
  index_group1 <- index_group1$id
  group1_expset <- expset[, index_group1]
  # print(group1_expset)
  index_group2 <- index[index[, 2] == group2, ]
  index_group2 <- index_group2$id
  group2_expset <- expset[, index_group2]
  # print(group2_expset)
  expset <- cbind(group1_expset, group2_expset)
  # print(expset)
  #  optional to choose calculate the difference method
  if (method == "t_test") {
    ### ————————————— t-test p-value  ——————————————————————————————-
    tstat.val <- apply(expset, 1, function(x) {
      if (AllEqual(x[1:ncol(expset)])) {
        return(1)
      } else if (AllEqual(x[1:sep1]) & AllEqual(x[sep2:ncol(expset)])) {
        return(0)
      } else {
        t.test(x[1:sep1], x[sep2:ncol(expset)], paired = F)$statistic
      }
    })
    tstat.pval <- apply(expset, 1, function(x) {
      if (AllEqual(x[1:ncol(expset)])) {
        return(1)
      } else if (AllEqual(x[1:sep1]) & AllEqual(x[sep2:ncol(expset)])) {
        return(0)
      } else {
        t.test(x[1:sep1], x[sep2:ncol(expset)], paired = F)$p.value
      }
    })
    p_val <- data.frame(tstat.val, tstat.pval)
    p_val[, "FDR_t-test"] <- NA
    p_val[, "FDR_t-test"] <- p.adjust(p_val[, "tstat.pval"], method = "BH", length(p_val[, "tstat.pval"]))
  } else if (method == "wilcoxn_test") {
    ### ————————————— wilcoxn-test p-value  ——————————————————————————————-
    wilcox_stat.val <- apply(expset, 1, function(x) {
      wilcox.test(x[1:sep1], x[sep2:ncol(expset)], paired = F)$statistic
    })
    wilcox_stat_pval <- apply(expset, 1, function(x) {
      wilcox.test(x[1:sep1], x[sep2:ncol(expset)], paired = F)$p.value
    })
    p_val <- data.frame(wilcox_stat.val, wilcox_stat_pval)
    p_val[, "FDR_wilcox"] <- NA
    p_val[, "FDR_wilcox"] <- p.adjust(p_val[, "wilcox_stat_pval"], method = "BH", length(p_val[, "wilcox_stat.val"]))
  } else {
    print("Please choose either `t-test` or `wilcoxn-test` method")
  }

  ### —————————————— Caculate group average and FC  ——————————————————————
  p_val[, index_1] <- apply(expset[, 1:sep1], 1, mean, na.rm = T)
  p_val[, index_2] <- apply(expset[, sep2:ncol(expset)], 1, mean, na.rm = T)
  p_val$count_1 <- rowSums(expset[, 1:sep1] >= 0, na.rm = T) ## 统计每组非 NA 值的个数
  p_val$count_2 <- rowSums(expset[, sep2:ncol(expset)] >= 0, na.rm = T) ## 统计每组非 NA 值的个数
  p_val$FC <- p_val[, index_2] / p_val[, index_1]
  p_val$log2_FC <- log2(p_val$FC)
  cat("index_1 =", index_1, "\n")
  cat("index_2 =", index_2, "\n")
  cat("FC = index_2/index_1:", index_2, "/", index_1)
  return(p_val)
}

data <- read.csv("~/DIANN_Data/merge_data/PXD025752.csv", header = T, row.names = 1, skip = 0)
ana_index <- read.csv("~/DIANN_Data/merge_data/RM_info.csv", header = T)
ana_index$id <- gsub("#", ".", ana_index$id)
table(ana_index$type)
## 原始数据处理
# data[is.na(data)] == 0
# data = log2(data+1)
res <- group_DEP(data, ana_index, group1 = "Plasma", group2 = "Serum", method = "wilcoxn_test")
# Export p-value results
write.csv(res, "~/results/clinical_EDA/continuous_clinical/significance_analysis/continuous_clinical_info_DEP.csv")


# =========================================================================
#     Multiple groups analysis of differentially expressed proteins
# =========================================================================
library(limma)
library(ComplexHeatmap)
library(Hmisc)
library(circlize)
library(BBmisc)

options(stringsAsFactors = F)
warnings("off")
rm(list = ls())
## --------- multi_group_DEP function ----------------------------
Z_Score <- function(x) {
  res <- (x - mean(x)) / sd(x)
  return(res)
}

Get_subtype_info <- function(meta) {
  f <- factor(meta$type)
  table(f)
  subtype_info <- data.frame(table(f))
  rownames(subtype_info) <- subtype_info[, 1]
  subtype_list <- rownames(subtype_info) ## 获取分组信息
  return(subtype_list)
}

## -------------- load file ----------------------
if (TRUE) {
  input_path <- "~/merge_data/data/processed_data/DEP/IPX0002924001/IPX0002924001_fot_normal.csv"
  anno_path <- "~/merge_data/data/processed_data/DEP/IPX0002924001/group_ana_normal.csv"
  input_data <- read.csv(input_path, header = T, row.names = 1, stringsAsFactors = F)
  meta <- read.csv(anno_path)
  exprSet <- as.matrix(input_data)
  ## 数据处理 ，磷酸化数据需要进行log转化
  exprSet[is.na(exprSet)] <- 0
  # exprSet<- log2(exprSet + 1)

  ## 数据存储路径
  setwd("~/merge_data/data/processed_data/DEP/IPX0002924001/DEP_normal_res")
  outpath <- getwd()
}
### 相关卡值
AdjPvalueCutoff <- 1
topNum <- 50
height <- 16
width <- 14
heatmap_breakNum <- 4
## the heatmap_annotation color
Tar_group <- "#2080C3"
The_remaining <- "#89BCDF" ## Tar_group ="#2080C3";The_remaining ="#89BCDF" 常用颜色
## the heatmap color
heatmap_col <- c("#2166ac", "#f7f7f7", "red")

### 获取分组信息
subtype_list <- Get_subtype_info(meta)

## ----------- using limma to conduct DE analysis ------------
for (i in subtype_list)
{
  print("-----------------------------------------------------")
  cat(i, "differential protein analysis is starting", sep = " ")
  print(i)
  # print("-----------------------------------------------------")
  Group <- meta
  Group$type <- ifelse(Group$type == i, i, "The_remaining") ## 表示亚型—1，和剩余其他亚型的比较
  group <- factor(Group$type)
  design <- model.matrix(~ 0 + group)
  # gsub("group","", colnames(design))
  colnames(design) <- c(i, "The_remaining") ## 更改名字
  rownames(design) <- colnames(exprSet)

  fit <- lmFit(exprSet, design)
  Contrasts <- paste(i, "The_remaining", sep = "-")
  cont.matrix <- makeContrasts(Contrasts, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  diff_gene <- topTable(fit2, adjust = "BH", number = Inf, p.value = AdjPvalueCutoff)
  diff_gene$cluster <- ifelse(diff_gene$t > 0, i, "The_remaining")
  # return(diff_gene)

  ## ------ to plot we take top 20 pathways by logFC -------------------------------

  diff_gene <- sortByCol(diff_gene, c("logFC")) # ,asc = FALSE 升序排列
  topGene <- c(rev(rownames(diff_gene))[1:topNum], rownames(diff_gene)[1:topNum])
  mat <- exprSet[topGene, ]

  top_Gene <- data.frame(topGene)
  topGeneset <- data.frame(mat)

  #----------- crete dir  --------------------------
  outpath <- paste(getwd(), i, sep = "/")
  if (file.exists(outpath)) {
    print("Path already exists")
  } else {
    dir.create(file.path(outpath))
  }

  #----------- Print results objects to workspace  --------------------------

  write.csv(diff_gene, paste(getwd(), i, paste("diff_gene_results_pvalue", ".csv", sep = ""), sep = "/"))
  write.csv(top_Gene, paste(getwd(), i, paste("top_Gene", ".csv", sep = ""), sep = "/"))
  write.csv(topGeneset, paste(getwd(), i, paste("topGeneset_results", ".csv", sep = ""), sep = "/"))

  ## -----------------------format pathway names ----------------------------------
  # rownames(mat) <- gsub("_", " ", rownames(mat))
  rownames(mat) <- gsub("KEGG ", "", rownames(mat))

  ## -----------------------format pathway names 大小写转换 -----------------------
  pathway_name <- rownames(mat)
  # pathway_name = tolower(pathway_name) ## 将x中的字符全部转换为小写字母
  # pathway_name = capitalize(pathway_name) ## 将y中的字符全部转换为大写字母
  rownames(mat) <- pathway_name
  if (TRUE) {
    Group$type <- ifelse(Group$type != "The_remaining", "Tar_group", "The_remaining")
    topAnn <- HeatmapAnnotation(
      df = Group[, "type", drop = F],
      col = list(type = c("Tar_group" = Tar_group, "The_remaining" = The_remaining)),
      # col = c("#2080C3", "#89BCDF"),
      annotation_height = unit(1, "cm")
    )

    heat.col <- colorRamp2(c(-heatmap_breakNum, 0, heatmap_breakNum), heatmap_col) # c('#2166ac', '#f7f7f7', 'red'))  ## 原始图例c('#2166ac', '#f7f7f7', '#b2182b'))
    mat <- apply(mat, 1, Z_Score)
    mat <- t(mat)
    ht <- Heatmap(mat,
      name = "Z-score", col = heat.col, top_annotation = topAnn,
      cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = TRUE,
      row_names_side = "right"
    )
    ## pdf oupput pathway
    pdf(paste(getwd(), i, paste("top_20_pathway", ".pdf", sep = ""), sep = "/"), height = height, width = width)
    maxChar <- rownames(mat)[nchar(rownames(mat)) == max(nchar(rownames(mat)))]

    padding <- unit.c(
      unit(2, "mm"),
      grobWidth(textGrob(maxChar)) - unit(50, "mm"),
      unit(c(2, 2), "mm")
    )
    draw(ht, padding = padding, merge_legends = TRUE)
    dev.off()
  }
  cat(i, "differential protein analysis is Done", sep = " ", "\n")
}
