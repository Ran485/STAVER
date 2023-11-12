# =========================================================================
#       Panel=Fig4 - Heatmap of the reported biomarkers
# =========================================================================

library(pheatmap)
library(RColorBrewer)
matrix <- read.csv("~/DIANN_Data/merge_data/data/processed_data/Barplot_Normal/heatmap_biomarker.csv",
    header = T,
    skip = 0, row.names = 1
)
# Preprocessing the data
# matrix[is.na(matrix)] = 1
# matrix <- log2(matrix+1)
# matrix[matrix == 1] <- NA

# color for subtype
annotation_col <- read.csv("~/DIANN_Data/merge_data/data/processed_data/Barplot_Normal/heatmap_ana.csv", header = T, row.names = 1)
rownames(annotation_col) <- colnames(matrix)
# ann_colors = list(type = c(iCC_N= "#48cae4", eCC_N = "#f48c06", iCC_T = "#0077b6" , eCC_T = "#dc2f02")) ## #0772BC 蓝色
# ann_colors = list(Subtype = c(S_I = "#F94F21",S_II = "#EFC000",S_III = "#00AEBA"))
chending_color <- c(
    colorRampPalette(c("#0772BC", "white"))(50), ## #184D87  常用#1E90FF ##14417C ## 0772BC ## #395F90黑蓝   ## 淡蓝#1E90FF
    colorRampPalette(c("white", "red"))(50)
)

breaks <- unique(c(seq(-6, 0, length = 50), 0, seq(0, 6, length = 50)))

## 选取基因进行可视化
# df = read.csv("~/results/biomarkers.csv")
# x <- read.delim(pipe("pbpaste"),row.names = 1) ## 读取剪切板的数据
# gene = x$topGene
# # gene = c("ROCK1","MAPK3")
# data = matrix[gene,]
# # data[is.na(data)] = 0
# write.csv(data,"~/results/heatmap_biomarkers.csv")
# matrix = log2(matrix+1)
pheatmap(matrix,
    cluster_rows = 5, cluster_cols = 0,
    clustering_distance_cols = "correlation", fill = T, breaks = breaks,
    clustering_distance_rows = "correlation", border_color = "white", na_col = "white",
    col = chending_color, show_rownames = T, show_colnames = F, display_numbers = F,
    width = 4.85, height = 5, fontsize_number = 12, number_color = "black", number_format = "%.1f",
    annotation_col = annotation_col
) # , annotation_colors = ann_colors) #, annotation_row = annotation_row)
