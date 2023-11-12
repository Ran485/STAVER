## --------------- PCA  ---------------------------------
library("factoextra")
library("FactoMineR")
library("ggthemes")
library("imputePCA")
library("missMDA")
library("data.table")


matrix <- read.csv("~/merge_data/PCA/PXD024707/PCA_QC.csv",
    row.names = 1, header = 1
)
## drop NA row,  return the rows that have at least TWO non-NA value.
matrix[rowSums(is.na(matrix)) < (length(matrix) - 1), ]
matrix <- matrix[rowSums(is.na(matrix)) < (length(matrix) - 1), ]
# ## if not exist outliers, comment the follow two code
drop_outliers <- c("Exp068462", "Exp068465")
matrix <- matrix[, !names(matrix) %in% drop_outliers] # drop outliers
matrix <- t(matrix)
group_anatation <- read.csv("~/DIANN_data/merge_data/DEP_data/PXD024707/DEP_info.csv")
## if not exist outliers, comment the follow one code
group_anatation <- group_anatation[!group_anatation$id %in% drop_outliers, ] # drop outliers
group <- group_anatation$type
matrix[is.na(matrix)] <- 0.0001
# matrix[is.na(matrix)] = 0
matrix <- log2(matrix + 1)
# matrix = t(scale(t(matrix)))
matrix <- as.matrix(matrix)
res.pca <- PCA(matrix, # [,-c(1)],
    graph = TRUE, scale.unit = FALSE
)

# PCA plot
fviz_pca_ind(res.pca,
    label = "none", # hide individual labels "none""all"
    repel = TRUE,
    col.ind = group,
    legend.title = "Group",
    # palette = "npg",
    # habillage = annotation$BAP1, # color by groups dc3023
    axes = c(1, 2),
    palette = c("#dc3023", "#177cb0", "#EA686B", "#dc3023", "#EA686B", "#dc3023"),
    # palette = "npg", # ,"#EA686B" "#dc3023
    addEllipses = TRUE, # Concentration ellipses
    ellipse.level = 0.45,
    title = "PCA - The plasma raw data of COVID-19 cohort"
) + theme_base(base_size = 12)
# scale_shape_manual(values = c(16,17,15,19,19,19)) # Customize the shape of the points to 15, 19, and 17.
