library(circlize)
library(ComplexHeatmap)
library(Hmisc)

sx_syn <- read.csv("~/processed_data/Barplot_Normal/heatmap_pval.csv", header = TRUE)

syn_spearman <- Heatmap(sx_syn, circlize::colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c(
    "#d7191c", "#fdae61",
    "#FFFFFF", "#abd9e9",
    "#2c7bb6"
  )
),
column_names_gp = grid::gpar(
  fontsize = 20,
  col = c(rep("#4daf4a", 8), rep("#984ea3", 7))
),
row_names_gp = grid::gpar(
  fontsize = 20,
  col = c(rep("#4daf4a", 8), rep("#984ea3", 7))
),
heatmap_legend_param = list(
  title = "Spearman's rho",
  direction = "horizontal",
  at = c(-1, -0.5, 0, 0.5, 1)
),
cell_fun = function(j, i, x, y, w, h, fill) {
  if (sx_syn_p[i, j] < 0.05) {
    grid.text("*", x, y)
  }
}
)
sx_syn_spear <- draw(syn_spearman, heatmap_legend_side = "top")
