library(ggplot2)
library(scales)
library(ggrepel)
require(extrafont)
# need only do this once!
library(showtext)
font_add("Arial", "Arial.ttf")
showtext_auto()
# extrafont::loadfonts()
# font_import(pattern="[A/a]rial", prompt=T)
## --------------- scatter CV plot ------------
QC_data_cv <- read.csv("~/QC_repeat90_highCI_peptides_cv/90QC_data_cv.csv", check.names = F)
colnames(QC_data_cv)
genes <- QC_data_cv
genes["LFQ_intensity"] <- genes$Abundance
# genes$Significant <- ifelse(genes$p < 0.05 & genes$HR > 1 , "Bad prognosis",
#                             ifelse(genes$p < 0.05 & genes$HR < 1 , "good prognosis","Not Sig"))
genes$CV_threshold <- ifelse(genes$`Coefficient of Variation [%]` < 30, "CV% < 30",
      ifelse(genes$`Coefficient of Variation [%]` < 50, "CV% < 50",
            ifelse(genes$`Coefficient of Variation [%]` < 100, "CV% < 100",
                  ifelse(genes$`Coefficient of Variation [%]` < 250, "CV% < 250", "CV% >= 250")
            )
      )
)

# pdf(file = "./figure/CV_protein_abundance_scatter.pdf", width = 8, height = 5.5)
p1 <- ggplot(genes, aes(x = LFQ_intensity, y = `Coefficient of Variation [%]`)) +
      geom_point(aes(color = CV_threshold), shape = 1) + # changen the shape of the dot point
      scale_size_continuous(range = c(1, 4)) +
      # c("gray", "#7AA9CE") c( "#EA686B","#7AA9CE","gray"))
      scale_color_manual(values = c("#515151", "red", "gray", "#7AA9CE", "#EA686B")) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right") +
      geom_text_repel(
            data = subset(genes, CV_threshold == "select"),
            aes(label = index),
            size = 5,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
      ) +
      scale_x_continuous(
            trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
      ) +
      theme(
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
      ) +
      ggtitle("Quality control for the top3 peptides") +
      theme(
            axis.text.x = element_text(angle = 0, hjust = 0.6, colour = "black", family = "Arial", size = 16),
            axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", family = "Arial", size = 16),
            axis.text.y = element_text(family = "Arial", size = 16, face = "plain"),
            axis.title.y = element_text(family = "Arial", size = 20, face = "plain"),
            panel.border = element_blank(), axis.line = element_line(colour = "black", size = 0.8),
            legend.text = element_text(
                  face = "plain", family = "Arial", colour = "black",
                  size = 12
            ),
            legend.title = element_text(
                  face = "plain", family = "Arial", colour = "black",
                  size = 12
            ),
            panel.grid.major = element_blank(), # 不显示网格线
            panel.grid.minor = element_blank()
      ) +
      geom_hline(yintercept = 30, color = "gray", linetype = "dashed") +
      geom_hline(yintercept = 50, color = "gray", linetype = "dashed") +
      geom_hline(yintercept = 100, color = "gray", linetype = "dashed") +
      geom_hline(yintercept = 250, color = "gray", linetype = "dashed") +
      # xlim(min, max)+
      ylim(0, 30)
# geom_vline(xintercept = -0.15,color = 'gray', linetype="dashed") +
# geom_vline(xintercept = 0.15,color = 'gray', linetype="dashed")
p1
# ggsave(p1,file = "/Volumes/Samsung_T5/DIA-QC/results/Figure/CV_TOP1_peptide_abundance_scatter.pdf", width = 8, height = 5.5)


## ------ the idendification protein_numbers_of_diff_CV_thresholds ----------
protein_num <- as.data.frame(table(genes$CV_threshold))
colnames(protein_num) <- c("CV_thresholds", "Quantative protein numbers")
protein_num$`Quantative protein numbers [%]` <- round((protein_num$`Quantative protein numbers` / sum(protein_num$`Quantative protein numbers`) * 100), 1)
# specifically order ggplot 2x axis instead of alphabetical order
protein_num$CV_thresholds <- factor(protein_num$CV_thresholds, levels = c("CV% < 30", "CV% < 50", "CV% < 100", "CV% < 250", "CV% >= 250"))
## Inside bars
# pdf(file = "./figure/protein_num_barplot.pdf", width = 7, height = 5.5)
p2 <- ggplot(data = protein_num, aes(x = CV_thresholds, y = `Quantative protein numbers [%]`)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_text(aes(label = `Quantative protein numbers [%]`), vjust = 1.6, color = "white", size = 6) +
      theme(
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
      ) +
      ggtitle("The idendification peptide numbers of diff CV thresholds") +
      theme(
            axis.text.x = element_text(angle = 0, hjust = 0.6, colour = "black", family = "Arial", size = 16),
            axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", family = "Arial", size = 16),
            axis.text.y = element_text(family = "Arial", size = 16, face = "plain"),
            axis.title.y = element_text(family = "Arial", size = 20, face = "plain"),
            panel.border = element_blank(), axis.line = element_line(colour = "black", size = 0.8),
            legend.text = element_text(
                  face = "plain", family = "Arial", colour = "black",
                  size = 12
            ),
            legend.title = element_text(
                  face = "plain", family = "Arial", colour = "black",
                  size = 12
            ),
            panel.grid.major = element_blank(), # 不显示网格线
            panel.grid.minor = element_blank()
      ) +
      geom_hline(yintercept = 25, color = "gray", linetype = "dashed")

p2
ggsave(p2, file = "~/results/Figure/TOP3_peptide_num_barplot.pdf", width = 7, height = 5.5)



## --------------- cumulative CV plot ------------
genes <- read.csv("~/QC_repeat90_highCI_peptides_cv/90QC_data_cv.csv", check.names = F)
colnames(genes)
genes$CV_threshold <- ifelse(genes$`Coefficient of Variation [%]` < 30, "CV% < 30",
      ifelse(genes$`Coefficient of Variation [%]` < 50, "CV% < 50",
            ifelse(genes$`Coefficient of Variation [%]` < 100, "CV% < 100",
                  ifelse(genes$`Coefficient of Variation [%]` < 300, "CV% < 250", "CV% >= 250")
            )
      )
)

p3 <- ggplot(genes, aes(x = Abundance_Rank, y = `Coefficient of Variation [%]`)) +
      geom_point(aes(color = CV_threshold), shape = 1) + ## changen the shape of the dot point
      scale_size_continuous(range = c(1, 4)) +
      scale_color_manual(values = c("#515151", "red", "gray", "#7AA9CE", "#EA686B")) + # ("gray", "#7AA9CE") c( "#EA686B","#7AA9CE","gray"))
      theme_bw(base_size = 12) +
      theme(legend.position = "right") +
      geom_text_repel(
            data = subset(genes, CV_threshold == "select"),
            aes(label = index),
            size = 5,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.3, "lines")
      ) +
      # scale_x_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
      #                    labels = trans_format("log10", math_format(10^.x))) +
      theme(
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
      ) +
      ggtitle("Cumulative CV plot for the top3 peptides") +
      theme(
            axis.text.x = element_text(angle = 0, hjust = 0.6, colour = "black", family = "Arial", size = 16),
            axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", family = "Arial", size = 16),
            axis.text.y = element_text(family = "Arial", size = 16, face = "plain"),
            axis.title.y = element_text(family = "Arial", size = 20, face = "plain"),
            panel.border = element_blank(), axis.line = element_line(colour = "black", size = 0.8),
            legend.text = element_text(
                  face = "plain", family = "Arial", colour = "black",
                  size = 12
            ),
            legend.title = element_text(
                  face = "plain", family = "Arial", colour = "black",
                  size = 12
            ),
            panel.grid.major = element_blank(), # 不显示网格线
            panel.grid.minor = element_blank()
      ) +
      geom_hline(yintercept = 30, color = "gray", linetype = "dashed") +
      geom_hline(yintercept = 50, color = "gray", linetype = "dashed") +
      geom_hline(yintercept = 100, color = "gray", linetype = "dashed") +
      geom_hline(yintercept = 250, color = "gray", linetype = "dashed") +
      annotate("text", x = 5000, y = 15, label = "Nums = 16673", colour = "black", family = "Arial", size = 6) +
      annotate("text", x = 5000, y = 40, label = "Nums = 3823", colour = "black", family = "Arial", size = 6) +
      annotate("text", x = 5000, y = 70, label = "Nums = 5929", colour = "black", family = "Arial", size = 6) +
      annotate("text", x = 5000, y = 170, label = "Nums = 4271", colour = "black", family = "Arial", size = 6) +
      annotate("text", x = 5000, y = 260, label = "Nums =  4", colour = "black", family = "Arial", size = 6) +
      # xlim(min, max)
      ylim(0, 30)

table(genes$CV_threshold)
# geom_vline(xintercept = -0.15,color = 'gray', linetype="dashed") +
# geom_vline(xintercept = 0.15,color = 'gray', linetype="dashed")
p3
# ggsave(p1,file = "/Volumes/Samsung_T5/DIA-QC/results/Figure/CV_TOP3_peptide_cumulative_plot_with_nums.pdf", width = 8, height = 5.5)
