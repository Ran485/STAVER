library(ggplot2)
library(scales)


data = read.csv("~/results/plasma_data_QC_proteins.csv", row.names = 1)
data$LFQ = data$LFQ_intensity 
# protein_hihglight = c('ALB','APOA1','C3','APOA2','ITIH1','C1QA','C1QC','ARL2','FETUB') #'SERPINA6',
# data = data[data$QC_proteins %in% protein_hihglight,]
# Change color/shape by groups
# Remove confidence bands
p <- ggplot(data, aes(x=index, y=LFQ, color=QC_proteins)) + 
  geom_point(shape=1)+ ## # 使用空心圆
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  # xlim(10000, 1000000000000)+
  labs(title="Analytical Variability",
       x="Sample ID", 
       y = "LFQ intensity")

# Continuous colors
p + scale_color_brewer(palette="Paired") + 
  theme_classic() + 
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=0), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=16),
        axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                  size=14),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())

