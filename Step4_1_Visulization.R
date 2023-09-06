setwd('~/analysis/zika_new/')
# 根据注释好的细胞进行一些可视化
# 加载包
library(tidyverse)
library(reshape2)
library(patchwork)
library(Seurat)
library(SingleR)
library(celldex)
library(BiocParallel)
library(cowplot)
library(RColorBrewer)
library(ggplot2)
library(ggalluvial)
library(ggsci)
library(ggpubr)
library(plyr)
library(dplyr)
library(ggthemes)

#  先把数据load进来，load注释好的数据
load('20220820recluster_noerythroid.Rda')

# 画图
plot.data1 <- seurat.obj@meta.data[c("Group", "cell_type1")]
head(plot.data1)
plot.data1$Group <- factor(plot.data1$Group, levels = c("MOCK","V_FS","V_59"))
plot.data2 <- plot.data1 %>%
  group_by(Group, (plot.data1)) %>%
  dplyr::summarise(num = n())
head(plot.data2)

plot.data2$Group <- factor(plot.data2$Group, levels = c("MOCK","V_FS","V_59"))
ggplot(plot.data2, aes(x = cell_type1,
                       weight = num, 
                       fill = Group)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = groupcol) +
  theme_base() + theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
                    legend.position = 'top') +
  geom_text(aes(label=num, y = num), size = 3, position=position_dodge(0.9), vjust=-0.5) 

ggsave('cellnumber.pdf',
       ggplot2::last_plot(),
       width = 12,
       height = 6)

Cell.col <-c("#9479ad", 
             "#cc4a50",
             "#f8766d", 
             "#83ccd2", 
             "#47b9c9",
             "#e77e2c",
             "#9b9c9b",
             "#7fb687",
             "#1d933f", 
             "#d570a9", 
             "#caa980",
             "#00b6eb",
             "#986156",
             "#1271b4")


ggplot(plot.data2, aes(x = Group,
                       weight = num)) +
  geom_bar(position = "stack",aes( fill = cell_type1)) +
  scale_fill_manual(values = Cell.col) +
  theme_base() + 
  ylab("Cell Number") +
  xlab("")
ggsave('cellnumber_stack0822.pdf',
       ggplot2::last_plot(),
       width = 6,
       height = 6)

# 计算百分比
plot.data3 <- plot.data2 %>%
  ddply("Group", transform, percent = num/sum(num) * 100) %>%
  arrange(Group)
head(plot.data3)
# plot.data3$Group <- factor(plot.data3$group, levels = c("Befort_Treat", "After_Treat"))

ggplot(plot.data3, aes(x = Group,
                       weight = percent, 
                       fill = cell_type1)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = Cell.col)+
  theme_bw() + 
  # guides(fill = "none") + 
  ylab("Cell Percentage") +
  xlab("Treat") + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x =element_text(size = 12),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))  
  # theme(panel.grid =element_blank()) + 
  # theme(panel.border = element_blank())
ggsave('cellpercent_stack0822.pdf',
       ggplot2::last_plot(),
       width = 6,
       height = 6)

######新的比例图展示######
# 字体报错
library(extrafont)
extrafont::font_import() #这里用处不大，自己电脑应该可以
# Error in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  : 
# nvalid font type
library(ggstatsplot)  
  # (2)绘图
  figure = ggbarstats(data = plot.data1, 
                      x = Group, y = cell_type1,
                      package = 'ggsci',
                      palette = 'category20c_d3',
                      results.subtitle = FALSE,
                      bf.message = FALSE,
                      proportion.test = FALSE,
                      label.args = list(size = 4, 
                                        fill = 'white', 
                                        alpha = 0.85,
                                       #family = 'Arial',
                                        fontface = 'bold'),
                      perc.k = 2,
                      title = '',
                      xlab = '',
                      legend.title = 'Seurat Cluster',
                      ggtheme = ggpubr::theme_pubclean()) +
    scale_fill_manual(values = c("#e44349","#5797bc","#84bb7a")) +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = 'black', lineend = 'round'),
        legend.position = 'right',
        axis.text.x = element_text(angle = 45,hjust = 1,size = 12, color = 'black'),#, family = 'Arial'
        axis.text.y = element_text(size = 12, color = 'black'), #, family = 'Arial'
        legend.text = element_text( size = 10, color = 'black'), #family = 'Arial',
        legend.title = element_text( size = 13, color = 'black'))  #family = 'Arial',
  pdf('modified_percentage_0822.pdf',width = 12,height = 6)
  figure
  dev.off()

# (3)去除柱子下面的样本量标识：
gginnards::delete_layers(x = figure, match_type = 'GeomText')

#######分组各细胞类型批量做图##############
plot_b <- dplyr::filter(plot.data3, cell_type1 == 'B cell')
plot_b$percent <- round(plot_b$percent,3)

colp <- ggplot(plot_b, aes(x = Group,
                       weight = percent, 
                       fill = Group)) +
  geom_bar(position = "dodge") + #width = 0.9默认0.9调剂柱子宽度
  scale_fill_brewer(palette="Dark2") +
  theme_classic2() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = 'black'),
                     #legend.position = 'top', aspect.ratio = 0.4
        axis.text.y = element_text(size = 10, color = 'black')) +
  geom_text(aes(label=percent, y = percent), size = 4, position=position_dodge(0.9), vjust=-0.5)+
  coord_cartesian(ylim = c(0, 1.5))
pdf(paste0(names(table(plot.data3$cell_type1))[1],'.pdf'),width = 5,height = 5)
colp
dev.off()

# 批量将图做出来
for (i in 1:length(names(table(plot.data3$cell_type1)))) {
  plot_b <- dplyr::filter(plot.data3, cell_type1 == names(table(plot.data3$cell_type1))[i])
  plot_b$percent <- round(plot_b$percent,3)
  
  colp <- ggplot(plot_b, aes(x = Group,
                             weight = percent, 
                             fill = Group)) +
    geom_bar(position = "dodge") + #width = 0.9默认0.9调剂柱子宽度
    scale_fill_manual(values = groupcol) +
    theme_classic2() + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = 'black'),
          #legend.position = 'top', aspect.ratio = 0.4
          axis.text.y = element_text(size = 10, color = 'black')) +
    geom_text(aes(label=percent, y = percent), size = 4, position=position_dodge(0.9), vjust=-0.5)
    #coord_cartesian(ylim = c(0, 1.5))
  
  pdf(paste0(names(table(plot.data3$cell_type1))[i],'.pdf'),width = 5,height = 5)
  print(colp) # 最后打印图片要用print
  dev.off()
  
}


# 单独算一下细胞占比
num <- table(seurat.obj@meta.data$cell_type1) %>% as.data.frame() %>% 
 mutate(num, prop = Freq/sum(Freq))


