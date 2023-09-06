setwd('~/analysis/zika/')
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
library(plyr)
library(dplyr)

#  先把数据load进来，load注释好的数据
# load('zika_all_mthb_filter_SCT_cluster_PC20_23cluster_immune_marker.Rda')

######桑基图 不能画，没有细胞名称同时存在######
head(zika_all_mthb_filter_SCT@meta.data)
plot.data <- zika_all_mthb_filter_SCT@meta.data[,c("Group", "cell_type")]
head(plot.data)

plot.data1 <- plot.data %>%
  rownames_to_column("id") %>%
  separate(id, into = c("Treat", "Num", "cell"), sep = "_") 
for (i in c(25116:37796)) {
  plot.data1[i,3] <- plot.data1[i,2]
}

plot.data1 <- plot.data1 %>%
  dplyr::select(3:5)
plot.data2 <- pivot_wider(plot.data1, id_cols = cell,
              names_from = Group,
              values_from = cell_type) %>%
  na.omit()
head(plot.data2)

corLodes <- to_lodes_form(plot.data1, axes = 2:3)

display.brewer.all()
mycol <- rep(brewer.pal(9, "Set1"), 15)
ggplot(corLodes, 
       aes(x = x,
           stratum = stratum,
           alluvium = cell,
           fill = stratum,
           label = stratum)) +
  geom_flow(width = 2/10, aes.flow = "forward") +
  geom_stratum(alpha = 0.8) +
  scale_fill_manual(values = mycol) +
  geom_text(stat = "stratum",
            size = 2,
            color = "black") +
  xlab("") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank()) +
  guides(fill = "none") 
ggsave("cell_type.pdf")

###### 细胞比例图#########
plot.data <- zika_all_mthb_filter_SCT@meta.data[c("Group", "cell_type")]
head(plot.data)
plot.data2 <- plot.data %>%
  group_by(Group, (plot.data)) %>%
  dplyr::summarise(num = n())
head(plot.data2)
ggplot(plot.data2, aes(x = Group,
                       weight = num, 
                       fill = cell_type)) +
  geom_bar(position = "dodge") +
  scale_fill_lancet() +
  theme_bw() + 
  guides(fill = "none")

display.brewer.all()
length(table(alldata@meta.data$cell_type))
Cell.col <- c(brewer.pal(9,"Pastel1"),brewer.pal(9,"Set1"),brewer.pal(9,"Set3"))

ggplot(plot.data2, aes(x = Group,
                       weight = num, 
                       fill = cell_type)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = Cell.col) +
  theme_bw() + 
  ylab("Cell Number") +
  xlab("Treat")


# 计算百分比
plot.data3 <- plot.data2 %>%
  ddply("Group", transform, percent = num/sum(num) * 100) %>%
  arrange(Group)
head(plot.data3)
# plot.data3$Group <- factor(plot.data3$group, levels = c("Befort_Treat", "After_Treat"))
ggplot(plot.data3, aes(x = Group,
                       weight = percent, 
                       fill = cell_type)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = Cell.col)+
  theme_bw() + 
  # guides(fill = "none") + 
  ylab("Cell Number") +
  xlab("Treat") + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x =element_text(size = 20),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank())


# 巨噬细胞marker 

load('data/step6_7_Macroxifen/M.data_cluster_PC30_4cluster_marker.Rda')
VlnPlot(M.data,
        features = c('Cd68', 'Clec4a2','Cd14',
                     'Aoah', 'Apoe','Cd84'
        ),
        group.by = "seurat_clusters",pt.size = 0)+NoLegend()

# 神经细胞marker
VlnPlot(M.data,
        features = c('Rbfox3', 'Map2','Meg3',
                     'Gap43', 'Miat','Ccnd2'
        ),
        group.by = "seurat_clusters",pt.size = 0)+NoLegend()

