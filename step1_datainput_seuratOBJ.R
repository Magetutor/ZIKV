#读取数据并构建Seurat对象
setwd('~/analysis/zika/')
rm(list = ls())
# BiocManager::install('Seurat')
# BiocManager::install('scater')
# BiocManager::install('cowplot')
# install.packages("tidyverse")
library(Seurat)
library(tidyverse)

#######单个样本Seurat对象构建######
scRNA_9A59 <- Read10X(data.dir = "~/analysis/zika/data/210009A_59_filtered_feature_bc_matrix/")
scRNA_9AFS <- Read10X(data.dir = "~/analysis/zika/data/210009A_FS_filtered_feature_bc_matrix/")
scRNA_9AMOCK <- Read10X(data.dir = "~/analysis/zika/data/210009A_MOCK_filtered_feature_bc_matrix/")

#######创建Seurat对象#############
V_59 <- CreateSeuratObject(counts = scRNA_9A59,
                           min.cells = 3,
                           min.features = 200)
V_FS <- CreateSeuratObject(counts = scRNA_9AFS,
                           min.cells = 3,
                           min.features = 200)
MOCK <- CreateSeuratObject(counts = scRNA_9AMOCK,
                           min.cells = 3,
                           min.features = 200)
# 看一下结构
dim(V_59)
dim(V_FS)
dim(MOCK)

#把样本放在一起
zika_list <- list(V_59 = V_59,
                 V_FS = V_FS,
                 MOCK = MOCK)

#将样本merge在一起，进行数据指控
zika_all <- merge(x = zika_list[[1]], # 第⼀个
                    y = zika_list[c(2:length(zika_list))], # 除了第⼀个以外的Seurat对象
                    add.cell.ids = names(zika_list)) # cell id 添加前缀
dim(zika_all)

#metadata添加分组，标记不同株感染和对照组
# 添加分组信息
zika_all@meta.data <- zika_all@meta.data %>%
  mutate(Group = if_else(str_sub(rownames(zika_all@meta.data), 1, 4) == "V_59",
                         "V_59",
                         if_else(str_sub(rownames(zika_all@meta.data), 1, 4) == "V_FS",
                                 "V_FS",
                                 "MOCK")))
head(zika_all@meta.data)

save(zika_all, file = "zika_all_group_filter.Rda")
