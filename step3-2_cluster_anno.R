#####降维、聚类####
setwd('~/analysis/Covid_brain/data/')
rm(list = ls())
library(Seurat)
library(tidyverse)
library(stringr)
library(patchwork)
load("cortex_all_mtinf_SCT_final.Rda")
dim(cortex_all_mtinf_SCT)

############PCA分析###################
#进行PCA分析
cortex_all_mtinf_SCT <- RunPCA(object = cortex_all_mtinf_SCT,
                               npcs = 50,
                               rev.pca = FALSE,
                               weight.by.var = TRUE,
                               verbose = TRUE, # 输出特定PC值的特征基因
                               ndims.print = 1:5, # 输出前5个PC值的特征基因
                               nfeatures.print = 30, # 输出每个PC值的30个特征基因
                               reduction.key = "PC_")
#选择适合的PCA
ElbowPlot(cortex_all_mtinf_SCT,
          ndims = 50) # 寻找对细胞差异贡献度较大的主成分

# 每个主成分的标志物热图展示
DimHeatmap(object = cortex_all_mtinf_SCT,
           reduction = "pca",
           dims = c(1:5),
           nfeatures = 30,
           disp.min = -2.5,
           balanced = TRUE,
           projected = FALSE,
           fast = TRUE,
           raster = TRUE,
           slot = "scale.data",
           combine = TRUE)

#########从降维的图⽚观察有无批次效应#########
#
p1 <- DimPlot(object = cortex_all_mtinf_SCT,
              reduction = "pca",
              group.by = "Group",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p1
###看细胞周期是否影响
p2 <- DimPlot(object = cortex_all_mtinf_SCT,
              reduction = "pca",
              group.by = "Phase",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p2

############tSNE分析################
cortex_all_mtinf_SCT <- RunTSNE(cortex_all_mtinf_SCT,
                                dims = 1:30)# 根据上⾯的贡献度图选择了19个PC值进⾏tSNE分析，所以我们直接运⾏ RunTSNE 函数
#按照group画tSNE图
p3 <- DimPlot(object = cortex_all_mtinf_SCT,
              reduction = "tsne",
              group.by = "Group",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p3

#再按照细胞周期分组画tSNE图
p4 <- DimPlot(object = cortex_all_mtinf_SCT,
              reduction = "tsne",
              group.by = "Phase",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p4
##############UMAP分析##############
cortex_all_mtinf_SCT <- RunUMAP(cortex_all_mtinf_SCT,
                                dims = 1:30)
#按照group进行区分
p5 <- DimPlot(object = cortex_all_mtinf_SCT,
              reduction = "umap",
              group.by = "Group",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p5

#按照细胞周期区分
p6 <- DimPlot(object = cortex_all_mtinf_SCT,
              reduction = "umap",
              group.by = "Phase",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p6

#拼一下这几种方法的图
library(cowplot)
#拼一下按照group的
plot_grid(p1, p3, p5,
          labels = c("A", "B", "C"),
          nrow = 1,
          align = "h")

#拼一下按照细胞周期的
p2 + p4 + p6 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

#UMAP图展示特定基因
FeaturePlot(cortex_all_mtinf_SCT,
            reduction = "umap",
            features = "ST18")

#######2.单细胞数据的聚类分析######
#K- Nearest Neighbor法进行聚类
cortex_all_mtinf_SCT <- FindNeighbors(cortex_all_mtinf_SCT,
                                      k.param = 20,
                                      dims = 1:30)
cortex_all_mtinf_SCT <- FindClusters(cortex_all_mtinf_SCT,
                                     resolution = 0.1, # 最重要参数，该值越大，cluster 越多
                                     method = "igraph", # 根据结果调整
                                     algorithm = 1, # 根据结果调整'
                                     random.seed = 2022) 
#统计一下各个cluster
table(cortex_all_mtinf_SCT@meta.data$seurat_clusters)

#对细胞类型进行可视化
DimPlot(object = cortex_all_mtinf_SCT,
        group.by = "seurat_clusters",
        reduction = "pca",
        label = TRUE)
PCAPlot(object = cortex_all_mtinf_SCT,
        group.by = "seurat_clusters",
        pt.size = 0.5,
        label = TRUE)
#tSNE的聚类结果最好
TSNEPlot(object = cortex_all_mtinf_SCT,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)
UMAPPlot(object = cortex_all_mtinf_SCT,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)
#保存一下分群之后的细胞
save(cortex_all_mtinf_SCT, file = "cortex_all_mtinf_SCT_cluster_PC30_10cluster.Rda")

#############3.细胞注释###################
#一般是load进来上一步骤的数据
#load(file = "cortex_all_mtinf_SCT_cluster.Rda")
#在第一步尝试的cellmarker以及singler效果均不好，直接用文章中给出的marker

cell_marker <- read.csv("mouse_brain_experiment_cellmarker.csv")
View(cell_marker)

#如果有好几个标记基因，使用下面的代码分开
cell_marker <- cell_marker %>%
  separate_rows(Cell.Marker, sep = ",") %>%
  na.omit() %>%
  distinct() %>%
  arrange(Cell.Type)

#有的标志物前⾯有空格符号，因此我们要把这个空格符号去掉
cell_marker$Cell.Marker <- str_trim(cell_marker$Cell.Marker, "left")

#进行注释
p7 <- DimPlot(cortex_all_mtinf_SCT,
              reduction = "tsne", # pca, umap, tsne
              group.by = "seurat_clusters",
              label = T)
p8 <- DotPlot(cortex_all_mtinf_SCT,
              features = unique(cell_marker$Cell.Marker)) +
  theme(axis.text = element_text(size = 8,
                                 angle = 90,
                                 hjust = 1))
p7 + p8
#UMAP看一下
p9 <- DimPlot(cortex_all_mtinf_SCT,
              reduction = "umap", # pca, umap, tsne
              group.by = "seurat_clusters",
              label = T)
p9

#聚成10类 9 10细胞很少，再缩小一下聚类范围
#######2.单细胞数据的聚类分析######
#K- Nearest Neighbor法进行聚类 聚成了8类别
cortex_all_mtinf_SCT <- FindNeighbors(cortex_all_mtinf_SCT,
                                      k.param = 20,
                                      dims = 1:30)
cortex_all_mtinf_SCT <- FindClusters(cortex_all_mtinf_SCT,
                                     resolution = 0.06, # 最重要参数，该值越大，cluster 越多
                                     method = "igraph", # 根据结果调整
                                     algorithm = 1, # 根据结果调整'
                                     random.seed = 2022) 
#统计一下各个cluster
table(cortex_all_mtinf_SCT@meta.data$seurat_clusters)

#还是画一下tsne和UMAP
#tSNE的聚类结果最好
TSNEPlot(object = cortex_all_mtinf_SCT,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)
UMAPPlot(object = cortex_all_mtinf_SCT,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)

#按照covid文章的marker进行注释
p10 <- DimPlot(cortex_all_mtinf_SCT,
              reduction = "tsne", # pca, umap, tsne
              group.by = "seurat_clusters",
              label = T)
p11 <- DotPlot(cortex_all_mtinf_SCT,
              features = unique(cell_marker$Cell.Marker)) +
  theme(axis.text = element_text(size = 8,
                                 angle = 90,
                                 hjust = 1))
p12 <- DimPlot(cortex_all_mtinf_SCT,
               reduction = "umap", # pca, umap, tsne
               group.by = "seurat_clusters",
               label = T)
p10 + p11
p12 + p11

#写入细胞类型
Cell_type <- c("0" = "Oligodendrocyte",
               "1" = "Astrocyte",
               "2" = "L5/6 CC excitatory neuron", #皮质-皮质
               "3" = "Microglia",
               "4" = "OPC", #少突胶质细胞前体细胞
               "5" = "SST interneuron", #生长抑素中间神经元
               "6" = "L2/3 excitatory neuron",
               "7" = "Endothelial", #内皮细胞
               "8" = "PV interneuron") 

#将定义各个亚型对应的细胞名称，然后写⼊到meta.data⽂件
cortex_all_mtinf_SCT[['cell_type']] <- unname(Cell_type[cortex_all_mtinf_SCT@meta.data$seurat_clusters])

#看一下结果
View(srt3_all@meta.data)

#接着可以重新画⼀个UMAP图带上细胞名称
#UMAP
DimPlot(cortex_all_mtinf_SCT,
        reduction = "umap",
        group.by = "cell_type",
        label = TRUE,
        pt.size = 0.5) +
  NoLegend()
#tsne
DimPlot(cortex_all_mtinf_SCT,
        reduction = "tsne",
        group.by = "cell_type",
        label = TRUE,
        pt.size = 0.5) +
  NoLegend()

#保存一下注释好的数据
save(cortex_all_mtinf_SCT, file = "cortex_all_mtinf_SCT_cluster_PC30_8cluster_marker.Rda")
