#####降维、聚类####
setwd('~/analysis/zika_new/')
# rm(list = ls())
library(Seurat)
library(tidyverse)
library(stringr)
library(patchwork)
library(viridis)
load("zika_all_mthb_filter_SCT_scalemtcycle.Rda")
dim(zika_all_mthb_filter_SCT)

############PCA分析###################
#进行PCA分析
zika_all_mthb_filter_SCT <- RunPCA(object = zika_all_mthb_filter_SCT,
                   npcs = 50,
                   rev.pca = FALSE,
                   weight.by.var = TRUE,
                   verbose = TRUE, # 输出特定PC值的特征基因
                   ndims.print = 1:5, # 输出前5个PC值的特征基因
                   nfeatures.print = 30, # 输出每个PC值的30个特征基因
                   reduction.key = "PC_")
#这儿可以画个热图https://www.jianshu.com/p/21c1934e3061
View(zika_all_mthb_filter_SCT@reductions[["pca"]]@feature.loadings)#调去改变基因对应PC值
View(zika_all_mthb_filter_SCT@reductions[["pca"]]@cell.embeddings)#每个细胞在低维PCA轴上的映射坐标
DimHeatmap(object = zika_all_mthb_filter_SCT,dims = 1,cells = 1000,
           balanced = T, fast = F,nfeatures = 30)+ scale_fill_viridis()

#选择适合的PCA
ElbowPlot(zika_all_mthb_filter_SCT,
          ndims = 50) # 寻找对细胞差异贡献度较大的主成分

# 每个主成分的标志物热图展示
DimHeatmap(object = zika_all_mthb_filter_SCT,
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
#PCA可视化
#PCA可视化1: 一维可视化，对每个维度中，权重较大的基因进行可视化
#pdf(file="06.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = zika_all_mthb_filter_SCT, dims = 1:4, 
               reduction = "pca",nfeatures = 20)
#dev.off()

#PCA可视化2: 二维可视化
#pdf(file="06.PCA.pdf",width=6.5,height=6)
DimPlot(object = zika_all_mthb_filter_SCT, reduction = "pca")
#dev.off()


#########从降维的图⽚观察有无批次效应#########
#
p1 <- DimPlot(object = zika_all_mthb_filter_SCT,
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
p2 <- DimPlot(object = zika_all_mthb_filter_SCT,
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
zika_all_mthb_filter_SCT <- RunTSNE(zika_all_mthb_filter_SCT,
                    dims = 1:20)# 根据上⾯的贡献度图选择了20个PC值进⾏tSNE分析，所以我们直接运⾏ RunTSNE 函数
#按照group画tSNE图
p3 <- DimPlot(object = zika_all_mthb_filter_SCT,
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
p4 <- DimPlot(object = zika_all_mthb_filter_SCT,
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
zika_all_mthb_filter_SCT <- RunUMAP(zika_all_mthb_filter_SCT,
                    dims = 1:20)
#按照group进行区分
p5 <- DimPlot(object = zika_all_mthb_filter_SCT,
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
p6 <- DimPlot(object = zika_all_mthb_filter_SCT,
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
library(patchwork)
#拼一下按照group的
plot_grid(p1, p3, p5,
          labels = c("A", "B", "C"),
          nrow = 1,
          align = "h")
plot_grid(p2, p4, p6,
          labels = c("A", "B", "C"),
          nrow = 1,
          align = "h")
#拼一下按照细胞周期的
p2 + p4 + p6 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

#UMAP图展示特定基因
FeaturePlot(zika_all_mthb_filter_SCT,
            reduction = "tsne",
            features = "Mki67") #G2M期标志物
FeaturePlot(zika_all_mthb_filter_SCT,
            reduction = "tsne",
            features = "Top2a") #S期标志物 无

#######2.单细胞数据的聚类分析######
#K- Nearest Neighbor法进行聚类
zika_all_mthb_filter_SCT <- FindNeighbors(zika_all_mthb_filter_SCT,
                          k.param = 20,
                          dims = 1:20)
zika_all_mthb_filter_SCT <- FindClusters(zika_all_mthb_filter_SCT,
                         resolution = 0.5, # 最重要参数，该值越大，cluster 越多
                         method = "igraph", # 根据结果调整
                         algorithm = 1, # 根据结果调整'
                         random.seed = 2022) 
#统计一下各个cluster
table(zika_all_mthb_filter_SCT@meta.data$seurat_clusters)

#对细胞类型进行可视化
DimPlot(object = zika_all_mthb_filter_SCT,
        group.by = "seurat_clusters",
        reduction = "pca",
        label = TRUE)
PCAPlot(object = zika_all_mthb_filter_SCT,
        group.by = "seurat_clusters",
        pt.size = 0.5,
        label = TRUE)
#tSNE的聚类结果最好
TSNEPlot(object = zika_all_mthb_filter_SCT,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)
UMAPPlot(object = zika_all_mthb_filter_SCT,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)
#保存一下分群之后的细胞
save(zika_all_mthb_filter_SCT, file = "zika_all_mthb_filter_SCT_cluster23.Rda")
#细胞周期影响比较大，返回上一步去除细胞周期影响

#############3.细胞注释###################
#一般是load进来上一步骤的数据
load(file = "zika_all_mthb_filter_SCT_cluster23.Rda")
#读取总结好的marker信息 
cell_marker <- read.csv("mouse_brain_cell_atlas_paper.csv")
#提取文件中需要的信息
cell_marker <- cell_marker[c(3, 5)]
View(cell_marker)

#如果有好几个标记基因，使用下面的代码分开
cell_marker <- cell_marker %>%
  separate_rows(Cell.Marker, sep = ", ") %>%
  na.omit() %>%
  distinct() %>%
  arrange(Cell.Type)

#有的标志物前⾯有空格符号，因此我们要把这个空格符号去掉
cell_marker$Cell.Marker <- str_trim(cell_marker$Cell.Marker, "left")

#进行注释
p7 <- DimPlot(zika_all_mthb_filter_SCT,
              reduction = "tsne", # pca, umap, tsne
              group.by = "seurat_clusters",
              label = T)
p8 <- DotPlot(zika_all_mthb_filter_SCT,
              features = unique(cell_marker$Cell.Marker)) +
  theme(axis.text = element_text(size = 8,
                                 angle = 90,
                                 hjust = 1))
p7 + p8

#看一下表达范围
range(zika_all_mthb_filter_SCT@assays$integrated@scale.data)

# 看一下具体基因的表达
# 以B细胞为例子
p9 <- FeaturePlot(zika_all_mthb_filter_SCT,
                  reduction = "tsne",
                  features = c("Ms4a1", "Cd19", "Pxk", "Cd74", "Cd79a"),
                  label = TRUE)
p9
p7 + p9



#接下来使用singler试一试
BiocManager::install('SingleR')
BiocManager::install('celldex')
BiocManager::install('BiocParallel')
library(SingleR)
library(celldex)
library(BiocParallel)
#下载注释文件
ref <- celldex::MouseRNAseqData()
save(ref, file = "MouseRNAseqData.Rda")


#画一下每个cluster的差异基因
#load('zika_all_mthb_filter_SCT_cell_SingleR.Rda')
Idents(zika_all_mthb_filter_SCT) <- 'cell_type1'
Idents(zika_all_mthb_filter_SCT)

logFCfilter=0.5
adjPvalFilter=0.05
zika_all_mthb_filter_SCT.markers <- FindAllMarkers(object = zika_all_mthb_filter_SCT,
                                only.pos = TRUE, #TRUE表示只求高表达的基因
                                min.pct = 0.25, #表达的最小百分比
                                logfc.threshold = logFCfilter) # 这一步花了1-2h
# 把marker文件保存下来
save(zika_all_mthb_filter_SCT.markers, file = "zika_all_SCT.markers_0725.Rda")

sig.markers = zika_all_mthb_filter_SCT.markers[(abs(as.numeric(as.vector(zika_all_mthb_filter_SCT.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(zika_all_mthb_filter_SCT.markers$p_val_adj))<adjPvalFilter),]
write.table (sig.markers,file="PC20_09.markers.xls",sep="\t",row.names=F,quote=F)
top5 <- zika_all_mthb_filter_SCT.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- zika_all_mthb_filter_SCT.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- zika_all_mthb_filter_SCT.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table (top20 ,file="PC20_top20.markers.xls",sep="\t",row.names=F,quote=F)


#绘制marker在各个cluster的热图
pdf(file="09.tsneHeatmap.pdf",width=32,height=8)
DoHeatmap(object = zika_all_mthb_filter_SCT, features = top5$gene) + NoLegend()
dev.off()  # 绘制速度较慢

check.marker <- c('Ninj1','Cd68', # 0 Macrophage
                  'Ctss','C1qa', # 1 Microglia_Ctss high
                  'Lyz2','Ifitm3', # 2 Monocyte_Lyz2
                  'Cadm3','Syp', # 3 Neuron_Syp high
                  'Nkg7', 'Klrd1', # 4 NK cells 
                  'Hirip3', 'Mcm7', # 5 erythroid precursor cells
                  'C1ql1', # 6 OPC
                  'Fabp7','Gfap', # 7 Astrocyte_Fabp7
                  'Sox2','Sox4', # 8 Neural stem cell
                  'Sst','Dcc', # 9 Neuron_Sst
                  'Slc1a3','Gfap', # 10 Astrocyte_Slc1a3 high
                  'Gngt2','Csf1r', # 11 Microglia_Csf1r
                  'Nbea','Pax6', # 12 Neuron_Nbea
                  'Samhd1','Cd68', # 13 Monocytes & Neutrophils
                  'Camk2a','Ahi1', # 14 Excitatory_neuron
                  'Olig2','Mog', # 15 Oligodendrocyte
                  'Grrp1','Egfl7', # 16 Endothelial_cell_Egfl7
                  'Siglech','Cyb561a3', # 17 Microglia_Siglech
                  'Ttr','Clic6', # 18 Ependymal cell Ttr 
                  'Cd74', # 19 B cell_Cd74
                  'Rgs5','Mgp', # 20 Endothelial_cell_Mgp
                  'Ms4a1','Cd79a', # 21 B cell
                  'S100a9','S100a8', # 22 Neutrophils
                  'Cyp11a1','Hdc' # 23 Granulocyte_Cyp11a1
                  )

# 使用check.marker看一下
DotPlot(zika_all_mthb_filter_SCT,
        features = unique(check.marker)) +
  theme(axis.text = element_text(size = 8,
                                 angle = 90,
                                 hjust = 1))

#绘制marker的小提琴图
pdf(file="02.markerViolin.pdf",width=10,height=6)
VlnPlot(object = zika_all_mthb_filter_SCT, features = c("Cd83", "Cxcl2"))
dev.off()

#绘制marker在各个cluster的散点图
pdf(file="02.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = zika_all_mthb_filter_SCT, features = c("Cd83", "Cxcl2"))
dev.off()

#绘制marker在各个cluster的气泡图
pdf(file="06.markerBubble.pdf",width=12,height=6)
cluster0Marker=c("Cbr2", "Lyve1", "Selenop", "Folr2", "Ednrb", "F13a1", "Mrc1", "Igf1", "Slc40a1
", "Cd163")
DotPlot(object = zika_all_mthb_filter_SCT, features = cluster0Marker)
dev.off()

load('20220720recluster.Rda')
####写入细胞类型
Cell_type <- c(# "0" = "Macrophage",
               # "1" = "Microglia_Ctss high",
               # "2" = "Monocyte_Lyz2", #
               # "3" = "Neuron_Syp hig",
               # "4" = "NK cells", # Nkg7,Il2rb,Trbc1
               # "5" = "erythroid precursor cells", # Hirip3, Mcm7
               # "6" = "OPC",
               # "7" = "Astrocyte_Fabp7", 
               # "8" = "Neural stem cell",
               # "9" = "Neuron_Sst",
               # "10" = "Astrocyte_Slc1a3",
               # "11" = "Microglia_Csf1r",
               # "12" = "Neuron_Nbea",
               # "13" = "Monocytes & Neutrophils", #Dync1li2,
               # "14" = "Excitatory_neuron",
               # "15" = "Oligodendrocyte",
               # "16" = "Endothelial_cell_Egfl7", #Egfl7
               # "17" = "Microglia_Siglech",
               # "18" = "Ependymal_cell_Ttr",
               # "19" = "B cell_Cd74",
               # "20" = "Endothelial_cell_Mgp",
               # "21" = "B cell",
               # "22" = "Neutrophils",
               "23" = "Granulocyte_Cyp11a1") 

#将定义各个亚型对应的细胞名称，然后写⼊到meta.data⽂件
zika_all_mthb_filter_SCT[['cell_type']] <- unname(Cell_type[zika_all_mthb_filter_SCT@meta.data$seurat_clusters])

# 重新定义细胞类型
#将新矩阵cluster0，1，2，5，8分为Tcells。将cluster3分为Bcell，其他的为myeloid

zika_all_mthb_filter_SCT@meta.data$cell_type1  <-  ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(1,11,17),'Microglia',
                                 ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(2,13),'Monocyte',
                                        ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(3,9,12,14),'Neuron',
                                               ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(7,10),'Astrocyte',
                                                      ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(19,21),'Bcell',
                                                             ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(16,20),'Endothelial',
                                                                    ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(0),'Macrophage',
                                                                           ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(4),'NK',
                                                                                         ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(5),'erythroid',
                                                                                                ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(6),'OPC',
                                                                                                       ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(8),'NSC',
                                                                                                              ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(15),'Oligodendrocyte',
                                                                                                                     ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(18),'Ependymal',
                                                                                                                            ifelse(zika_all_mthb_filter_SCT@meta.data$seurat_clusters %in% c(22),'Neutrophils','Granulocyte'))))))))))))))
#看一下结果
table(zika_all_mthb_filter_SCT@meta.data$cell_type1)
zika_all_mthb_filter_SCT@meta.data$Group <- factor(zika_all_mthb_filter_SCT@meta.data$Group,
                                                      levels = c('MOCK','V_FS','V_59'))

#接着可以重新画⼀个UMAP图带上细胞名称
#UMAP
pdf('totalumap.pdf',width = 8,height = 8)
DimPlot(zika_all_mthb_filter_SCT,
        reduction = "umap",
        group.by = "cell_type1",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()
dev.off()

#tsne
pdf('totaltsne.pdf',width = 8,height = 8)
DimPlot(zika_all_mthb_filter_SCT,
        reduction = "tsne",
        group.by = "cell_type1",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()
dev.off()
#按照组别上色看看
pdf('group_umap.pdf',width = 20,height = 8)
DimPlot(zika_all_mthb_filter_SCT,
        reduction = "umap",
        split.by = "Group",
        group.by = "cell_type1",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()
dev.off()
#保存一下注释好的数据
save(zika_all_mthb_filter_SCT, file = "20220720recluster.Rda")


#### marker基因可视化 -----
# 散点图
FeaturePlot(zika_all_mthb_filter_SCT, features = c("C1qa", "Gfap", "Meg3", "Lyz2"),
            order = T)
# 小提琴图
VlnPlot(zika_all_mthb_filter_SCT, features = c("C1qa", "Gfap", "Meg3", "Lyz2"), ncol = 1)

# 气泡图
DotPlot(zika_all_mthb_filter_SCT, features = c("C1qa", "Gfap", "Meg3", "Lyz2"))

# 山脊图
RidgePlot(zika_all_mthb_filter_SCT, c("C1qa", "Gfap", "Meg3", "Lyz2"), ncol = 2)

# 热图
DoHeatmap(zika_all_mthb_filter_SCT, features = top10$gene)

load('~/analysis/zika/data/step6_7_Macroxifen/M.data_cluster_PC30_4cluster_marker.Rda')
# 气泡图 用来看巨噬细胞marker表达
DotPlot(M.data, features = c("Cd68", "Clec4a2", "Cd14",
                             'Aoah', 'Apoe', 'Arhgap15',
                             'Bank1', 'CD40', 'Cd84'))
FeaturePlot(M.data, features = c("Cd68", "Clec4a2", "Cd14"))
# 山脊图
RidgePlot(M.data, c("Cd68", "Clec4a2", "Cd14",
                                      'Aoah', 'Apoe', 'Arhgap15',
                                      'Bank1', 'CD40', 'Cd84'), group.by = 'cell_type',
           ncol = 3)
DoHeatmap(M.data, features = c("Cd68", "Clec4a2", "Cd14",
                                                 'Dock2', 'Abca9', 'Aoah', 'Apoe', 'Arhgap15',
                                                 'Bank1', 'CD40', 'Cd84', 'Cfh'))
# 神经细胞marker表达
# 山脊图 NeuN, Rph3a, B3gat1, Calb1, Calb2, Gap43, Pcdh20, Pcsk2, Actn2, Adamts5
RidgePlot(M.data, c("Rbfox3", "Map2", "Meg3",
                    'Rph3a', 'B3gat1', 
                    'Calb2', 'Gap43', 
                    'Snhg11', 'Miat','Ccnd2'), group.by = 'cell_type',
          ncol = 3)
