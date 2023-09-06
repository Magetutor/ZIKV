# 设置工作目录
setwd('~/analysis/zika/')
# 清空变量 加载包
# rm(list = ls())
library(Seurat)
library(DESeq2)
library(tidyverse)

# 加载数据
load('zika_all_mthb_filter_SCT_cluster_PC20_23cluster_immune_marker.Rda')
# 巨噬细胞亚分群
# 挑选T细胞进⾏进⼀步分群
target.cell <- row.names(zika_all_mthb_filter_SCT@meta.data)[which(zika_all_mthb_filter_SCT@meta.data$cell_type == "NK cells")]
length(target.cell)
NK.data <- subset(zika_all_mthb_filter_SCT, cells = target.cell)
NK.data
NK.data <- SplitObject(NK.data,
                      split.by = "Group")
NK.data

# 还是先看一下线粒体这些,整合之后才可以看
# VlnPlot(object = NK.data,
#         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#         group.by = "Group",
#         log = F,
#         pt.size = 0.1,
#         raster=FALSE)
# VlnPlot(object = NK.data,
#         features = c("percent.ribo","percent.hb"),
#         group.by = "Group",
#         log = F,
#         pt.size = 0.1,
#         raster=FALSE)

# 第一次做过线粒体、核糖体、细胞周期质控了所以直接从CT开始
# SCTransform
# 速度较慢，样本越多越慢
# 循环运⾏
for(i in 1:length(NK.data)){
  NK.data[[i]] <- SCTransform(
    NK.data[[i]],
    variable.features.n = 3000,
    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), # 去除线粒体和细胞周期的影响
    verbose = FALSE)
}

# 数据整合
# 选择⽤于整合的基因
# load("NK.data.SCT.Rda")
features <- SelectIntegrationFeatures(object.list = NK.data,
                                      nfeatures = 3000)
# 准备整合
NK.data <- PrepSCTIntegration(object.list = NK.data,
                             anchor.features = features)
names(NK.data)
# 寻找锚定基因(速度很慢)
AnchorSet <- FindIntegrationAnchors(object.list = NK.data,
                                    reference = 3,
                                    normalization.method = "SCT",
                                    anchor.features = features
)
# save(AnchorSet, file = "AnchorSet.Rda")
# 整合数据
# load("AnchorSet.Rda")
NK.data <- IntegrateData(anchorset = AnchorSet,
                        normalization.method = "SCT")
DefaultAssay(NK.data) <- "integrated"
save(NK.data, file = "NK.data.final.Rda")
############降维分析##########
rm(list = ls())
library(Seurat)
library(tidyverse)
library(stringr)
library(patchwork)
load("NK.data.final.Rda")
dim(NK.data)
# PCA分析
NK.data <- RunPCA(object = NK.data,
                 npcs = 50,
                 rev.pca = FALSE,
                 weight.by.var = TRUE,
                 verbose = TRUE, # 输出特定PC值的特征基因
                 ndims.print = 1:5, # 输出前5个PC值的特征基因
                 nfeatures.print = 30, # 输出每个PC值的30个特征基因
                 reduction.key = "PC_")
NK.data <- RunPCA(NK.data)
ElbowPlot(NK.data,
          ndims = 50) # 寻找对细胞差异贡献度较⼤的主成分
# 选择30个主成分进⾏后续分析
# t-SNE分析，基于PCA分析结果进⾏⾮线性降维分析
NK.data <- RunTSNE(NK.data,
                  dims = 1:40)
# UMAP分析
NK.data <- RunUMAP(NK.data,
                  dims = 1:40)
# 展示PCA降维特征及其特征基因的热图
DimHeatmap(object = NK.data,
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
# 按照分组观察是否具有批次效应
p1 <- DimPlot(object = NK.data,
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
# 按照细胞周期阶段观察是否具有批次效应
p2 <- DimPlot(object = NK.data,
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
p3 <- DimPlot(object = NK.data,
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
p4 <- DimPlot(object = NK.data,
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
p5 <- DimPlot(object = NK.data,
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
p6 <- DimPlot(object = NK.data,
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
# 合并三种降维⽅法的结果
library(cowplot)
plot_grid(p1, p3, p5,
          labels = c("A", "B", "C"),
          nrow = 1,
          align = "h")
library(patchwork)
p2 + p4 + p6 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p1 + p3 + p5 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# 从上⾯结果看，未⻅明显批次效应
# 检查线粒体⽐例是否为协变量
FeaturePlot(NK.data,
            reduction = "umap",
            features = "percent.mt")
# 展示特定基因
FeaturePlot(NK.data,
            reduction = "umap",
            features = "Cxcl2")
##########聚类分析##########
# KNN聚类法
NK.data <- FindNeighbors(NK.data,
                        k.param = 20, # 不要⼩于5
                        dims = 1:40)
NK.data <- FindClusters(NK.data,
                       resolution = 0.4,
                       method = "igraph",
                       algorithm = 1,
                       random.seed = 2021)
table(NK.data@meta.data$seurat_clusters)

PCAPlot(object = NK.data,
        group.by = "seurat_clusters",
        pt.size = 0.5,
        label = TRUE)
TSNEPlot(object = NK.data,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)
UMAPPlot(object = NK.data,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)

# 保存一下分好群的数据
save(NK.data, file = "NK.data_SCT_cluster5.Rda")

# 计算每个cluster的差异基因
# 没有数据的话需要load('NK.data_SCT_cluster5.Rda')

logFCfilter=0.5
adjPvalFilter=0.05
NK.data.markers <- FindAllMarkers(object = NK.data,
                                 only.pos = TRUE, #TRUE表示只求高表达的基因
                                 min.pct = 0.25, #表达的最小百分比
                                 logfc.threshold = logFCfilter) # 这一步花了1-2h
# 把marker文件保存下来
save(NK.data.markers, file = "NK.data.markers_cluster5.Rda")

sig.markers = NK.data.markers[(abs(as.numeric(as.vector(NK.data.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(NK.data.markers$p_val_adj))<adjPvalFilter),]
write.table (sig.markers,file="PC40_NK.markers.xls",sep="\t",row.names=F,quote=F)
top5 <- NK.data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- NK.data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- NK.data.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table (top20 ,file="PC30_top20_Macro.markers.xls",sep="\t",row.names=F,quote=F)


#绘制marker在各个cluster的热图
pdf(file="09.tsneHeatmap.pdf",width=32,height=8)
DoHeatmap(object = NK.data, features = top10$gene,size = 3,angle = 0) + NoLegend()
dev.off()  # 绘制速度较慢

FeaturePlot(object = NK.data, features = c("Irak3", "Bcl2a1b"))
FeaturePlot(object = NK.data, features = c( "Irak3", "Bcl2a1b",
                                           'Sell', 'Lef1',
                                           "Cd7", "Tyrobp", 
                                           'Ltb4r1', 'Blk',
                                           'Nfib', 'Mir99ahg'),
            reduction = "umap")
VlnPlot(NK.data,
        features = c('Sell', 'Lef1'
                     ),
        group.by = "seurat_clusters")

# 先给每个群命名一下
####写入细胞类型
Cell_type <- c("0" = "Irak3_NK",
               "1" = "Sell_NK",
               "2" = "Cd7_NK", #
               "3" = "Ltb4r1_NK",
               "4" = "Nfib_NK"
) 

#将定义各个亚型对应的细胞名称，然后写⼊到meta.data⽂件
NK.data[['cell_type']] <- unname(Cell_type[NK.data@meta.data$seurat_clusters])

#看一下结果
View(NK.data@meta.data)

#接着可以重新画⼀个UMAP图带上细胞名称
#UMAP
DimPlot(NK.data,
        reduction = "umap",
        group.by = "cell_type",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()
#tsne
DimPlot(NK.data,
        reduction = "tsne",
        group.by = "cell_type",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()

#按照组别上色看看
DimPlot(NK.data,
        reduction = "umap",
        split.by = "Group",
        group.by = "cell_type",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()
#保存一下注释好的数据
save(NK.data, file = "NK.data_cluster_PC40_5cluster_marker.Rda")
#### 保存细胞亚型 ----
saveRDS(NK.data, file = "NK.rds")

###找一下各个亚群的特异性marker，看看是否存在########
# 这一部分效果不好，不用管
# install.packages('BiocManager')
# # BiocManager::install('multtest')
# library(multtest)
# # install.packages('metap')
# library(metap)
# 
# # 先给定细胞身份
# Idents(NK.data) <- NK.data@meta.data$cell_type
# levels(NK.data)
# table(NK.data@meta.data$cell_type)
# 
# View(NK.data@meta.data)
# Macro_Apoe_high2 <- FindConservedMarkers(NK.data,
#                                       ident.1 = "Macro_Apoe_high",
#                                       grouping.var = "Group")
# # 查看原来注释的B细胞的marker
# 
# Macro_Apoe_high2.marker.res <- Macro_Apoe_high2[c("Cd68", "Ninj1"), ]
# Macro_Apoe_high2.marker.res
# 
# FeaturePlot(NK.data,
#             features = c("Cd68", "Ninj1"))
# DotPlot(NK.data, 
#         assay = "SCT",
#         features = c("Cd68", "Ninj1"))
# 
# # 可视化对比
# dput(rownames(head(Macro_Apoe_high2)))
# dput(rownames(tail(Macro_Apoe_high2)))
# DotPlot(NK.data, 
#               features = c("Gadd45b", "Atf3", "Ccl4", "Ier5", "Basp1", "Nfkbia"))
# DotPlot(NK.data, 
#               features = c("Brk1", "Dek", "Ndufs4", "Cd74", "Pdia6", "Ldha"))
# p3 + p4
# 
# # 筛选
# head(Macro_Apoe_high2)
# Macro_Apoe_high2.marker.marker <- Macro_Apoe_high2 %>%
#   dplyr::filter(After_Treat_avg_log2FC > 20 & 
#                   Befort_Treat_avg_log2FC > 20 &
#                   After_Treat_p_val_adj < 0.05 &
#                   Befort_Treat_p_val_adj < 0.05)
# dput(rownames(head(Bcell.marker.marker)))
# DotPlot(NK.data, 
#         features = c(c("IGHM", "IGLC2", "IGLC3", "PTMA", "NPM1", "CCL4")))
# 
# # 分组切割
# DotPlot(NK.data, 
#         features = c(c("IGHM", "IGLC2", "IGLC3", "PTMA", "NPM1", "CCL4")),
#         cols = c("blue", "darkred"),
#         dot.scale = 8,
#         split.by = "group") + 
#   RotatedAxis()
# 
# # 展示目标marker
# FeaturePlot(NK.data,
#             features = c("IGHM"))
# VlnPlot(NK.data, 
#         features = c("IGHM"),
#         split.by = "group",
#         pt.size = 0.2)

load('/home/zhangc/analysis/zika/data/NK/NK.data_cluster_PC40_5cluster_marker.Rda')
########## 巨噬细胞亚群的差异基因分析#############
# 批量找各个细胞在V59中的差异基因
dput(unique(NK.data@meta.data$cell_type))
table(NK.data@meta.data$cell_type)
Cells <- c("Irak3_NK", "Ltb4r1_NK", "Cd7_NK", "Sell_NK", "Nfib_NK")
diff_V59_NK <- list()
for (i in 1:length(Cells)) {
  cells <- subset(NK.data,
                  cell_type == Cells[i])
  Idents(cells) <- cells@meta.data$Group
  cells.diff <- FindMarkers(cells,
                            ident.1 = "V_59",
                            ident.2 = "MOCK",
                            test.use = "MAST",
                            min.pct = 0.5, # 过滤掉那些在50%以下细胞中检测到的基因
                            logfc.threshold = log(2) # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
  )
  diff_V59_NK[[i]] <- cells.diff
}
names(diff_V59_NK) <- Cells
save(diff_V59_NK, file = "diff_V59_NK_alldiff.Rda")

# 每个细分亚群的功能富集结果
# ######富集分析#######
# 基因ID转换
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

length(diff_V59_NK)
rownames(diff_V59_NK[[1]])
ENTREZID <- bitr(rownames(diff_V59_NK[[1]]), 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = "org.Mm.eg.db")
Irak3_NK.diff <- list(ENTREZID)

# 批量转换
scRNA.diff_V59_NK <- list()
for (i in 1:length(diff_V59_NK)) {
  ENTREZID <- bitr(rownames(diff_V59_NK[[i]]), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
  ENTREZID <- ENTREZID$ENTREZID
  scRNA.diff_V59_NK[[i]] <- ENTREZID
}
names(scRNA.diff_V59_NK) <- names(diff_V59_NK)

# 多个cluster比较
# KEGG分析
xx <- compareCluster(scRNA.diff_V59_NK, 
                     fun="enrichKEGG",
                     organism="mmu",
                     pvalueCutoff=0.05,
                     qvalueCutoff=0.05)
xx <- pairwise_termsim(xx) # 获取相似性矩阵
dotplot(xx, showCategory = 5) #,label_format = 60 可以修改通路一行还是两行

p1 <- emapplot(xx)
p2 <- emapplot(xx, 
               legend_n = 2) 
p3 <- emapplot(xx, 
               pie = "count")
p4 <- emapplot(xx, 
               pie = "count", 
               cex_category = 1.5, 
               layout="kk")


# GO分析
# all 整体分析会有NA出现，单独分析看看
xx_BP <- compareCluster(scRNA.diff_V59_NK, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "BP",
                        pvalueCutoff=0.0001,
                        qvalueCutoff=0.000001)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_BP, showCategory = 5) +
  facet_grid(.~ONTOLOGY) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15,  vjust = 1, 
                                   hjust = 1, angle = 45), # color = "navy",
        axis.text.y = element_text(size = 12,  vjust = 1, 
                                   hjust = 1, angle = -8)) # color = "darkred",
ggsave("compareCluster.GO.pdf")

xx_MF <- compareCluster(scRNA.diff_V59_NK, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "MF",
                        pvalueCutoff=0.001,
                        qvalueCutoff=0.001)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_MF, showCategory = 5) +
  facet_grid(.~ONTOLOGY) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15,  vjust = 1, 
                                   hjust = 1, angle = 45), # color = "navy",
        axis.text.y = element_text(size = 12,  vjust = 1, 
                                   hjust = 1, angle = -8)) # color = "darkred",
# ggsave("compareCluster.GO.pdf")

xx_CC <- compareCluster(scRNA.diff_V59_NK, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "CC",
                        pvalueCutoff=0.001,
                        qvalueCutoff=0.01)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_CC, showCategory = 5) +
  facet_grid(.~ONTOLOGY) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15,  vjust = 1, 
                                   hjust = 1, angle = 45), # color = "navy",
        axis.text.y = element_text(size = 12,  vjust = 1, 
                                   hjust = 1, angle = -8)) # color = "darkred",
# ggsave("compareCluster.GO.pdf")

# 循环做一下Gsea吧
## 解释原理
library(clusterProfiler)

#获得基因列表
gene <- rownames(diff_V59_NK[[5]])
## 转换
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=diff_V59_NK[[5]]$avg_log2FC,
                     SYMBOL = rownames(diff_V59_NK[[5]]))
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneListgse <- gene_df$logFC
## 2.命名
names(geneListgse) = as.character(gene_df$ENTREZID)
## 3.排序很重要
geneListgse = sort(geneListgse, decreasing = TRUE)
head(geneListgse)
##真正进行的gsea分析
gsea <- gseKEGG(geneListgse, organism = "mmu", 
                pvalueCutoff = 0.999,
                scoreType = "pos") #物种查询https://www.genome.jp/kegg/catalog/org_list.html
head(gsea)
a <- gsea
dim(gsea)
gsea <- setReadable(gsea, 'org.Mm.eg.db', 'ENTREZID')
gsea_result<-as.data.frame(gsea@result)
names(diff_V59_NK)

write.csv(gsea_result,'Nfib_NK_gsea_result.csv')

#我们先把结果提取出来，方便查看
library(enrichplot)
gseaplot2(gsea, "mmu04142",color = 'steelblue',  title = "Lysosome")
gseaplot2(gsea, "mmu04015",color = 'steelblue',  title = "Rap1 signaling pathway")
gseaplot2(gsea, "hsa04727",color = 'steelblue',  title = "GABAergic synapse")

############VFS组各个细胞亚型差异功能#############
# 批量找各个细胞在VFS中的差异基因
dput(unique(NK.data@meta.data$cell_type))
table(NK.data@meta.data$cell_type)
Cells <- c("Irak3_NK", "Ltb4r1_NK", "Cd7_NK", "Sell_NK", "Nfib_NK")
diff_VFS_NK <- list()
for (i in 1:length(Cells)) {
  cells <- subset(NK.data,
                  cell_type == Cells[i])
  Idents(cells) <- cells@meta.data$Group
  cells.diff <- FindMarkers(cells,
                            ident.1 = "V_FS",
                            ident.2 = "MOCK",
                            test.use = "MAST",
                            min.pct = 0.5, # 过滤掉那些在50%以下细胞中检测到的基因
                            logfc.threshold = log(2) # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
  )
  diff_VFS_NK[[i]] <- cells.diff
}
names(diff_VFS_NK) <- Cells
save(diff_VFS_NK, file = "diff_VFS_NK_alldiff.Rda")


# ######富集分析
# 基因ID转换
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

length(diff_VFS_NK)
rownames(diff_VFS_NK[[1]])
ENTREZID <- bitr(rownames(diff_VFS_NK[[1]]), 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = "org.Mm.eg.db")
Macrophage.diff <- list(ENTREZID)

# 批量转换
scRNA.diff_VFS_NK <- list()
for (i in 1:length(diff_VFS_NK)) {
  ENTREZID <- bitr(rownames(diff_VFS_NK[[i]]), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
  ENTREZID <- ENTREZID$ENTREZID
  scRNA.diff_VFS_NK[[i]] <- ENTREZID
}
names(scRNA.diff_VFS_NK) <- names(diff_VFS_NK)

# 多个cluster比较
# KEGG分析
xx <- compareCluster(scRNA.diff_VFS_NK, 
                     fun="enrichKEGG",
                     organism="mmu",
                     pvalueCutoff=0.05,
                     qvalueCutoff=0.05)
xx <- pairwise_termsim(xx) # 获取相似性矩阵
dotplot(xx, showCategory = 5)

p1 <- emapplot(xx)
p2 <- emapplot(xx, 
               legend_n = 2) 
p3 <- emapplot(xx, 
               pie = "count")
p4 <- emapplot(xx, 
               pie = "count", 
               cex_category = 1.5, 
               layout="kk")

cowplot::plot_grid(p1,
                   p2, 
                   p3,
                   p4,
                   ncol = 2, 
                   labels = LETTERS[1:4])
ggsave("compareCluster.kegg.pdf")

# GO分析
# all 整体分析会有NA出现，单独分析看看
xx_BP <- compareCluster(scRNA.diff_VFS_NK, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "BP",
                        pvalueCutoff=0.0001,
                        qvalueCutoff=0.00001)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_BP, showCategory = 5)

xx_MF <- compareCluster(scRNA.diff_VFS_NK, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "MF",
                        pvalueCutoff=0.001,
                        qvalueCutoff=0.001)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_MF, showCategory = 5) 

xx_CC <- compareCluster(scRNA.diff_VFS_NK, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "CC",
                        pvalueCutoff=0.01,
                        qvalueCutoff=0.01)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_CC, showCategory = 5)


################细胞亚型数量及比例##############
plot.data <- NK.data@meta.data[c("Group", "cell_type")]
head(plot.data)
plot.data1 <- plot.data %>%
  group_by(Group, (plot.data)) %>%
  dplyr::summarise(num = n())
head(plot.data1)
plot.data1$Group <- factor(plot.data1$Group, levels = c("MOCK","V_FS","V_59"))
ggplot(plot.data1, aes(x = cell_type,
                       weight = num, 
                       fill = Group)) +
  geom_bar(position = "dodge") +
  scale_fill_lancet() +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'),
                     legend.position = 'top',
                     aspect.ratio = 0.4) +
  geom_text(aes(label=num, y = num), size = 3, position=position_dodge(0.9), vjust=-0.5) 

display.brewer.all()
Cell.col2 <- c(brewer.pal(8,"Accent"),brewer.pal(8,"Set2"))
Cell.col <- brewer.pal(8,"Set2")
ggplot(plot.data1, aes(x = Group,
                       weight = num, 
                       fill = cell_type)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = Cell.col2) +
  theme_bw() + 
  ylab("Cell Number") +
  xlab("Treat")


# 计算百分比
plot.data2 <- plot.data1 %>%
  ddply("Group", transform, percent = num/sum(num) * 100) %>%
  arrange(Group)
head(plot.data2)
# plot.data3$Group <- factor(plot.data3$group, levels = c("Befort_Treat", "After_Treat"))
ggplot(plot.data2, aes(x = Group,
                       weight = percent, 
                       fill = cell_type)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = Cell.col)+
  theme_bw() + 
  # guides(fill = "none") + 
  ylab("Cell Percentage") +
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

####批量画图百分比
for (i in 1:length(names(table(plot.data2$cell_type)))) {
  plot_n <- dplyr::filter(plot.data2, cell_type == names(table(plot.data2$cell_type))[i])
  plot_n$percent <- round(plot_n$percent,3)
  
  colp <- ggplot(plot_n, aes(x = Group,
                             weight = percent, 
                             fill = Group)) +
    geom_bar(position = "dodge") + #width = 0.9默认0.9调剂柱子宽度
    scale_fill_brewer(palette="Dark2") +
    theme_classic2() + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = 'black'),
          #legend.position = 'top', aspect.ratio = 0.4
          axis.text.y = element_text(size = 10, color = 'black')) +
    geom_text(aes(label=percent, y = percent), size = 4, position=position_dodge(0.9), vjust=-0.5)
  #coord_cartesian(ylim = c(0, 1.5))
  
  pdf(paste0(names(table(plot.data2$cell_type))[i],'.pdf'),width = 5,height = 5)
  print(colp) # 最后打印图片要用print
  dev.off()
  
}
