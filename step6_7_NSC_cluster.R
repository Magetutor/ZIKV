# 设置工作目录
setwd('~/analysis/zika/')
# 清空变量 加载包
# rm(list = ls())
library(Seurat)
library(DESeq2)
library(tidyverse)

# 加载数据
load('zika_all_mthb_filter_SCT_cluster_PC20_23cluster_immune_marker.Rda')
# 单核细胞亚分群
# 挑选T细胞进⾏进⼀步分群
nsc <- row.names(zika_all_mthb_filter_SCT@meta.data)[which(zika_all_mthb_filter_SCT@meta.data$cell_type == "Neural stem cell")]

target.cell <- nsc
length(target.cell)
nsc.data <- subset(zika_all_mthb_filter_SCT, cells = target.cell)
nsc.data
DefaultAssay(nsc.data) <- "RNA"
nsc.data <- SplitObject(nsc.data,
                          split.by = "Group")
nsc.data

# 第一次做过线粒体、核糖体、细胞周期质控了所以直接从CT开始
# SCTransform
# 速度较慢，样本越多越慢
# 循环运⾏
for(i in 1:length(nsc.data)){
  nsc.data[[i]] <- SCTransform(
    nsc.data[[i]],
    variable.features.n = 3000,
    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), # 去除线粒体和细胞周期的影响
    verbose = FALSE)
}

# 数据整合
# 选择⽤于整合的基因
features <- SelectIntegrationFeatures(object.list = nsc.data,
                                      nfeatures = 3000)
# 准备整合
nsc.data <- PrepSCTIntegration(object.list = nsc.data,
                                 anchor.features = features)
names(nsc.data)
# 寻找锚定基因(速度很慢)
AnchorSet <- FindIntegrationAnchors(object.list = nsc.data,
                                    reference = 1,
                                    normalization.method = "SCT",
                                    anchor.features = features
)
# save(AnchorSet, file = "AnchorSet.Rda")
# 整合数据
# load("AnchorSet.Rda")
nsc.data <- IntegrateData(anchorset = AnchorSet,
                            normalization.method = "SCT")
DefaultAssay(nsc.data) <- "integrated"
save(nsc.data, file = "nsc.data.final.Rda")
############降维分析##########
rm(list = ls())
library(Seurat)
library(tidyverse)
library(stringr)
library(patchwork)
load("nsc.data.final.Rda")
dim(nsc.data)
# PCA分析
nsc.data <- RunPCA(object = nsc.data,
                     npcs = 50,
                     rev.pca = FALSE,
                     weight.by.var = TRUE,
                     verbose = TRUE, # 输出特定PC值的特征基因
                     ndims.print = 1:5, # 输出前5个PC值的特征基因
                     nfeatures.print = 30, # 输出每个PC值的30个特征基因
                     reduction.key = "PC_")

ElbowPlot(nsc.data,
          ndims = 50) # 寻找对细胞差异贡献度较⼤的主成分
# 选择30个主成分进⾏后续分析
# t-SNE分析，基于PCA分析结果进⾏⾮线性降维分析
nsc.data <- RunTSNE(nsc.data,
                      dims = 1:30)
# UMAP分析
nsc.data <- RunUMAP(nsc.data,
                      dims = 1:30)
# 展示PCA降维特征及其特征基因的热图
DimHeatmap(object = nsc.data,
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
p1 <- DimPlot(object = nsc.data,
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
p2 <- DimPlot(object = nsc.data,
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
p3 <- DimPlot(object = nsc.data,
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
p4 <- DimPlot(object = nsc.data,
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
p5 <- DimPlot(object = nsc.data,
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
p6 <- DimPlot(object = nsc.data,
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
FeaturePlot(nsc.data,
            reduction = "umap",
            features = "percent.mt")
# 展示特定基因
FeaturePlot(nsc.data,
            reduction = "umap",
            features = "Cxcl2")
save(nsc.data, file = "nsc.data.before_cluster.Rda")
##########聚类分析##########
# KNN聚类法
nsc.data <- FindNeighbors(nsc.data,
                            k.param = 20, # 不要⼩于5
                            dims = 1:30)
nsc.data <- FindClusters(nsc.data,
                           resolution = 0.4,
                           method = "igraph",
                           algorithm = 1,
                           random.seed = 2021)
table(nsc.data@meta.data$seurat_clusters)

PCAPlot(object = nsc.data,
        group.by = "seurat_clusters",
        pt.size = 0.5,
        label = TRUE)
TSNEPlot(object = nsc.data,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)
UMAPPlot(object = nsc.data,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)

DimPlot(nsc.data,
        reduction = "umap",
        split.by = "Group",
        group.by = "seurat_clusters",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()
# 保存一下分好群的数据,可以下次都进来重新聚类
save(nsc.data, file = "nsc.data_SCT_cluster3.Rda")


# 没有数据的话需要load('nsc.data_SCT_cluster8.Rda')
logFCfilter=0.5
adjPvalFilter=0.05
nsc.data.markers <- FindAllMarkers(object = nsc.data,
                                     only.pos = TRUE, #TRUE表示只求高表达的基因
                                     min.pct = 0.25, #表达的最小百分比
                                     logfc.threshold = logFCfilter) # 这一步花了1-2h

# 把marker文件
保存下来
save(nsc.data.markers, file = "nsc.data.markers_cluster3.Rda")

sig.markers = nsc.data.markers[(abs(as.numeric(as.vector(nsc.data.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(nsc.data.markers$p_val_adj))<adjPvalFilter),]
write.table (sig.markers,file="PC30_astro.markers_cluster6.xls",sep="\t",row.names=F,quote=F)
top5 <- nsc.data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- nsc.data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- nsc.data.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table (top20 ,file="PC30_top20_astro.markers_cluster6.xls",sep="\t",row.names=F,quote=F)


#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap.pdf",width=32,height=8)
DoHeatmap(object = nsc.data, features = top10$gene,size = 3,angle = 0) + NoLegend()
dev.off()  # 绘制速度较慢

FeaturePlot(object = nsc.data, features = c("Irak3", "Bcl2a1b"))
FeaturePlot(object = nsc.data, features = c( 'Ier3', 'Fos'),
            reduction = "umap")
VlnPlot(nsc.data,
        features = c('Gap43', 'Fabp7','Selenom'
        ),
        group.by = "seurat_clusters",pt.size = 0)+NoLegend()
VlnPlot(nsc.data,
        features = c('Nrxn1', 'Slc6a11','Zbtb20',
                     'Sfrp1', 'Nid1','Gria4','Elof1', 'Dynlrb2','Tyrobp',
                     'Lgmn'
        ),
        group.by = "seurat_clusters",pt.size = 0)+NoLegend()


# 先给每个群命名一下
####写入细胞类型
Cell_type1 <- c("0" = "Meg",
                "1" = "Ccnd2",
                "2" = "Fcer1g") 

#将定义各个亚型对应的细胞名称，然后写⼊到meta.data⽂件
nsc.data[['cell_type1']] <- unname(Cell_type1[nsc.data@meta.data$seurat_clusters])

#看一下结果
View(nsc.data@meta.data)

#接着可以重新画⼀个UMAP图带上细胞名称
#UMAP
DimPlot(nsc.data,
        reduction = "umap",
        group.by = "cell_type1",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()
#tsne
DimPlot(nsc.data,
        reduction = "tsne",
        group.by = "cell_type1",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()

#按照组别上色看看
DimPlot(nsc.data,
        reduction = "umap",
        split.by = "Group",
        group.by = "cell_type1",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()
#保存一下注释好的数据
save(nsc.data, file = "nsc.data_3cluster_ann.Rda")
load('/home/zhangc/analysis/zika/data/step6_7_astro/nsc.data_3cluster_ann.Rda')

#### 保存细胞亚型 ----
saveRDS(nsc.data, file = "nsc.data_3cluster_ann.rds")

########## 单核细胞亚群的差异基因分析#############
# 批量找各个细胞在V59中的差异基因
dput(unique(nsc.data@meta.data$cell_type1))
table(nsc.data@meta.data$cell_type1)
Cells <-c("Ccnd2", "Meg", "Fcer1g")

diff_V59_nsc <- list()
for (i in 1:length(Cells)) {
  cells <- subset(nsc.data,
                  cell_type1 == Cells[i])
  Idents(cells) <- cells@meta.data$Group
  cells.diff <- FindMarkers(cells,
                            ident.1 = "V_59",
                            ident.2 = "MOCK",
                            test.use = "MAST",
                            min.pct = 0.3, # 过滤掉那些在50%以下细胞中检测到的基因
                            logfc.threshold = log(2) # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
  )
  diff_V59_nsc [[i]] <- cells.diff
}
names(diff_V59_nsc ) <- Cells
save(diff_V59_nsc , file = "diff_V59_nsc_alldiff.Rda")

# 每个细分亚群的功能富集结果
# ######富集分析#######
# 基因ID转换
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

length(diff_V59_nsc)
rownames(diff_V59_nsc[[1]])
ENTREZID <- bitr(rownames(diff_V59_nsc[[1]][diff_V59_nsc[[1]]$p_val_adj < 0.05,]), 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = "org.Mm.eg.db")
Irak3_NK.diff <- list(ENTREZID)

# 批量转换
scRNA.diff_V59_nsc <- list()
for (i in 1:length(diff_V59_nsc)) {
  ENTREZID <- bitr(rownames(diff_V59_nsc[[i]][diff_V59_nsc[[i]]$p_val_adj < 0.05,]), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
  ENTREZID <- ENTREZID$ENTREZID
  scRNA.diff_V59_nsc[[i]] <- ENTREZID
}
names(scRNA.diff_V59_nsc) <- names(diff_V59_nsc)

############## 多个cluster比较 ###################
# KEGG分析
xx <- compareCluster(scRNA.diff_V59_nsc, 
                     fun="enrichKEGG",
                     organism="mmu",
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.05)
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

cowplot::plot_grid(p1,
                   p2, 
                   p3,
                   p4,
                   ncol = 2, 
                   labels = LETTERS[1:4])
ggsave("compareCluster.kegg.pdf")

# GO分析
# all 整体分析会有NA出现，单独分析看看
xx_BP <- compareCluster(scRNA.diff_V59_nsc, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "BP",
                        pvalueCutoff = 0.05,
                        qvalueCutoff= 0.05)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_BP, showCategory = 5) 
#ggsave("compareCluster.GO.pdf")

xx_MF <- compareCluster(scRNA.diff_V59_nsc, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "MF",
                        pvalueCutoff = 0.05,
                        qvalueCutoff= 0.05)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_MF, showCategory = 5) #,label_format = 60
# ggsave("compareCluster.GO.pdf")

xx_CC <- compareCluster(scRNA.diff_V59_nsc, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "CC",
                        pvalueCutoff = 0.05,
                        qvalueCutoff= 0.05) #astroglia_compare_BP
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_CC, showCategory = 5) 
# ggsave("compareCluster.GO.pdf")

###################单独做每个组中差异的KEGG和GO#################
##clusterprofiler加载包
library("clusterProfiler")
library("enrichplot")
library("org.Hs.eg.db")
library("tidyverse")
library('viridis')
library('ggpubr')
#########富集分析##########################
####VFS_Mock_diffgene_分析######
diffgene_flna <- scRNA.diff_V59_nsc[[2]]
diffgene_flna_BP=enrichGO(diffgene_flna,
                          OrgDb="org.Mm.eg.db",   ##物种来源
                          keyType= 'ENTREZID',
                          ont = "BP",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          minGSSize = 2,
                          maxGSSize = 500)
diffgene_flna_BP=setReadable(diffgene_flna_BP, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
V59_Mock_diffgene_flna_BP <- as.data.frame(diffgene_flna_BP@result)
write.csv(V59_Mock_diffgene_flna_BP,'V59_Mock_diffgene_flna_BP.csv')
barplot(diffgene_flna_BP,categorySize="p.adjust", showCategory=10,color = "p.adjust")+
  scale_fill_gradient(low='#2c728e', high='#85d44a')+ #barplot改颜色要用scale_fill_gradient
  theme_classic2() #层级图和柱形图

#########KEGG
diffgene_flna_KE=enrichKEGG(diffgene_flna,
                            organism = "mmu",
                            keyType = "kegg",
                            pvalueCutoff =0.05, 
                            qvalueCutoff =0.05)

#没有结果就把ENTREZ ID写出来 用KOBAS
diffgene_flna_KE=setReadable(diffgene_flna_KE, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
View(diffgene_flna_KE@result)
V59_Mock_diffgene_flna_KE <- as.data.frame(diffgene_flna_KE@result)
write.csv(V59_Mock_diffgene_flna_KE,'V59_Mock_diffgene_unknow_KE.csv')

#柱形图
barplot(diffgene_flna_KE,categorySize="p.adjust", showCategory=10,color = "p.adjust")+
  scale_fill_gradient(low='#2c728e', high='#85d44a')+ #barplot改颜色要用scale_fill_gradient
  theme_classic2() #层级图和柱形图

#V59_Mock_diffgene_Tcf4分析
diffgene_Tcf4 <- scRNA.diff_V59_nsc[[3]]
diffgene_Tcf4_BP=enrichGO(diffgene_Tcf4,
                          OrgDb="org.Mm.eg.db",   ##物种来源
                          keyType= 'ENTREZID',
                          ont = "BP",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          minGSSize = 2,
                          maxGSSize = 500)
diffgene_Tcf4_BP=setReadable(diffgene_Tcf4_BP, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
V59_Mock_diffgene_Tcf4_BP <- as.data.frame(diffgene_Tcf4_BP@result)
write.csv(V59_Mock_diffgene_Tcf4_BP,'V59_Mock_diffgene_Tcf4_BP.csv')
barplot(diffgene_Tcf4_BP,categorySize="p.adjust", showCategory=10,color = "p.adjust")+
  scale_fill_gradient(low='#2c728e', high='#85d44a')+ #barplot改颜色要用scale_fill_gradient
  theme_classic2() #层级图和柱形图

#########KEGG
diffgene_Tcf4_KE=enrichKEGG(diffgene_Tcf4,
                            organism = "mmu",
                            keyType = "kegg",
                            pvalueCutoff =0.99, 
                            qvalueCutoff =0.99)

#没有结果就把ENTREZ ID写出来 用KOBAS
diffgene_Tcf4_KE=setReadable(diffgene_Tcf4_KE, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
View(diffgene_Tcf4_KE@result)
V59_Mock_diffgene_Tcf4_KE <- as.data.frame(diffgene_Tcf4_KE@result)
write.csv(V59_Mock_diffgene_Tcf4_KE,'V59_Mock_diffgene_Tcf4_KE.csv')

#柱形图
barplot(diffgene_Tcf4_KE,categorySize="pvalue", showCategory=10,color = "pvalue")+
  scale_fill_gradient(low='#2c728e', high='#85d44a')+ #barplot改颜色要用scale_fill_gradient
  theme_classic2() #层级图和柱形图


####diff_VFS_diffgene_unkonw分析######
diffgene_unknow <- scRNA.diff_V59_nsc[[1]]
diffgene_unknow_BP=enrichGO(diffgene_unknow,
                            OrgDb="org.Mm.eg.db",   ##物种来源
                            keyType= 'ENTREZID',
                            ont = "BP",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            minGSSize = 2,
                            maxGSSize = 500)
diffgene_unknow_BP=setReadable(diffgene_unknow_BP, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
VFS_Mock_diffgene_unknow_BP <- as.data.frame(diffgene_unknow_BP@result)
write.csv(VFS_Mock_diffgene_unknow_BP,'VFS_Mock_diffgene_unknow_BP.csv')
barplot(diffgene_unknow_BP,categorySize="p.adjust", showCategory=10,color = "p.adjust")+
  scale_fill_gradient(low='#2c728e', high='#85d44a')+ #barplot改颜色要用scale_fill_gradient
  theme_classic2() #层级图和柱形图

#########KEGG
diffgene_unknow_KE=enrichKEGG(diffgene_unknow,
                              organism = "mmu",
                              keyType = "kegg",
                              pvalueCutoff =0.99, 
                              qvalueCutoff =0.99)

#没有结果就把ENTREZ ID写出来 用KOBAS
diffgene_unknow_KE=setReadable(diffgene_unknow_KE, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
View(diffgene_unknow_KE@result)
VFS_Mock_diffgene_unknow_KE <- as.data.frame(diffgene_unknow_KE@result)
write.csv(VFS_Mock_diffgene_unknow_KE,'VFS_Mock_diffgene_unknow_KE.csv')

#柱形图
barplot(diffgene_unknow_KE,categorySize="p.adjust", showCategory=10,color = "p.adjust")+
  scale_fill_gradient(low='#2c728e', high='#85d44a')+ #barplot改颜色要用scale_fill_gradient
  theme_classic2() #层级图和柱形图

####diff_VFS_diffgene_flna分析######
diffgene_flna <- scRNA.diff_V59_nsc[[2]]
diffgene_flna_BP=enrichGO(diffgene_flna,
                          OrgDb="org.Mm.eg.db",   ##物种来源
                          keyType= 'ENTREZID',
                          ont = "BP",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          minGSSize = 2,
                          maxGSSize = 500)
diffgene_flna_BP=setReadable(diffgene_flna_BP, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
VFS_Mock_diffgene_flna_BP <- as.data.frame(diffgene_flna_BP@result)
write.csv(VFS_Mock_diffgene_flna_BP,'VFS_Mock_diffgene_flna_BP.csv')
barplot(diffgene_flna_BP,categorySize="p.adjust", showCategory=10,color = "p.adjust")+
  scale_fill_gradient(low='#2c728e', high='#85d44a')+ #barplot改颜色要用scale_fill_gradient
  theme_classic2() #层级图和柱形图

#########KEGG
diffgene_flna_KE=enrichKEGG(diffgene_flna,
                            organism = "mmu",
                            keyType = "kegg",
                            pvalueCutoff =0.99, 
                            qvalueCutoff =0.99)

#没有结果就把ENTREZ ID写出来 用KOBAS
diffgene_flna_KE=setReadable(diffgene_flna_KE, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
View(diffgene_flna_KE@result)
VFS_Mock_diffgene_flna_KE <- as.data.frame(diffgene_flna_KE@result)
write.csv(VFS_Mock_diffgene_flna_KE,'VFS_Mock_diffgene_flna_KE.csv')

#柱形图
barplot(diffgene_flna_KE,categorySize="pvalue", showCategory=10,color = "pvalue")+
  scale_fill_gradient(low='#2c728e', high='#85d44a')+ #barplot改颜色要用scale_fill_gradient
  theme_classic2() #层级图和柱形图

#V59_Mock_diffgene_Tcf4分析
diffgene_Tcf4 <- scRNA.diff_V59_nsc[[3]]
diffgene_Tcf4_BP=enrichGO(diffgene_Tcf4,
                          OrgDb="org.Mm.eg.db",   ##物种来源
                          keyType= 'ENTREZID',
                          ont = "BP",
                          pvalueCutoff  = 0.99,
                          qvalueCutoff  = 0.99,
                          minGSSize = 10,
                          maxGSSize = 500)
diffgene_Tcf4_BP=setReadable(diffgene_Tcf4_BP, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
VFS_Mock_diffgene_Tcf4_BP <- as.data.frame(diffgene_Tcf4_BP@result)
write.csv(VFS_Mock_diffgene_Tcf4_BP,'VFS_Mock_diffgene_Tcf4_BP.csv')
barplot(diffgene_Tcf4_BP,categorySize="p.adjust", showCategory=10,color = "p.adjust")+
  scale_fill_gradient(low='#2c728e', high='#85d44a')+ #barplot改颜色要用scale_fill_gradient
  theme_classic2() #层级图和柱形图

#########KEGG
diffgene_Tcf4_KE=enrichKEGG(diffgene_Tcf4,
                            organism = "mmu",
                            keyType = "kegg",
                            pvalueCutoff =0.99, 
                            qvalueCutoff =0.99)

#没有结果就把ENTREZ ID写出来 用KOBAS
diffgene_Tcf4_KE=setReadable(diffgene_Tcf4_KE, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
View(diffgene_Tcf4_KE@result)
VFS_Mock_diffgene_Tcf4_KE <- as.data.frame(diffgene_Tcf4_KE@result)
write.csv(VFS_Mock_diffgene_Tcf4_KE,'VFS_Mock_diffgene_Tcf4_KE.csv')

#柱形图
barplot(diffgene_Tcf4_KE,categorySize="pvalue", showCategory=10,color = "pvalue")+
  scale_fill_gradient(low='#2c728e', high='#85d44a')+ #barplot改颜色要用scale_fill_gradient
  theme_classic2() #层级图和柱形图


##########循环做一下Gsea吧###############
## 解释原理
library(clusterProfiler)

#获得基因列表
gene <- rownames(diff_V59_nsc[[5]])
## 转换
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=diff_V59_nsc[[5]]$avg_log2FC,
                      SYMBOL = rownames(diff_V59_nsc[[5]]))
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
names(diff_V59_nsc)

write.csv(gsea_result,'Nfib_NK_gsea_result.csv')

#我们先把结果提取出来，方便查看
library(enrichplot)
gseaplot2(gsea, "mmu04142",color = 'steelblue',  title = "Lysosome")
gseaplot2(gsea, "mmu04015",color = 'steelblue',  title = "Rap1 signaling pathway")
gseaplot2(gsea, "hsa04727",color = 'steelblue',  title = "GABAergic synapse")


load('nsc.data_3cluster_ann.Rda')

# 批量找各个细胞在VFS中的差异基因
dput(unique(nsc.data@meta.data$cell_type1))
table(nsc.data@meta.data$cell_type1)
Cells <- c("Ccnd2", "Meg", "Fcer1g")
diff_VFS_nsc <- list()
for (i in 1:length(Cells)) {
  cells <- subset(nsc.data,
                  cell_type1 == Cells[i])
  Idents(cells) <- cells@meta.data$Group
  cells.diff <- FindMarkers(cells,
                            ident.1 = "V_FS",
                            ident.2 = "MOCK",
                            test.use = "MAST",
                            min.pct = 0.3, # 过滤掉那些在50%以下细胞中检测到的基因
                            logfc.threshold = log(2) # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
  )
  diff_VFS_nsc [[i]] <- cells.diff
}
names(diff_VFS_nsc ) <- Cells
save(diff_VFS_nsc , file = "diff_VFS_nsc_alldiff.Rda")

# 每个细分亚群的功能富集结果
# ######富集分析#######
# 基因ID转换
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

length(diff_VFS_nsc)
rownames(diff_VFS_nsc[[1]])
ENTREZID <- bitr(rownames(diff_VFS_nsc[[1]][diff_VFS_nsc[[1]]$p_val_adj < 0.05,]), 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = "org.Mm.eg.db")
Irak3_NK.diff <- list(ENTREZID)

# 批量转换
scRNA.diff_VFS_nsc <- list()
for (i in 1:length(diff_VFS_nsc)) {
  ENTREZID <- bitr(rownames(diff_VFS_nsc[[i]][diff_VFS_nsc[[i]]$p_val_adj < 0.05,]), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
  ENTREZID <- ENTREZID$ENTREZID
  scRNA.diff_VFS_nsc[[i]] <- ENTREZID
}
names(scRNA.diff_VFS_nsc) <- names(diff_VFS_nsc)

############## 多个cluster比较 ###################
# KEGG分析
xx <- compareCluster(scRNA.diff_VFS_nsc, 
                     fun="enrichKEGG",
                     organism="mmu",
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.05)
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

cowplot::plot_grid(p1,
                   p2, 
                   p3,
                   p4,
                   ncol = 2, 
                   labels = LETTERS[1:4])
ggsave("compareCluster.kegg.pdf")

# GO分析
# all 整体分析会有NA出现，单独分析看看
xx_BP <- compareCluster(scRNA.diff_VFS_nsc, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "BP",
                        pvalueCutoff = 0.05,
                        qvalueCutoff= 0.05)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_BP, showCategory = 5) #,label_format = 60
#ggsave("compareCluster.GO.pdf")

xx_MF <- compareCluster(scRNA.diff_VFS_nsc, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "MF",
                        pvalueCutoff = 0.05,
                        qvalueCutoff= 0.05)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_MF, showCategory = 5) 
# ggsave("compareCluster.GO.pdf")

xx_CC <- compareCluster(scRNA.diff_VFS_nsc, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "CC",
                        pvalueCutoff = 0.05,
                        qvalueCutoff= 0.05) #astroglia_compare_BP
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_CC, showCategory = 5) 
# ggsave("compareCluster.GO.pdf")

################细胞亚型数量及比例##############
plot.data <- nsc.data@meta.data[c("Group", "cell_type1")]
head(plot.data)
plot.data1 <- plot.data %>%
  group_by(Group, (plot.data)) %>%
  dplyr::summarise(num = n())
head(plot.data1)
plot.data1$Group <- factor(plot.data1$Group, levels = c("MOCK","V_FS","V_59"))
ggplot(plot.data1, aes(x = cell_type1,
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
for (i in 1:length(names(table(plot.data2$cell_type1)))) {
  plot_n <- dplyr::filter(plot.data2, cell_type1 == names(table(plot.data2$cell_type1))[i])
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
