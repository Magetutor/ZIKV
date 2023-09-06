# 加载数据 按照五类聚类
load('micro.data.before_cluster.Rda')

##########聚类分析##########
# KNN聚类法
micro.data <- FindNeighbors(micro.data,
                            k.param = 20, # 不要⼩于5
                            dims = 1:30)
micro.data <- FindClusters(micro.data,
                           resolution = 0.1,
                           method = "igraph",
                           algorithm = 1,
                           random.seed = 2021)
table(micro.data@meta.data$seurat_clusters)

PCAPlot(object = micro.data,
        group.by = "seurat_clusters",
        pt.size = 0.5,
        label = TRUE)
TSNEPlot(object = micro.data,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)
UMAPPlot(object = micro.data,
         group.by = "seurat_clusters",
         pt.size = 0.5,
         label = TRUE)


# 保存一下分好群的数据,可以下次都进来重新聚类
save(micro.data, file = "micro.data_SCT_cluster8.Rda")


# 没有数据的话需要load('micro.data_SCT_cluster8.Rda')
logFCfilter=0.5
adjPvalFilter=0.05
micro.data.markers <- FindAllMarkers(object = micro.data,
                                     only.pos = TRUE, #TRUE表示只求高表达的基因
                                     min.pct = 0.25, #表达的最小百分比
                                     logfc.threshold = logFCfilter) # 这一步花了1-2h

# 把marker文件
# 保存下来
save(micro.data.markers, file = "micro.data.markers_cluster5.Rda")

sig.markers = micro.data.markers[(abs(as.numeric(as.vector(micro.data.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(micro.data.markers$p_val_adj))<adjPvalFilter),]
write.table (sig.markers,file="PC30_micro.markers_cluster5.xls",sep="\t",row.names=F,quote=F)
top5 <- micro.data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- micro.data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- micro.data.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table (top20 ,file="PC30_top20_micro.markers_cluster5.xls",sep="\t",row.names=F,quote=F)


#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap.pdf",width=32,height=8)
DoHeatmap(object = micro.data, features = top10$gene,size = 3,angle = 0) + NoLegend()
dev.off()  # 绘制速度较慢

FeaturePlot(object = micro.data, features = c("Irak3", "Bcl2a1b"))
FeaturePlot(object = micro.data, features = c( 'Ier3', 'Fos'),
            reduction = "umap")
VlnPlot(micro.data,
        features = c('Slc40a1', 'Dab2','Adam33'
        ),
        group.by = "seurat_clusters",pt.size = 0)+NoLegend()
VlnPlot(micro.data,
        features = c('Rbpms', 'Flna','Adgre5',
                     'S100a4', 'Ifitm6','Smpdl3a'
        ),
        group.by = "seurat_clusters",pt.size = 0)+NoLegend()
VlnPlot(micro.data,
        features = c('Gpr27', 'Itgb8','Mpped2',
                     'Tcf4', 'Gm4258','Pax6'
        ),
        group.by = "seurat_clusters",pt.size = 0)+NoLegend()

# 先给每个群命名一下
####写入细胞类型
Cell_type1 <- c("0" = "Syngr1",
                "1" = "Cd44",
                "2" = "Ccr2",
                "3" = "Hist1h2ai",
                "4" = "Ptn")

#将定义各个亚型对应的细胞名称，然后写⼊到meta.data⽂件
micro.data[['cell_type1']] <- unname(Cell_type1[micro.data@meta.data$seurat_clusters])

#看一下结果
View(micro.data@meta.data)

#接着可以重新画⼀个UMAP图带上细胞名称
#UMAP
DimPlot(micro.data,
        reduction = "umap",
        group.by = "cell_type1",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()
#tsne
DimPlot(micro.data,
        reduction = "tsne",
        group.by = "cell_type1",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()

#按照组别上色看看
DimPlot(micro.data,
        reduction = "umap",
        split.by = "Group",
        group.by = "cell_type",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()

DimPlot(micro.data,
        reduction = "umap",
        split.by = "Group",
        group.by = "cell_type1",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()

#保存一下注释好的数据
save(micro.data, file = "micro.data_5cluster_ann.Rda")
#### 保存细胞亚型 ----
saveRDS(micro.data, file = "micro_5cluster_ann.rds")

########## 单核细胞亚群的差异基因分析#############
# 批量找各个细胞在V59中的差异基因
dput(unique(micro.data@meta.data$cell_type1))
table(micro.data@meta.data$cell_type1)
Cells <- c("Syngr1", "Cd44", "Ptn", "Hist1h2ai", "Ccr2")
diff_V59_micro <- list()
for (i in 1:length(Cells)) {
  cells <- subset(micro.data,
                  cell_type1 == Cells[i])
  Idents(cells) <- cells@meta.data$Group
  cells.diff <- FindMarkers(cells,
                            ident.1 = "V_59",
                            ident.2 = "MOCK",
                            test.use = "MAST",
                            min.pct = 0.3, # 过滤掉那些在50%以下细胞中检测到的基因
                            logfc.threshold = log(2) # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
  )
  diff_V59_micro [[i]] <- cells.diff
}
names(diff_V59_micro ) <- Cells
save(diff_V59_micro , file = "diff_V59_micro_alldiff_cluster5.Rda")

# 每个细分亚群的功能富集结果
# ######富集分析#######
# 基因ID转换
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

length(diff_V59_micro)
rownames(diff_V59_micro[[1]])
ENTREZID <- bitr(rownames(diff_V59_micro[[1]][diff_V59_micro[[1]]$p_val_adj < 0.05,]), 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = "org.Mm.eg.db")
Irak3_NK.diff <- list(ENTREZID)

# 批量转换
scRNA.diff_V59_micro <- list()
for (i in 1:length(diff_V59_micro)) {
  ENTREZID <- bitr(rownames(diff_V59_micro[[i]][diff_V59_micro[[i]]$p_val_adj < 0.05,]), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
  ENTREZID <- ENTREZID$ENTREZID
  scRNA.diff_V59_micro[[i]] <- ENTREZID
}
names(scRNA.diff_V59_micro) <- names(diff_V59_micro)

############## 多个cluster比较 ###################
# KEGG分析
xx <- compareCluster(scRNA.diff_V59_micro, 
                     fun="enrichKEGG",
                     organism="mmu",
                     pvalueCutoff= 0.01,
                     qvalueCutoff= 0.01)
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
xx_BP <- compareCluster(scRNA.diff_V59_micro, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "BP",
                        pvalueCutoff = 0.00001,
                        qvalueCutoff= 0.00001)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_BP, showCategory = 5) 
ggsave("compareCluster.GO.pdf")

xx_MF <- compareCluster(scRNA.diff_V59_micro, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "MF",
                        pvalueCutoff = 0.01,
                        qvalueCutoff= 0.01)
#save(xx, file = "compareCluster.GO.Rda")

dotplot(xx_MF, showCategory = 5) 
# ggsave("compareCluster.GO.pdf")

xx_CC <- compareCluster(scRNA.diff_V59_micro, 
                        fun = "enrichGO",
                        OrgDb = "org.Mm.eg.db",
                        ont = "CC",
                        pvalueCutoff = 0.001,
                        qvalueCutoff= 0.001) #microglia_compare_BP
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
####V59_Mock_diffgene_flna分析######
diffgene_flna <- scRNA.diff_V59_micro[[2]]
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
diffgene_Tcf4 <- scRNA.diff_V59_micro[[3]]
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
diffgene_unknow <- scRNA.diff_V59_micro[[1]]
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
diffgene_flna <- scRNA.diff_V59_micro[[2]]
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
diffgene_Tcf4 <- scRNA.diff_V59_micro[[3]]
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
gene <- rownames(diff_V59_micro[[5]])
## 转换
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=diff_V59_micro[[5]]$avg_log2FC,
                      SYMBOL = rownames(diff_V59_micro[[5]]))
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
names(diff_V59_micro)

write.csv(gsea_result,'Nfib_NK_gsea_result.csv')

#我们先把结果提取出来，方便查看
library(enrichplot)
gseaplot2(gsea, "mmu04142",color = 'steelblue',  title = "Lysosome")
gseaplot2(gsea, "mmu04015",color = 'steelblue',  title = "Rap1 signaling pathway")
gseaplot2(gsea, "hsa04727",color = 'steelblue',  title = "GABAergic synapse")
