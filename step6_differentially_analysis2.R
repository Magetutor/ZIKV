# 设置工作目录
setwd('~/analysis/zika/')
# 加载包
library(Seurat)
library(DESeq2)
library(tidyverse)
####差异分析######
## load('zika_all_mthb_filter_SCT_cluster_PC20_23cluster_immune_marker.Rda')

# 单个细胞在组之间差异基因，可能是很多细胞在组之间都存在，所以批量找差异基因
# 找某个细胞类型中特已存在的基因可以作为研究对象

# 写成循环
dput(unique(zika_all_mthb_filter_SCT@meta.data$cell_type))
table(zika_all_mthb_filter_SCT@meta.data$cell_type)
Cells <- c("Macrophage", "erythroid precursor cells", "Neural stem cell", 
         "NK cells", "Microglia_Ctss high", "Monocyte_Lyz2", "Neuron_Syp hig", 
         "Monocytes & Neutrophils", "Endothelial_cell_Egfl7", "Neuron_Nbea", 
         "Oligodendrocyte", "Neutrophils", "B cell", "B cell_Cd74", "Microglia_Csf1r", 
         "Astrocyte_Slc1a3", "OPC", "Astrocyte_Fabp7", "Neuron_Sst", "Excitatory_neuron", 
         "Endothelial_cell_Mgp", "Microglia_Siglech", "Ependymal_cell_Ttr", 
         "Granulocyte_Cyp11a1")

# 批量找各个细胞在V59中的差异基因
diff_V59 <- list()
for (i in 1:length(Cells)) {
  cells <- subset(zika_all_mthb_filter_SCT,
                  cell_type == Cells[i])
  Idents(cells) <- cells@meta.data$Group
  cells.diff <- FindMarkers(cells,
                            ident.1 = "V_59",
                            ident.2 = "MOCK",
                            test.use = "MAST",
                            min.pct = 0.5, # 过滤掉那些在50%以下细胞中检测到的基因
                            logfc.threshold = log(2) # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
  )
  diff_V59[[i]] <- cells.diff
}
names(diff_V59) <- Cells
save(diff_V59, file = "diff_V59_alldiff.Rda")
# 批量取交集
# load(file = "diff_V59_alldiff.Rda")
intersect.Macrophage <- intersect(rownames(diff_V59[["Macrophage"]]), rownames(diff_V59[["erythroid precursor cells"]]))
intersect.mat <- diff_V59[["Macrophage"]][intersect.Macrophage, ]

Cells <- c("Macrophage", "erythroid precursor cells", "Neural stem cell", 
           "NK cells", "Microglia_Ctss high", "Monocyte_Lyz2", "Neuron_Syp hig", 
           "Monocytes & Neutrophils", "Endothelial_cell_Egfl7", "Neuron_Nbea", 
           "Oligodendrocyte", "Neutrophils", "B cell", "B cell_Cd74", "Microglia_Csf1r", 
           "Astrocyte_Slc1a3", "OPC", "Astrocyte_Fabp7", "Neuron_Sst", "Excitatory_neuron", 
           "Endothelial_cell_Mgp", "Microglia_Siglech", "Ependymal_cell_Ttr", 
           "Granulocyte_Cyp11a1")

intersect.mat.all <- data.frame()
for (i in 1:length(Cells)) {
  intersect.Macrophage <- intersect(rownames(diff_V59[["Macrophage"]]), rownames(diff_V59[[i]]))
  intersect.mat <- diff_V59[["Macrophage"]][intersect.Macrophage, ]
  
  # 数据合并
  if(nrow(intersect.mat.all) == 0){
    intersect.mat.all <- intersect.mat
  }else{
    intersect.mat.all <- rbind(intersect.mat.all, intersect.mat)
  }
}

# 去除重复行
intersect.mat.all1 <- distinct(intersect.mat.all)

# 所以这种办法不适用于所有数据，只适合用来找关键分子
save(diff, file = "alldiff.Rda")

######富集分析#######
# 基因ID转换
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

length(diff_V59)
rownames(diff_V59[[1]])
ENTREZID <- bitr(rownames(diff_V59[[1]]), 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = "org.Mm.eg.db")
Macrophage.diff <- list(ENTREZID)

# 批量转换
scRNA.diff_V59 <- list()
for (i in 1:length(diff_V59)) {
  ENTREZID <- bitr(rownames(diff_V59[[i]]), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = "org.Mm.eg.db")
  ENTREZID <- ENTREZID$ENTREZID
  scRNA.diff_V59[[i]] <- ENTREZID
}
names(scRNA.diff_V59) <- names(diff_V59)

# 多个cluster比较
# KEGG分析
xx <- compareCluster(scRNA.diff_V59, 
                     fun="enrichKEGG",
                     organism="mmu",
                     pvalueCutoff=0.999)
xx <- pairwise_termsim(xx) # 获取相似性矩阵

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
xx <- compareCluster(scRNA.diff_V59, 
                     fun = "enrichGO",
                     OrgDb = "org.Mm.eg.db",
                     ont = "ALL",
                     pvalueCutoff = 0.999)
save(xx, file = "compareCluster.GO.Rda")

dotplot(xx, showCategory = 5) +
  facet_grid(.~ONTOLOGY = 'BP') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15,  vjust = 1, hjust = 1, angle = 45), # color = "navy",
        axis.text.y = element_text(size = 12,  vjust = 1, hjust = 1, angle = -8)) # color = "darkred",
ggsave("compareCluster.GO.pdf")

dotplot(xx, showCategory = 5) +
  facet_grid(.~Cluster) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, color = "navy", vjust = 1, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 12, color = "darkred", vjust = 1, hjust = 1, angle = -8))
ggsave("compareCluster.GO2.pdf", width = 35, height = 10)

# 单个差异基因富集分析
# 以Macrophage细胞为例
gsym.fc <- diff_V59$'Macrophage' %>%
  dplyr::select(2) %>%
  rownames_to_column("gene")
gsym.id <- bitr(gsym.fc$gene, 
                fromType = "SYMBOL", 
                toType = "ENTREZID",
                OrgDb = "org.Mm.eg.db") 
head(gsym.id)

gsym.dat <- gsym.fc %>%
  dplyr::select(SYMBOL = gene,
                logFC = avg_log2FC) %>%
  inner_join(gsym.id, "SYMBOL") %>%
  na.omit() 
write_csv(gsym.dat, "Macrophage.annotation.csv")

##############GO/KEGG/GSEA富集分析##############
# 分组GO和KEGG
input <- read.csv("Macrophage.annotation.csv", sep = ",", header = T, check.names = F)
head(input)
geneList <- input[,2] # 构建list
names(geneList) = as.character(input[,3])
head(geneList)
gene <- names(geneList)[abs(geneList) > 1]
head(gene)
geneList1 <- sort(geneList, decreasing = TRUE)

mydf <- data.frame(Entrez = names(geneList1), FC = geneList1)

mydf$group <- "Up"
mydf$group[mydf$FC < 0] <- "Down"
# mydf$othergroup <- "A"
# mydf$othergroup[abs(mydf$FC) > 1] <- "B"
table(mydf$group)
# table(mydf$othergroup)
# KEGG分析
formula_res <- compareCluster(Entrez ~ group,
                              data = mydf, 
                              fun = "enrichKEGG",
                              organism = "mmu",
                              pvalueCutoff = 0.99)

head(formula_res)
dotplot(formula_res, x = "group", showCategory = 10)
dotplot(formula_res, showCategory = 10) + facet_grid(~ group)

formula_res <- pairwise_termsim(formula_res) # 获取相似性矩阵

p1 <- emapplot(formula_res)
p2 <- emapplot(formula_res, 
               legend_n = 2) 
p3 <- emapplot(formula_res, 
               pie = "count")
p4 <- emapplot(formula_res, 
               pie = "count", 
               cex_category = 1.5, 
               layout="kk")

cowplot::plot_grid(p1,
                   p2, 
                   p3,
                   p4,
                   ncol = 2, 
                   labels = LETTERS[1:4])
# 与前面的多细胞网络图不同，这里没有共享通路

# GO分析
formula_res <- compareCluster(Entrez ~ group,
                              data = mydf, 
                              fun = "enrichGO",
                              OrgDb = "org.Mm.eg.db",
                              ont = "ALL",
                              pvalueCutoff = 0.999)
head(formula_res)
dotplot(formula_res, x = "group", showCategory = 10)
dotplot(formula_res, showCategory = 30) + facet_grid(~ group)
dotplot(formula_res, x = "group", showCategory = 10) + facet_grid(~ ONTOLOGY)
dotplot(formula_res, x = "ONTOLOGY", showCategory = 10) + facet_grid(~ group)

# GO/KEGG/GSEA富集分析
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(GOplot)
library(ggplot2)
library(ggnewscale)
library(enrichplot)
library(ggupset)
library(org.Mm.eg.db)

input <- read.csv("Macrophage.annotation.csv", sep = ",", header = T, check.names = F)
head(input)
geneList <- input[,2] # 构建list
names(geneList) = as.character(input[,3])
head(geneList)
gene <- names(geneList)[abs(geneList) > 1]
head(gene)
geneList1 <- sort(geneList, decreasing = TRUE)

# GO分析
ego <- enrichGO(gene = gene,
                OrgDb = "org.Mm.eg.db",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.99,
                qvalueCutoff = 0.99,
                readable = TRUE)
dim(ego)
# 也可以不用readable = TRUE，用下面这个代码
ego2 <- setReadable(ego, OrgDb = org.Hs.eg.db)

# simplify函数不支持all的富集分析结果
ego_BP <- enrichGO(gene = gene,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.99,
                   qvalueCutoff = 0.99,
                   readable = TRUE)
dim(ego_BP)
# 横轴为GeneRatio，代表该GO term下富集到的基因个数占列表基因总数的比例
# 纵轴为富集到的GO Terms的描述信息，showCategory指定展示的GO Terms的个数
p1 <- dotplot(ego_BP, showCategory = 10) + 
  ggtitle("dotplot for GO")

library(GOSemSim)
egosimp <- simplify(ego_BP,
                    cutoff = 0.5,
                    by = "p.adjust",
                    select_fun = min)
dim(egosimp)
p2 <- dotplot(egosimp, showCategory = 10) + 
  ggtitle("dotplot for GO Simplify")

# GSEA_GO分析
geneList1 <- sort(geneList, decreasing = TRUE)
ego3 <- gseGO(geneList     = geneList1,
              OrgDb        = org.Mm.eg.db,
              ont          = "ALL",
              pvalueCutoff = 0.999,
              verbose      = FALSE)
dim(ego3)

p3 <- dotplot(ego3, showCategory = 10, color = 'pvalue') + 
  ggtitle("dotplot for GSEA")

plot_grid(p1, p2, p3, nrow = 3,
          rel_heights = c(2,2),
          rel_widths = c(4,2))

# KEGG富集分析
kk <- enrichKEGG(gene         = gene,
                 organism     = 'mmu',
                 pvalueCutoff = 0.999,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.999)
head(kk)
kk <- setReadable(kk, 'org.Mm.eg.db', 'ENTREZID')
dotplot(kk , showCategory = 10)
# GSEA_KEGG分析
kk2 <- gseKEGG(geneList     = geneList1,
               organism     = 'mmu',
               pvalueCutoff = 0.99,
               verbose      = FALSE,
               pAdjustMethod = "BH")
head(kk2)
kk2 <- setReadable(kk2, 'org.Mm.eg.db', 'ENTREZID')
dotplot(kk2 , showCategory = 10)
gseaplot2(kk2,1:4,  #ES_geom = "dot",  base_size = 20, #按照第一个作图         
          pvalue_table = T) 
# 可视化
library(enrichplot)
# 横轴为基因个数，纵轴为富集到的GO Terms的描述信息
# 颜色对应p.adjust值，红色p值小，蓝色p值大
# showCategory指定展示的GO Terms的个数，默认为10，即p.adjust最小的10个
barplot(ego, showCategory=20) 

# 使用all的GO分析
dotplot(ego,
        title = "Top5 GO terms",
        showCategory = 5, 
        split = "ONTOLOGY")+ 
  facet_grid(ONTOLOGY~.,
             scale="free")

# 网络图
p1 <- cnetplot(kk, 
               foldChange = geneList)
p2 <- cnetplot(kk, 
               categorySize = "pvalue", 
               foldChange = geneList)
p3 <- cnetplot(kk,
               foldChange = geneList,
               circular = TRUE, 
               colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3,
                   nrow = 3, 
                   labels = LETTERS[1:3], 
                   rel_widths = c(.8, .8, 1.2))
#对于基因和富集的GO terms之间的对应关系进行展示
#图中灰色的点代表基因，黄色的点代表富集到的GO terms
#如果一个基因位于一个GO Terms下，则将该基因与GO连线
#黄色节点的大小对应富集到的基因个数，默认画top5富集到的GO terms
p1 <- cnetplot(kk, 
               node_label = "category") 
p2 <- cnetplot(kk,
               node_label = "gene") 
p3 <- cnetplot(kk, 
               node_label = "all") 
p4 <- cnetplot(kk, 
               node_label = "none") 
cowplot::plot_grid(p1, p2, p3, p4, 
                   ncol = 2, 
                   labels = LETTERS[1:4])

p1 <- heatplot(kk2)
p2 <- heatplot(kk2,
               foldChange = geneList)
cowplot::plot_grid(p1, p2,
                   ncol = 1, 
                   labels = LETTERS[1:2])

upsetplot(ego)
upsetplot(ego3) 

ridgeplot(ego3)

gseaplot(kk2, geneSetID = 1, title = kk2$Description[1])

gseaplot2(kk2, geneSetID = 1, title = kk2$Description[1])

gseaplot2(ego3, geneSetID = 1:3)
gseaplot2(ego3, geneSetID = 1:3, pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"),
          ES_geom = "dot")

p1 <- gseaplot2(ego3, geneSetID = 1:3, subplots = 1)
p2 <- gseaplot2(ego3, geneSetID = 1:3, subplots = 1:2)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])


# 富集分析的问题之一是寻找进一步研究的方向。
# 我们可以使用pmcplot函数，根据PubMed Central的查询结果绘制出版物趋势的数量/比例
library(europepmc)
terms <- ego3$Description[1:5]
p <- pmcplot(terms,
             2010:2021)
p2 <- pmcplot(terms,
              2010:2021,
              proportion = FALSE)
plot_grid(p, p2, nrow = 2)

