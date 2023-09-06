# 设置工作目录
setwd('~/analysis/zika/')
# 加载包
library(tidyverse)
library(reshape2)
library(patchwork)
library(Seurat)
# 导入数据
## load('zika_all_mthb_filter_SCT_cluster_PC20_23cluster_immune_marker.Rda')

##########差异分析########
rm(list = ls())
# 提取特定细胞单细胞数据
# B.cells <- subset(alldata, 
#                   idents = "B Cells")
Macrophage.cells <- subset(zika_all_mthb_filter_SCT,
                  cell_type == "Macrophage")
Idents(Macrophage.cells)
Idents(Macrophage.cells) <- Macrophage.cells@meta.data$Group
Idents(Macrophage.cells)

Macrophage.diff.V59 <- FindMarkers(Macrophage.cells,
                            ident.1 = "V_59",
                            ident.2 = "MOCK")

dput(rownames(head(Macrophage.diff.V59_1)))

VlnPlot(Macrophage.cells, 
        features = c("Gm5150", "Arg1", "Cebpe", "Trem1", "Tnip3", "Ptgs2"),
        split.by = "Group")

# 直接筛选结果
Macrophage.diff.V59_1 <- FindMarkers(Macrophage.cells,
                             ident.1 = "V_59",
                             ident.2 = "MOCK",
                             min.pct = 0.5, # 过滤掉那些在50%以下细胞中检测到的基因
                             logfc.threshold = log(2), # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
)
write.table (Macrophage.diff.V59_1 ,file="Macrophage.diff.V59_1.xls",sep="\t",row.names=F,quote=F)

Macrophage.diff.V_FS <- FindMarkers(Macrophage.cells,
                                     ident.1 = "V_FS",
                                     ident.2 = "MOCK",
                                     min.pct = 0.5, # 过滤掉那些在50%以下细胞中检测到的基因
                                     logfc.threshold = log(2), # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
)
dput(rownames(head(Macrophage.diff.V_FS)))

DotPlot(Macrophage.cells, 
        features = c("Gm15987", "Ltb4r1", "Gm5150", "Trem1", "Ifitm6", "Ptgs2"))
        # split.by = "Group") 
write.table (Macrophage.diff.V_FS ,file="Macrophage.diff.V_FS.xls",sep="\t",row.names=F,quote=F)

# 其他算法进行差异分析(需要提前安装对应R包)
# BiocManager::install("MAST")

B.cells.diff2 <- FindMarkers(B.cells,
                             ident.1 = "Befort_Treat",
                             ident.2 = "After_Treat",
                             test.use = "MAST",
                             min.pct = 0.5, # 过滤掉那些在50%以下细胞中检测到的基因
                             logfc.threshold = log(2) # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
)

B.cells.diff3 <- FindMarkers(B.cells,
                             ident.1 = "Befort_Treat",
                             ident.2 = "After_Treat",
                             assay = "RNA",
                             test.use = "DESeq2",
                             min.pct = 0.5, # 过滤掉那些在50%以下细胞中检测到的基因
                             logfc.threshold = log(2) # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
)