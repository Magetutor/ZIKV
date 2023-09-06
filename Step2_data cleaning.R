###数据的预处理
setwd('~/analysis/zika/')
rm(list = ls())
#如果数据重新加载，进行Rda的load即可
load(file = 'zika_all_group.Rda')

#### 1、数据的质控 1：去除线粒体基因和极端值##
##寻找线粒体基因
library(tidyverse)
library(Seurat)
#查看线粒体基因
mito.genes <- str_subset(string = rownames(zika_all),
                         pattern = "mt-") #mouse mitochondrial genes are flagged with ^mt-, rather than ^Mt-
mito.genes 
#查看核糖体基因
ribo.genes <-  rownames(zika_all)[grep("^Rp[sl]", rownames(zika_all),ignore.case = T)]
ribo.genes
#查看红细胞基因
hb.genes <- rownames(zika_all)[grep("^Hb[^(p)]", rownames(zika_all),ignore.case = T)]
hb.genes
#这样取出的红细胞基因有些并不是
# "Hbb-bt" "Hbb-bs" "Hbs1l"非红  "Hba-a1" "Hbq1b"  "Hba-a2" "Hbq1a"  "Hbegf"非红  "Hbb-y"
# 定义一下红细胞基因集合
hb.genes <- hb.genes[-c(3,8)]
hb.genes

#添加线粒体百分比

zika_all[["percent.mt"]] <- PercentageFeatureSet(zika_all,
                                             pattern = "^mt-")
head(zika_all[["percent.mt"]] )

#添加核糖体比例
zika_all[["percent.ribo"]] <- PercentageFeatureSet(zika_all,
                                                 pattern = "^Rp[sl]")
head(zika_all[["percent.ribo"]] )

#添加红细胞比例
zika_all[["percent.hb"]] <- PercentageFeatureSet(zika_all,
                                                   features = hb.genes)
head(zika_all[["percent.hb"]] )


#画图看两组中的比例
#1.小提琴图
VlnPlot(object = zika_all,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "Group",
        log = F,
        pt.size = 0.1,
        raster=FALSE)

VlnPlot(object = zika_all,
        features = c("percent.ribo","percent.hb"),
        group.by = "Group",
        log = F,
        pt.size = 0.1,
        raster=FALSE)

#2.山峦图
RidgePlot(object = zika_all,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"),
          log = T, #对数据是否取log
          ncol = 1,
          group.by = "Group")
#3.散点图
p1 <- FeatureScatter(object = zika_all,
                     feature1 = "nCount_RNA",
                     feature2 = "nFeature_RNA",
                     group.by = "Group")
p2 <- FeatureScatter(object = zika_all,
                     feature1 = "nFeature_RNA",
                     feature2 = "percent.mt",
                     group.by = "Group")
p1+p2

#去除线粒体基因
# 去除线粒体基因和极端值
# 小鼠线粒体基因综合网站https://www.broadinstitute.org/files/shared/metabolism/mitocarta/mouse.mitocarta3.0.html
# 关于小鼠线粒体基因讨论https://www.biostars.org/p/347428/
# 关于小鼠脑的线粒体去除比例，参考文献综合判断 adv sci 10%
zika_all_mthb_filter <- subset(zika_all,
                   subset = nFeature_RNA > 200 &
                     nFeature_RNA < 2500 &
                     percent.mt >  -Inf &  # 极端值
                     percent.mt < 10 & #线粒体选择小于10%
                     percent.ribo > 3 & #核糖体选择大于3%
                     percent.hb < 1 ) #红细胞基因小于1%

#做一下细胞筛选之后的可视化
VlnPlot1 <- VlnPlot(object =zika_all_mthb_filter,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "Group",
        log = F,
        pt.size = 0.1,
        raster=FALSE)
tiff("real_filter_mito_ribo_hb_after_Vlnplot.tiff",width = 750,height = 400)
VlnPlot1
dev.off()

VlnPlot2 <- VlnPlot(object = zika_all_mthb_filter,
        features = c("percent.ribo","percent.hb"),
        group.by = "Group",
        log = F,
        pt.size = 0.1,
        raster=FALSE)
tiff("real_filter_mito_ribo_hb_after_Vlnplot2.tiff",width = 750,height = 400)
VlnPlot2
dev.off()

#2.山峦图
RidgePlot2 <- RidgePlot(object = zika_all_mthb_filter,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"),
          log = T, #对数据是否取log
          ncol = 1,
          group.by = "Group")
tiff("real_filter_mito_ribo_hb_after_RidgePlot.tiff",width = 840,height = 1200)
RidgePlot2
dev.off()

#3.散点图
p3 <- FeatureScatter(object = zika_all_mthb_filter,
                     feature1 = "nCount_RNA",
                     feature2 = "nFeature_RNA",
                     group.by = "Group")
p4 <- FeatureScatter(object = zika_all_mthb_filter,
                     feature1 = "nFeature_RNA",
                     feature2 = "percent.mt",
                     group.by = "Group")
tiff("real_filter_mito_ribo_hb_after_ScattorPlot.tiff",width = 1000,height = 400)
p3 + p4
dev.off()

#### 2、数据的质控2: 消除细胞周期#####
# 细胞周期打分
#细胞周期打分因为是人的物种，将大写全部转变为小写，再把第一个字母大写
library(Hmisc)
cc.genes$s.genes <- cc.genes$s.genes %>% tolower() %>% capitalize()
cc.genes$g2m.genes<- cc.genes$g2m.genes %>% tolower() %>% capitalize()

zika_all_mthb_filter <- CellCycleScoring(zika_all_mthb_filter,
                             s.features = cc.genes$s.genes,
                             g2m.features = cc.genes$g2m.genes)
head(zika_all_mthb_filter@meta.data[ , 8:10])

# ggplo2对细胞周期结果可视化看看,后期聚类后再处理看是否要去除细胞周期
library(ggplot2)
library(ggsci)
library(ggpubr)
filter_cycle <- ggplot(data = zika_all_mthb_filter@meta.data,
                  aes(x = zika_all_mthb_filter$S.Score,
                  y = zika_all_mthb_filter$G2M.Score,
                 fill = factor(zika_all_mthb_filter$Phase))) +
  geom_point(shape = 21,
             size = 5.5,
             alpha = 0.5) +
  theme_classic2()+
  ggtitle("Cell Cycle Score of Single Cell RNA") +
  theme(plot.title = element_text(hjust = 0.50,
                                  color="black",
                                  size = 15))
tiff("real_filter_mito_ribo_hb_after_Cellcycle",width = 670,height = 433)
filter_cycle
dev.off()
#其实这里可以保存一下去除线粒体\核糖体、血红蛋白和极值之后的数据，不保存也可以
save(zika_all_mthb_filter, file = "zika_all_mthb_filter.Rda")
load("zika_all_mthb_filter.Rda")

#### 3、数据的标准化、批次矫正####
zika_all_mthb_filter

#添加样本信息，分别那每个样本分开做标准化
# 添加分组信息
zika_all_mthb_filter@meta.data <- zika_all_mthb_filter@meta.data %>%
  mutate(Sample = if_else(str_sub(rownames(zika_all_mthb_filter@meta.data), 1, 4) == "V_59",
                         "V_59",
                         if_else(str_sub(rownames(zika_all_mthb_filter@meta.data), 1, 4) == "V_FS",
                                 "V_FS",
                                 'MOCK')))

#将seurat对象按照样本拆分
zika_all_mthb_filter_list <- SplitObject(zika_all_mthb_filter,
                         split.by = "Sample")
zika_all_mthb_filter_list

#先对一个样本进行SCT标准化看看
test_Cont <- zika_all_mthb_filter_list[[1]]
test_Cont <- SCTransform(test_Cont,
                         variable.features.n = 2000,
                         # vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                         verbose = FALSE)

#查看一下计算结果
HVFInfo(object = test_Cont)[1:15,1:3]
top10 <- head(VariableFeatures(test_Cont), 10)

#可视化结果
plot1 <- VariableFeaturePlot(test_Cont)
LabelPoints(plot = plot1,
            points = top10,
            repel = TRUE)
#多个样品的话抽出一个进行测试



#多个seurat对象进行批量标准化、计算高变基因
for(i in 1:length(zika_all_mthb_filter_list)){
  zika_all_mthb_filter_list[[i]] <- SCTransform(
    zika_all_mthb_filter_list[[i]], 
    variable.features.n = 3000, #选择3000高变基因
    vars.to.regress = c("nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"),#第一次跑没有去除细胞周期
    verbose = FALSE)
}

# ⽤ SelectIntegrationFeatures 来筛选⽤于整合的基因
features <- SelectIntegrationFeatures(object.list =  zika_all_mthb_filter_list,
                                      nfeatures = 3000)

# 准备要整合的数据，⽤ PrepSCTIntegration
zika_all_mthb_filter_list <- PrepSCTIntegration(object.list =  zika_all_mthb_filter_list,
                                anchor.features = features)

#选择锚定基因
#这一步运行起来相当慢，因此考虑使用并行计算
# linux查看总的线程 grep 'processor' /proc/cpuinfo |  wc -l 数https://www.jianshu.com/p/b9bd94ef2703 https://www.cnblogs.com/52-qq/p/9772680.html
library(future)
# check the current active plan
plan()
#增加线程数 用100个线程跑一下
plan("multiprocess", workers = 100)
plan()
#并行之后会报错，增大内存数量
options(future.globals.maxSize = 10000 * 1024^2)

AnchorSet <- FindIntegrationAnchors(object.list = zika_all_mthb_filter_list,
                                    reference = 1,
                                    normalization.method = "SCT",
                                    anchor.features = features)

#整合数据
zika_all_mthb_filter_SCT <- IntegrateData(anchorset = AnchorSet,
                          normalization.method = "SCT")

DefaultAssay(zika_all_mthb_filter_SCT) <- "integrated"

#将分析好的数据保存起来
save(zika_all_mthb_filter_SCT, file = "zika_all_mthb_filter_SCT_scalemtcycle.Rda")








