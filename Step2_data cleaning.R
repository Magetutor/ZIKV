#Preprocessing of Data
setwd('~/analysis/zika/')
rm(list = ls())

#If the data is reloaded, use the load function for Rda
load(file = 'zika_all_group.Rda')

#1. Data Quality Control 1: Removal of Mitochondrial Genes and Outliers##
Finding mitochondrial genes
library(tidyverse)
library(Seurat)

#View mitochondrial genes
mito.genes <- str_subset(string = rownames(zika_all),
pattern = "mt-") # mouse mitochondrial genes are flagged with ^mt-, rather than ^Mt-
mito.genes

#View ribosomal genes
ribo.genes <-  rownames(zika_all)[grep("^Rp[sl]", rownames(zika_all),ignore.case = T)]
ribo.genes

#View hemoglobin genes
hb.genes <- rownames(zika_all)[grep("^Hb[^(p)]", rownames(zika_all),ignore.case = T)]
hb.genes

#Some of the genes obtained for red blood cells are not valid:
#"Hbb-bt" "Hbb-bs" "Hbs1l" (not red) "Hba-a1" "Hbq1b" "Hba-a2" "Hbq1a" "Hbegf" (not red) "Hbb-y"
#Define the set of valid red blood cell genes
hb.genes <- hb.genes[-c(3,8)]
hb.genes

#Add percentage of mitochondrial genes
zika_all[["percent.mt"]] <- PercentageFeatureSet(zika_all,
pattern = "^mt-")
head(zika_all[["percent.mt"]] )

#Add ribosomal gene percentage
zika_all[["percent.ribo"]] <- PercentageFeatureSet(zika_all,
pattern = "^Rp[sl]")
head(zika_all[["percent.ribo"]] )

#Add red blood cell gene percentage
zika_all[["percent.hb"]] <- PercentageFeatureSet(zika_all,
features = hb.genes)
head(zika_all[["percent.hb"]] )

#Plot the proportions in the two groups
#1. Violin plots
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

#2. Ridge plots
RidgePlot(object = zika_all,
features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"),
log = T,
ncol = 1,
group.by = "Group")

#3. Scatter plots
p1 <- FeatureScatter(object = zika_all,
feature1 = "nCount_RNA",
feature2 = "nFeature_RNA",
group.by = "Group")
p2 <- FeatureScatter(object = zika_all,
feature1 = "nFeature_RNA",
feature2 = "percent.mt",
group.by = "Group")
p1 + p2

# Remove mitochondrial genes
zika_all_mthb_filter <- subset(zika_all,
                   subset = nFeature_RNA > 200 &
                     nFeature_RNA < 2500 &
                     percent.mt >  -Inf &  # 极端值
                     percent.mt < 10 & #线粒体选择小于10%
                     percent.ribo > 3 & #核糖体选择大于3%
                     percent.hb < 1 ) #红细胞基因小于1%

#Visualize after cell filtering
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

# Ridge plot
RidgePlot2 <- RidgePlot(object = zika_all_mthb_filter,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"),
          log = T, #对数据是否取log
          ncol = 1,
          group.by = "Group")
tiff("real_filter_mito_ribo_hb_after_RidgePlot.tiff",width = 840,height = 1200)
RidgePlot2
dev.off()

# Scatter plot
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

#### 2、Cell cycle scoring#####
library(Hmisc)
cc.genes$s.genes <- cc.genes$s.genes %>% tolower() %>% capitalize()
cc.genes$g2m.genes<- cc.genes$g2m.genes %>% tolower() %>% capitalize()

zika_all_mthb_filter <- CellCycleScoring(zika_all_mthb_filter,
                             s.features = cc.genes$s.genes,
                             g2m.features = cc.genes$g2m.genes)
head(zika_all_mthb_filter@meta.data[ , 8:10])

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
save(zika_all_mthb_filter, file = "zika_all_mthb_filter.Rda")
load("zika_all_mthb_filter.Rda")

# 3. Data normalization and batch correction
zika_all_mthb_filter

zika_all_mthb_filter@meta.data <- zika_all_mthb_filter@meta.data %>%
  mutate(Sample = if_else(str_sub(rownames(zika_all_mthb_filter@meta.data), 1, 4) == "V_59",
                         "V_59",
                         if_else(str_sub(rownames(zika_all_mthb_filter@meta.data), 1, 4) == "V_FS",
                                 "V_FS",
                                 'MOCK')))


zika_all_mthb_filter_list <- SplitObject(zika_all_mthb_filter,
                         split.by = "Sample")
zika_all_mthb_filter_list

for(i in 1:length(zika_all_mthb_filter_list)){
  zika_all_mthb_filter_list[[i]] <- SCTransform(
    zika_all_mthb_filter_list[[i]], 
    variable.features.n = 3000, 
    vars.to.regress = c("nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"), # The first run without removing cell cycle
    verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list =  zika_all_mthb_filter_list,
                                      nfeatures = 3000)

zika_all_mthb_filter_list <- PrepSCTIntegration(object.list =  zika_all_mthb_filter_list,
                                anchor.features = features)


AnchorSet <- FindIntegrationAnchors(object.list = zika_all_mthb_filter_list,
                                    reference = 1,
                                    normalization.method = "SCT",
                                    anchor.features = features)

# Integrate data
zika_all_mthb_filter_SCT <- IntegrateData(anchorset = AnchorSet,
                          normalization.method = "SCT")

DefaultAssay(zika_all_mthb_filter_SCT) <- "integrated"

# Save the analyzed data








