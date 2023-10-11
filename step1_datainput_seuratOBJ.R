setwd('~/analysis/zika/')
rm(list = ls())

BiocManager::install('Seurat')
BiocManager::install('scater')
BiocManager::install('cowplot')
install.packages("tidyverse")
library(Seurat)
library(tidyverse)

####### Construction of individual sample Seurat objects ######
scRNA_9A59 <- Read10X(data.dir = "~/analysis/zika/data/210009A_59_filtered_feature_bc_matrix/")
scRNA_9AFS <- Read10X(data.dir = "~/analysis/zika/data/210009A_FS_filtered_feature_bc_matrix/")
scRNA_9AMOCK <- Read10X(data.dir = "~/analysis/zika/data/210009A_MOCK_filtered_feature_bc_matrix/")

####### Create Seurat objects #############
V_59 <- CreateSeuratObject(counts = scRNA_9A59,
min.cells = 3,
min.features = 200)
V_FS <- CreateSeuratObject(counts = scRNA_9AFS,
min.cells = 3,
min.features = 200)
MOCK <- CreateSeuratObject(counts = scRNA_9AMOCK,
min.cells = 3,
min.features = 200)

#Check the structure
dim(V_59)
dim(V_FS)
dim(MOCK)

#Combine samples together
zika_list <- list(V_59 = V_59,
V_FS = V_FS,
MOCK = MOCK)

#Merge the samples together for data integration
zika_all <- merge(x = zika_list[[1]], # The first one
y = zika_list[c(2:length(zika_list))], # Other Seurat objects except the first one
add.cell.ids = names(zika_list)) # Add cell ID prefix
dim(zika_all)

#Add grouping metadata to indicate different strains of infection and control group
#Add grouping information
zika_all@meta.data <- zika_all@meta.data %>%
mutate(Group = if_else(str_sub(rownames(zika_all@meta.data), 1, 4) == "V_59",
"V_59",
if_else(str_sub(rownames(zika_all@meta.data), 1, 4) == "V_FS",
"V_FS",
"MOCK")))
head(zika_all@meta.data)

save(zika_all, file = "zika_all_group_filter.Rda")
