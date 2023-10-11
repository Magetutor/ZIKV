setwd('~/analysis/zika_new/')
# rm(list = ls())
library(Seurat)
library(tidyverse)
library(stringr)
library(patchwork)
library(viridis)
load("zika_all_mthb_filter_SCT_scalemtcycle.Rda")
dim(zika_all_mthb_filter_SCT)

############ PCA Analysis ###################
# Perform PCA analysis
zika_all_mthb_filter_SCT <- RunPCA(object = zika_all_mthb_filter_SCT,
                   npcs = 50,
                   rev.pca = FALSE,
                   weight.by.var = TRUE,
                   verbose = TRUE, # Output feature genes for specific PC values
                   ndims.print = 1:5, # Output feature genes for the first 5 PC values
                   nfeatures.print = 30, # Output 30 feature genes for each PC value
                   reduction.key = "PC_")
# You can plot a heatmap here: https://www.jianshu.com/p/21c1934e3061
View(zika_all_mthb_filter_SCT@reductions[["pca"]]@feature.loadings) # Modify gene-PC value mapping
View(zika_all_mthb_filter_SCT@reductions[["pca"]]@cell.embeddings) # Mapping coordinates of each cell on low-dimensional PCA axes
DimHeatmap(object = zika_all_mthb_filter_SCT, dims = 1, cells = 1000,
           balanced = TRUE, fast = FALSE, nfeatures = 30) + scale_fill_viridis()

# Select appropriate PCs
ElbowPlot(zika_all_mthb_filter_SCT,
          ndims = 50) # Find principal components contributing significantly to cell differences

# Heatmap display of markers for each principal component
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

# PCA visualization
# PCA visualization 1: One-dimensional visualization, visualize genes with higher weights in each dimension
# pdf(file="06.pcaGene.pdf", width=10, height=8)
VizDimLoadings(object = zika_all_mthb_filter_SCT, dims = 1:4, 
               reduction = "pca", nfeatures = 20)
# dev.off()

# PCA visualization 2: Two-dimensional visualization
# pdf(file="06.PCA.pdf", width=6.5, height=6)
DimPlot(object = zika_all_mthb_filter_SCT, reduction = "pca")
# dev.off()


######### Identify Batch Effects from Dimensionality Reduction Plots #########

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
### Check if cell cycle affects the analysis
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

############ tSNE Analysis ################
zika_all_mthb_filter_SCT <- RunTSNE(zika_all_mthb_filter_SCT,
                    dims = 1:20) # Perform tSNE analysis with the selected 20 PC values from the contribution plot
# Plot tSNE based on group
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
#Plot tSNE graph grouped by cell cycle phase
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

############## UMAP Analysis ##############
zika_all_mthb_filter_SCT <- RunUMAP(zika_all_mthb_filter_SCT,
dims = 1:20)

#Grouping by 'Group'
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

#Grouping by cell cycle phase
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

#Combining the plots for different methods
library(cowplot)
library(patchwork)

#Grouping by 'Group'
plot_grid(p1, p3, p5,
labels = c("A", "B", "C"),
nrow = 1,
align = "h")
plot_grid(p2, p4, p6,
labels = c("A", "B", "C"),
nrow = 1,
align = "h")

#Grouping by cell cycle phase
p2 + p4 + p6 +
plot_layout(guides = "collect") &
theme(legend.position = "bottom")

#Displaying specific genes in UMAP plot
FeaturePlot(zika_all_mthb_filter_SCT,
reduction = "tsne",
features = "Mki67") # G2M phase marker
FeaturePlot(zika_all_mthb_filter_SCT,
reduction = "tsne",
features = "Top2a") # S phase marker (none)

####### 2. Single Cell Clustering Analysis ######

#Perform clustering using K-Nearest Neighbor algorithm
zika_all_mthb_filter_SCT <- FindNeighbors(zika_all_mthb_filter_SCT,
k.param = 20,
dims = 1:20)
zika_all_mthb_filter_SCT <- FindClusters(zika_all_mthb_filter_SCT,
resolution = 0.5, # The higher the value, the more clusters formed (most important parameter)
method = "igraph", # Adjust based on results
algorithm = 1, # Adjust based on results'
random.seed = 2022)

##Count the cells in each cluster
table(zika_all_mthb_filter_SCT@meta.data$seurat_clusters)

#Visualize cell types
DimPlot(object = zika_all_mthb_filter_SCT,
group.by = "seurat_clusters",
reduction = "pca",
label = TRUE)
PCAPlot(object = zika_all_mthb_filter_SCT,
group.by = "seurat_clusters",
pt.size = 0.5,
label = TRUE)

#tSNE clustering results are better
TSNEPlot(object = zika_all_mthb_filter_SCT,
group.by = "seurat_clusters",
pt.size = 0.5,
label = TRUE)
UMAPPlot(object = zika_all_mthb_filter_SCT,
group.by = "seurat_clusters",
pt.size = 0.5,
label = TRUE)

Save clustered cells
save(zika_all_mthb_filter_SCT, file = "zika_all_mthb_filter_SCT_cluster23.Rda")

#Cell cycle has a significant impact, so remove the effect in the previous step
############# 3. Cell Annotation ###################

#Generally, load the data from the previous step
load(file = "zika_all_mthb_filter_SCT_cluster23.Rda")

#Read pre-annotated marker information
cell_marker <- read.csv("mouse_brain_cell_atlas_paper.csv")

#Extract relevant information from the file
cell_marker <- cell_marker[c(3, 5)]
View(cell_marker)

# If there are multiple marker genes, separate them using the following code
cell_marker <- cell_marker %>%
  separate_rows(Cell.Marker, sep = ", ") %>%
  na.omit() %>%
  distinct() %>%
  arrange(Cell.Type)

# Remove leading whitespace from some markers
cell_marker$Cell.Marker <- str_trim(cell_marker$Cell.Marker, "left")

# Perform annotation
p7 <- DimPlot(zika_all_mthb_filter_SCT,
              reduction = "tsne",
              group.by = "seurat_clusters",
              label = TRUE)
p8 <- DotPlot(zika_all_mthb_filter_SCT,
              features = unique(cell_marker$Cell.Marker)) +
  theme(axis.text = element_text(size = 8,
                                 angle = 90,
                                 hjust = 1))
p7 + p8

# Check expression range
range(zika_all_mthb_filter_SCT@assays$integrated@scale.data)

# Check specific gene expression
# Using B cells as an example
p9 <- FeaturePlot(zika_all_mthb_filter_SCT,
                  reduction = "tsne",
                  features = c("Ms4a1", "Cd19", "Pxk", "Cd74", "Cd79a"),
                  label = TRUE)
p9
p7 + p9

# Next, try using SingleR
BiocManager::install('SingleR')
BiocManager::install('celldex')
BiocManager::install('BiocParallel')
library(SingleR)
library(celldex)
library(BiocParallel)

# Download annotation file
ref <- celldex::MouseRNAseqData()
save(ref, file = "MouseRNAseqData.Rda")

# Plot differential genes for each cluster
# load('zika_all_mthb_filter_SCT_cell_SingleR.Rda')
Idents(zika_all_mthb_filter_SCT) <- 'cell_type1'
Idents(zika_all_mthb_filter_SCT)

logFCfilter=0.5
adjPvalFilter=0.05
zika_all_mthb_filter_SCT.markers <- FindAllMarkers(object = zika_all_mthb_filter_SCT,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = logFCfilter)
# Save the marker file
save(zika_all_mthb_filter_SCT.markers, file = "zika_all_SCT.markers_0725.Rda")

sig.markers = zika_all_mthb_filter_SCT.markers[(abs(as.numeric(as.vector(zika_all_mthb_filter_SCT.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(zika_all_mthb_filter_SCT.markers$p_val_adj))<adjPvalFilter),]
write.table (sig.markers,file="PC20_09.markers.xls",sep="\t",row.names=F,quote=F)
top5 <- zika_all_mthb_filter_SCT.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- zika_all_mthb_filter_SCT.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- zika_all_mthb_filter_SCT.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table (top20 ,file="PC20_top20.markers.xls",sep="\t",row.names=F,quote=F)


# Draw heatmap of markers in each cluster
pdf(file="09.tsneHeatmap.pdf",width=32,height=8)
DoHeatmap(object = zika_all_mthb_filter_SCT, features = top5$gene) + NoLegend()
dev.off()  # Slow to generate

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

# Use check.marker to visualize
DotPlot(zika_all_mthb_filter_SCT,
        features = unique(check.marker)) +
  theme(axis.text = element_text(size = 8,
                                 angle = 90,
                                 hjust = 1))

# Draw violin plot of markers
pdf(file="02.markerViolin.pdf",width=10,height=6)
VlnPlot(object = zika_all_mthb_filter_SCT, features = c("Cd83", "Cxcl2"))
dev.off()

# Draw scatter plot of markers in each cluster
pdf(file="02.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = zika_all_mthb_filter_SCT, features = c("Cd83", "Cxcl2"))
dev.off()

# Draw bubble plot of markers in each cluster
pdf(file="06.markerBubble.pdf",width=12,height=6)
cluster0Marker=c("Cbr2", "Lyve1", "Selenop", "Folr2", "Ednrb", "F13a1", "Mrc1", "Igf1", "Slc40a1", "Cd163")
DotPlot(object = zika_all_mthb_filter_SCT, features = cluster0Marker)
dev.off()

load('20220720recluster.Rda')
# Define cell types
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

#The cell names corresponding to each subtype are defined and then written to the meta-.data file
zika_all_mthb_filter_SCT[['cell_type']] <- unname(Cell_type[zika_all_mthb_filter_SCT@meta.data$seurat_clusters])


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

table(zika_all_mthb_filter_SCT@meta.data$cell_type1)
zika_all_mthb_filter_SCT@meta.data$Group <- factor(zika_all_mthb_filter_SCT@meta.data$Group,
                                                      levels = c('MOCK','V_FS','V_59'))


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

pdf('group_umap.pdf',width = 20,height = 8)
DimPlot(zika_all_mthb_filter_SCT,
        reduction = "umap",
        split.by = "Group",
        group.by = "cell_type1",
        label = TRUE,
        pt.size = 0.3) +
  NoLegend()
dev.off()
# savedata
save(zika_all_mthb_filter_SCT, file = "20220720recluster.Rda")



