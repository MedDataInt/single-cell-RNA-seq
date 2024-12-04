# scRNA-seq
This repository is designed to help beginners learn the fundamental aspects of scRNA-seq, with a focus on both theoretical knowledge and practical applications.  
This is based on one of our published project.
```r
# loading library
library(Seurat)
library(dplyr)
library(magrittr)
library(data.table)
library(Matrix)
library(devtools)
library(RcppArmadillo)
library(Rcpp)
library(scales)
library(pheatmap)
library(gplots)
library(ggplot2)
library(cowplot)
library(foreach)
library(doParallel)
library(stringr)
library(here)
library(purrr)
library(readr)
library(harmony)
library(scCustomize)
library(SeuratDisk)
library(patchwork)
library(viridis)
library(qs)
library(DropletUtils)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(base)
```
## preparing data for scRNA-seq
h5 file: a type of file format used to store large amounts of data in a hierarchical manner. This makes it easy to store and organize complex datasets. In scRNA-seq, you might have groups for "gene expression", "cell metadata", and "quality control" within the same .h5 file.

```r
# convert to h5 file
scramble_filter_matrix <- Read10X("filtered_feature_bc_matrix/")
write10xCounts("filtered_feature_bc_matrix.h5", scramble_filter_matrix, type = "HDF5",
               genome = "mm10", version = "3", overwrite = TRUE,
               gene.id = rownames(scramble_filter_matrix),
               gene.symbol = rownames(scramble_filter_matrix))

# create seurat object
scramble <- Read10X_h5('filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)
scramble <- CreateSeuratObject(counts = scramble, project = "scramble", min.cells = 3, min.features = 200)
```


## scRNA-seq quality control
QC is to ensure that the data you are working with is reliable and that any low-quality cells or technical artifacts are filtered out before performing downstream analyses.  
* Number of Genes Detected per Cell
* Number of Unique Molecules  
* Mitochondrial RNA Content  
* Cell Doublets  
```r
# filtering cells
scramble[["percent.mt"]] <- PercentageFeatureSet(scramble, pattern = "^mt-") # add percent of mitochondrial
scramble <- subset(scramble, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

## Normalize and scale the Data
Normalize gene expression levels across cells to account for differences in sequencing depth. Sequencing depth refers to the total number of reads in a cell, which can vary across cells.The most common approach is total-count normalization, where each cell is scaled so that it has the same total count of reads. This makes it easier to compare cells with different sequencing depths.  While scaling is a standard pre-processing step prior to dimensional reduction, it allow the mean expression across cells is 0, and variance across cells is 1, in this way, the highly expressed genes do not dominate.In this case, only the selected features are scaled.
```r
# normalization
scramble <- NormalizeData(scramble)

#Identify highly variable feaures (they are highly expressed in some cells, and lowly expressed in others)
scramble <- FindVariableFeatures(scramble, selection.method = "vst", nfeatures = 2000) # default

# scaling the data
scramble_all.genes <- rownames(scramble)
scramble <- ScaleData(scramble, features = scramble_all.genes)
```
## Perform linear dimensional reduction 
PCA will be performed on the scaled data.for the first principal components, the output is a list of genes with the most positive and negative loadings.
```r
# linear dimensional reduction
scramble <- RunPCA(scramble, features = VariableFeatures(object = scramble))
```

## Cluster the cells
First need to decide how many components should we choose, then construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (it takes previous defined dimensionality of dataset, first 20 PCs); Last, we need to apply modularity optimization tech (Louvain algorithm- default or SLM) to iteratively group cells togeter. The resolution can decide the granularity of the downstream clustering. 
```r
# determine the dimensionality of the dataset using ElbowPlot
pdf('Elbow_scramble.pdf', width = 8, height = 8)
ElbowPlot(scramble)
dev.off()

scramble <- FindNeighbors(scramble, dims = 1:20) # say 20,  
scramble <- FindClusters(scramble, resolution = 0.5)
head(Idents(scramble), 5) # check on the cluster ID of the first 5 cells.
```
## Run non-linear dimensional reduction (UMAP/tSNE)
encourage to leverage techniques like UMAP for visualization, cells that are grouped together within graph-based clusters determined above should co-localize on these dimension reduction plots.
```r
scramble <- RunUMAP(scramble, dims = 1:20)
pdf('UMAP_scramble.pdf', width = 8, height = 8)
DimPlot_scCustom(seurat_object = scramble)
dev.off()

# save seurat object 
saveRDS(scramble, file = "scramble.rds")
```

## Integrating scRNA-seq
In this section, we aim to integrate three samples in different conditions together for subsequent analysis--scramble, mut1, and mut2. 
Merge based on the normalized data.By default, merge() will combine the Seurat objects based on the raw count matrices, erasing any previously normalized and scaled data matrices. If you want to merge the normalized data matrices as well as the raw count matrices, 
simply pass merge.data = TRUE. This should be done if the same normalization approach was applied to all objects. Without integration, cells are grouping both cel type and by underlying method, and cluster this dataset will return predominatly batch-specific clusters.
While integrative analysis using IntegrateLayers can return a dimensional reduction that aims to co-embed shared cell types across batches. Five integration methods can be used:

* Anchor-based CCA integration (method=CCAIntegration)
* Anchor-based RPCA integration (method=RPCAIntegration) --faster and more conservative 
* Harmony (method=HarmonyIntegration)
* FastMNN (method= FastMNNIntegration)
* scVI (method=scVIIntegration)
Once integrative analysis is complete, you can rejoin the layers - which collapses the individual datasets together and recreates the original counts and data layers. You will need to do this before performing any differential expression analysis. However, you can always resplit the layers in case you would like to reperform integrative analysis.
```r
# merge two mut samples first
Mut_normalized <- merge(mut1, y= mut2, add.cell.ids = c('mut1', 'mut2'), project = 'mut', merge.data = TRUE)
LayerData(Mut_normalized)[1:10, 1:15]

# merge with scramble
normalized <- merge(Mut_normalized, y= scramble, add.cell.ids = c('mut', 'scramble'), project = 'combined', merge.data = TRUE)
LayerData(Marchf8_normalized)[1:10, 1:15]

# Perform analysis without integration
# run standard anlaysis workflow
normalized <- NormalizeData(normalized)
normalized <- FindVariableFeatures(normalized)
normalized <- ScaleData(normalized)
normalized <- RunPCA(normalized)
print('done')

# find clusters without integration
normalized <- FindNeighbors(normalized, dims = 1:30, reduction = "pca")
normalized <- FindClusters(normalized, resolution = 0.5, cluster.name = "unintegrated_clusters")
normalized <- RunUMAP(normalized, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

pdf('combined_umap.unintegrated.pdf', width = 16, height = 8)
DimPlot(normalized, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))
dev.off()

# Perform integration
normalized <- IntegrateLayers(object = normalized, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                   verbose = FALSE)
print('done')


# find clusters in integration
normalized <- FindNeighbors(normalized, reduction = "integrated.cca", dims = 1:30)
normalized <- FindClusters(normalized, resolution =0.5)
normalized <- RunUMAP(normalized, dims = 1:30, reduction = "integrated.cca")
print('done')

pdf('Combined_umap.integrated_05.pdf', width = 16, height = 8)
DimPlot_scCustom(Marchf8_normalized, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
dev.off()

# calculate the cell number in each of cluster 
table(Idents(normalized))
md <- normalized@meta.data %>% as.data.table
md[, .N, by = c("orig.ident", "seurat_clusters")]
cellcount = md[, .N, by = c("orig.ident", "seurat_clusters")]
file_name <- "Combine_integrated_metadata_cellcount05.csv"
fwrite(cellcount, file = file_name)

# re-join layers after integration
# normalized[["RNA"]] <- JoinLayers(normalized[["RNA"]]) 

```
## Finding differentially expressed expression (DEs)
```r
# find differentially expression
normalized <- JoinLayers(normalized, verbose = TRUE)
de = FindAllMarkers(normalized, only.pos = TRUE)

csv_file <- "Combined_differential_genes_05.csv"
write.csv(de, file = csv_file, row.names = FALSE)

# Find markers and limit to those expressed in greater than 25% of target population and visulize
all_markers <- de %>%
  Add_Pct_Diff() %>%
  filter(pct_diff > 0.25)
top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7, named_vector = FALSE,
                                   make_unique = TRUE)
pdf('bubble05.pdf', width = 8.5, height = 14)
Clustered_DotPlot(seurat_object = normalized, features = top_markers, k=16)
dev.off()

# check on specific gene expression in umap
pdf('combined_CD4.pdf', width = 5, height = 5)
FeaturePlot_scCustom(seurat_object = normalized,  reduction = "umap",features = "Cd4",alpha_exp = 0.75, pt.size = 0.8)
dev.off()
saveRDS(normalized, file = "normalized.rds")
```
## subcluster analysis and identify the differential expressed genes between groups
In this section, if you are interested in a certain specific cluster, and want to do the subsequent cluster analysis on this cluster, then you can go continue the cluster just like what have done above. otherwise if you are interested on the DEGs of this cluster between the groups, we can follow this chunk of script.
```r
path_to_rds <-'normalized.rds'
normalized<-readRDS(path_to_rds)
C7 <- subset(normalized, subset = seurat_clusters %in% c(7))

Idents(C7) <- CD8@meta.data$orig.ident

# volin plot on the comparation of PD1 expression 
pdf('C7_Pdcd1.pdf', width =5, height = 6)
VlnPlot(C7, split.by = 'orig.ident', features = 'Pdcd1', pt.size = 0.1,) + NoLegend()
dev.off()

# Find markers, export the full different expressed genes list 
mut2_C7_DE <- FindMarkers(C7, ident.1 = "mut2", ident.2 = "scramble", verbose = TRUE, subset.ident = "0", logfc.threshold = log(1))
csv_file <- "mut2_C7_DE.csv"
write.csv(mut2_C7_DE, file = csv_file, row.names = TRUE)
```
## scRNA-seq gene signature/sets analysis
Calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin. detail can be found [here](https://satijalab.org/seurat/reference/addmodulescore)
```r
# antigen_processing_and_presentation, gene list can be download from KEGG database or other dataset 
data <- read.csv("antigen_processing_and_presentation_mus.csv")
gene_list <- data$Gene
print(gene_list)

antigen_marker <-list(gene_list)
normalized <- AddModuleScore(normalized, features = antigen_marker, ctrl = 5, name = 'antigen_marker')
print(head(normalized@meta.data)) # check the meta.data 

# vizualize the gene set in the s 
pdf('antigen_Split1d.pdf', width = 6, height = 3)
FeaturePlot(normalized,features = "antigen_marker1", label = FALSE, split.by = "orig.ident", cols = viridis(11), pt.size = 0.8, repel = TRUE) 
dev.off()

# Extract orig.ident and autophage_marker1 from Seurat object
data <- data.frame(orig.ident = normalized@meta.data$orig.ident, 
                   antigen_marker1 = normalized@meta.data$antigen_marker1)

# Write the data to a CSV file
write.csv(data, file = "normalized_antigen_marker_data.csv", row.names = FALSE)
```
## scRNA-seq cell-cell interaction 
there are many approaches allow us to do the cell-cell interaction, and [this paper](https://www.nature.com/articles/s41467-022-30755-0) is good review regarding the CCI, and in this section, I will show how to perform cell-cell interaction in single cells using CellChat.
```r
# Label the clusters 
normalized@meta.data$labels <- NA
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(5,14)] <- "B"
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(1, 2, 3)] <- "MDSCs"
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(4,13)] <- "macrophages"
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(10)] <- "NK"
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(7)] <- "CD8_T"
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(6,9)] <- "CD4_T"
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(11)] <- "DCs"
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(12)] <- "epithelial"
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(15)] <- "fibroblasts"
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(0)] <- "tumor"
normalized@meta.data$labels[integrated@meta.data$seurat_clusters %in% c(8)] <- "eosinophils"

# prepare the input data, we work on mut2 sample only here
data.input <- normalized@assays$RNA$data
meta <- normalized@meta.data 
cell.use<- rownames(meta)[meta$orig.ident == 'mut2'] # extract the cell names 
data.input <- data.input[, cell.use]
meta <- meta[cell.use, ]
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# check on the database
CellChatDB <- CellChatDB.mouse # use  if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat) # returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors 
head(df.net)
write.csv(df.net, file = "mut2_df_net.csv", row.names = FALSE)

df.net1 <- subsetCommunication(cellchat, slot.name = "netP") # to access the the inferred communications at the level of signaling pathways
head(df.net1)
write.csv(df.net1, file = "mut2_df_net1.csv", row.names = FALSE)

cellchat <- computeCommunProbPathway(cellchat)
str(cellchat@net)
str(cellchat@netP)
cellchat <- aggregateNet(cellchat)

pdf("mut2_CCI_default.pdf")
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

# we can examine the signaling sent from each cell group
pdf("mut2_CCI_default_individal.pdf")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

```
## scRNA-seq interactive visulization
