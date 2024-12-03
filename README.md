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

## integrating scRNA-seq
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

# re-join layers after integration
normalized[["RNA"]] <- JoinLayers(normalized[["RNA"]]) 

```


## scRNA-seq cell-cell interaction 

## scRNA-seq gene signature/sets analysis

## scRNA-seq gene differential expression 

## scRNA-seq interactive visulization
