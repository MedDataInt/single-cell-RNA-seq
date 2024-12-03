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

## Normalization of the Data
Normalize gene expression levels across cells to account for differences in sequencing depth. Sequencing depth refers to the total number of reads in a cell, which can vary across cells.The most common approach is total-count normalization, where each cell is scaled so that it has the same total count of reads. This makes it easier to compare cells with different sequencing depths.  
```r
# normalization
scramble <- NormalizeData(scramble)

#Identify highly variable feaures (features selection)
scramble <- FindVariableFeatures(scramble, selection.method = "vst", nfeatures = 2000)

# scaling the data
scramble_all.genes <- rownames(scramble)
scramble <- ScaleData(scramble, features = scramble_all.genes)
```
## Perform linear dimensional reduction 
```r
# linear dimensional reduction
scramble <- RunPCA(scramble, features = VariableFeatures(object = scramble))
```

## Cluster the cells
```r
scramble <- FindNeighbors(scramble, dims = 1:20)
scramble <- FindClusters(scramble, resolution = 0.5)
head(Idents(scramble), 5)
```
## Run non-linear dimensional reduction (UMAP/tSNE)
```r
scramble <- RunUMAP(scramble, dims = 1:20)
pdf('UMAP_scramble.pdf', width = 8, height = 8)
DimPlot_scCustom(seurat_object = scramble)
dev.off()

# save seurat object 
saveRDS(scramble, file = "scramble.rds")
```


## integrating scRNA-seq

## scRNA-seq clustering

## scRNA-seq cell type annotation

## scRNA-seq cell-cell interaction 

## scRNA-seq gene signature/sets analysis

## scRNA-seq gene differential expression 

## scRNA-seq interactive visulization
