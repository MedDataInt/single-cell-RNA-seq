# scRNA-seq
This repository is designed to help beginners learn the fundamental aspects of scRNA-seq, with a focus on both theoretical knowledge and practical applications.

## preparing data for scRNA-seq
h5 file: a type of file format used to store large amounts of data in a hierarchical manner. This makes it easy to store and organize complex datasets. In scRNA-seq, you might have groups for "gene expression", "cell metadata", and "quality control" within the same .h5 file.
convert to h5 file
```r
library(Seurat)
library(here)
# read data as seruat object
```

## scRNA-seq quality control
QC is to ensure that the data you are working with is reliable and that any low-quality cells or technical artifacts are filtered out before performing downstream analyses.  
a. Number of Genes Detected per Cell  
b. Number of Unique Molecules  
c. Mitochondrial RNA Content  
d. Cell Doublets  
```r
# filtering cells

```

## Normalization of the Data
Normalize gene expression levels across cells to account for differences in sequencing depth. Sequencing depth refers to the total number of reads in a cell, which can vary across cells.The most common approach is total-count normalization, where each cell is scaled so that it has the same total count of reads. This makes it easier to compare cells with different sequencing depths.  
```r
# normalization

```



## integrating scRNA-seq

## scRNA-seq clustering

## scRNA-seq cell type annotation

## scRNA-seq cell-cell interaction 

## scRNA-seq gene signature/sets analysis

## scRNA-seq gene differential expression 

## scRNA-seq interactive visulization
