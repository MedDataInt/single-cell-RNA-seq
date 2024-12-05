###### Project:  single cell RNA-seq Analysis
###### Date: 12-2024

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


# Set working directory
# setwd("/home/jiew/H5")

# Load helper functions
source("../R/helpers.R")

# Convert data to HDF5 format
write_h5("scramble/filtered_feature_bc_matrix/",
         "scramble_filtered_feature_bc_matrix.h5")
write_h5("mut1/filtered_feature_bc_matrix/",
         "mut1_filtered_feature_bc_matrix.h5")
write_h5("mut2/filtered_feature_bc_matrix/",
         "mut2_filtered_feature_bc_matrix.h5")

# Load HDF5 files into Seurat objects
scramble <- Read10X_h5("scramble_filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject(project = "scramble", min.cells = 3, min.features = 200)

mut1 <- Read10X_h5("mut1_filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject(project = "mut1", min.cells = 3, min.features = 200)

mut2 <- Read10X_h5("mut2_filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject(project = "mut2", min.cells = 3, min.features = 200)

# add mitochondrial percentage
scramble <- add_mito_percentage(scramble)
mut1 <- add_mito_percentage(mut1)
mut2 <- add_mito_percentage(mut2)

# Generate QC plots
generate_QC_plot(scramble, "scramble_QC.pdf")
generate_QC_plot(mut1, "mut1_QC.pdf")
generate_QC_plot(mut2, "mut2_QC.pdf")

# Filter cells based on quality control metrics
scramble <- filter_cells(scramble)
mut1 <- filter_cells(mut1)
mut2 <- filter_cells(mut2)

# Normalize, find variable features, and scale data
scramble <- process_data(scramble)
mut1 <- process_data(mut1)
mut2 <- process_data(mut2)

# Cluster cells and perform UMAP dimensionality reduction
scramble <- cluster_and_visualize(scramble, "UMAP_scramble.pdf")
mut1 <- cluster_and_visualize(mut1, "UMAP_mut1.pdf")
mut2 <- cluster_and_visualize(mut2, "UMAP_mut2.pdf")

# Save processed Seurat objects
saveRDS(scramble, "scramble.rds")
saveRDS(mut1, "mut1.rds")
saveRDS(mut2, "mut2.rds")
