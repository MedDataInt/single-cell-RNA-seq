# Project: MSU_xx_Single-cell RNA-seq Analysis
# Description: This script performs primiarly single-cell RNA-seq analysis on MSU xx project.

# Load required libraries
library(Seurat)
library(dplyr)
library(magrittr)
library(Matrix)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(DropletUtils)
library(here)


# Set working directory dynamically
project_dir <- here::here()
setwd(project_dir)

# Define helper functions-----------------------------------------------------------------

#' Convert 10x data to HDF5 format
#' @param input_path Path to the 10x data directory.
#' @param output_path Path to save the HDF5 file.
#' @param genome Genome assembly (default: "mm10").
#' @examples
#' write_h5("/path/to/data", "/path/to/output.h5", genome = "hg38")
write_h5 <- function(input_path, output_path, genome = "mm10") {
  if (!dir.exists(input_path)) stop(paste0("Input path does not exist: ", input_path))
  tryCatch({
    matrix_data <- Read10X(input_path)
    write10xCounts(output_path, matrix_data, type = "HDF5", genome = genome,
                   version = "3", overwrite = TRUE,
                   gene.id = rownames(matrix_data),
                   gene.symbol = rownames(matrix_data))
    message(paste0("Successfully wrote HDF5 to ", output_path))
  }, error = function(e) {
    stop(paste0("Error writing HDF5: ", e$message))
  })
}

#' Add mitochondrial gene percentage as metadata.
#' @param obj Seurat object.
add_mito_percentage <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  return(obj)
}

#' Generate QC violin plots.
#' @param obj Seurat object.
#' @param filename Output PDF filename.
generate_QC_plot <- function(obj, filename) {
  pdf(filename, width = 12, height = 8) 
  print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1) + NoLegend())
  dev.off()
}

#' Filter cells based on QC metrics.
#' @param obj Seurat object.
filter_cells <- function(obj, min_features = 200, max_features = 2500, max_percent_mt = 5) {
  subset(obj, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_percent_mt)
}

#' Process data (normalization, feature selection, scaling, PCA).
#' @param obj Seurat object.
process_data <- function(obj) {
  obj <- NormalizeData(obj) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(features = rownames(obj)) %>%
    RunPCA(features = VariableFeatures(obj))
  return(obj)
}



# process analysis -----------------------------------------------------------------
# Data loading and preprocessing
data_dirs <- list(scramble = file.path(project_dir, "scramble/filtered_feature_bc_matrix/"),
                  mut1 = file.path(project_dir, "mut1/filtered_feature_bc_matrix/"),
                  mut2 = file.path(project_dir, "mut2/filtered_feature_bc_matrix/"))

h5_files <- file.path(project_dir, paste0(names(data_dirs), "_filtered_feature_bc_matrix.h5"))

mapply(write_h5, data_dirs, h5_files)

seurat_objects <- lapply(h5_files, Read10X_h5)
names(seurat_objects) <- names(data_dirs)

seurat_objects <- lapply(names(seurat_objects), function(name) {
  obj <- CreateSeuratObject(seurat_objects[[name]], project = name, min.cells = 3, min.features = 200)
  obj <- add_mito_percentage(obj)
  generate_QC_plot(obj, paste0(name, "_QC.pdf"))
  obj <- filter_cells(obj)
  process_data(obj)
  obj
})
names(seurat_objects) <- names(data_dirs)

# Save individual Seurat objects
save_rds_files <- file.path(project_dir, paste0(names(seurat_objects), ".rds"))
mapply(saveRDS, seurat_objects, save_rds_files)

# Merge datasets and perform integrated analysis
mut_combined <- merge(seurat_objects$mut1, y = seurat_objects$mut2, add.cell.ids = c("mut1", "mut2"), project = "mut")
all_normalized <- merge(mut_combined, y = seurat_objects$scramble, add.cell.ids = c("mut", "scramble"), project = "combined")

# Process and visualize unintegrated data
all_normalized <- all_normalized %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:30, reduction = "pca")

pdf(file.path(project_dir,"UMAP_all_unintegrated.pdf"), width = 16, height = 8)
print(DimPlot(all_normalized, reduction = "umap", group.by = c("orig.ident", "seurat_clusters")))
dev.off()

# Integration using CCA
all_normalized <- IntegrateLayers(object = all_normalized,
                                  method = CCAIntegration,
                                  orig.reduction = "pca",
                                  new.reduction = "integrated.cca",
                                  verbose = FALSE)

# Clustering and visualization for integrated data
all_normalized <- all_normalized %>%
  FindNeighbors(reduction = "integrated.cca", dims = 1:30) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:30, reduction = "integrated.cca")

pdf(file.path(project_dir,"UMAP_all_integrated.pdf"), width = 16, height = 8)
print(DimPlot_scCustom(all_normalized, reduction = "umap", group.by = "seurat_clusters"))
dev.off()

# Save final integrated object
saveRDS(all_normalized, file.path(project_dir,"all_normalized_integrated.rds"))

message("Clustering analysis complete.")