## Define helper functions for sample processing 
# write h5 file function
write_h5 <- function(input_path, output_path, genome = "mm10") {
  # Convert 10x data to HDF5 format
  if (!file.exists(input_path)) stop("Input path does not exist: ", input_path)
  matrix_data <- Read10X(input_path)
  write10xCounts(output_path, matrix_data, type = "HDF5", genome = genome, 
                 version = "3", overwrite = TRUE, 
                 gene.id = rownames(matrix_data), 
                 gene.symbol = rownames(matrix_data))
}


# add mitochondrial percentage function
add_mito_percentage <- function(obj) {
  # Add mitochondrial gene percentage as metadata
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  return(obj)
}


# generate QC plot function
generate_QC_plot <- function(obj, filename) {
  # Generate QC violin plots
  pdf(filename, width = 8, height = 8)
  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()
}


# filter cells function
filter_cells <- function(obj) {
  # Filter cells based on QC metrics
  subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
}


# process data function
process_data <- function(obj) {
  # Normalize, identify variable features, scale data, and perform PCA
  obj <- NormalizeData(obj) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(features = rownames(obj)) %>%
    RunPCA(features = VariableFeatures(obj))
  return(obj)
}


# Cluster cells and perform UMAP dimensionality reduction
cluster_and_visualize <- function(obj, filename, dims = 1:20, resolution = 0.5) {
  obj <- FindNeighbors(obj, dims = dims)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, dims = dims)
  
  pdf(filename, width = 8, height = 8)
  DimPlot_scCustom(obj)
  dev.off()
  
  return(obj)
}