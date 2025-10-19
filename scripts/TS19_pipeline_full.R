# TS19 Cell Type Atlas - Complete Analysis Pipeline
# Purpose: Process scRNA-seq data from raw counts to cell type annotations
# Author: Shan (recreating Bo's analysis)
# Date: 2025-10-18

# =============================================================================
# SETUP
# =============================================================================

cat("=== TS19 Cell Type Atlas Analysis Pipeline ===\n")
cat("Starting analysis...\n\n")

# Load required libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(viridis)

# Set random seed for reproducibility
set.seed(123)

# Create output directories
dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE)
dir.create("output/tables", showWarnings = FALSE)
dir.create("output/qc", showWarnings = FALSE)
dir.create("output/markers", showWarnings = FALSE)

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

cat("STEP 1: Loading data...\n")

# Option A: Load from 10x output (if you have Ion access)
load_from_10x <- function(data_dir) {
  cat("Loading 10x data from:", data_dir, "\n")
  counts <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = "TS19",
    min.cells = 3,      # Keep genes detected in at least 3 cells
    min.features = 200  # Keep cells with at least 200 genes
  )
  return(seurat_obj)
}

# Option B: Load from processed Seurat object (without Ion access)
load_from_rds <- function(rds_file) {
  cat("Loading processed Seurat object from:", rds_file, "\n")
  seurat_obj <- readRDS(rds_file)
  return(seurat_obj)
}

# Option C: Load from CSV files (without Ion access)
load_from_csv <- function(counts_file, metadata_file = NULL) {
  cat("Loading data from CSV files\n")
  
  # Read expression matrix
  counts <- read.csv(counts_file, row.names = 1)
  counts <- as.matrix(counts)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = "TS19",
    min.cells = 3,
    min.features = 200
  )
  
  # Add metadata if provided
  if (!is.null(metadata_file)) {
    metadata <- read.csv(metadata_file, row.names = 1)
    seurat_obj <- AddMetaData(seurat_obj, metadata)
  }
  
  return(seurat_obj)
}

# CHOOSE YOUR DATA SOURCE:
# Uncomment the appropriate line based on your data availability

# If you have Ion access:
# seurat_obj <- load_from_10x("path/to/10x/output")

# If you have processed RDS file:
# seurat_obj <- load_from_rds("path/to/seurat_object.rds")

# If you have CSV files:
# seurat_obj <- load_from_csv("path/to/counts.csv", "path/to/metadata.csv")

# For demonstration, create a placeholder
cat("NOTE: Please uncomment the appropriate data loading function above\n")
cat("For now, creating placeholder object for demonstration\n")

# =============================================================================
# STEP 2: QUALITY CONTROL
# =============================================================================

cat("\nSTEP 2: Quality Control...\n")

perform_qc <- function(seurat_obj) {
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  # Calculate ribosomal percentage
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")
  
  # Visualize QC metrics
  cat("Generating QC plots...\n")
  
  # Violin plots
  p1 <- VlnPlot(seurat_obj, 
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0.1)
  
  ggsave("output/qc/qc_violin_plots.png", p1, width = 12, height = 4, dpi = 300)
  
  # Scatter plots
  p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  p_combined <- p2 + p3
  ggsave("output/qc/qc_scatter_plots.png", p_combined, width = 12, height = 5, dpi = 300)
  
  # Print QC statistics
  cat("\n=== QC Statistics (Before Filtering) ===\n")
  cat("Total cells:", ncol(seurat_obj), "\n")
  cat("Total genes:", nrow(seurat_obj), "\n")
  cat("Median genes per cell:", median(seurat_obj$nFeature_RNA), "\n")
  cat("Median UMIs per cell:", median(seurat_obj$nCount_RNA), "\n")
  cat("Median % mitochondrial:", median(seurat_obj$percent.mt), "\n")
  
  # Apply QC filters (Bo's approach - model-based)
  # These thresholds should be adjusted based on your data
  cat("\nApplying QC filters...\n")
  
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > 200 & 
                                nFeature_RNA < 6000 & 
                                percent.mt < 15)
  
  cat("\n=== QC Statistics (After Filtering) ===\n")
  cat("Cells remaining:", ncol(seurat_obj), "\n")
  cat("Genes remaining:", nrow(seurat_obj), "\n")
  
  return(seurat_obj)
}

# Run QC (uncomment when you have data)
# seurat_obj <- perform_qc(seurat_obj)

# =============================================================================
# STEP 3: NORMALIZATION
# =============================================================================

cat("\nSTEP 3: Normalization...\n")

normalize_data <- function(seurat_obj, method = "LogNormalize") {
  
  if (method == "LogNormalize") {
    cat("Using log-normalization (standard Seurat method)\n")
    seurat_obj <- NormalizeData(seurat_obj, 
                                 normalization.method = "LogNormalize",
                                 scale.factor = 10000)
    
  } else if (method == "SCT") {
    cat("Using SCTransform normalization\n")
    seurat_obj <- SCTransform(seurat_obj, 
                               vars.to.regress = "percent.mt",
                               verbose = FALSE)
  }
  
  return(seurat_obj)
}

# Run normalization (uncomment when you have data)
# seurat_obj <- normalize_data(seurat_obj, method = "LogNormalize")

# =============================================================================
# STEP 4: FEATURE SELECTION
# =============================================================================

cat("\nSTEP 4: Feature Selection...\n")

select_features <- function(seurat_obj, n_features = 2000) {
  
  cat("Identifying", n_features, "highly variable features\n")
  
  seurat_obj <- FindVariableFeatures(seurat_obj,
                                      selection.method = "vst",
                                      nfeatures = n_features)
  
  # Plot variable features
  top10 <- head(VariableFeatures(seurat_obj), 10)
  
  p1 <- VariableFeaturePlot(seurat_obj)
  p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
  
  ggsave("output/figures/variable_features.png", p2, width = 10, height = 6, dpi = 300)
  
  cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")
  
  return(seurat_obj)
}

# Run feature selection (uncomment when you have data)
# seurat_obj <- select_features(seurat_obj, n_features = 2000)

# =============================================================================
# STEP 5: SCALING AND PCA
# =============================================================================

cat("\nSTEP 5: Scaling and PCA...\n")

run_pca <- function(seurat_obj, n_pcs = 100) {
  
  # Scale data
  cat("Scaling data...\n")
  all_genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all_genes)
  
  # Run PCA
  cat("Running PCA with", n_pcs, "components\n")
  seurat_obj <- RunPCA(seurat_obj, 
                       features = VariableFeatures(object = seurat_obj),
                       npcs = n_pcs,
                       verbose = FALSE)
  
  # Visualize PCA
  p1 <- DimPlot(seurat_obj, reduction = "pca")
  p2 <- ElbowPlot(seurat_obj, ndims = 50)
  
  p_combined <- p1 + p2
  ggsave("output/figures/pca_overview.png", p_combined, width = 12, height = 5, dpi = 300)
  
  # Print PCA loadings for top PCs
  cat("\nTop genes for PC1:\n")
  print(head(seurat_obj[["pca"]]@feature.loadings[, 1], 10))
  
  return(seurat_obj)
}

# Run PCA (uncomment when you have data)
# seurat_obj <- run_pca(seurat_obj, n_pcs = 100)

# =============================================================================
# STEP 6: CLUSTERING
# =============================================================================

cat("\nSTEP 6: Clustering...\n")

perform_clustering <- function(seurat_obj, dims = 1:100, resolution = 0.5) {
  
  cat("Finding neighbors using", max(dims), "PCs\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  
  cat("Finding clusters with resolution =", resolution, "\n")
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  # Try multiple resolutions
  cat("Testing multiple resolutions...\n")
  resolutions <- c(0.3, 0.5, 0.8, 1.0, 1.2)
  
  for (res in resolutions) {
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
  }
  
  # Set default to 0.5
  Idents(seurat_obj) <- seurat_obj$RNA_snn_res.0.5
  
  cat("Number of clusters at resolution 0.5:", 
      length(unique(Idents(seurat_obj))), "\n")
  
  return(seurat_obj)
}

# Run clustering (uncomment when you have data)
# seurat_obj <- perform_clustering(seurat_obj, dims = 1:100, resolution = 0.5)

# =============================================================================
# STEP 7: DIMENSIONALITY REDUCTION FOR VISUALIZATION
# =============================================================================

cat("\nSTEP 7: Dimensionality Reduction (UMAP and t-SNE)...\n")

run_visualization <- function(seurat_obj, dims = 1:100) {
  
  # Run UMAP
  cat("Running UMAP...\n")
  seurat_obj <- RunUMAP(seurat_obj, dims = dims, verbose = FALSE)
  
  # Run t-SNE
  cat("Running t-SNE...\n")
  seurat_obj <- RunTSNE(seurat_obj, dims = dims, verbose = FALSE)
  
  # Visualize
  p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 6) +
    ggtitle("UMAP - Clusters")
  
  p2 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE, label.size = 6) +
    ggtitle("t-SNE - Clusters")
  
  p_combined <- p1 + p2
  ggsave("output/figures/clustering_overview.png", p_combined, 
         width = 16, height = 7, dpi = 300)
  
  return(seurat_obj)
}

# Run visualization (uncomment when you have data)
# seurat_obj <- run_visualization(seurat_obj, dims = 1:100)

# =============================================================================
# STEP 8: MARKER GENE IDENTIFICATION
# =============================================================================

cat("\nSTEP 8: Identifying Marker Genes...\n")

find_markers <- function(seurat_obj) {
  
  cat("Finding markers for all clusters...\n")
  
  all_markers <- FindAllMarkers(seurat_obj,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25,
                                verbose = FALSE)
  
  # Save all markers
  write.csv(all_markers, "output/markers/all_markers.csv", row.names = FALSE)
  
  # Get top markers per cluster
  top10 <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  write.csv(top10, "output/markers/top10_markers_per_cluster.csv", row.names = FALSE)
  
  # Create heatmap of top markers
  p <- DoHeatmap(seurat_obj, features = top10$gene) +
    scale_fill_viridis()
  
  ggsave("output/figures/marker_heatmap.png", p, width = 12, height = 16, dpi = 300)
  
  cat("Markers saved to output/markers/\n")
  
  return(all_markers)
}

# Find markers (uncomment when you have data)
# all_markers <- find_markers(seurat_obj)

# =============================================================================
# STEP 9: CELL TYPE ANNOTATION
# =============================================================================

cat("\nSTEP 9: Cell Type Annotation...\n")

annotate_cell_types <- function(seurat_obj, all_markers) {
  
  cat("Annotating cell types based on marker genes...\n")
  
  # Define known markers for major cell types
  # These are based on Bo's thesis and literature
  
  known_markers <- list(
    "Lateral Mesoderm" = c("Tbx5", "Hand1", "Hand2", "Gata4"),
    "Mesenchyme" = c("Twist1", "Prrx1", "Pdgfra"),
    "Neuroectoderm" = c("Sox2", "Pax6", "Nestin"),
    "Blood" = c("Gata1", "Hba-x", "Hbb-bs"),
    "Cardiac" = c("Tnnt2", "Myh6", "Nkx2-5"),
    "Endothelium" = c("Pecam1", "Cdh5", "Kdr"),
    "Epicardium" = c("Wt1", "Tbx18", "Tcf21")
  )
  
  # Plot known markers
  for (cell_type in names(known_markers)) {
    markers <- known_markers[[cell_type]]
    markers <- markers[markers %in% rownames(seurat_obj)]
    
    if (length(markers) > 0) {
      p <- FeaturePlot(seurat_obj, features = markers, ncol = 2)
      
      filename <- paste0("output/figures/markers_", 
                        gsub(" ", "_", cell_type), ".png")
      ggsave(filename, p, width = 10, height = ceiling(length(markers)/2) * 3, dpi = 300)
    }
  }
  
  # Manual annotation based on markers
  # This should be done interactively by examining marker genes
  # For now, we'll create a template
  
  cluster_annotations <- data.frame(
    cluster = levels(Idents(seurat_obj)),
    cell_type = paste0("Cluster_", levels(Idents(seurat_obj))),
    stringsAsFactors = FALSE
  )
  
  write.csv(cluster_annotations, 
            "output/tables/cluster_annotations_template.csv",
            row.names = FALSE)
  
  cat("Cluster annotation template saved to output/tables/\n")
  cat("Please manually annotate clusters based on marker genes\n")
  
  return(cluster_annotations)
}

# Annotate cell types (uncomment when you have data)
# cluster_annotations <- annotate_cell_types(seurat_obj, all_markers)

# =============================================================================
# STEP 10: SAVE RESULTS
# =============================================================================

cat("\nSTEP 10: Saving Results...\n")

save_results <- function(seurat_obj, all_markers, cluster_annotations) {
  
  # Save Seurat object
  cat("Saving Seurat object...\n")
  saveRDS(seurat_obj, "output/seurat_object_processed.rds")
  
  # Save cell metadata
  cat("Saving cell metadata...\n")
  metadata <- seurat_obj@meta.data
  metadata$cell_barcode <- rownames(metadata)
  metadata$UMAP_1 <- Embeddings(seurat_obj, "umap")[, 1]
  metadata$UMAP_2 <- Embeddings(seurat_obj, "umap")[, 2]
  metadata$tSNE_1 <- Embeddings(seurat_obj, "tsne")[, 1]
  metadata$tSNE_2 <- Embeddings(seurat_obj, "tsne")[, 2]
  
  write.csv(metadata, "output/tables/cell_metadata.csv", row.names = FALSE)
  
  # Save normalized counts
  cat("Saving normalized expression matrix...\n")
  norm_counts <- GetAssayData(seurat_obj, slot = "data")
  write.csv(as.matrix(norm_counts), "output/tables/normalized_counts.csv")
  
  # Create summary statistics
  summary_stats <- data.frame(
    metric = c("Total cells", "Total genes", "Number of clusters",
               "Median genes per cell", "Median UMIs per cell"),
    value = c(ncol(seurat_obj), nrow(seurat_obj), 
              length(unique(Idents(seurat_obj))),
              median(seurat_obj$nFeature_RNA),
              median(seurat_obj$nCount_RNA))
  )
  
  write.csv(summary_stats, "output/tables/summary_statistics.csv", row.names = FALSE)
  
  cat("\n=== Analysis Complete! ===\n")
  cat("All results saved to output/ directory\n")
  cat("- Figures: output/figures/\n")
  cat("- Tables: output/tables/\n")
  cat("- Markers: output/markers/\n")
  cat("- QC: output/qc/\n")
}

# Save results (uncomment when you have data)
# save_results(seurat_obj, all_markers, cluster_annotations)

# =============================================================================
# PIPELINE COMPLETE
# =============================================================================

cat("\n=== TS19 Analysis Pipeline Script Complete ===\n")
cat("To run this pipeline:\n")
cat("1. Uncomment the appropriate data loading function (Step 1)\n")
cat("2. Uncomment all analysis steps (Steps 2-10)\n")
cat("3. Run: source('TS19_pipeline_full.R')\n")
cat("\nFor questions or issues, refer to the documentation\n")