# TS19 Paper - Environment Setup Script
# Purpose: Install all required R packages for scRNA-seq analysis
# Author: Shan (with AI assistance)
# Date: 2025-10-18

cat("=== TS19 Paper Environment Setup ===\n")
cat("This script will install all required packages for scRNA-seq analysis\n\n")

# Function to install packages if not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  } else {
    cat(paste(pkg, "already installed\n"))
  }
}

# Function to install Bioconductor packages
install_bioc_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "from Bioconductor...\n"))
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(paste(pkg, "already installed\n"))
  }
}

cat("\n=== Installing CRAN Packages ===\n")

# Core packages
cran_packages <- c(
  "Seurat",           # Main scRNA-seq analysis
  "tidyverse",        # Data manipulation and visualization
  "dplyr",            # Data manipulation
  "ggplot2",          # Plotting
  "patchwork",        # Combine plots
  "cowplot",          # Publication-ready plots
  "RColorBrewer",     # Color palettes
  "viridis",          # Color scales
  "pheatmap",         # Heatmaps
  "Matrix",           # Sparse matrices
  "data.table",       # Fast data manipulation
  "readr",            # Read CSV files
  "writexl",          # Write Excel files
  "openxlsx"          # Read/write Excel files
)

for (pkg in cran_packages) {
  install_if_missing(pkg)
}

cat("\n=== Installing Bioconductor Packages ===\n")

# Bioconductor packages
bioc_packages <- c(
  "SingleCellExperiment",  # Single-cell data structure
  "scater",                # QC and visualization
  "scran",                 # Normalization
  "limma",                 # Differential expression
  "edgeR",                 # RNA-seq statistics
  "clusterProfiler",       # GO/KEGG enrichment
  "org.Mm.eg.db",          # Mouse gene annotations
  "fgsea",                 # Fast GSEA
  "ComplexHeatmap"         # Advanced heatmaps
)

for (pkg in bioc_packages) {
  install_bioc_if_missing(pkg)
}

cat("\n=== Verifying Installation ===\n")

# Test loading key packages
test_packages <- c("Seurat", "tidyverse", "SingleCellExperiment")
all_loaded <- TRUE

for (pkg in test_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("✓", pkg, "loaded successfully\n"))
  } else {
    cat(paste("✗", pkg, "failed to load\n"))
    all_loaded <- FALSE
  }
}

if (all_loaded) {
  cat("\n=== SUCCESS! All packages installed and loaded ===\n")
  cat("You can now run the TS19 analysis pipeline\n")
} else {
  cat("\n=== WARNING: Some packages failed to install ===\n")
  cat("Please check error messages above\n")
}

# Print session info
cat("\n=== Session Info ===\n")
sessionInfo()