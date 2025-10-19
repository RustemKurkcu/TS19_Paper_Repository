# TS19 Paper - Figure and Data Compilation Script
# Purpose: Compile all important figures and data for the manuscript
# Author: Shan
# Date: 2025-10-18

# =============================================================================
# SETUP
# =============================================================================

library(tidyverse)
library(Seurat)
library(patchwork)
library(cowplot)

cat("=== TS19 Figure Compilation Script ===\n\n")

# Create manuscript figure directory
dir.create("manuscript_figures", showWarnings = FALSE)
dir.create("manuscript_figures/main", showWarnings = FALSE)
dir.create("manuscript_figures/supplemental", showWarnings = FALSE)
dir.create("manuscript_tables", showWarnings = FALSE)

# =============================================================================
# FIGURE 1: EXPERIMENTAL DESIGN AND PIPELINE
# =============================================================================

cat("Compiling Figure 1: Experimental Design and Pipeline\n")

compile_figure1 <- function() {
  # This figure is typically created manually or with BioRender
  # Components needed:
  # A) Schematic of E11.5 embryo dissection
  # B) 10x Genomics workflow
  # C) Computational pipeline overview
  # D) Quality control metrics
  
  cat("Figure 1 components:\n")
  cat("  A) Embryo dissection schematic (from Bo's slides)\n")
  cat("  B) 10x workflow diagram\n")
  cat("  C) Analysis pipeline flowchart\n")
  cat("  D) QC metrics (violin plots)\n")
  
  # Copy QC plots if they exist
  if (file.exists("output/qc/qc_violin_plots.png")) {
    file.copy("output/qc/qc_violin_plots.png",
              "manuscript_figures/main/Figure1D_QC_metrics.png",
              overwrite = TRUE)
  }
  
  # Look for existing pipeline figures in Bo's files
  bo_fig_paths <- c(
    "ts19_project/TS19 manuscript/Manuscript figures/Fig2 Dissection of embryos",
    "ch2_ts19_celltype/Figure pile"
  )
  
  for (path in bo_fig_paths) {
    if (dir.exists(path)) {
      cat("  Found Bo's figures in:", path, "\n")
      # Copy relevant files
      fig_files <- list.files(path, pattern = "\\.(png|pdf|tiff)$", 
                             full.names = TRUE)
      if (length(fig_files) > 0) {
        cat("  Copying", length(fig_files), "files\n")
        for (f in fig_files) {
          file.copy(f, "manuscript_figures/main/", overwrite = TRUE)
        }
      }
    }
  }
}

compile_figure1()

# =============================================================================
# FIGURE 2: CLUSTERING AND CELL TYPE IDENTIFICATION
# =============================================================================

cat("\nCompiling Figure 2: Clustering and Cell Type Identification\n")

compile_figure2 <- function(seurat_obj = NULL) {
  # Components:
  # A) UMAP with clusters
  # B) t-SNE with clusters
  # C) Cluster dendrogram
  # D) Number of cells per cluster
  
  if (!is.null(seurat_obj)) {
    # Generate fresh plots
    p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, 
                  label.size = 6, pt.size = 0.5) +
      ggtitle("UMAP - Cell Clusters") +
      theme_cowplot()
    
    p2 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE,
                  label.size = 6, pt.size = 0.5) +
      ggtitle("t-SNE - Cell Clusters") +
      theme_cowplot()
    
    # Cell counts per cluster
    cluster_counts <- table(Idents(seurat_obj))
    p3 <- ggplot(data.frame(cluster = names(cluster_counts),
                           count = as.numeric(cluster_counts)),
                aes(x = cluster, y = count, fill = cluster)) +
      geom_bar(stat = "identity") +
      theme_cowplot() +
      theme(legend.position = "none") +
      labs(x = "Cluster", y = "Number of Cells",
           title = "Cells per Cluster")
    
    # Combine and save
    p_combined <- (p1 | p2) / p3
    ggsave("manuscript_figures/main/Figure2_Clustering.png",
           p_combined, width = 16, height = 12, dpi = 300)
    
  } else {
    # Copy from Bo's existing figures
    bo_fig_paths <- c(
      "ts19_project/TS19 manuscript/Manuscript figures/Fig6 Clustering results",
      "ch2_ts19_celltype/Figure pile"
    )
    
    for (path in bo_fig_paths) {
      if (dir.exists(path)) {
        cat("  Found Bo's clustering figures in:", path, "\n")
        fig_files <- list.files(path, pattern = "\\.(png|pdf|tiff)$",
                               full.names = TRUE)
        for (f in fig_files) {
          file.copy(f, "manuscript_figures/main/", overwrite = TRUE)
        }
      }
    }
  }
}

compile_figure2()

# =============================================================================
# FIGURE 3: MAJOR LINEAGES
# =============================================================================

cat("\nCompiling Figure 3: Major Lineages\n")

compile_figure3 <- function(seurat_obj = NULL) {
  # Components:
  # A) Lateral Mesoderm lineage
  # B) Mesenchyme subclusters
  # C) Neuroectoderm subclusters
  # D) Blood lineage
  
  # Copy from Bo's existing figures
  bo_fig_paths <- c(
    "ts19_project/TS19 manuscript/Manuscript figures/Fig5 LM lineage tree",
    "ts19_project/TS19 manuscript/Manuscript figures/Fig7 ME subclustering",
    "ts19_project/TS19 manuscript/Manuscript figures/Fig8 NE subclustering",
    "ts19_project/TS19 manuscript/Manuscript figures/Fig9 Blood lineage",
    "ch2_ts19_celltype/LM lineage tree",
    "ch2_ts19_celltype/blood lineage"
  )
  
  for (path in bo_fig_paths) {
    if (dir.exists(path)) {
      cat("  Found Bo's lineage figures in:", path, "\n")
      fig_files <- list.files(path, pattern = "\\.(png|pdf|tiff)$",
                             full.names = TRUE, recursive = TRUE)
      for (f in fig_files) {
        basename_f <- basename(f)
        file.copy(f, paste0("manuscript_figures/main/Figure3_", basename_f),
                 overwrite = TRUE)
      }
    }
  }
}

compile_figure3()

# =============================================================================
# FIGURE 4: MARKER GENES
# =============================================================================

cat("\nCompiling Figure 4: Marker Genes\n")

compile_figure4 <- function(seurat_obj = NULL) {
  # Components:
  # A) Heatmap of top marker genes
  # B) Dot plot of key markers
  # C) Feature plots of selected markers
  # D) Violin plots of key markers
  
  if (!is.null(seurat_obj)) {
    # Read marker genes
    if (file.exists("output/markers/top10_markers_per_cluster.csv")) {
      top_markers <- read.csv("output/markers/top10_markers_per_cluster.csv")
      
      # Create heatmap
      p1 <- DoHeatmap(seurat_obj, features = top_markers$gene) +
        scale_fill_viridis_c()
      
      ggsave("manuscript_figures/main/Figure4A_Marker_Heatmap.png",
             p1, width = 12, height = 16, dpi = 300)
      
      # Create dot plot
      p2 <- DotPlot(seurat_obj, features = unique(top_markers$gene[1:50])) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      
      ggsave("manuscript_figures/main/Figure4B_Marker_DotPlot.png",
             p2, width = 14, height = 8, dpi = 300)
    }
  } else {
    # Copy from Bo's existing figures
    bo_fig_paths <- c(
      "ts19_project/TS19 manuscript/Manuscript figures/Sup figure Marker projection",
      "ch2_ts19_celltype/Figure pile"
    )
    
    for (path in bo_fig_paths) {
      if (dir.exists(path)) {
        cat("  Found Bo's marker figures in:", path, "\n")
        fig_files <- list.files(path, pattern = "\\.(png|pdf|tiff)$",
                               full.names = TRUE, recursive = TRUE)
        for (f in fig_files) {
          basename_f <- basename(f)
          file.copy(f, paste0("manuscript_figures/main/Figure4_", basename_f),
                   overwrite = TRUE)
        }
      }
    }
  }
}

compile_figure4()

# =============================================================================
# FIGURE 5: CARDIAC DEVELOPMENT
# =============================================================================

cat("\nCompiling Figure 5: Cardiac Development\n")

compile_figure5 <- function() {
  # Components:
  # A) Cardiac dumbbell structure
  # B) Cardiac marker expression
  # C) Developmental trajectory
  
  bo_fig_paths <- c(
    "ts19_project/TS19 manuscript/Manuscript figures/Fig10 Cardiac dumbbell",
    "ch2_ts19_celltype/89.Cardiac_dumbell_individual_PCA"
  )
  
  for (path in bo_fig_paths) {
    if (dir.exists(path)) {
      cat("  Found Bo's cardiac figures in:", path, "\n")
      fig_files <- list.files(path, pattern = "\\.(png|pdf|tiff)$",
                             full.names = TRUE, recursive = TRUE)
      for (f in fig_files) {
        basename_f <- basename(f)
        file.copy(f, paste0("manuscript_figures/main/Figure5_", basename_f),
                 overwrite = TRUE)
      }
    }
  }
}

compile_figure5()

# =============================================================================
# FIGURE 6: HOX GENES AND SPATIAL PATTERNING
# =============================================================================

cat("\nCompiling Figure 6: Hox Genes and Spatial Patterning\n")

compile_figure6 <- function() {
  # Components:
  # A) Hox gene expression patterns
  # B) Spatial organization
  
  bo_fig_paths <- c(
    "ts19_project/TS19 manuscript/Manuscript figures/Fig11 Hox genes"
  )
  
  for (path in bo_fig_paths) {
    if (dir.exists(path)) {
      cat("  Found Bo's Hox gene figures in:", path, "\n")
      fig_files <- list.files(path, pattern = "\\.(png|pdf|tiff)$",
                             full.names = TRUE, recursive = TRUE)
      for (f in fig_files) {
        basename_f <- basename(f)
        file.copy(f, paste0("manuscript_figures/main/Figure6_", basename_f),
                 overwrite = TRUE)
      }
    }
  }
}

compile_figure6()

# =============================================================================
# SUPPLEMENTAL FIGURES
# =============================================================================

cat("\nCompiling Supplemental Figures\n")

compile_supplemental_figures <- function() {
  # Supplemental Figure 1: GSEA results
  # Supplemental Figure 2: Additional marker projections
  # Supplemental Figure 3: Comparison with other atlases
  
  bo_fig_paths <- c(
    "ts19_project/TS19 manuscript/Manuscript figures/Sup figure GSEA results",
    "ts19_project/TS19 manuscript/Manuscript figures/Sup figure Marker projection",
    "ch2_ts19_celltype/GSEA output"
  )
  
  for (path in bo_fig_paths) {
    if (dir.exists(path)) {
      cat("  Found Bo's supplemental figures in:", path, "\n")
      fig_files <- list.files(path, pattern = "\\.(png|pdf|tiff)$",
                             full.names = TRUE, recursive = TRUE)
      for (f in fig_files) {
        basename_f <- basename(f)
        file.copy(f, paste0("manuscript_figures/supplemental/SupFig_", basename_f),
                 overwrite = TRUE)
      }
    }
  }
}

compile_supplemental_figures()

# =============================================================================
# COMPILE TABLES
# =============================================================================

cat("\nCompiling Tables\n")

compile_tables <- function() {
  # Table 1: Summary statistics
  # Table 2: Cluster annotations
  # Table 3: Top marker genes per cluster
  
  # Copy existing tables if available
  if (file.exists("output/tables/summary_statistics.csv")) {
    file.copy("output/tables/summary_statistics.csv",
              "manuscript_tables/Table1_Summary_Statistics.csv",
              overwrite = TRUE)
  }
  
  if (file.exists("output/tables/cluster_annotations_template.csv")) {
    file.copy("output/tables/cluster_annotations_template.csv",
              "manuscript_tables/Table2_Cluster_Annotations.csv",
              overwrite = TRUE)
  }
  
  if (file.exists("output/markers/top10_markers_per_cluster.csv")) {
    file.copy("output/markers/top10_markers_per_cluster.csv",
              "manuscript_tables/Table3_Top_Markers.csv",
              overwrite = TRUE)
  }
  
  # Look for Bo's existing tables
  bo_table_paths <- c(
    "ch2_ts19_celltype/DGE",
    "ch2_ts19_celltype/cell cluster annotation"
  )
  
  for (path in bo_table_paths) {
    if (dir.exists(path)) {
      cat("  Found Bo's tables in:", path, "\n")
      table_files <- list.files(path, pattern = "\\.(csv|xlsx|txt)$",
                               full.names = TRUE, recursive = TRUE)
      for (f in table_files) {
        basename_f <- basename(f)
        file.copy(f, paste0("manuscript_tables/Bo_", basename_f),
                 overwrite = TRUE)
      }
    }
  }
}

compile_tables()

# =============================================================================
# CREATE FIGURE MANIFEST
# =============================================================================

cat("\nCreating Figure Manifest\n")

create_manifest <- function() {
  # List all compiled figures
  main_figs <- list.files("manuscript_figures/main", full.names = FALSE)
  supp_figs <- list.files("manuscript_figures/supplemental", full.names = FALSE)
  tables <- list.files("manuscript_tables", full.names = FALSE)
  
  manifest <- data.frame(
    Type = c(rep("Main Figure", length(main_figs)),
             rep("Supplemental Figure", length(supp_figs)),
             rep("Table", length(tables))),
    Filename = c(main_figs, supp_figs, tables),
    Location = c(rep("manuscript_figures/main", length(main_figs)),
                 rep("manuscript_figures/supplemental", length(supp_figs)),
                 rep("manuscript_tables", length(tables)))
  )
  
  write.csv(manifest, "manuscript_figure_manifest.csv", row.names = FALSE)
  
  cat("\n=== Figure Compilation Summary ===\n")
  cat("Main Figures:", length(main_figs), "\n")
  cat("Supplemental Figures:", length(supp_figs), "\n")
  cat("Tables:", length(tables), "\n")
  cat("\nManifest saved to: manuscript_figure_manifest.csv\n")
  
  return(manifest)
}

manifest <- create_manifest()

# =============================================================================
# PRINT USAGE GUIDE
# =============================================================================

cat("\n=== Figure Usage Guide ===\n")
cat("\nMAIN FIGURES (for manuscript):\n")
cat("Figure 1: Experimental Design and Pipeline\n")
cat("  - Use for Introduction/Methods\n")
cat("  - Shows experimental workflow\n\n")

cat("Figure 2: Clustering and Cell Type Identification\n")
cat("  - Use for Results (main clustering)\n")
cat("  - Shows UMAP/t-SNE with clusters\n\n")

cat("Figure 3: Major Lineages\n")
cat("  - Use for Results (lineage analysis)\n")
cat("  - Shows LM, ME, NE, Blood lineages\n\n")

cat("Figure 4: Marker Genes\n")
cat("  - Use for Results (cell type validation)\n")
cat("  - Shows heatmaps and marker expression\n\n")

cat("Figure 5: Cardiac Development\n")
cat("  - Use for Results (cardiac analysis)\n")
cat("  - Shows cardiac dumbbell structure\n\n")

cat("Figure 6: Hox Genes and Spatial Patterning\n")
cat("  - Use for Results (spatial organization)\n")
cat("  - Shows Hox gene expression patterns\n\n")

cat("\nSUPPLEMENTAL FIGURES:\n")
cat("- GSEA enrichment results\n")
cat("- Additional marker projections\n")
cat("- Quality control details\n")
cat("- Comparison with other atlases\n\n")

cat("\nTABLES:\n")
cat("Table 1: Summary Statistics\n")
cat("Table 2: Cluster Annotations\n")
cat("Table 3: Top Marker Genes per Cluster\n\n")

cat("=== Compilation Complete! ===\n")
cat("All figures and tables are in:\n")
cat("  - manuscript_figures/main/\n")
cat("  - manuscript_figures/supplemental/\n")
cat("  - manuscript_tables/\n")