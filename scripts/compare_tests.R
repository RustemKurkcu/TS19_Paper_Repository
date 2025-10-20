suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/compare_tests.R --heart=<path.rds> --pbmc=<path.rds>")
}

parse_arg <- function(flag) {
  x <- args[grepl(paste0("^", flag, "="), args)]
  if (length(x) == 0) return(NA_character_)
  sub(paste0("^", flag, "="), "", x)
}

heart_rds <- parse_arg("--heart")
pbmc_rds  <- parse_arg("--pbmc")

if (!file.exists(heart_rds)) stop("Heart RDS not found: ", heart_rds)
if (!file.exists(pbmc_rds))  stop("PBMC RDS not found: ", pbmc_rds)

outdir <- "comparisons"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Helper: ensure UMAP exists (compute minimal workflow if needed)
ensure_umap <- function(obj) {
  if (!"umap" %in% names(Reductions(obj))) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
    obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  }
  obj
}

# Load objects
heart <- readRDS(heart_rds)
pbmc  <- readRDS(pbmc_rds)

heart <- ensure_umap(heart)
pbmc  <- ensure_umap(pbmc)

# UMAPs
p_heart <- DimPlot(heart, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Heart (test)")
p_pbmc  <- DimPlot(pbmc,  reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("PBMC3k (test)")
umap_panel <- p_heart + p_pbmc
ggsave(file.path(outdir, "comparison_umap.png"), umap_panel, width = 14, height = 6, dpi = 300)

# Cluster proportion barplots
prop_df <- function(obj, tag) {
  tibble(cluster = as.character(Idents(obj))) |>
    count(cluster) |>
    mutate(prop = n / sum(n), dataset = tag)
}
df <- bind_rows(prop_df(heart, "Heart"), prop_df(pbmc, "PBMC3k"))

p_bar <- ggplot(df, aes(x = cluster, y = prop, fill = dataset)) +
  geom_col(position = "dodge") +
  labs(title = "Cluster proportions", x = "Cluster", y = "Fraction") +
  theme_bw(base_size = 12)
ggsave(file.path(outdir, "comparison_cluster_proportions.png"), p_bar, width = 12, height = 5, dpi = 300)

cat("Wrote:\n  - comparisons/comparison_umap.png\n  - comparisons/comparison_cluster_proportions.png\n")