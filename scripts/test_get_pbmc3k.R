# scripts/test_get_pbmc3k.R
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)  # provides pbmc3k dataset
})

dir.create("data/test", recursive = TRUE, showWarnings = FALSE)
rds_out <- "data/test/pbmc3k_seurat.rds"

# Install the dataset if not already present
if (!"pbmc3k.SeuratData" %in% rownames(installed.packages())) {
  SeuratData::InstallData("pbmc3k")
}

# Load as a Seurat object
obj <- SeuratData::LoadData("pbmc3k")  # returns a Seurat object

# Save for our pipeline
saveRDS(obj, rds_out)
message("PBMC3k Seurat RDS saved to: ", normalizePath(rds_out))