# scripts/test_get_heart10x.R
suppressPackageStartupMessages({
  library(Seurat)
})

dir.create("data/test", recursive = TRUE, showWarnings = FALSE)
h5_url  <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_1k_v2/heart_1k_v2_filtered_feature_bc_matrix.h5"
h5_path <- "data/test/heart_1k_v2_filtered_feature_bc_matrix.h5"
rds_out <- "data/test/heart_1k_v2_seurat.rds"

if (!file.exists(h5_path)) {
  message("Downloading 10x heart_1k_v2 h5...")
  download.file(h5_url, destfile = h5_path, mode = "wb", quiet = TRUE)
  message("Saved: ", normalizePath(h5_path))
}

# Read 10x h5 and build a Seurat object
message("Reading 10x h5 -> Seurat...")
m <- Read10X_h5(h5_path)  # Will return a list if multiple modalities
if (is.list(m)) m <- m[["Gene Expression"]]
obj <- CreateSeuratObject(counts = m, project = "heart_1k_v2", min.cells = 3, min.features = 200)

saveRDS(obj, rds_out)
message("Test Seurat RDS saved to: ", normalizePath(rds_out))