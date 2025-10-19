# TS19 Paper Repository

This repository contains data, scripts, and manuscript assets for the TS19 project.
- scripts/01_setup_environment.R  install/load R packages, set options/paths
- scripts/02_pipeline_full.R  end-to-end pipeline (from raw or processed inputs)
- scripts/03_compile_figures.R  compiles/collects figures/tables for the manuscript

## Getting Started
1. Run Rscript scripts/01_setup_environment.R
2. Run Rscript scripts/03_compile_figures.R (uses existing processed objects, if available)
3. (Optional) Run Rscript scripts/02_pipeline_full.R to reproduce from scratch
