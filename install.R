# Package installation for the ARIA scRNA-seq book chapters.
# Run once in a fresh R session (R >= 4.5) before rendering the .qmd chapters.
# Scope: the 9 committed book chapters listed in _quarto.yml. The gitignored
# trajectory and scanpy steps are out of scope and need extra packages.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

installed <- rownames(installed.packages())

# CRAN
cran <- c(
  "tidyverse", "Seurat", "patchwork", "kableExtra", "corrplot", "future",
  "compareGroups", "openxlsx", "ggtext", "gtable", "jsonlite", "pheatmap",
  "lme4", "scCustomize", "clustree", "HGNChelper"
)
install.packages(setdiff(cran, installed))

# Bioconductor
bioc <- c("glmGamPoi", "EnhancedVolcano")
BiocManager::install(setdiff(bioc, installed), update = FALSE, ask = FALSE)

# GitHub (not on CRAN/Bioconductor)
if (!requireNamespace("presto", quietly = TRUE))         remotes::install_github("immunogenomics/presto")
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) remotes::install_github("satijalab/seurat-wrappers")
if (!requireNamespace("monocle3", quietly = TRUE))       remotes::install_github("cole-trapnell-lab/monocle3")

# Note: 20250620-ScType.qmd source()s the ScType helper from a GitHub raw URL at
# render time, so that chapter additionally needs internet access (no install).
