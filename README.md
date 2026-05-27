# ARIA single-cell RNA-seq

## Summary

Amyloid-related imaging abnormalities (ARIA) are the principal dose-limiting toxicity of
anti-amyloid immunotherapy in Alzheimer's disease, yet the cellular response in the brain is
poorly resolved. 
This project profiles the brain of aged female APOE4 (E4/E4) amyloid (FAD)
mice treated with the anti-amyloid antibody aducanumab versus an IgG isotype control (a saline
arm was collected but excluded from the integrated analysis), using droplet-based single-cell
RNA sequencing (10x Genomics, three mice per arm).

The rendered analysis is published at https://ljohnsonlab.github.io/ARIA_scRNAseq/.

## Setup and reproduction

To re-run the analysis after cloning:

1. **Install R packages.** From the repo root, in R: `source("install.R")`. This covers the
   nine book chapters (CRAN, Bioconductor, and three GitHub packages). It does not cover the
   trajectory or scanpy steps, which are out of scope and gitignored.
2. **Download the data** into `data/` (see below). All chapters load from these saved objects;
   the raw 10x ingestion chunks in `20250615-cellRanger.qmd` and `20250620-seuratObject.qmd`
   are gated `eval: false` and point to inputs that are not part of this repo, so the pipeline
   is reproduced from the saved objects rather than rebuilt from the CellRanger matrices.
3. **Render** with `quarto render` (chapter order is defined in `_quarto.yml`), or open
   individual `.qmd` files in RStudio. `20250620-ScType.qmd` needs internet access because it
   `source()`s the ScType helper from a URL at render time.

See `data/README.md` for the file-by-file lineage (which chapter produces and consumes each object).

### Tested with

R 4.5.2 on macOS. Key package versions: Seurat 5.4.0, SeuratObject 5.3.0, tidyverse 2.0.0,
glmGamPoi 1.22.0, EnhancedVolcano 1.28.2, presto 1.0.0, SeuratWrappers 0.4.0, scCustomize 3.3.0,
compareGroups 4.10.2, lme4 1.1.38. Other versions may produce slightly different numbers.

## Data Folder with seurat objects
The data folder contains all the relevant datasets required for this project. Due to file size constraints, the actual data files are hosted externally on Dropbox.

__How to Access the Data__
Download the full dataset from the following [Dropbox link](https://www.dropbox.com/scl/fo/m6naiy85qayf6y3rtfvz6/AJ1G9P1TkVL_aQTn7urEj_o?rlkey=f5x824ae2vdd3v6nspslsz0y7&dl=0)
and place the files inside the data folder _before running or analyzing the project_.

