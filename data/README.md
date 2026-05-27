# `data/` contents and lineage

These objects are **not** stored in git (the folder is gitignored); download them from the
[Dropbox link](../README.md) and drop them here before rendering the notebooks.

The table maps each file to the notebook that **produces** it and the notebook(s) that
**ingest** it. Notebook titles come from each `.qmd` YAML header.

| File | Produced by | Producer title | Ingested by | Ingester title |
|------|-------------|----------------|-------------|----------------|
| `20200620-list_of_Seurat_objects.Rdata` | `20250620-seuratObject.qmd` | Seurat Object | `20250620-SCT&intehra.qmd` | SCTransform & integration |
| `20250620-Seurat_merged.Rdata` | `20250620-seuratObject.qmd` | Seurat Object | `20250620-seuratObject.qmd` | Seurat Object |
| `20250619-seurat_integrated_nFeat>2000.Rdata` | `20250620-SCT&intehra.qmd` | SCTransform & integration | `20250620-ymas.qmd` | Unsupervised Plots |
| `20250620-seurat_integrated_UMAP_nFeat>600.Rds` | `20250620-ymas.qmd` | Unsupervised Plots | `20250620-ScType.qmd` | ScType (Automatic annotation) |
| `20250620-seurat_integrated_UMAP_nFeat-ScTypeannotated.Rdata` | `20250620-ScType.qmd` | ScType (Automatic annotation) | `20250620-ScType.qmd` | ScType (Automatic annotation) |
| `20250812-alldata_sctype_annotated_CHLOE.Rdata` | external (Chloe), not in repo | – | `20250728-markers.qmd` | Find Markers |
| `20250812-all_markers.csv` | `20250728-markers.qmd` † | Find Markers | `20250728-markers.qmd` | Find Markers |
| `20250812-alldata_sctype_annotated_CHLOE_Lance.Rdata` | external (Chloe/Lance), not in repo | – | `20250818-microglia1.1.qmd`, `20250825-pseudobulk.qmd` | Microglial subpopulations; Pseudobulk |
| `20250818-microglia_reclustered_reintegrated.Rdata` | `20250818-microglia1.1.qmd` † | Microglial subpopulations | `20250818-microglia1.1.qmd`, `20250818-microglia2.qmd`, `20250829-trajectory.qmd` ‡, `20250903-scanpy.qmd` ‡ | Microglial subpopulations; Differences in microglia by sample & treatment; Trajectory analysis; Scanpy |
| `20250818-microglia_markers.json` | `20250818-microglia1.1.qmd` † | Microglial subpopulations | – (pasted into an external LLM annotation prompt) | – |
| `20250826-alldata_sctype_annotated_CHLOE_Lance_reintegrated.Rdata` | `20250825-pseudobulk.qmd` † | Pseudobulk | `20250825-pseudobulk.qmd` | Pseudobulk |
| `microglia.h5Seurat` | `20250903-scanpy.qmd` ‡ | Scanpy | `20250903-scanpy.qmd` ‡ | Scanpy |
| `microglia.h5ad` | `20250903-scanpy.qmd` ‡ | Scanpy | external (scanpy / Python) | – |

## Legend

- `–` no producer/consumer in the repo.
- `†` the `save()`/`write` line is currently commented out in the notebook; the file exists from a prior run.
- `‡` the notebook is gitignored (`20250829-trajectory.qmd`, `20250903-scanpy.qmd`), so it is not a committed file.
- External objects (`CHLOE`, `Lance`) were produced by a collaborator outside this repository.
