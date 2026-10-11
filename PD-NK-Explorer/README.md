# PD-NK Explorer

**An Interactive Platform for Exploring Natural Killer Cell Transcriptomic Biomarkers in Parkinson's Disease**

This is a **GUI prototype** built with R Shiny, using **simulated placeholder data**.
It accompanies the research project *"Cross-Dataset Validation of Natural Killer Cell
Transcriptomic Biomarkers in Parkinson's Disease Using Single-Cell and Bulk RNA Sequencing."*

## 1. Install required packages

Run this once in R (or RStudio console):

```r
install.packages(c(
  "shiny", "bs4Dash", "shinyWidgets", "DT", "plotly",
  "ggplot2", "pheatmap", "dplyr", "tidyr", "pROC"
))
```

You do **not** need Seurat, clusterProfiler, fgsea, or ComplexHeatmap to run this
prototype — those are only needed later when you plug in real analysis results
(see "Replacing placeholder data" below). Install them then with:

```r
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "enrichplot", "fgsea"))
# Seurat:
install.packages("Seurat")
# Optional heatmap upgrade:
BiocManager::install("ComplexHeatmap")
```

## 2. Folder structure

```
PD-NK-Explorer/
├── app.R                    # Main UI + server assembly
├── global.R                 # Package loading, sourcing, constants
├── modules/                 # One module per page (UI + server)
│   ├── mod_dashboard.R
│   ├── mod_dataset.R
│   ├── mod_singlecell.R
│   ├── mod_nkcell.R
│   ├── mod_dge.R
│   ├── mod_pathway.R
│   ├── mod_validation.R
│   ├── mod_biomarker.R
│   ├── mod_download.R
│   └── mod_about.R
├── data/                    # Placeholder CSVs
│   ├── example_metadata.csv
│   ├── example_deg_results.csv
│   ├── example_validation_results.csv
│   └── example_biomarkers.csv
├── R/
│   └── helper_functions.R   # All data-loading / simulation functions
├── www/
│   ├── custom.css
│   └── images/
└── README.md
```

## 3. Run the application

Open `PD-NK-Explorer/` as your working directory in R/RStudio, then:

```r
shiny::runApp()
```

or open `app.R` in RStudio and click **Run App**.

## 4. Replacing placeholder data with real results

Every place that generates or loads placeholder data lives in **`R/helper_functions.R`**,
in functions named `load_*()` (read real files later) or `simulate_*()` (currently
generate fake data). To connect real results:

| Function | Currently | Replace with |
|---|---|---|
| `load_dataset_metadata()` | reads `data/example_metadata.csv` | your real GEO series metadata table |
| `simulate_singlecell_data()` | simulated UMAP + expression | `Seurat::Embeddings()` + `Seurat::GetAssayData()` from your real Seurat object |
| `load_deg_results()` | reads `data/example_deg_results.csv` | your real DESeq2/edgeR/limma output |
| `simulate_go_enrichment()`, `simulate_kegg_pathways()`, `simulate_gsea_curve()` | simulated | real `clusterProfiler::enrichGO()`/`enrichKEGG()` and `fgsea::fgsea()` output |
| `load_validation_results()` | reads `data/example_validation_results.csv` | your real cross-dataset validation table |
| `load_biomarker_list()` / `simulate_biomarker_expression()` | reads CSV / simulated | your real quantitative expression data (the Biomarker page already computes genuine `pROC::roc()` results from whatever data frame it's given) |

No UI/module code changes should be required — each module calls these functions and
renders whatever they return, as long as the column names stay the same.

## 5. Notes

- Every page shows a small blue "This uses placeholder data" note where relevant —
  remove these once real results are wired in.
- `NK_MARKER_GENES` in `global.R` is the example marker gene panel used throughout
  the app; update it once your real candidate gene list is finalized.
