# =============================================================
# global.R
# PD-NK Explorer
# Loads packages, sources helper functions & modules, sets
# global options. Runs once when the app starts.
# =============================================================

# ---- Core packages ----
library(shiny)
library(bs4Dash)       # modern admin-style dashboard UI
library(shinyWidgets)  # extra UI widgets (pickers, sliders, etc.)
library(DT)             # interactive tables
library(plotly)         # interactive plots
library(ggplot2)        # static/base plots (wrapped by plotly where needed)
library(pheatmap)       # heatmaps (swap for ComplexHeatmap later if desired)
library(dplyr)
library(tidyr)
library(pROC)           # ROC / AUC — lightweight, safe to use now with simulated data

# ---- FUTURE / OPTIONAL packages (NOT required for this prototype) ----
# Uncomment and install these ONLY when you plug in real analysis results:
# library(Seurat)            # for real scRNA-seq objects (Single-Cell / NK Explorer pages)
# library(clusterProfiler)   # for real GO/KEGG/GSEA results (Pathway Analysis page)
# library(enrichplot)        # for real enrichment plots
# library(fgsea)             # for real GSEA
# library(ComplexHeatmap)    # optional upgrade over pheatmap

# ---- Global options ----
options(shiny.maxRequestSize = 200 * 1024^2)  # allow larger file uploads later (e.g. RDS/Seurat objects)
set.seed(42)  # reproducible simulated placeholder data across app restarts

# ---- Source helper functions ----
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\R\\helper_functions.R")

# ---- Source all modules ----
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\modules\\mod_dashboard.R")
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\modules\\mod_dataset.R")
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\modules\\mod_singlecell.R")
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\modules\\mod_nkcell.R")
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\modules\\mod_dge.R")
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\modules\\mod_pathway.R")
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\modules\\mod_validation.R")
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\modules\\mod_biomarker.R")
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\modules\\mod_download.R")
source("D:\\academics\\L4S1\\r\\PD-NK-Explorer\\PD-NK-Explorer\\modules\\mod_about.R")

# ---- App-wide constants ----
APP_TITLE    <- "PD-NK Explorer"
APP_SUBTITLE <- "An Interactive Platform for Exploring Natural Killer Cell Transcriptomic Biomarkers in Parkinson's Disease"

# Example NK-cell marker genes used across several pages.
# TODO: Replace/extend with the final marker gene panel from your real analysis.
NK_MARKER_GENES <- c("NKG7", "KLRD1", "GNLY", "PRF1", "GZMB", "XCL1", "XCL2")
