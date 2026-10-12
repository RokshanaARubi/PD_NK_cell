library(shiny)
library(bs4Dash)       
library(shinyWidgets)  #
library(DT)             # interactive tables
library(plotly)         # interactive plots
library(ggplot2)        # static/base plots (wrapped by plotly where needed)
library(pheatmap)       
library(dplyr)
library(tidyr)
library(pROC)           # ROC / AUC — lightweight, safe to use now with simulated data


# ---- Global options ----
options(shiny.maxRequestSize = 200 * 1024^2) 
set.seed(42)  

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


APP_TITLE    <- "PD-NK Explorer"
APP_SUBTITLE <- "An Interactive Platform for Exploring Natural Killer Cell Transcriptomic Biomarkers in Parkinson's Disease"


NK_MARKER_GENES <- c("NKG7", "KLRD1", "GNLY", "PRF1", "GZMB", "XCL1", "XCL2")
