# =============================================================
# modules/mod_about.R
# Page 10: About
# =============================================================

mod_about_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      bs4Card(
        title = "Project Description", width = 12, status = "primary", solidHeader = TRUE,
        p(strong("Title: "), "Cross-Dataset Validation of Natural Killer Cell Transcriptomic ",
          "Biomarkers in Parkinson's Disease Using Single-Cell and Bulk RNA Sequencing"),
        p("PD-NK Explorer is an interactive companion tool that accompanies this research, ",
          "allowing readers and reviewers to explore the datasets, single-cell populations, ",
          "differential expression, functional enrichment, and cross-dataset validation ",
          "results underlying the study.")
      )
    ),
    fluidRow(
      bs4Card(
        title = "Research Objective & Question", width = 6, status = "info", solidHeader = TRUE,
        p(strong("Objective: "), "Identify NK-cell-specific transcriptomic biomarkers of ",
          "Parkinson's disease from single-cell RNA-seq and evaluate their reproducibility ",
          "across independent bulk RNA-seq datasets."),
        p(strong("Question: "), "Can NK-cell transcriptomic signatures serve as reproducible ",
          "peripheral biomarkers for PD?")
      ),
      bs4Card(
        title = "Proposed Methodology", width = 6, status = "info", solidHeader = TRUE,
        tags$ol(
          tags$li("Quality control and preprocessing of public scRNA-seq PD datasets"),
          tags$li("Cell clustering and annotation"),
          tags$li("NK cell identification using canonical markers"),
          tags$li("Differential gene expression analysis (PD vs Control)"),
          tags$li("Candidate biomarker gene selection"),
          tags$li("Validation across independent bulk RNA-seq cohorts"),
          tags$li("Functional enrichment (GO / KEGG / GSEA)"),
          tags$li("ROC analysis of candidate biomarkers")
        )
      )
    ),
    fluidRow(
      bs4Card(
        title = "Technologies & R Packages Used", width = 6, status = "secondary", solidHeader = TRUE,
        tags$ul(
          tags$li("R / Shiny, bs4Dash, shinyWidgets"),
          tags$li("DT, plotly, ggplot2, pheatmap"),
          tags$li("pROC (ROC / AUC analysis)"),
          tags$li("Seurat (single-cell QC, clustering, annotation)"),
          tags$li("edgeR / limma (pseudobulk and bulk differential expression), fgsea + msigdbr (GSEA)")
        )
      ),
      bs4Card(
        title = "Limitations", width = 6, status = "secondary", solidHeader = TRUE,
        tags$ul(
          tags$li("Discovery has only 2 Control and 4 PD donors (Xiong GSE223138): low power, results are exploratory"),
          tags$li("NK-specific signatures did not replicate in the independent bulk dataset GSE100054 (microarray, mixed cell types)"),
          tags$li("No second scRNA-seq replication yet (Moquin-Beaudry raw data not processed) and no PPMI"),
          tags$li("Real analysis results depend on availability and quality of public datasets"),
          tags$li("Cross-dataset validation is limited by batch effects and technology differences"),
          tags$li("Candidate biomarkers require independent clinical validation before use")
        )
      )
    ),
    fluidRow(
      bs4Card(
        title = "Repository", width = 12, status = "primary", solidHeader = TRUE,
        p(strong("GitHub Repository: "), a(href = "#", "https://github.com/RokshanaARubi/PD_NK_cell.git "))
      )
    )
  )
}

mod_about_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # No server logic needed for this static page
  })
}