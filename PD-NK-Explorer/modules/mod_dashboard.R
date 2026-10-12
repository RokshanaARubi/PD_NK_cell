
mod_dashboard_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        width = 12,
        div(class = "pdnk-hero",
            h2(APP_TITLE),
            h4(APP_SUBTITLE),
        )
      )
    ),
    
    # ---- Summary stat cards (REAL VALUES from the analysis pipeline) ----
    fluidRow(
      bs4ValueBoxOutput(ns("box_datasets"), width = 2),
      bs4ValueBoxOutput(ns("box_pd"), width = 2),
      bs4ValueBoxOutput(ns("box_control"), width = 2),
      bs4ValueBoxOutput(ns("box_cells"), width = 2),
      bs4ValueBoxOutput(ns("box_nk"), width = 2),
      bs4ValueBoxOutput(ns("box_biomarkers"), width = 2)
    ),
    
    fluidRow(
      bs4Card(
        title = "About Parkinson's Disease", width = 6, status = "primary",
        solidHeader = TRUE, collapsible = TRUE,
        p("Parkinson's disease (PD) is a progressive neurodegenerative disorder ",
          "characterized primarily by the loss of dopaminergic neurons in the ",
          "substantia nigra, leading to motor symptoms such as tremor, rigidity, ",
          "and bradykinesia. Increasing evidence points to a role for the immune ",
          "system, including innate immune cells, in PD pathogenesis and progression.")
      ),
      bs4Card(
        title = "About Natural Killer (NK) Cells", width = 6, status = "primary",
        solidHeader = TRUE, collapsible = TRUE,
        p("Natural Killer (NK) cells are innate lymphocytes capable of recognizing ",
          "and eliminating stressed, infected, or abnormal cells without prior ",
          "sensitization. Emerging evidence suggests altered NK cell frequency and ",
          "transcriptomic profiles in Parkinson's disease, making them a promising ",
          "source of peripheral biomarkers.")
      )
    ),
    
    fluidRow(
      bs4Card(
        title = "Research Question", width = 12, status = "info", solidHeader = TRUE,
        p(strong("Can NK-cell-specific transcriptomic signatures, identified from single-cell ",
                 "RNA sequencing and validated across independent bulk RNA-seq datasets, serve as ",
                 "reproducible peripheral biomarkers for Parkinson's disease?"))
      )
    ),
    
    fluidRow(
      bs4Card(
        title = "Project Workflow", width = 12, status = "secondary",
        solidHeader = TRUE, collapsible = TRUE,
        div(class = "workflow-diagram",
            uiOutput(ns("workflow_ui"))
        )
      )
    )
  )
}

mod_dashboard_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # REAL VALUES from the downloaded metadata and the analysis results.
    # If the results file has not been built yet, boxes show "Pending".
    stats <- load_dashboard_stats()
    fmt <- function(x) if (is.na(x)) "Pending" else format(x, big.mark = ",")
    
    output$box_datasets <- renderbs4ValueBox({
      bs4ValueBox(value = stats$n_datasets, subtitle = "Datasets", color = "primary")
    })
    output$box_pd <- renderbs4ValueBox({
      bs4ValueBox(value = stats$n_pd, subtitle = "PD Samples", color = "danger")
    })
    output$box_control <- renderbs4ValueBox({
      bs4ValueBox(value = stats$n_control, subtitle = "Control Samples", color = "success")
    })
    output$box_cells <- renderbs4ValueBox({
      bs4ValueBox(value = fmt(stats$cells_profiled), subtitle = "Cells Profiled (after QC)", color = "info")
    })
    output$box_nk <- renderbs4ValueBox({
      bs4ValueBox(value = fmt(stats$nk_cells), subtitle = "NK Cells", color = "warning")
    })
    output$box_biomarkers <- renderbs4ValueBox({
      bs4ValueBox(value = fmt(stats$n_biomarkers), subtitle = if (is.na(stats$n_reproducible)) "Candidate Biomarkers" else paste0("Exploratory Candidates (", stats$n_reproducible, " reproducible)"), color = "secondary")
    })
    
    output$workflow_ui <- renderUI({
      steps <- c(
        "Public scRNA-seq Dataset", "Quality Control", "Cell Clustering & Annotation",
        "NK Cell Identification", "Differential Gene Expression",
        "Candidate Biomarker Selection", "Independent Bulk Dataset Validation",
        "Functional Analysis", "ROC Evaluation", "Reproducible Biomarkers"
      )
      tagList(
        lapply(seq_along(steps), function(i) {
          tagList(
            div(class = "workflow-step", steps[i]),
            if (i < length(steps)) div(class = "workflow-arrow", icon("arrow-down"))
          )
        })
      )
    })
  })
}