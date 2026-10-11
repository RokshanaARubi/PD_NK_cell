# =============================================================
# modules/mod_dataset.R
# Page 2: Dataset Explorer
# =============================================================

mod_dataset_ui <- function(id) {
  ns <- NS(id)
  
  meta <- load_dataset_metadata()
  
  tagList(
    fluidRow(
      bs4Card(
        title = "Select a Dataset", width = 12, status = "primary", solidHeader = TRUE,
        pickerInput(ns("dataset"), "Dataset:",
                    choices = setNames(meta$dataset_id, meta$name),
                    options = list(style = "btn-outline-primary", container = "body"))
      )
    ),
    fluidRow(
      bs4Card(title = "Dataset Summary", width = 4, status = "info", solidHeader = TRUE,
              uiOutput(ns("summary_ui"))),
      bs4Card(title = "Sample Distribution", width = 8, status = "info", solidHeader = TRUE,
              plotlyOutput(ns("dist_plot"), height = "320px"))
    ),
    fluidRow(
      bs4Card(title = "Sample Metadata", width = 12, status = "secondary", solidHeader = TRUE,
              DTOutput(ns("meta_table")))
    )
  )
}

mod_dataset_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    meta <- load_dataset_metadata()
    
    selected_meta <- reactive({
      req(input$dataset)
      meta[meta$dataset_id == input$dataset, ]
    })
    
    sample_meta <- reactive({
      req(input$dataset)
      row <- selected_meta()
      simulate_sample_metadata(input$dataset, n_pd = row$n_pd, n_control = row$n_control)
    })
    
    output$summary_ui <- renderUI({
      row <- selected_meta()
      tagList(
        tags$table(class = "table table-sm",
                   tags$tr(tags$td(strong("Tissue:")), tags$td(row$tissue)),
                   tags$tr(tags$td(strong("Technology:")), tags$td(row$technology)),
                   tags$tr(tags$td(strong("PD samples:")), tags$td(row$n_pd)),
                   tags$tr(tags$td(strong("Control samples:")), tags$td(row$n_control))
        ),
        hr(),
        p(row$description)
      )
    })
    
    output$dist_plot <- renderPlotly({
      row <- selected_meta()
      df <- data.frame(group = c("PD", "Control"), n = c(row$n_pd, row$n_control))
      p <- ggplot(df, aes(x = group, y = n, fill = group)) +
        geom_col(width = 0.5) +
        scale_fill_manual(values = c("PD" = "#d9534f", "Control" = "#5bc0de")) +
        labs(x = NULL, y = "Number of samples", title = "PD vs Control sample counts") +
        theme_pdnk()
      ggplotly(p)
    })
    
    output$meta_table <- renderDT({
      datatable(sample_meta(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    })
  })
}