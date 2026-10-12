

mod_dge_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      bs4Card(
        title = "Filters", width = 12, status = "primary", solidHeader = TRUE,
        fluidRow(
          column(3, pickerInput(ns("dataset"), "Dataset:",
                                choices = c(DS_XIONG, DS_BULK), options = list(container = "body"))),
          column(3, pickerInput(ns("comparison"), "Comparison:", choices = "PD vs Control", options = list(container = "body")),
                 checkboxInput(ns("use_nominal"), "Use nominal P-value (exploratory)", value = TRUE)),
          column(3, sliderInput(ns("padj_cut"), "P-value / FDR threshold:",
                                min = 0.001, max = 1, value = 0.05, step = 0.001)),
          column(3, sliderInput(ns("lfc_cut"), "|log2FC| threshold:",
                                min = 0, max = 3, value = 1, step = 0.1))
        )
      )
    ),
    fluidRow(
      bs4Card(title = "Volcano Plot", width = 6, status = "info", solidHeader = TRUE,
              plotlyOutput(ns("volcano"), height = "380px")),
      bs4Card(title = "MA Plot", width = 6, status = "info", solidHeader = TRUE,
              plotlyOutput(ns("ma_plot"), height = "380px"))
    ),
    fluidRow(
      bs4Card(title = "Top DEGs (Bar Plot)", width = 12, status = "info", solidHeader = TRUE,
              plotlyOutput(ns("top_bar"), height = "380px"))
    ),
    fluidRow(
      bs4Card(title = "DEG Table", width = 12, status = "secondary", solidHeader = TRUE,
              DTOutput(ns("deg_table")))
    )
  )
}

mod_dge_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(input$dataset, {
      ch <- if (identical(input$dataset, DS_XIONG))
        c("PD vs Control", "Early PD vs Control", "Late PD vs Control", "Late PD vs Early PD")
      else "PD vs Control"
      updatePickerInput(session, "comparison", choices = ch, selected = "PD vs Control")
    }, ignoreInit = TRUE)
    
    deg_filtered <- reactive({
      req(input$dataset, input$comparison)
      deg <- load_deg_results()
      deg <- deg[deg$dataset == input$dataset & deg$comparison == input$comparison, ]
      validate(need(nrow(deg) > 0, "No results for this dataset / comparison."))
      deg$p_used <- if (isTRUE(input$use_nominal)) deg$pvalue else deg$padj
      deg$significant <- deg$p_used < input$padj_cut & abs(deg$log2FC) > input$lfc_cut
      deg
    })
    
    p_label <- reactive(if (isTRUE(input$use_nominal)) "nominal P-value" else "adjusted p-value (FDR)")
    
    output$volcano <- renderPlotly({
      d <- deg_filtered()
      p <- ggplot(d, aes(x = log2FC, y = -log10(p_used), color = significant,
                         text = paste0("Gene: ", gene))) +
        geom_point(alpha = 0.7) +
        geom_vline(xintercept = c(-input$lfc_cut, input$lfc_cut), linetype = "dashed") +
        geom_hline(yintercept = -log10(input$padj_cut), linetype = "dashed") +
        scale_color_manual(values = c("TRUE" = "#c0392b", "FALSE" = "grey70")) +
        labs(title = "Volcano plot", y = paste0("-log10(", p_label(), ")"), color = "Significant") +
        theme_pdnk()
      ggplotly(p, tooltip = "text")
    })
    
    output$ma_plot <- renderPlotly({
      d <- deg_filtered()
      p <- ggplot(d, aes(x = avg_expr, y = log2FC, color = significant,
                         text = paste0("Gene: ", gene))) +
        geom_point(alpha = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("TRUE" = "#c0392b", "FALSE" = "grey70")) +
        labs(x = "Mean expression (logCPM / AveExpr)", y = "log2FC", title = "MA plot", color = "Significant") +
        theme_pdnk()
      ggplotly(p, tooltip = "text")
    })
    
    output$top_bar <- renderPlotly({
      d <- deg_filtered()
      top <- d[order(d$p_used), ][1:min(15, nrow(d)), ]
      p <- ggplot(top, aes(x = reorder(gene, -log10(p_used)), y = -log10(p_used), fill = direction)) +
        geom_col() + coord_flip() +
        scale_fill_manual(values = c("Upregulated" = "#c0392b", "Downregulated" = "#2980b9")) +
        labs(x = NULL, y = paste0("-log10(", p_label(), ")"), title = "Top 15 DEGs by significance") +
        theme_pdnk()
      ggplotly(p)
    })
    
    output$deg_table <- renderDT({
      d <- deg_filtered()
      d <- d[order(d$p_used), ]
      datatable(d[, c("gene", "log2FC", "pvalue", "padj", "direction", "dataset", "validation_status")],
                filter = "top", options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) %>%
        formatSignif(columns = c("log2FC", "pvalue", "padj"), digits = 3)
    })
  })
}