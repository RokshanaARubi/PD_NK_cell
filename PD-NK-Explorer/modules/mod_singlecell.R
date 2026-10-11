# =============================================================
# modules/mod_singlecell.R
# Page 3: Single-Cell Explorer
# REAL DATA: Xiong et al. GSE223138 PBMC scRNA-seq (after QC), exported by
# scripts/11_export_ui_data.R. UMAP/expression are shown for a random
# stratified subsample of cells (speed); the distribution plot uses ALL cells.
# =============================================================

mod_singlecell_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      bs4Card(
        title = "Controls", width = 3, status = "primary", solidHeader = TRUE,
        pickerInput(ns("meta_var"), "Color UMAP by:",
                    choices = c("cell_type", "group"), selected = "cell_type",
                    options = list(container = "body")),
        pickerInput(ns("cell_type_filter"), "Highlight cell type:",
                    choices = c("All", "NK cells", "T cells", "B cells",
                                "Monocytes", "Other / unassigned"),
                    selected = "All", options = list(container = "body")),
        pickerInput(ns("gene"), "Gene for Feature/Violin plot:",
                    choices = c(NK_MARKER_GENES, "CD3D", "CD3E", "MS4A1", "CD14"),
                    selected = "NKG7", options = list(container = "body"))
      ),
      bs4Card(
        title = "UMAP — Cell Clusters", width = 9, status = "primary", solidHeader = TRUE,
        plotlyOutput(ns("umap_plot"), height = "420px")
      )
    ),
    fluidRow(
      bs4Card(title = "Feature Plot (Selected Gene)", width = 6, status = "info",
              solidHeader = TRUE, plotlyOutput(ns("feature_plot"), height = "350px")),
      bs4Card(title = "Violin Plot (Selected Gene)", width = 6, status = "info",
              solidHeader = TRUE, plotlyOutput(ns("violin_plot"), height = "350px"))
    ),
    fluidRow(
      bs4Card(title = "Cell Type Distribution", width = 12, status = "secondary",
              solidHeader = TRUE, plotlyOutput(ns("celltype_dist"), height = "320px"))
    )
  )
}

mod_singlecell_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    sc_data <- reactive({ load_singlecell_data() })
    
    output$umap_plot <- renderPlotly({
      sc <- sc_data(); d <- sc$cells
      color_var <- input$meta_var
      d$highlight <- if (input$cell_type_filter != "All") {
        ifelse(d$cell_type == input$cell_type_filter, d$cell_type, "Other")
      } else {
        d[[color_var]]
      }
      p <- ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = highlight,
                         text = paste0("Cell type: ", cell_type, "<br>Group: ", group,
                                       "<br>Sample: ", sample_id))) +
        geom_point(size = 0.8, alpha = 0.7) +
        labs(color = if (input$cell_type_filter != "All") "Highlight" else color_var,
             title = sprintf("UMAP of PBMC cells (%s of %s cells shown)",
                             format(sc$n_cells_shown, big.mark = ","),
                             format(sc$n_cells_total, big.mark = ","))) +
        theme_pdnk()
      ggplotly(p, tooltip = "text")
    })
    
    output$feature_plot <- renderPlotly({
      sc <- sc_data()
      df <- cbind(sc$cells, expression = sc_gene(sc, input$gene))
      p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = expression)) +
        geom_point(size = 0.9) +
        scale_color_gradient(low = "grey85", high = "#b30000") +
        labs(title = paste("Expression of", input$gene), color = "log-norm") +
        theme_pdnk()
      ggplotly(p)
    })
    
    output$violin_plot <- renderPlotly({
      sc <- sc_data()
      df <- cbind(sc$cells, expression = sc_gene(sc, input$gene))
      p <- ggplot(df, aes(x = cell_type, y = expression, fill = cell_type)) +
        geom_violin(scale = "width") +
        labs(x = NULL, y = "Expression (log-normalised)",
             title = paste(input$gene, "expression by cell type")) +
        theme_pdnk() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
      ggplotly(p)
    })
    
    output$celltype_dist <- renderPlotly({
      df <- sc_data()$counts_full %>% group_by(cell_type, group) %>%
        summarise(n = sum(n_cells), .groups = "drop")
      p <- ggplot(df, aes(x = cell_type, y = n, fill = group)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = c("PD" = "#d9534f", "Control" = "#5bc0de")) +
        labs(x = NULL, y = "Number of cells (all cells)",
             title = "Cell type distribution: PD (4 donors) vs Control (2 donors)") +
        theme_pdnk() +
        theme(axis.text.x = element_text(angle = 20, hjust = 1))
      ggplotly(p)
    })
  })
}