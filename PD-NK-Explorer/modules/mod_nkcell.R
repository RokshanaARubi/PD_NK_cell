
mod_nkcell_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      bs4Card(
        title = "NK Cell Highlight — UMAP", width = 6, status = "primary",
        solidHeader = TRUE, plotlyOutput(ns("umap_nk"), height = "380px")
      ),
      bs4Card(
        title = "Marker Gene Dot Plot", width = 6, status = "primary",
        solidHeader = TRUE, plotOutput(ns("dot_plot"), height = "380px")
      )
    ),
    fluidRow(
      bs4Card(
        title = "Controls", width = 3, status = "info", solidHeader = TRUE,
        pickerInput(ns("marker_gene"), "NK marker gene:",
                    choices = NK_MARKER_GENES, selected = "NKG7", options = list(container = "body"))
      ),
      bs4Card(
        title = "Marker Expression (Violin)", width = 4, status = "info",
        solidHeader = TRUE, plotlyOutput(ns("marker_violin"), height = "340px")
      ),
      bs4Card(
        title = "PD vs Control in NK cells", width = 5, status = "info",
        solidHeader = TRUE, plotlyOutput(ns("pd_vs_control"), height = "340px")
      )
    ),
    fluidRow(
      bs4Card(
        title = "NK-Cell Exploratory Candidate Genes (Tier 2: P < 0.01 and |log2FC| >= 1)",
        width = 12, status = "secondary",
        solidHeader = TRUE, DTOutput(ns("candidate_table"))
      )
    )
  )
}

mod_nkcell_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    sc_data <- reactive({ load_singlecell_data() })
    
    output$umap_nk <- renderPlotly({
      d <- sc_data()$cells
      d$is_nk <- ifelse(d$cell_type == "NK cells", "NK cells", "Other")
      p <- ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = is_nk, alpha = is_nk)) +
        geom_point(size = 0.9) +
        scale_color_manual(values = c("NK cells" = "#c0392b", "Other" = "grey80")) +
        scale_alpha_manual(values = c("NK cells" = 0.9, "Other" = 0.35)) +
        labs(title = sprintf("NK cells highlighted on UMAP (%s NK cells in full data)",
                             format(sc_data()$n_nk_total, big.mark = ",")), color = NULL) +
        guides(alpha = "none") +
        theme_pdnk()
      ggplotly(p)
    })
    
    output$dot_plot <- renderPlot({
      sc <- sc_data(); cells <- sc$cells
      genes <- sc_genes_available(NK_MARKER_GENES)
      df <- do.call(rbind, lapply(genes, function(g) {
        vals <- sc_gene(sc, g)
        do.call(rbind, lapply(split(seq_along(vals), cells$cell_type), function(ix) {
          data.frame(gene = g, cell_type = cells$cell_type[ix[1]],
                     mean_expr = mean(vals[ix]), pct_expressing = mean(vals[ix] > 0) * 100)
        }))
      }))
      ggplot(df, aes(x = cell_type, y = gene, size = pct_expressing, color = mean_expr)) +
        geom_point() +
        scale_color_gradient(low = "grey80", high = "#8e0000") +
        scale_size(range = c(1, 10)) +
        labs(x = NULL, y = NULL, title = "NK marker genes across cell types",
             size = "% expressing", color = "Mean\nexpression") +
        theme_pdnk() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    })
    
    output$marker_violin <- renderPlotly({
      sc <- sc_data()
      df <- cbind(sc$cells, expression = sc_gene(sc, input$marker_gene))
      df <- df[df$cell_type == "NK cells", ]
      p <- ggplot(df, aes(x = group, y = expression, fill = group)) +
        geom_violin() + geom_jitter(width = 0.1, size = 0.5, alpha = 0.4) +
        scale_fill_manual(values = c("PD" = "#d9534f", "Control" = "#5bc0de")) +
        labs(x = NULL, y = "Expression (log-normalised)",
             title = paste(input$marker_gene, "in NK cells")) +
        theme_pdnk() + theme(legend.position = "none")
      ggplotly(p)
    })
    
    output$pd_vs_control <- renderPlotly({
      sc <- sc_data(); cells <- sc$cells
      nk <- cells$cell_type == "NK cells"
      genes <- sc_genes_available(NK_MARKER_GENES)
      df <- do.call(rbind, lapply(genes, function(g) {
        vals <- sc_gene(sc, g)[nk]
        data.frame(gene = g, group = cells$group[nk], expression = vals)
      }))
      summary_df <- df %>% group_by(gene, group) %>% summarise(mean_expr = mean(expression), .groups = "drop")
      p <- ggplot(summary_df, aes(x = gene, y = mean_expr, fill = group)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = c("PD" = "#d9534f", "Control" = "#5bc0de")) +
        labs(x = NULL, y = "Mean expression (NK cells)", title = "NK marker means: PD vs Control") +
        theme_pdnk() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
      ggplotly(p)
    })
    
    output$candidate_table <- renderDT({
      df <- load_nk_candidates()
      datatable(df, filter = "top", options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) %>%
        formatSignif(columns = c("logFC_xiong", "P_xiong", "FDR_xiong", "logFC_bulk", "P_bulk"), digits = 3)
    })
  })
}