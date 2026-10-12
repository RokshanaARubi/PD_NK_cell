

mod_pathway_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    bs4TabCard(
      id = ns("pathway_tabs"), width = 12, status = "primary", solidHeader = TRUE,
      title = "Functional Enrichment", side = "left", type = "tabs",
      
      tabPanel(
        title = "GO Enrichment",
        pickerInput(ns("go_category"), "Gene-set collection:",
                    choices = c("GO Biological Process" = "GO_BP", "Reactome" = "Reactome",
                                "Hallmark" = "Hallmark"), selected = "GO_BP", options = list(container = "body")),
        fluidRow(
          column(6, plotlyOutput(ns("go_dot"), height = "380px")),
          column(6, plotlyOutput(ns("go_bar"), height = "380px"))
        ),
        DTOutput(ns("go_table"))
      ),
      
      tabPanel(
        title = "KEGG Pathways",
        fluidRow(
          column(6, plotlyOutput(ns("kegg_dot"), height = "380px")),
          column(6, plotlyOutput(ns("kegg_bar"), height = "380px"))
        ),
        DTOutput(ns("kegg_table"))
      ),
      
      tabPanel(
        title = "GSEA",
        pickerInput(ns("gene_set"), "Gene set (top Hallmark / KEGG by P):",
                    choices = c("Loading..." = ""), options = list(container = "body")),
        fluidRow(
          column(8, plotlyOutput(ns("gsea_curve"), height = "380px")),
          column(4,
                 bs4Card(title = "GSEA Statistics", width = 12, status = "info",
                         uiOutput(ns("gsea_stats"))))
        )
      )
    )
  )
}

mod_pathway_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    go_data   <- reactive({ req(input$go_category); load_pathway_results(input$go_category) })
    kegg_data <- reactive({ load_pathway_results("KEGG") })
    
    observe({
      ch <- load_gsea_curve_pathways()
      if (length(ch) > 0) updatePickerInput(session, "gene_set", choices = ch, selected = ch[1])
    })
    
    dot_plot <- function(d, ttl) {
      d <- head(d, 20)
      ggplot(d, aes(x = NES, y = reorder(term, NES), size = size, color = -log10(padj))) +
        geom_point() +
        scale_color_gradient(low = "grey60", high = "#c0392b") +
        labs(x = "Normalized enrichment score (NES)", y = NULL, title = ttl, color = "-log10(padj)") +
        theme_pdnk()
    }
    bar_plot <- function(d, fillc, ttl) {
      d <- head(d, 10)
      ggplot(d, aes(x = reorder(term, -log10(padj)), y = -log10(padj), fill = direction)) +
        geom_col() + coord_flip() +
        scale_fill_manual(values = c("Up in PD" = "#c0392b", "Down in PD" = "#2980b9")) +
        labs(x = NULL, y = "-log10(adjusted p-value)", title = ttl) +
        theme_pdnk()
    }
    tbl <- function(d) {
      datatable(d[, c("term", "NES", "pval", "padj", "size", "direction", "NES_bulk", "P_bulk")],
                options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) %>%
        formatSignif(columns = c("NES", "pval", "padj", "NES_bulk", "P_bulk"), digits = 3)
    }
    
    output$go_dot   <- renderPlotly({ ggplotly(dot_plot(go_data(), paste("Top 20:", names(PATHWAY_COLLECTIONS)[PATHWAY_COLLECTIONS == input$go_category]))) })
    output$go_bar   <- renderPlotly({ ggplotly(bar_plot(go_data(), "#2980b9", "Top 10 gene sets by adjusted P")) })
    output$go_table <- renderDT({ tbl(go_data()) })
    
    output$kegg_dot   <- renderPlotly({ ggplotly(dot_plot(kegg_data(), "KEGG pathways (Xiong NK GSEA)")) })
    output$kegg_bar   <- renderPlotly({ ggplotly(bar_plot(kegg_data(), "#8e44ad", "Top KEGG pathways")) })
    output$kegg_table <- renderDT({ tbl(kegg_data()) })
    
    output$gsea_curve <- renderPlotly({
      req(nzchar(input$gene_set))
      d <- load_gsea_curve(input$gene_set)
      p <- ggplot(d, aes(x = rank_position, y = enrichment_score)) +
        geom_line(color = "#c0392b", linewidth = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        labs(x = "Rank position (normalized; left = up in PD)", y = "Running enrichment score",
             title = paste("GSEA enrichment curve:", d$gene_set[1])) +
        theme_pdnk()
      ggplotly(p)
    })
    
    output$gsea_stats <- renderUI({
      req(nzchar(input$gene_set))
      s <- load_gsea_stats(input$gene_set)
      r <- .ui_data()$pathway_replication
      r <- r[r$pathway == input$gene_set, ][1, ]
      cam <- .ui_data()$camera_robustness
      cp <- if (!is.null(cam)) cam$P_camera_bulk[match(input$gene_set, cam$pathway)] else NA
      row <- function(a, b) tags$tr(tags$td(strong(a)), tags$td(b))
      tagList(tags$table(class = "table table-sm",
                         row("Normalized Enrichment Score:", round(s$NES, 2)),
                         row("P-value (fgsea):", signif(s$pval, 3)),
                         row("Adjusted p-value:", signif(s$padj, 3)),
                         row("Genes in set:", s$size),
                         row("NES in GSE100054:", if (is.na(r$NES_bulk)) "n/a" else round(r$NES_bulk, 2)),
                         row("camera P in GSE100054:", if (is.na(cp)) "n/a" else signif(cp, 3))
      ))
    })
  })
}