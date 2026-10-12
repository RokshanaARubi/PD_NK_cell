

STATUS_COLORS <- c("REPRODUCIBLE NK BIOMARKER" = "#27ae60",
                   "Suggestive (nominal in both datasets, not FDR-significant)" = "#f39c12",
                   "Not replicated" = "grey60",
                   "Contradicted by bulk data" = "#c0392b",
                   "Not measurable in bulk" = "grey85")

mod_validation_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      bs4Card(title = "Discovery (Xiong NK) vs Validation (GSE100054) Status by Gene", width = 12,
              status = "primary", solidHeader = TRUE, DTOutput(ns("validation_table")))
    ),
    fluidRow(
      bs4Card(title = "log2FC Comparison Across Datasets", width = 6, status = "info",
              solidHeader = TRUE, plotlyOutput(ns("lfc_compare"), height = "380px")),
      bs4Card(title = "Reproducibility Score", width = 6, status = "info",
              solidHeader = TRUE, plotlyOutput(ns("repro_score"), height = "380px"))
    ),
    fluidRow(
      bs4Card(title = "Overlap of Nominally Significant Genes (Venn-style)", width = 12, status = "secondary",
              solidHeader = TRUE, plotOutput(ns("venn_plot"), height = "380px"))
    )
  )
}

mod_validation_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    top_candidates <- function(n = 20) {
      d <- load_validation_results()
      d <- d[!is.na(d$P_bulk), ]
      head(d[order(d$P_xiong), ], n)
    }
    
    output$validation_table <- renderDT({
      d <- load_validation_results()
      d <- d[order(d$P_xiong), ]
      datatable(d, filter = "top", options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) %>%
        formatSignif(columns = c("logFC_xiong", "P_xiong", "FDR_xiong", "logFC_bulk", "P_bulk", "FDR_candset"), digits = 3)
    })
    
    output$lfc_compare <- renderPlotly({
      d <- top_candidates(20)
      long <- rbind(data.frame(gene = d$gene, dataset = "Xiong NK (discovery)", log2FC = d$logFC_xiong),
                    data.frame(gene = d$gene, dataset = "GSE100054 (validation)", log2FC = d$logFC_bulk))
      long$gene <- factor(long$gene, levels = d$gene)
      p <- ggplot(long, aes(x = gene, y = log2FC, fill = dataset)) +
        geom_col(position = "dodge") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(x = NULL, y = "log2FC (PD vs Control)", title = "Top 20 Xiong NK candidates: log2FC in both datasets") +
        theme_pdnk() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggplotly(p)
    })
    
    output$repro_score <- renderPlotly({
      d <- top_candidates(25)
      d$score <- ifelse(d$final_status == "REPRODUCIBLE NK BIOMARKER", 1,
                        ifelse(grepl("^Suggestive", d$final_status), 0.5, 0))
      p <- ggplot(d, aes(x = reorder(gene, score), y = score, fill = final_status)) +
        geom_col() + coord_flip() +
        scale_fill_manual(values = STATUS_COLORS) +
        labs(x = NULL, y = "Reproducibility score",
             title = "Reproducibility of top 25 candidates", fill = NULL) +
        theme_pdnk()
      ggplotly(p)
    })
    
    output$venn_plot <- renderPlot({
      ui <- .ui_data()
      x <- ui$xiong_de; b <- ui$bulk_de
      x$g <- toupper(x$gene); b$g <- toupper(b$gene)
      shared <- intersect(x$g, b$g)
      x <- x[match(shared, x$g), ]; b <- b[match(shared, b$g), ]
      A <- x$PValue < 0.01; B <- b$P.Value < 0.05
      both <- A & B; same <- both & sign(x$logFC) == sign(b$logFC)
      expected <- sum(A) * mean(B)
      # circles drawn as polygons on a fixed 1:1 axis, so they stay round at any card width
      circle <- function(cx, cy, r) { t <- seq(0, 2 * pi, length.out = 200); cbind(cx + r * cos(t), cy + r * sin(t)) }
      par(mar = c(1, 1, 3, 1))
      plot.new(); plot.window(xlim = c(0, 16), ylim = c(0, 10), asp = 1)
      polygon(circle(6, 5.2, 3.6), col = adjustcolor("#c0392b", 0.25), border = "#c0392b", lwd = 2)
      polygon(circle(10, 5.2, 3.6), col = adjustcolor("#2980b9", 0.25), border = "#2980b9", lwd = 2)
      text(4.4, 5.2, sum(A & !B), cex = 1.6); text(11.6, 5.2, sum(B & !A), cex = 1.6)
      text(8, 5.2, sprintf("%d\n(%d same dir.)", sum(both), sum(same)), cex = 1.2)
      text(4.6, 9.4, "Xiong NK\nP < 0.01", font = 2); text(11.4, 9.4, "GSE100054\nP < 0.05", font = 2)
      text(8, 0.2, sprintf("Overlap expected by chance: %.1f  |  observed: %d  |  genes tested in both: %s",
                           expected, sum(both), format(length(shared), big.mark = ",")), cex = 0.9)
      title("Overlap of nominally significant genes")
    })
  })
}