
mod_download_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      bs4Card(
        title = "Download Data Tables", width = 6, status = "primary", solidHeader = TRUE,
        p("Download the real results tables as CSV files."),
        downloadButton(ns("dl_deg"), "Download DEG Results (CSV)", class = "btn-outline-primary mb-2"), br(),
        downloadButton(ns("dl_biomarkers"), "Download Candidate Biomarkers (CSV)", class = "btn-outline-primary mb-2"), br(),
        downloadButton(ns("dl_validation"), "Download Validation Results (CSV)", class = "btn-outline-primary mb-2"), br(),
        downloadButton(ns("dl_enrichment"), "Download Enrichment Results (CSV)", class = "btn-outline-primary mb-2")
      ),
      bs4Card(
        title = "Download Figures", width = 6, status = "primary", solidHeader = TRUE,
        p("Export a figure (Xiong NK vs GSE100054 log2FC of top candidates) as PNG or PDF."),
        downloadButton(ns("dl_png"), "Download Figure (PNG)", class = "btn-outline-secondary mb-2"), br(),
        downloadButton(ns("dl_pdf"), "Download Figure (PDF)", class = "btn-outline-secondary mb-2")
      )
    )
  )
}

mod_download_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    output$dl_deg <- downloadHandler(
      filename = function() "deg_results.csv",
      content = function(file) write.csv(load_deg_results(), file, row.names = FALSE)
    )
    output$dl_biomarkers <- downloadHandler(
      filename = function() "candidate_biomarkers_auc_gse100054.csv",
      content = function(file) write.csv(load_biomarker_list(), file, row.names = FALSE)
    )
    output$dl_validation <- downloadHandler(
      filename = function() "validation_results.csv",
      content = function(file) write.csv(load_validation_results(), file, row.names = FALSE)
    )
    output$dl_enrichment <- downloadHandler(
      filename = function() "gsea_pathway_results_xiong_nk.csv",
      content = function(file) write.csv(do.call(rbind, lapply(PATHWAY_COLLECTIONS, load_pathway_results)), file, row.names = FALSE)
    )
    
    example_plot <- function() {
      d <- load_validation_results()
      d <- d[!is.na(d$P_bulk), ]
      d <- head(d[order(d$P_xiong), ], 20)
      long <- rbind(data.frame(gene = d$gene, dataset = "Xiong NK (discovery)", log2FC = d$logFC_xiong),
                    data.frame(gene = d$gene, dataset = "GSE100054 (validation)", log2FC = d$logFC_bulk))
      long$gene <- factor(long$gene, levels = d$gene)
      ggplot(long, aes(x = gene, y = log2FC, fill = dataset)) +
        geom_col(position = "dodge") +
        labs(title = "Top 20 Xiong NK candidates: log2FC (PD vs Control)", x = NULL) +
        theme_pdnk() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    
    output$dl_png <- downloadHandler(
      filename = function() "candidate_log2fc_figure.png",
      content = function(file) ggsave(file, plot = example_plot(), width = 8, height = 5, dpi = 300)
    )
    output$dl_pdf <- downloadHandler(
      filename = function() "candidate_log2fc_figure.pdf",
      content = function(file) ggsave(file, plot = example_plot(), width = 8, height = 5)
    )
  })
}