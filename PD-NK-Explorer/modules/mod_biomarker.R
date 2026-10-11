# =============================================================
# modules/mod_biomarker.R
# Page 8: Biomarker Analysis (ROC / AUC via pROC)
# REAL DATA: GSE100054 bulk PBMC (9 Control vs 10 PD). The expected direction
# of each gene is fixed from the Xiong NK result (not from the GSE100054
# labels), so AUC < 0.5 means "wrong direction".
# =============================================================

mod_biomarker_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      bs4Card(
        title = "Select Candidate Biomarker", width = 3, status = "primary", solidHeader = TRUE,
        pickerInput(ns("gene"), "Gene:", choices = NK_MARKER_GENES, selected = "NKG7",
                    options = list(container = "body", `live-search` = TRUE)),
        uiOutput(ns("gene_info"))
      ),
      bs4Card(title = "Expression in GSE100054: PD vs Control", width = 4, status = "info",
              solidHeader = TRUE, plotlyOutput(ns("boxplot"), height = "340px")),
      bs4Card(title = "ROC Curve", width = 5, status = "info", solidHeader = TRUE,
              plotOutput(ns("roc_plot"), height = "340px"))
    ),
    fluidRow(
      bs4ValueBoxOutput(ns("auc_box"), width = 4),
      bs4ValueBoxOutput(ns("sens_box"), width = 4),
      bs4ValueBoxOutput(ns("spec_box"), width = 4)
    ),
    fluidRow(
      bs4Card(title = "AUC of Xiong-derived Signatures in GSE100054 (95% CI)", width = 12,
              status = "secondary", solidHeader = TRUE, plotlyOutput(ns("auc_compare"), height = "360px"))
    )
  )
}

mod_biomarker_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    bm <- reactive({ load_biomarker_list() })
    
    observe({
      b <- bm()
      updatePickerInput(session, "gene", choices = b$gene,
                        selected = if ("NKG7" %in% b$gene) "NKG7" else b$gene[1])
    })
    
    biomarker_expr <- reactive({
      req(input$gene)
      load_biomarker_expression(input$gene)
    })
    
    # direction fixed a priori from the Xiong NK logFC
    gene_dir <- reactive({ xiong_direction(input$gene) })
    
    roc_obj <- reactive({
      d <- biomarker_expr()
      pROC::roc(response = d$group, predictor = gene_dir() * d$expression,
                levels = c("Control", "PD"), direction = "<", quiet = TRUE)
    })
    
    output$gene_info <- renderUI({
      b <- bm(); row <- b[b$gene == input$gene, ]
      req(nrow(row) == 1)
      tags$small(paste0("Expected direction: ", row$xiong_direction,
                        ifelse(is.na(row$P_xiong), "", paste0(" | Xiong NK P = ", signif(row$P_xiong, 3)))))
    })
    
    output$boxplot <- renderPlotly({
      d <- biomarker_expr()
      p <- ggplot(d, aes(x = group, y = expression, fill = group)) +
        geom_boxplot(outlier.alpha = 0.4) + geom_jitter(width = 0.12, size = 1.4, alpha = 0.7) +
        scale_fill_manual(values = c("PD" = "#d9534f", "Control" = "#5bc0de")) +
        labs(x = NULL, y = "log2 expression (GSE100054)", title = paste(input$gene, "expression")) +
        theme_pdnk() + theme(legend.position = "none")
      ggplotly(p)
    })
    
    output$roc_plot <- renderPlot({
      r <- roc_obj()
      plot(r, main = paste("ROC curve —", input$gene), col = "#c0392b", lwd = 2,
           print.auc = TRUE, legacy.axes = TRUE)
    })
    
    output$auc_box <- renderbs4ValueBox({
      r <- roc_obj()
      bs4ValueBox(value = round(as.numeric(pROC::auc(r)), 3), subtitle = "AUC (direction fixed from Xiong)",
                  icon = icon("chart-line"), color = "primary")
    })
    output$sens_box <- renderbs4ValueBox({
      r <- roc_obj()
      co <- pROC::coords(r, "best", ret = "sensitivity", best.method = "youden")
      bs4ValueBox(value = round(as.numeric(co[[1]][1]), 3), subtitle = "Sensitivity (best cutoff, in-sample)",
                  icon = icon("check"), color = "success")
    })
    output$spec_box <- renderbs4ValueBox({
      r <- roc_obj()
      co <- pROC::coords(r, "best", ret = "specificity", best.method = "youden")
      bs4ValueBox(value = round(as.numeric(co[[1]][1]), 3), subtitle = "Specificity (best cutoff, in-sample)",
                  icon = icon("shield"), color = "warning")
    })
    
    output$auc_compare <- renderPlotly({
      s <- load_signature_roc()
      s$signature <- gsub("_", " ", s$signature)
      p <- ggplot(s, aes(x = reorder(signature, AUC), y = AUC,
                         text = paste0("Permutation P: ", perm_P_one_sided))) +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
        geom_pointrange(aes(ymin = CI_low, ymax = CI_high), color = "#2980b9") +
        coord_flip(ylim = c(0, 1)) +
        labs(x = NULL, y = "AUC (0.5 = chance)", title = "Signatures defined from Xiong NK data, tested in GSE100054") +
        theme_pdnk()
      ggplotly(p, tooltip = c("y", "text"))
    })
  })
}