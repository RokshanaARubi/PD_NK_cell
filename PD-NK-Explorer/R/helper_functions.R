
.read_real_metadata_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) NULL)
}

load_dataset_metadata <- function() {
  
  xiong_meta <- .read_real_metadata_csv("data/xiong/GSE223138_sample_metadata.csv")
  gse_meta   <- .read_real_metadata_csv("data/gse100054/GSE100054_sample_metadata.csv")
  moquin_sra <- .read_real_metadata_csv("data/moquin_beaudry/moquin_beaudry_SRA_run_metadata.csv")
  
  xiong_n_pd <- xiong_n_control <- NA_integer_
  if (!is.null(xiong_meta) && "disease state:ch1" %in% names(xiong_meta)) {
    ds <- xiong_meta[["disease state:ch1"]]
    xiong_n_pd      <- sum(grepl("Parkinson", ds, ignore.case = TRUE))
    xiong_n_control <- sum(grepl("Normal|Control|Healthy", ds, ignore.case = TRUE))
  }
  
  # --- GSE100054: real diagnosis counts from GEO sample metadata ---
  gse_n_pd <- gse_n_control <- NA_integer_
  if (!is.null(gse_meta) && "diagnosis:ch1" %in% names(gse_meta)) {
    dx <- gse_meta[["diagnosis:ch1"]]
    gse_n_pd      <- sum(grepl("Parkinson", dx, ignore.case = TRUE))
    gse_n_control <- sum(grepl("healthy|control", dx, ignore.case = TRUE))
  }
  

  moquin_n_runs <- if (!is.null(moquin_sra)) nrow(moquin_sra) else NA_integer_
  
  data.frame(
    dataset_id = c("xiong", "moquin_beaudry", "gse100054"),
    name = c("Xiong et al. 2024 (GSE223138)",
             "Moquin-Beaudry et al. 2024/25 (PRJNA1145007)",
             "GSE100054"),
    tissue = c("PBMC", "PBMC", "PBMC"),
    technology = c("scRNA-seq", "scRNA-seq", "Bulk transcriptomics (microarray)"),
    n_pd = c(xiong_n_pd, NA_integer_, gse_n_pd),
    n_control = c(xiong_n_control, NA_integer_, gse_n_control),
    description = c(
      sprintf(paste0("PRIMARY DISCOVERY dataset. Real 10x PBMC scRNA-seq for %s donors ",
                     "(%s PD [Early/Late], %s Control/Normal): QC, clustering, NK-cell ",
                     "identification and donor-level pseudobulk differential expression are ",
                     "complete (see the Single-Cell, NK Cell and DGE pages)."),
              if (!is.null(xiong_meta)) nrow(xiong_meta) else "NOT VERIFIED",
              ifelse(is.na(xiong_n_pd), "NOT VERIFIED", xiong_n_pd),
              ifelse(is.na(xiong_n_control), "NOT VERIFIED", xiong_n_control)),
      sprintf(paste0("INDEPENDENT scRNA-seq REPLICATION dataset. Only SRA run-level metadata ",
                     "is available locally (%s runs); no donor/disease-group labels or processed ",
                     "expression matrix have been acquired yet, so PD/Control counts are not available."),
              ifelse(is.na(moquin_n_runs), "NOT VERIFIED", moquin_n_runs)),
      sprintf(paste0("INDEPENDENT PBMC VALIDATION dataset. READY: real expression matrix (%s ",
                     "features x %s samples) and real sample metadata are both loaded, with ",
                     "%s PD and %s Control samples confirmed by diagnosis labels."),
              "138,745", "19",
              ifelse(is.na(gse_n_pd), "NOT VERIFIED", gse_n_pd),
              ifelse(is.na(gse_n_control), "NOT VERIFIED", gse_n_control))
    ),
    stringsAsFactors = FALSE
  )
}


simulate_sample_metadata <- function(dataset_id, n_pd = 12, n_control = 10) {
  
  if (identical(dataset_id, "xiong")) {
    xiong_meta <- .read_real_metadata_csv("data/xiong/GSE223138_sample_metadata.csv")
    if (!is.null(xiong_meta)) {
      return(data.frame(
        sample_id = xiong_meta[["geo_accession"]],
        group = dplyr::case_when(
          grepl("Parkinson", xiong_meta[["disease state:ch1"]], ignore.case = TRUE) ~ "PD",
          grepl("Normal|Control|Healthy", xiong_meta[["disease state:ch1"]], ignore.case = TRUE) ~ "Control",
          TRUE ~ NA_character_
        ),
        age = NA_integer_,   
        sex = NA_character_, 
        disease_stage = xiong_meta[["disease state:ch1"]],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (identical(dataset_id, "gse100054")) {
    gse_meta <- .read_real_metadata_csv("data/gse100054/GSE100054_sample_metadata.csv")
    if (!is.null(gse_meta)) {
      return(data.frame(
        sample_id = gse_meta[["geo_accession"]],
        group = dplyr::case_when(
          grepl("Parkinson", gse_meta[["diagnosis:ch1"]], ignore.case = TRUE) ~ "PD",
          grepl("healthy|control", gse_meta[["diagnosis:ch1"]], ignore.case = TRUE) ~ "Control",
          TRUE ~ NA_character_
        ),
        age = suppressWarnings(as.integer(gsub("[^0-9]", "", gse_meta[["age:ch1"]]))),
        sex = tools::toTitleCase(gse_meta[["gender:ch1"]]),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (identical(dataset_id, "moquin_beaudry")) {
    moquin_sra <- .read_real_metadata_csv("data/moquin_beaudry/moquin_beaudry_SRA_run_metadata.csv")
    if (!is.null(moquin_sra) && "run_accession" %in% names(moquin_sra)) {
      run_id <- stringr::str_extract(moquin_sra[["run_accession"]], "SRR[0-9]+")
      return(data.frame(
        sample_id = ifelse(is.na(run_id), paste0("run_", seq_len(nrow(moquin_sra))), run_id),
        group = "Not available (no disease labels in SRA run metadata)",
        age = NA_integer_,
        sex = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Fallback only if a real file genuinely can't be found/read for a
  # recognized dataset_id (keeps the app from crashing, but makes the
  # fallback obvious rather than passing it off as real).
  n_pd <- ifelse(is.na(n_pd), 0, n_pd)
  n_control <- ifelse(is.na(n_control), 0, n_control)
  n_fallback <- max(1, n_pd + n_control)
  data.frame(
    sample_id = paste0(dataset_id, "_sample_", seq_len(n_fallback)),
    group = "UNKNOWN - real metadata file not found",
    age = NA_integer_,
    sex = NA_character_,
    stringsAsFactors = FALSE
  )
}




.ui_env <- new.env(parent = emptyenv())

.ui_data <- function() {
  if (is.null(.ui_env$data)) {
    f <- file.path("results", "validated_data", "ui_app_data.rds")
    if (!file.exists(f))
      stop("Results file not found: ", f, ". Run scripts/10_final_summary.R and ",
           "scripts/11_export_ui_data.R from the project root first.", call. = FALSE)
    .ui_env$data <- readRDS(f)
  }
  .ui_env$data
}

has_real_results <- function() file.exists(file.path("results", "validated_data", "ui_app_data.rds"))


load_singlecell_data <- function() .ui_data()$singlecell

# numeric expression vector for one gene, aligned to sc$cells
sc_gene <- function(sc, gene) {
  if (!gene %in% rownames(sc$expr)) return(rep(NA_real_, nrow(sc$cells)))
  as.numeric(sc$expr[gene, sc$cells$cell_id])
}

sc_cell_types <- function() {
  ct <- unique(load_singlecell_data()$cells$cell_type)
  ct[order(ct == "Other / unassigned", ct)]
}

sc_genes_available <- function(genes) intersect(genes, rownames(load_singlecell_data()$expr))

DS_XIONG <- "Xiong NK cells (scRNA-seq pseudobulk)"
DS_BULK  <- "GSE100054 (bulk PBMC)"

load_deg_results <- function() {
  if (!is.null(.ui_env$deg)) return(.ui_env$deg)
  ui <- .ui_data()
  status <- setNames(ui$verdicts$final_status, toupper(ui$verdicts$gene))
  
  fmt_edger <- function(d, comp) {
    data.frame(gene = d$gene, log2FC = d$logFC, pvalue = d$PValue, padj = d$FDR,
               avg_expr = d$logCPM, comparison = comp, stringsAsFactors = FALSE)
  }
  x <- list(fmt_edger(ui$xiong_de, "PD vs Control"))
  nm_map <- c(Early_vs_Control = "Early PD vs Control", Late_vs_Control = "Late PD vs Control",
              Late_vs_Early = "Late PD vs Early PD")
  for (n in names(ui$stage_de)) if (!is.null(ui$stage_de[[n]]))
    x[[length(x) + 1]] <- fmt_edger(ui$stage_de[[n]], nm_map[[n]])
  x <- do.call(rbind, x); x$dataset <- DS_XIONG
  x$validation_status <- ifelse(is.na(status[toupper(x$gene)]), "Not a candidate",
                                status[toupper(x$gene)])
  x$validation_status[x$comparison != "PD vs Control"] <- "n/a (stage contrast)"
  
  b <- ui$bulk_de
  b <- data.frame(gene = b$gene, log2FC = b$logFC, pvalue = b$P.Value, padj = b$adj.P.Val,
                  avg_expr = b$AveExpr, comparison = "PD vs Control", dataset = DS_BULK,
                  validation_status = "n/a (validation dataset)", stringsAsFactors = FALSE)
  d <- rbind(x, b)
  d$direction <- ifelse(d$log2FC >= 0, "Upregulated", "Downregulated")
  .ui_env$deg <- d[, c("gene", "log2FC", "pvalue", "padj", "avg_expr", "direction",
                       "dataset", "comparison", "validation_status")]
  .ui_env$deg
}

# NK pseudobulk log2(CPM+1) per donor: list(logcpm = genes x donors, info = sample_id/group/stage)
load_pseudobulk <- function() .ui_data()$pseudobulk

# GSE100054 gene-level expression subset: list(expr = genes x samples, group = factor)
load_bulk_expression <- function() .ui_data()$bulk_expr


pathway_label <- function(x) {
  x <- sub("^(HALLMARK|KEGG|REACTOME|GOBP)_", "", x)
  x <- tolower(gsub("_", " ", x))
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}

PATHWAY_COLLECTIONS <- c("GO Biological Process" = "GO_BP", "Reactome" = "Reactome",
                         "Hallmark" = "Hallmark", "KEGG" = "KEGG")

load_pathway_results <- function(collection = "KEGG") {
  ui <- .ui_data()
  g <- ui$gsea_xiong
  g <- g[g$collection == collection, ]
  r <- ui$pathway_replication
  r <- r[r$collection == collection, c("pathway", "NES_bulk", "P_bulk")]
  g <- merge(g, r, by = "pathway", all.x = TRUE)
  g$term <- pathway_label(g$pathway)
  g$direction <- ifelse(g$NES > 0, "Up in PD", "Down in PD")
  g <- g[order(g$pval), c("term", "pathway", "collection", "NES", "pval", "padj", "size",
                          "direction", "NES_bulk", "P_bulk", "leadingEdge")]
  rownames(g) <- NULL
  g
}

load_gsea_curve_pathways <- function() {
  cv <- .ui_data()$gsea_curves
  if (is.null(cv)) return(character(0))
  g <- .ui_data()$gsea_xiong
  p <- unique(cv$pathway)
  p <- p[order(g$pval[match(p, g$pathway)])]
  setNames(p, pathway_label(p))
}

load_gsea_curve <- function(pathway) {
  cv <- .ui_data()$gsea_curves
  d <- cv[cv$pathway == pathway, c("rank_position", "enrichment_score")]
  d$gene_set <- pathway_label(pathway)
  d
}

load_gsea_stats <- function(pathway) {
  g <- .ui_data()$gsea_xiong
  g[g$pathway == pathway, ][1, ]
}


load_validation_results <- function() {
  d <- .ui_data()$verdicts
  d[, c("gene", "set", "tier", "logFC_xiong", "P_xiong", "FDR_xiong", "logFC_bulk", "P_bulk",
        "FDR_candset", "concordant", "robust_to_NK_adjustment", "final_status")]
}

load_set_level_tests <- function() .ui_data()$set_level_tests


.auc_fast <- function(score, y) {
  r <- rank(score); n1 <- sum(y); n0 <- sum(!y)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

xiong_direction <- function(gene) {
  de <- .ui_data()$xiong_de
  lf <- de$logFC[match(toupper(gene), toupper(de$gene))]
  ifelse(is.na(lf) || lf >= 0, 1, -1)
}

# Genes that can be evaluated: candidate genes measurable in bulk + NK markers present in bulk
load_biomarker_list <- function() {
  if (!is.null(.ui_env$bm)) return(.ui_env$bm)
  ui <- .ui_data(); be <- ui$bulk_expr
  y <- be$group == "PD"
  v <- ui$verdicts
  genes <- unique(c(v$gene[!is.na(v$P_bulk)], intersect(c("NKG7", "KLRD1", "GNLY", "PRF1", "GZMB", "XCL1", "XCL2"), rownames(be$expr))))
  genes <- intersect(genes, rownames(be$expr))
  out <- do.call(rbind, lapply(genes, function(g) {
    dir <- xiong_direction(g)
    data.frame(gene = g,
               AUC = round(.auc_fast(dir * be$expr[g, ], y), 3),
               xiong_direction = ifelse(dir > 0, "Up in PD (Xiong)", "Down in PD (Xiong)"),
               stringsAsFactors = FALSE)
  }))
  p <- v$P_xiong[match(out$gene, v$gene)]
  out$P_xiong <- p
  out <- out[order(ifelse(is.na(out$P_xiong), 1, out$P_xiong)), ]
  rownames(out) <- NULL
  .ui_env$bm <- out
  out
}

# Per-sample bulk expression for one gene: sample_id, group, expression, gene
load_biomarker_expression <- function(gene) {
  be <- .ui_data()$bulk_expr
  data.frame(sample_id = colnames(be$expr), group = as.character(be$group),
             expression = as.numeric(be$expr[gene, ]), gene = gene,
             stringsAsFactors = FALSE)
}

load_signature_roc <- function() .ui_data()$roc_summary

# ---------------------------------------------------------------
# 7. NK CANDIDATE GENES (NK Cell Explorer table / downloads)
# ---------------------------------------------------------------
load_nk_candidates <- function() {
  d <- .ui_data()$verdicts
  d <- d[d$set == "tier2_primary", ]
  d[, c("gene", "tier", "logFC_xiong", "P_xiong", "FDR_xiong", "logFC_bulk", "P_bulk", "final_status")]
}

# ---------------------------------------------------------------
# 8. DASHBOARD SUMMARY STATS (Home / Dashboard page)
# ---------------------------------------------------------------
load_dashboard_stats <- function() {
  meta <- load_dataset_metadata()
  base <- list(
    n_datasets = nrow(meta),
    n_pd = sum(meta$n_pd, na.rm = TRUE),
    n_control = sum(meta$n_control, na.rm = TRUE),
    cells_profiled = NA, nk_cells = NA, n_biomarkers = NA, n_reproducible = NA)
  if (!has_real_results()) return(base)
  ui <- .ui_data()
  base$cells_profiled <- ui$singlecell$n_cells_total
  base$nk_cells <- ui$singlecell$n_nk_total
  base$n_biomarkers <- sum(grepl("^Tier", ui$verdicts$tier) & ui$verdicts$set == "tier2_primary")
  base$n_reproducible <- sum(ui$verdicts$final_status == "REPRODUCIBLE NK BIOMARKER")
  base
}

# ---------------------------------------------------------------
# Shared plotting helpers
# ---------------------------------------------------------------

# Simple themed ggplot base for consistency across the app
theme_pdnk <- function() {
  theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right"
    )
}