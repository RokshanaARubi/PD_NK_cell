#!/usr/bin/env Rscript
# ==============================================================================
# PD-NK Explorer
# scripts/11_export_ui_data.R
#
# Builds ONE light file for the Shiny app so it never has to load the 55k-cell
# Seurat object:   results/validated_data/ui_app_data.rds
#
# Contains (all REAL, nothing simulated):
#   singlecell   : UMAP + metadata for a stratified subsample of cells, sparse
#                  log-normalised expression for ~200 genes, and FULL-data
#                  cell-type counts
#   pseudobulk   : NK pseudobulk log2(CPM+1) per donor (for the DGE heatmap)
#   bulk_expr    : GSE100054 gene-level expression for candidate / top genes
#   gsea_curves  : running enrichment-score curves for top Xiong pathways
#   + everything from script 10 (ui_results.rds)
#
# Coarse cell types are assigned by explicit marker RULES applied to the cluster
# marker audit (results/qc/xiong_scrna/cluster_marker_audit_v2.csv), not by hand:
#   NK       : clusters called NK in script 06
#   B cells  : mean MS4A1 > 1
#   Monocytes: mean LYZ > 2 and mean NKG7 < 1
#   T cells  : mean CD3E > 1
#   otherwise: "Other / unassigned"  (mixed or unresolved clusters)
#
# Run AFTER scripts 06-10.  Safe to re-run.  Needs: Seurat, Matrix, dplyr, readr,
# stringr, tibble, (msigdbr optional, for GSEA curves)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(readr)
  library(stringr); library(tibble)
})
set.seed(42)

root <- getwd()
P <- function(...) file.path(root, ...)
nk_rds   <- P("results", "validated_data", "xiong_scrna_nk_annotated.rds")
ui_rds   <- P("results", "validated_data", "ui_results.rds")
val_rds  <- P("results", "validated_data", "gse100054_validation.rds")
fun_rds  <- P("results", "validated_data", "functional_analysis.rds")
audit_f  <- P("results", "qc", "xiong_scrna", "cluster_marker_audit_v2.csv")
out_rds  <- P("results", "validated_data", "ui_app_data.rds")
for (f in c(nk_rds, ui_rds, val_rds, fun_rds, audit_f))
  if (!file.exists(f)) stop("REQUIRED FILE MISSING: ", f, "  (run scripts 06-10 first)")

NK_MARKER_GENES <- c("NKG7", "KLRD1", "GNLY", "PRF1", "GZMB", "XCL1", "XCL2")
EXTRA_GENES <- c("NCAM1", "FCGR3A", "KLRF1", "NCR1", "KLRB1", "CD3D", "CD3E", "CD3G", "CD8A",
                 "IL7R", "MS4A1", "CD14", "LYZ", "CD79A", "CST3", "S100A8", "S100A9")
N_CELLS_UI <- 12000

cat("Loading annotated Seurat object (1-2 min)...\n")
obj <- readRDS(nk_rds)
ui  <- readRDS(ui_rds)
cat("Cells: ", ncol(obj), "\n")

# ------------------------------------------------------------------------
# 1. Cell types from marker rules
# ------------------------------------------------------------------------
audit <- read_csv(audit_f, show_col_types = FALSE) %>% mutate(cluster = as.character(cluster))
g <- function(a, nm) if (nm %in% names(a)) a[[nm]] else rep(0, nrow(a))
audit$cell_type <- case_when(
  audit$called_NK ~ "NK cells",
  g(audit, "MS4A1") > 1 ~ "B cells",
  g(audit, "LYZ") > 2 & g(audit, "NKG7") < 1 ~ "Monocytes",
  g(audit, "CD3E") > 1 ~ "T cells",
  TRUE ~ "Other / unassigned")
cat("\nCluster -> cell type:\n"); print(as.data.frame(audit[, c("cluster", "cell_type")]))

md <- obj@meta.data
md$cell_id <- rownames(md)
md$cluster <- as.character(md$seurat_clusters)
md$cell_type <- audit$cell_type[match(md$cluster, audit$cluster)]
md$group <- ifelse(str_detect(md$disease_state, "Parkinson"), "PD", "Control")
md$stage <- case_when(str_detect(md$disease_state, "Early") ~ "Early PD",
                      str_detect(md$disease_state, "Late") ~ "Late PD", TRUE ~ "Control")
# sanity: NK call must agree with script 06
stopifnot(all((md$cell_type == "NK cells") == (md$cell_type_call_v2 == "NK cell")))

umap_name <- grep("umap", names(obj@reductions), ignore.case = TRUE, value = TRUE)[1]
if (is.na(umap_name)) stop("No UMAP reduction found in object")
emb <- Embeddings(obj, umap_name)
md$UMAP_1 <- emb[md$cell_id, 1]; md$UMAP_2 <- emb[md$cell_id, 2]

counts_full <- md %>% count(cell_type, group, stage, sample_id, name = "n_cells")

# stratified subsample (by cluster, proportional, min 50 per cluster)
frac <- N_CELLS_UI / nrow(md)
keep_ids <- unlist(lapply(split(md$cell_id, md$cluster), function(ids)
  sample(ids, min(length(ids), max(50, round(length(ids) * frac))))))
sub <- md[md$cell_id %in% keep_ids, c("cell_id", "UMAP_1", "UMAP_2", "cell_type", "group", "stage",
                                      "sample_id", "cluster")]
cat("UI subsample: ", nrow(sub), " of ", nrow(md), " cells\n")

# ------------------------------------------------------------------------
# 2. Expression (sparse) for marker + candidate genes
# ------------------------------------------------------------------------
cand_genes <- unique(unlist(lapply(ui$verdicts$gene, toupper)))
want <- unique(c(NK_MARKER_GENES, EXTRA_GENES, cand_genes))
norm <- LayerData(obj, assay = "RNA", layer = "data")
genes_present <- intersect(want, rownames(norm))
cat("Expression genes exported: ", length(genes_present), " of ", length(want), " requested\n")
expr <- norm[genes_present, sub$cell_id, drop = FALSE]

# ------------------------------------------------------------------------
# 3. NK pseudobulk log2(CPM+1) per donor (same samples/cells as script 06)
# ------------------------------------------------------------------------
nk_ids <- md$cell_id[md$cell_type == "NK cells"]
cnt <- LayerData(obj, assay = "RNA", layer = "counts")[, nk_ids, drop = FALSE]
samp <- md$sample_id[match(nk_ids, md$cell_id)]
pb <- sapply(sort(unique(samp)), function(s) Matrix::rowSums(cnt[, samp == s, drop = FALSE]))
colnames(pb) <- sort(unique(samp))
pb <- pb[rownames(pb) %in% ui$xiong_de$gene, , drop = FALSE]
logcpm <- log2(t(t(pb) / colSums(pb) * 1e6) + 1)
pb_info <- md %>% distinct(sample_id, group, stage) %>% filter(sample_id %in% colnames(logcpm)) %>%
  arrange(match(sample_id, colnames(logcpm)))
cat("Pseudobulk matrix: ", nrow(logcpm), " genes x ", ncol(logcpm), " donors\n")
rm(cnt, norm, obj); invisible(gc())

# ------------------------------------------------------------------------
# 4. Bulk (GSE100054) expression subset
# ------------------------------------------------------------------------
v <- readRDS(val_rds)
Gb <- v$gene_expression
top_bulk <- head(v$bulk_de_modelA$gene[order(v$bulk_de_modelA$P.Value)], 1500)
bulk_genes <- intersect(unique(c(top_bulk, cand_genes, NK_MARKER_GENES, EXTRA_GENES, v$nk_proxy_genes)), rownames(Gb))
bulk_expr <- Gb[bulk_genes, , drop = FALSE]
bulk_group <- v$group
cat("Bulk expression exported: ", nrow(bulk_expr), " genes x ", ncol(bulk_expr), " samples\n")

# ------------------------------------------------------------------------
# 5. GSEA running-ES curves for top Xiong Hallmark / KEGG pathways
# ------------------------------------------------------------------------
fun <- readRDS(fun_rds)
curves <- NULL
if (requireNamespace("msigdbr", quietly = TRUE)) {
  library(msigdbr)
  new_api <- "collection" %in% names(formals(msigdbr))
  H <- tryCatch(if (new_api) msigdbr(species = "Homo sapiens", collection = "H") else msigdbr(species = "Homo sapiens", category = "H"), error = function(e) NULL)
  K <- tryCatch(if (new_api) msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") else msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"), error = function(e) NULL)
  sets_df <- bind_rows(Filter(Negate(is.null), list(H, K)))
  if (nrow(sets_df) > 0) {
    sets <- lapply(split(toupper(sets_df$gene_symbol), sets_df$gs_name), unique)
    de <- ui$xiong_de
    st <- sign(de$logFC) * -log10(pmax(de$PValue, 1e-300)); names(st) <- toupper(de$gene)
    st <- sort(st[!duplicated(names(st))], decreasing = TRUE)
    running_es <- function(stats, hits) {
      n <- length(stats); hit <- names(stats) %in% hits; nh <- sum(hit)
      if (nh < 2 || nh == n) return(NULL)
      w <- abs(stats); es <- cumsum(ifelse(hit, w, 0)) / sum(w[hit]) - cumsum(!hit) / (n - nh)
      idx <- sort(unique(c(seq(1, n, by = max(1, n %/% 400)), which(hit))))
      data.frame(rank_position = idx / n, enrichment_score = es[idx])
    }
    top_paths <- fun$gsea_xiong_primary %>% filter(collection %in% c("Hallmark", "KEGG")) %>%
      group_by(collection) %>% slice_min(pval, n = 30) %>% ungroup() %>% pull(pathway)
    curves <- bind_rows(lapply(top_paths, function(p) {
      r <- running_es(st, sets[[p]]); if (is.null(r)) NULL else cbind(pathway = p, r)
    }))
    cat("GSEA curves: ", length(unique(curves$pathway)), " pathways\n")
  }
} else cat("[skip] msigdbr not installed -> no GSEA curves\n")

# ------------------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------------------
ui$singlecell <- list(cells = as.data.frame(sub), expr = expr, counts_full = counts_full,
                      n_cells_total = nrow(md), n_cells_shown = nrow(sub),
                      n_nk_total = sum(md$cell_type == "NK cells"),
                      cluster_map = as.data.frame(audit[, c("cluster", "cell_type")]))
ui$pseudobulk <- list(logcpm = logcpm, info = as.data.frame(pb_info))
ui$bulk_expr  <- list(expr = bulk_expr, group = bulk_group)
ui$gsea_curves <- curves
saveRDS(ui, out_rds)
cat("\nSaved: results/validated_data/ui_app_data.rds  (",
    round(file.size(out_rds) / 1e6, 1), " MB)\n", sep = "")
cat("Now replace R/helper_functions.R and modules/*.R with the updated files and run the app.\n")