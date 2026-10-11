#!/usr/bin/env Rscript
# ==============================================================================
# PD-NK Explorer
# scripts/08_functional_analysis_nk.R
#
# METHODOLOGY STEP 8: Functional (pathway) analysis
#
# WHY RANKED GSEA (not a gene list + enrichment test):
#   No Xiong NK gene passed FDR < 0.05, so any "significant gene list" would be
#   an arbitrary cut. GSEA uses ALL ~12k tested genes ranked by evidence, so it
#   needs no cut-off and is more stable with 2 vs 4 donors.
#
# WHAT IT DOES
#   1. Ranks Xiong NK genes by signed -log10(P)  (PD vs Control, plus stage contrasts)
#   2. fgsea on Hallmark, KEGG, Reactome, GO:BP (whatever msigdbr can supply)
#   3. Pathway-level replication: runs the SAME gene sets on GSE100054 (bulk,
#      ranked by limma t) and checks whether the Xiong direction (NES sign)
#      is reproduced. Pathway-level is far more powerful than single genes.
#
# READS
#   results/de/xiong_nk/de_primary_PD_vs_Control.csv (+ de_stage_*.csv if present)
#   results/validated_data/gse100054_validation.rds   (from script 07)
# WRITES
#   results/functional/xiong_nk/*.csv, plots/*.png
#   results/functional/gse100054/gsea_bulk_*.csv
#   results/functional/pathway_replication.csv
#   results/validated_data/functional_analysis.rds
#
# Install once if needed:
#   options(timeout = 600)
#   install.packages("msigdbr")
#   BiocManager::install("fgsea")
# ==============================================================================

pkgs <- c("dplyr", "readr", "stringr", "tibble", "ggplot2", "fgsea", "msigdbr")
miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(miss) > 0) stop("Missing packages: ", paste(miss, collapse = ", "),
                           "\n  install.packages('msigdbr'); BiocManager::install('fgsea')")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tibble)
  library(ggplot2); library(fgsea); library(msigdbr)
})
set.seed(42)

# ---- SETTINGS ----------------------------------------------------------------
MIN_SIZE <- 10
MAX_SIZE <- 500
PADJ_CUT <- 0.05
P_EXPLORATORY <- 0.01
# ------------------------------------------------------------------------------

root    <- getwd()
de_dir  <- file.path(root, "results", "de", "xiong_nk")
val_rds <- file.path(root, "results", "validated_data", "gse100054_validation.rds")
out_x   <- file.path(root, "results", "functional", "xiong_nk")
out_b   <- file.path(root, "results", "functional", "gse100054")
plot_x  <- file.path(out_x, "plots")
for (d in c(out_x, out_b, plot_x)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

prim_f <- file.path(de_dir, "de_primary_PD_vs_Control.csv")
if (!file.exists(prim_f)) stop("REQUIRED FILE MISSING: ", prim_f, "  (run script 06)")
if (!file.exists(val_rds)) stop("REQUIRED FILE MISSING: ", val_rds, "  (run script 07)")

# ------------------------------------------------------------------------
# 1. Gene sets
# ------------------------------------------------------------------------
get_sets <- function(label, ...) {
  df <- tryCatch(msigdbr(species = "Homo sapiens", ...), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) { cat("  [skip] ", label, ": not available\n"); return(NULL) }
  sets <- split(toupper(df$gene_symbol), df$gs_name)
  sets <- lapply(sets, unique)
  cat("  ", label, ": ", length(sets), " gene sets\n", sep = "")
  sets
}
# msigdbr >= 10 uses collection/subcollection; older uses category/subcategory
new_api <- "collection" %in% names(formals(msigdbr))
cat("Loading gene sets (msigdbr ", as.character(packageVersion("msigdbr")), ")\n", sep = "")
gs <- list()
if (new_api) {
  gs$Hallmark <- get_sets("Hallmark", collection = "H")
  gs$KEGG     <- get_sets("KEGG", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  gs$Reactome <- get_sets("Reactome", collection = "C2", subcollection = "CP:REACTOME")
  gs$GO_BP    <- get_sets("GO:BP", collection = "C5", subcollection = "GO:BP")
} else {
  gs$Hallmark <- get_sets("Hallmark", category = "H")
  gs$KEGG     <- get_sets("KEGG", category = "C2", subcategory = "CP:KEGG")
  gs$Reactome <- get_sets("Reactome", category = "C2", subcategory = "CP:REACTOME")
  gs$GO_BP    <- get_sets("GO:BP", category = "C5", subcategory = "GO:BP")
}
gs <- Filter(Negate(is.null), gs)
if (length(gs) == 0) stop("No gene sets could be loaded from msigdbr (needs internet/data package).")

# ------------------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------------------
make_rank <- function(df, gene_col, p_col, fc_col) {
  s <- sign(df[[fc_col]]) * -log10(pmax(df[[p_col]], 1e-300))
  names(s) <- toupper(df[[gene_col]])
  s <- s[!is.na(s) & !is.na(names(s)) & nzchar(names(s))]
  s <- s[!duplicated(names(s))]
  sort(s, decreasing = TRUE)
}

run_gsea <- function(rank, tag) {
  res_all <- list()
  for (nm in names(gs)) {
    r <- suppressWarnings(fgsea(pathways = gs[[nm]], stats = rank,
                                minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = 0))
    if (nrow(r) == 0) next
    r <- as_tibble(r) %>% mutate(collection = nm, .before = 1,
                                 leadingEdge = vapply(leadingEdge, function(x) paste(head(x, 25), collapse = ";"), ""))
    # BH within collection (fgsea default); also across ALL tested sets below
    res_all[[nm]] <- r
  }
  out <- bind_rows(res_all)
  if (nrow(out) == 0) return(out)
  out$padj_all <- p.adjust(out$pval, "BH")
  out <- out %>% arrange(pval)
  write_csv(out, file.path(out_x, paste0("gsea_", tag, ".csv")))
  out
}

report <- function(res, tag) {
  cat(sprintf("\n[%s] sets tested: %d | padj(within collection)<%.2f: %d | padj(all sets)<%.2f: %d | nominal P<%.2f: %d\n",
              tag, nrow(res), PADJ_CUT, sum(res$padj < PADJ_CUT, na.rm = TRUE),
              PADJ_CUT, sum(res$padj_all < PADJ_CUT, na.rm = TRUE),
              P_EXPLORATORY, sum(res$pval < P_EXPLORATORY, na.rm = TRUE)))
}

# ------------------------------------------------------------------------
# 3. GSEA on Xiong NK cells
# ------------------------------------------------------------------------
cat("\n=== GSEA: Xiong NK pseudobulk ===\n")
prim <- read_csv(prim_f, show_col_types = FALSE)
rank_primary <- make_rank(prim, "gene", "PValue", "logFC")
cat("Ranked genes: ", length(rank_primary), "\n")
g_primary <- run_gsea(rank_primary, "PD_vs_Control")
report(g_primary, "PD_vs_Control")

stage_res <- list()
for (st in c("Early_vs_Control", "Late_vs_Control", "Late_vs_Early")) {
  f <- file.path(de_dir, paste0("de_stage_", st, ".csv"))
  if (file.exists(f)) {
    d <- read_csv(f, show_col_types = FALSE)
    stage_res[[st]] <- run_gsea(make_rank(d, "gene", "PValue", "logFC"), st)
    report(stage_res[[st]], st)
  }
}

cat("\nTop 15 pathways, PD vs Control (by P):\n")
print(as.data.frame(g_primary %>% select(collection, pathway, NES, pval, padj, padj_all, size) %>%
                      head(15) %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))

# ------------------------------------------------------------------------
# 4. Same gene sets in GSE100054 (pathway-level validation)
# ------------------------------------------------------------------------
cat("\n=== GSEA: GSE100054 bulk PBMC (independent validation) ===\n")
v <- readRDS(val_rds)
bulk <- v$bulk_de_modelA
rank_bulk <- make_rank(bulk, "gene", "P.Value", "logFC")
# rank by limma t instead of P*sign(logFC) is equivalent in order; keep consistent with Xiong
cat("Ranked genes: ", length(rank_bulk), "\n")
res_b <- list()
for (nm in names(gs)) {
  r <- suppressWarnings(fgsea(pathways = gs[[nm]], stats = rank_bulk,
                              minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = 0))
  if (nrow(r) == 0) next
  res_b[[nm]] <- as_tibble(r) %>%
    mutate(collection = nm, .before = 1,
           leadingEdge = vapply(leadingEdge, function(x) paste(head(x, 25), collapse = ";"), ""))
}
g_bulk <- bind_rows(res_b) %>% arrange(pval)
write_csv(g_bulk, file.path(out_b, "gsea_bulk_PD_vs_Control.csv"))
cat("Bulk sets tested: ", nrow(g_bulk), " | nominal P<0.05: ", sum(g_bulk$pval < 0.05), "\n")

# ------------------------------------------------------------------------
# 5. Pathway-level replication
# ------------------------------------------------------------------------
rep_df <- g_primary %>%
  select(collection, pathway, size_x = size, NES_xiong = NES, P_xiong = pval, padj_xiong = padj) %>%
  inner_join(g_bulk %>% select(collection, pathway, NES_bulk = NES, P_bulk = pval, padj_bulk = padj),
             by = c("collection", "pathway")) %>%
  mutate(concordant = sign(NES_xiong) == sign(NES_bulk),
         status = case_when(
           concordant & P_xiong < 0.05 & P_bulk < 0.05 ~ "Replicated (both P<0.05, same direction)",
           !concordant & P_xiong < 0.05 & P_bulk < 0.05 ~ "CONTRADICTED (both P<0.05, opposite)",
           TRUE ~ "Not replicated / not significant in both")) %>%
  arrange(P_xiong)
write_csv(rep_df, file.path(root, "results", "functional", "pathway_replication.csv"))

sub <- rep_df %>% filter(P_xiong < 0.05)
k <- sum(sub$concordant); n <- nrow(sub)
bt <- if (n >= 5) binom.test(k, n, 0.5, alternative = "greater")$p.value else NA_real_
sp <- suppressWarnings(cor(rep_df$NES_xiong, rep_df$NES_bulk, method = "spearman"))
cat(sprintf("\nPathway-level replication:\n  all shared sets: %d | Spearman(NES) = %.3f\n", nrow(rep_df), sp))
cat(sprintf("  Xiong-nominal-significant sets (P<0.05): %d | same direction in bulk: %d (%.1f%%) | binomial P = %s\n",
            n, k, 100 * k / max(n, 1), signif(bt, 3)))
cat("  Replicated in both: ", sum(str_detect(rep_df$status, "^Replicated")),
    " | contradicted: ", sum(str_detect(rep_df$status, "^CONTRADICTED")), "\n", sep = "")
cat("\nPathways replicated in BOTH datasets:\n")
print(as.data.frame(rep_df %>% filter(str_detect(status, "^Replicated")) %>%
                      select(collection, pathway, NES_xiong, P_xiong, NES_bulk, P_bulk) %>%
                      head(20) %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))

# ------------------------------------------------------------------------
# 6. Plots
# ------------------------------------------------------------------------
tidy_name <- function(x) str_to_sentence(str_replace_all(str_remove(x, "^(HALLMARK|KEGG|REACTOME|GOBP)_"), "_", " "))
topn <- g_primary %>% filter(collection == "Hallmark") %>% arrange(pval) %>% head(15) %>%
  mutate(label = tidy_name(pathway), label = factor(label, levels = rev(label)))
if (nrow(topn) > 0) {
  p1 <- ggplot(topn, aes(NES, label, fill = NES > 0)) + geom_col() +
    scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "steelblue"),
                      labels = c(`TRUE` = "Up in PD", `FALSE` = "Down in PD"), name = NULL) +
    labs(title = "Hallmark pathways in NK cells: PD vs Control",
         subtitle = "Xiong pseudobulk, top 15 by P (exploratory - see padj in CSV)", x = "NES", y = NULL) +
    theme_minimal()
  ggsave(file.path(plot_x, "hallmark_top15_nes.png"), p1, width = 8, height = 6, dpi = 150)
}
p2 <- ggplot(rep_df %>% filter(collection %in% c("Hallmark", "KEGG", "Reactome")),
             aes(NES_xiong, NES_bulk, colour = status)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(alpha = 0.6, size = 1.6) +
  labs(title = "Pathway NES: Xiong NK vs GSE100054 bulk",
       subtitle = "Top-right / bottom-left = same direction in both datasets",
       x = "NES, Xiong NK pseudobulk", y = "NES, GSE100054 bulk") +
  theme_minimal() + theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(plot_x, "pathway_nes_xiong_vs_bulk.png"), p2, width = 8, height = 7.5, dpi = 150)

# ------------------------------------------------------------------------
# 7. Save for Shiny
# ------------------------------------------------------------------------
saveRDS(list(gsea_xiong_primary = g_primary, gsea_xiong_stage = stage_res,
             gsea_bulk = g_bulk, pathway_replication = rep_df,
             settings = list(MIN_SIZE = MIN_SIZE, MAX_SIZE = MAX_SIZE),
             run_on = Sys.time()),
        file.path(root, "results", "validated_data", "functional_analysis.rds"))

cat("\n", strrep("=", 70), "\nSTEP 8 SUMMARY\n", strrep("=", 70), "\n", sep = "")
cat("Gene-set collections used: ", paste(names(gs), collapse = ", "), "\n")
cat("Xiong pathways tested:     ", nrow(g_primary), "\n")
cat("Significant (padj<0.05):   ", sum(g_primary$padj < PADJ_CUT, na.rm = TRUE), "\n")
cat("Replicated in GSE100054:   ", sum(str_detect(rep_df$status, "^Replicated")), "\n")
cat("Outputs: results/functional/  and  results/validated_data/functional_analysis.rds\n")
cat("NEXT: script 09 (ROC)\n")
cat(strrep("=", 70), "\n")
