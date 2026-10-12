
pkgs <- c("Seurat", "Matrix", "dplyr", "stringr", "ggplot2", "tibble", "readr", "edgeR", "limma")
miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(miss) > 0) {
  stop("Missing packages: ", paste(miss, collapse = ", "), "\n",
       "CRAN ones: install.packages(...). For edgeR/limma run:\n",
       "  options(timeout = 600)\n",
       "  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
       "  BiocManager::install(c('edgeR','limma'))")
}
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(stringr)
  library(ggplot2); library(tibble); library(readr); library(edgeR); library(limma)
})

set.seed(42)
root    <- getwd()
rds_in  <- file.path(root, "results", "validated_data", "xiong_scrna_processed.rds")
rds_out <- file.path(root, "results", "validated_data", "xiong_scrna_nk_annotated.rds")
qc_dir  <- file.path(root, "results", "qc", "xiong_scrna")
de_dir  <- file.path(root, "results", "de", "xiong_nk")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(rds_in)) stop("REQUIRED FILE MISSING: ", rds_in, "\nRun scripts/05 first.")

cat("Loading processed Seurat object (may take a minute)...\n")
merged <- readRDS(rds_in)
cat("Loaded: ", ncol(merged), " cells x ", nrow(merged), " genes\n\n")

# ------------------------------------------------------------------------
# A. Corrected NK-cluster call
# ------------------------------------------------------------------------
md <- merged@meta.data
md$cluster <- as.character(md$seurat_clusters)

cs <- md %>%
  group_by(cluster) %>%
  summarise(n_cells = n(), mean_NK = mean(NK_score1), mean_T = mean(Tcell_score1), .groups = "drop") %>%
  arrange(desc(mean_NK))

# Data-driven tier: cut at the largest gap in NK score among the top half of clusters.
gaps <- -diff(cs$mean_NK)
k <- which.max(gaps[seq_len(floor(nrow(cs) / 2))])
top_tier <- cs[seq_len(k), ]
nk_clusters <- top_tier$cluster[top_tier$mean_T < 0 & top_tier$mean_NK > 0.5]

cat("Cluster NK/T scores (sorted by NK score):\n"); print(as.data.frame(cs))
cat("\nLargest NK-score gap is after rank ", k, " -> top tier: ", paste(top_tier$cluster, collapse = ", "), "\n", sep = "")
cat("NK clusters (top tier AND mean T-cell score < 0 AND NK score > 0.5): ",
    paste(nk_clusters, collapse = ", "), "\n\n", sep = "")
if (length(nk_clusters) == 0) stop("No NK clusters identified - inspect cluster scores above.")

merged$cell_type_call_v2 <- ifelse(as.character(merged$seurat_clusters) %in% nk_clusters,
                                   "NK cell", "Other / not NK")

# Marker audit: mean log-normalised expression per cluster
audit_genes <- intersect(c("NKG7", "GNLY", "NCAM1", "KLRD1", "FCGR3A", "CD3E", "CD3D", "CD8A",
                           "CD14", "LYZ", "MS4A1", "IL7R"), rownames(merged))
dat <- as.matrix(LayerData(merged, assay = "RNA", layer = "data")[audit_genes, , drop = FALSE])
audit <- t(rowsum(t(dat), group = as.character(merged$seurat_clusters))) 
audit <- sweep(audit, 2, as.numeric(table(as.character(merged$seurat_clusters))[colnames(audit)]), "/")
audit_df <- as.data.frame(t(audit)) %>% rownames_to_column("cluster") %>%
  mutate(called_NK = cluster %in% nk_clusters)
write_csv(audit_df, file.path(qc_dir, "cluster_marker_audit_v2.csv"))
cat("Marker audit (mean log-norm expression; NK clusters should be NKG7/GNLY/NCAM1-high, CD3E/CD14/MS4A1-low):\n")
print(audit_df %>% mutate(across(where(is.numeric), ~ round(.x, 2))))

nk_by_sample <- merged@meta.data %>%
  group_by(sample_id, disease_state) %>%
  summarise(n_cells = n(), n_NK = sum(cell_type_call_v2 == "NK cell"),
            pct_NK = round(100 * n_NK / n_cells, 2), .groups = "drop")
write_csv(nk_by_sample, file.path(qc_dir, "nk_cell_counts_by_sample_v2.csv"))
cat("\nCorrected NK counts per sample:\n"); print(as.data.frame(nk_by_sample))

ggsave(file.path(qc_dir, "plots", "umap_nk_calls_v2.png"),
       DimPlot(merged, group.by = "cell_type_call_v2") + labs(title = "NK call v2 (corrected)"),
       width = 7, height = 6, dpi = 150)

saveRDS(merged, rds_out)
cat("\nSaved annotated object: results/validated_data/xiong_scrna_nk_annotated.rds\n\n")

# ------------------------------------------------------------------------
# B. Pseudobulk DE in NK cells
# ------------------------------------------------------------------------
nk <- subset(merged, subset = cell_type_call_v2 == "NK cell")
counts <- LayerData(nk, assay = "RNA", layer = "counts")
samples <- sort(unique(nk$sample_id))
n_per <- table(nk$sample_id)
MIN_CELLS <- 50
keep_samples <- samples[n_per[samples] >= MIN_CELLS]
cat("NK cells per sample: ", paste(names(n_per), n_per, sep = "=", collapse = ", "), "\n")
if (length(keep_samples) < length(samples))
  cat("Dropping samples with < ", MIN_CELLS, " NK cells: ", paste(setdiff(samples, keep_samples), collapse = ", "), "\n")

pb <- sapply(keep_samples, function(s) Matrix::rowSums(counts[, nk$sample_id == s, drop = FALSE]))
colnames(pb) <- keep_samples

sinfo <- merged@meta.data %>% distinct(sample_id, disease_state) %>%
  filter(sample_id %in% keep_samples) %>% arrange(match(sample_id, keep_samples)) %>%
  mutate(group = case_when(str_detect(disease_state, "Normal|Control|Healthy") ~ "Control",
                           TRUE ~ "PD"),
         stage = case_when(str_detect(disease_state, "Early") ~ "Early",
                           str_detect(disease_state, "Late") ~ "Late",
                           TRUE ~ "Control"))
print(sinfo)
if (length(unique(sinfo$group)) < 2) stop("Need both Control and PD samples after filtering.")

run_edger <- function(pb, group, contrasts_txt, tag) {
  y <- DGEList(counts = pb, group = group)
  y <- y[filterByExpr(y, group = group), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ 0 + group); colnames(design) <- levels(group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  out <- list()
  for (nm in names(contrasts_txt)) {
    con <- makeContrasts(contrasts = contrasts_txt[[nm]], levels = design)
    tt <- topTags(glmQLFTest(fit, contrast = con), n = Inf)$table %>%
      rownames_to_column("gene") %>% mutate(contrast = nm, .before = 1)
    write_csv(tt, file.path(de_dir, paste0("de_", tag, "_", nm, ".csv")))
    out[[nm]] <- tt
    cat(sprintf("  %-18s genes tested: %d | FDR<0.05: %d | nominal P<0.01: %d\n",
                nm, nrow(tt), sum(tt$FDR < 0.05), sum(tt$PValue < 0.01)))
  }
  list(y = y, res = out)
}

cat("\nPrimary contrast: PD vs Control (pseudobulk, edgeR QL)\n")
g1 <- factor(sinfo$group, levels = c("Control", "PD"))
primary <- run_edger(pb, g1, list(PD_vs_Control = "PD - Control"), "primary")

if (all(c("Control", "Early", "Late") %in% sinfo$stage)) {
  cat("\nStage contrasts (exploratory, 2 donors per group):\n")
  g2 <- factor(sinfo$stage, levels = c("Control", "Early", "Late"))
  stage <- run_edger(pb, g2, list(Early_vs_Control = "Early - Control",
                                  Late_vs_Control = "Late - Control",
                                  Late_vs_Early = "Late - Early"), "stage")
}

png(file.path(de_dir, "pseudobulk_mds.png"), width = 800, height = 700, res = 120)
plotMDS(primary$y, col = ifelse(g1 == "PD", "firebrick", "steelblue"),
        main = "NK pseudobulk MDS (red = PD, blue = Control)")
dev.off()

res <- primary$res$PD_vs_Control
vol <- ggplot(res, aes(logFC, -log10(PValue), colour = FDR < 0.05)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text(data = head(res[order(res$PValue), ], 12), aes(label = gene),
            size = 3, vjust = -0.6, check_overlap = TRUE, colour = "black") +
  scale_colour_manual(values = c("grey60", "firebrick")) +
  labs(title = "NK cells: PD vs Control (pseudobulk)", x = "log2 fold change", y = "-log10 P") +
  theme_minimal()
ggsave(file.path(de_dir, "volcano_PD_vs_Control.png"), vol, width = 7, height = 6, dpi = 150)

# ------------------------------------------------------------------------
# C. Candidate shortlist (step 6)
# ------------------------------------------------------------------------
cand <- res %>%
  mutate(direction = ifelse(logFC > 0, "Up in PD", "Down in PD"),
         tier = case_when(
           FDR < 0.05 & abs(logFC) >= 0.58 ~ "Tier 1: FDR<0.05 & |log2FC|>=0.58",
           PValue < 0.01 & abs(logFC) >= 1 ~ "Tier 2 (exploratory): P<0.01 & |log2FC|>=1",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(tier)) %>% arrange(tier, PValue)
write_csv(cand, file.path(de_dir, "nk_candidate_biomarkers.csv"))

cat("\n====================================================\n")
cat("NK DE SUMMARY (Xiong GSE223138, pseudobulk)\n")
cat("====================================================\n")
cat("NK clusters called:     ", paste(nk_clusters, collapse = ", "), "\n")
cat("NK cells total:         ", sum(n_per), "\n")
cat("Donors used:            ", paste(sinfo$sample_id, "(", sinfo$stage, ")", collapse = ", "), "\n")
cat("Genes tested:           ", nrow(res), "\n")
cat("Tier 1 candidates:      ", sum(str_detect(cand$tier, "^Tier 1")), "\n")
cat("Tier 2 (exploratory):   ", sum(str_detect(cand$tier, "^Tier 2")), "\n")
cat("\nTop 15 by P-value:\n"); print(head(as.data.frame(res[order(res$PValue), c("gene", "logFC", "PValue", "FDR")]), 15))
cat("\nFiles in results/de/xiong_nk/ and results/qc/xiong_scrna/ (v2 tables)\n")
cat("====================================================\n")

