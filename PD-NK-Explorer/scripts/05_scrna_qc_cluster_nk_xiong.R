#!/usr/bin/env Rscript
# ==============================================================================
# PD-NK Explorer
# scripts/05_scrna_qc_cluster_nk_xiong.R
#
# PURPOSE:
#   Methodology steps 2-4 on the REAL Xiong et al. 2024 (GSE223138) data,
#   now that scripts/04_acquire_xiong_scrna_supplementary.R has downloaded
#   all 6 samples' barcodes/features/matrix.mtx files:
#     Step 2: Quality Control        (per-cell filtering, real thresholds)
#     Step 3: Cell Clustering        (PCA -> graph clustering -> UMAP)
#     Step 4: NK Cell Identification (canonical marker-based, not fabricated)
#
#   This script does NOT perform differential expression, pathway analysis,
#   cross-dataset validation, or biomarker selection (methodology steps
#   5-10). Those remain separate, later scripts.
#
# IMPORTANT ON THRESHOLDS:
#   QC cutoffs (min genes/cell, max mito%) and clustering resolution are
#   inherently judgment calls in single-cell analysis. This script uses
#   adaptive, MAD-based (median absolute deviation) thresholds rather than
#   arbitrary fixed numbers, and labels every such choice explicitly in the
#   log. ALWAYS inspect the QC violin plots and UMAP plots this produces
#   before treating the resulting NK-cell calls as final - they are a
#   principled first pass, not a validated ground truth.
#
# SAFETY:
#   - Reads only from data/xiong/raw/<GSM>/ (already downloaded, untouched)
#     and data/xiong/GSE223138_sample_metadata.csv (already validated).
#   - Does not modify or delete any raw file.
#   - Does not perform differential expression or any hypothesis testing.
#   - NK-cell calls require BOTH positive NK markers AND negative T-cell
#     markers - a cluster is never called "NK" from a single marker alone.
# ==============================================================================

required_packages <- c("Seurat", "Matrix", "dplyr", "stringr", "ggplot2", "tibble", "readr")
missing_pkgs <- required_packages[!vapply(required_packages, requireNamespace,
                                          quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install with: install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(readr)
})

set.seed(42)  # reproducibility for clustering/UMAP stochastic steps

# ------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------

project_root  <- getwd()
xiong_dir     <- file.path(project_root, "data", "xiong")
raw_dir       <- file.path(xiong_dir, "raw")
meta_path     <- file.path(xiong_dir, "GSE223138_sample_metadata.csv")
validated_dir <- file.path(project_root, "results", "validated_data")
qc_dir        <- file.path(project_root, "results", "qc", "xiong_scrna")
qc_plot_dir   <- file.path(qc_dir, "plots")

dir.create(validated_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_plot_dir, recursive = TRUE, showWarnings = FALSE)

qc_log <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  qc_log <<- c(qc_log, msg)
}

if (!dir.exists(raw_dir)) {
  stop("REQUIRED DIRECTORY MISSING: ", raw_dir, "\n",
       "Run scripts/04_acquire_xiong_scrna_supplementary.R first.")
}
if (!file.exists(meta_path)) {
  stop("REQUIRED FILE MISSING: ", meta_path)
}

sample_meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)

log_msg("=================================================================")
log_msg("XIONG GSE223138 - scRNA-seq QC, CLUSTERING, NK-CELL IDENTIFICATION")
log_msg("Generated: ", as.character(Sys.time()))
log_msg("=================================================================\n")

# ------------------------------------------------------------------------
# 1. Load each sample's real 10x files (filenames carry a GEO-specific
#    numeric prefix, e.g. GSM6940192_16_matrix.mtx.gz, so Seurat::ReadMtx()
#    is used with explicit file paths rather than Read10X(), which expects
#    exact "matrix.mtx.gz" / "barcodes.tsv.gz" / "features.tsv.gz" names).
# ------------------------------------------------------------------------

sample_dirs <- list.dirs(raw_dir, recursive = FALSE)
if (length(sample_dirs) == 0) {
  stop("No per-sample subdirectories found under ", raw_dir,
       ". Expected one subfolder per GSM accession.")
}

find_one <- function(dir, pattern) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

seurat_list <- list()
load_log <- data.frame(gsm = character(0), n_cells = integer(0), n_genes = integer(0),
                       status = character(0), stringsAsFactors = FALSE)

for (sdir in sample_dirs) {
  gsm <- basename(sdir)
  mtx_path   <- find_one(sdir, "matrix\\.mtx\\.gz$")
  bc_path    <- find_one(sdir, "barcodes\\.tsv\\.gz$")
  feat_path  <- find_one(sdir, "features\\.tsv\\.gz$|genes\\.tsv\\.gz$")
  
  if (any(is.na(c(mtx_path, bc_path, feat_path)))) {
    log_msg(gsm, ": SKIPPED - missing one or more of matrix/barcodes/features files.")
    load_log <- rbind(load_log, data.frame(gsm = gsm, n_cells = NA, n_genes = NA,
                                           status = "SKIPPED - incomplete files", stringsAsFactors = FALSE))
    next
  }
  
  counts <- tryCatch(
    Seurat::ReadMtx(mtx = mtx_path, cells = bc_path, features = feat_path, feature.column = 2),
    error = function(e) {
      log_msg(gsm, ": FAILED to read 10x files - ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(counts)) {
    load_log <- rbind(load_log, data.frame(gsm = gsm, n_cells = NA, n_genes = NA,
                                           status = "FAILED to parse", stringsAsFactors = FALSE))
    next
  }
  
  so <- CreateSeuratObject(counts = counts, project = gsm, min.cells = 0, min.features = 0)
  so$sample_id <- gsm
  seurat_list[[gsm]] <- so
  
  log_msg(gsm, ": loaded ", ncol(so), " cells x ", nrow(so), " genes.")
  load_log <- rbind(load_log, data.frame(gsm = gsm, n_cells = ncol(so), n_genes = nrow(so),
                                         status = "Loaded", stringsAsFactors = FALSE))
}

if (length(seurat_list) == 0) {
  stop("No samples could be loaded. Check results/qc/xiong_scrna/qc_report.txt for details.")
}

write_csv(load_log, file.path(qc_dir, "sample_load_log.csv"))

# ------------------------------------------------------------------------
# 2. Merge samples and attach REAL disease-state metadata
#    (from the already-validated GSE223138_sample_metadata.csv - never
#    fabricated or inferred).
# ------------------------------------------------------------------------

log_msg("\n--- Merging ", length(seurat_list), " samples ---")

merged <- if (length(seurat_list) > 1) {
  merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
} else {
  seurat_list[[1]]
}

# Seurat v5 keeps each merged sample's counts as a separate "layer" until
# explicitly joined. Several downstream steps (notably FindAllMarkers,
# which needs to compare cells across all samples at once) silently
# produce empty/wrong results without this - join immediately after merge.
merged <- JoinLayers(merged)
log_msg("Joined per-sample count layers into one unified assay (required for ",
        "cross-sample marker detection in Seurat v5).")

# NOTE: disease_lookup[merged$sample_id] returns a vector *named* by GSM
# accession (from disease_lookup's own names), not by cell barcode. Seurat's
# `$<-` metadata assignment matches by name when names are present, finds no
# overlap with the object's cell barcodes, and fails ("No cell overlap...").
# unname() makes this a plain positional vector, which is correct here since
# it's already aligned to merged's cell order via merged$sample_id.
disease_lookup <- setNames(sample_meta[["disease state:ch1"]], sample_meta[["geo_accession"]])
merged$disease_state <- unname(disease_lookup[merged$sample_id])
merged$pd_status <- dplyr::case_when(
  str_detect(merged$disease_state, "Parkinson") ~ "PD",
  str_detect(merged$disease_state, "Normal|Control|Healthy") ~ "Control",
  TRUE ~ NA_character_
)

log_msg("Merged object: ", ncol(merged), " cells x ", nrow(merged), " genes, ", length(seurat_list), " samples.")
log_msg("Real disease-state breakdown (from GSE223138_sample_metadata.csv): ")
print(table(merged$disease_state, useNA = "ifany"))

# ------------------------------------------------------------------------
# 3. QC metrics + adaptive (MAD-based) filtering
# ------------------------------------------------------------------------

log_msg("\n--- QC metrics ---")

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

qc_pre <- tibble(
  cell_id = colnames(merged),
  sample_id = merged$sample_id,
  nCount_RNA = merged$nCount_RNA,
  nFeature_RNA = merged$nFeature_RNA,
  percent_mt = merged$percent.mt
)
write_csv(qc_pre, file.path(qc_dir, "qc_metrics_pre_filter.csv"))

# Adaptive thresholds: median +/- N*MAD, computed PER SAMPLE (different
# samples/runs can have different depth), rather than one fixed number
# applied to the whole cohort. This is a standard, defensible heuristic
# (e.g. used by scater::isOutlier) - NOT an arbitrary fixed cutoff.
MAD_MULTIPLIER <- 3

compute_adaptive_bounds <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)
  c(lower = med - MAD_MULTIPLIER * mad_val, upper = med + MAD_MULTIPLIER * mad_val)
}

keep_cell <- rep(FALSE, ncol(merged))
names(keep_cell) <- colnames(merged)

per_sample_thresholds <- list()

for (gsm in unique(merged$sample_id)) {
  idx <- which(merged$sample_id == gsm)
  nf <- merged$nFeature_RNA[idx]
  nc <- merged$nCount_RNA[idx]
  pmt <- merged$percent.mt[idx]
  
  nf_bounds <- compute_adaptive_bounds(nf)
  nc_bounds <- compute_adaptive_bounds(nc)
  # Mitochondrial % only gets an UPPER bound (high mito = likely dying cell;
  # low mito is not itself a quality problem), floored at a sane minimum of 5%
  # so healthy low-mito samples aren't over-filtered.
  pmt_upper <- max(median(pmt, na.rm = TRUE) + MAD_MULTIPLIER * mad(pmt, na.rm = TRUE), 5)
  
  keep <- nf >= max(200, nf_bounds["lower"]) & nf <= nf_bounds["upper"] &
    nc >= nc_bounds["lower"] & nc <= nc_bounds["upper"] &
    pmt <= pmt_upper
  
  keep_cell[idx] <- keep
  
  per_sample_thresholds[[gsm]] <- data.frame(
    gsm = gsm, n_cells_total = length(idx), n_cells_kept = sum(keep),
    nFeature_lower = max(200, nf_bounds["lower"]), nFeature_upper = nf_bounds["upper"],
    nCount_lower = nc_bounds["lower"], nCount_upper = nc_bounds["upper"],
    percent_mt_upper = pmt_upper
  )
  
  log_msg(gsm, ": ", sum(keep), " / ", length(idx), " cells pass QC ",
          "(nFeature in [", round(max(200, nf_bounds["lower"])), ", ", round(nf_bounds["upper"]), "], ",
          "nCount in [", round(nc_bounds["lower"]), ", ", round(nc_bounds["upper"]), "], ",
          "percent.mt <= ", round(pmt_upper, 1), "%)")
}

thresholds_df <- bind_rows(per_sample_thresholds)
write_csv(thresholds_df, file.path(qc_dir, "qc_adaptive_thresholds_per_sample.csv"))

# QC violin plots (pre-filter, all cells, for visual inspection)
vln_plot <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    group.by = "sample_id", ncol = 1, pt.size = 0)
ggsave(file.path(qc_plot_dir, "qc_violin_pre_filter.png"), vln_plot, width = 8, height = 10, dpi = 150)
log_msg("\nWrote plots/qc_violin_pre_filter.png - INSPECT THIS before trusting the filtering below.")

n_before <- ncol(merged)
merged <- subset(merged, cells = names(keep_cell)[keep_cell])
n_after <- ncol(merged)

log_msg("\nTotal cells before QC filtering: ", n_before)
log_msg("Total cells after QC filtering:  ", n_after, " (",
        sprintf("%.1f%%", 100 * n_after / n_before), " retained)")

# ------------------------------------------------------------------------
# 4. Normalization, variable features, scaling, PCA
# ------------------------------------------------------------------------

log_msg("\n--- Normalization and dimensionality reduction ---")

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)

log_msg("PCA complete (30 PCs computed on 2000 variable features).")


# ------------------------------------------------------------------------
# 5. Clustering + UMAP
#    Resolution = 0.5 is a commonly used Seurat default, explicitly
#    labeled here as a tunable heuristic - not a fabricated "correct"
#    value. Re-run with a different resolution if clusters look too
#    coarse/fine on the UMAP.
# ------------------------------------------------------------------------

CLUSTER_RESOLUTION <- 0.5
N_PCS_FOR_CLUSTERING <- 20

merged <- FindNeighbors(merged, dims = 1:N_PCS_FOR_CLUSTERING, verbose = FALSE)
merged <- FindClusters(merged, resolution = CLUSTER_RESOLUTION, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:N_PCS_FOR_CLUSTERING, verbose = FALSE)

n_clusters <- length(unique(Idents(merged)))
log_msg("\nClustering complete: ", n_clusters, " clusters found ",
        "(resolution = ", CLUSTER_RESOLUTION, ", ", N_PCS_FOR_CLUSTERING, " PCs - heuristic parameters).")

umap_by_cluster <- DimPlot(merged, reduction = "umap", label = TRUE) +
  labs(title = "Xiong GSE223138 - UMAP by cluster")
ggsave(file.path(qc_plot_dir, "umap_by_cluster.png"), umap_by_cluster, width = 7, height = 6, dpi = 150)

umap_by_sample <- DimPlot(merged, reduction = "umap", group.by = "sample_id") +
  labs(title = "Xiong GSE223138 - UMAP by sample (check for batch effects)")
ggsave(file.path(qc_plot_dir, "umap_by_sample.png"), umap_by_sample, width = 7, height = 6, dpi = 150)

umap_by_disease <- DimPlot(merged, reduction = "umap", group.by = "pd_status") +
  labs(title = "Xiong GSE223138 - UMAP by disease status")
ggsave(file.path(qc_plot_dir, "umap_by_disease_status.png"), umap_by_disease, width = 7, height = 6, dpi = 150)

log_msg("Wrote plots/umap_by_cluster.png, umap_by_sample.png, umap_by_disease_status.png.")
log_msg("IMPORTANT: inspect umap_by_sample.png for batch effects - no batch-correction/integration ",
        "(e.g. Harmony, CCA) has been applied in this script. If clusters separate mainly by sample ",
        "rather than by biology, integration should be added before trusting cluster identities.")

# ------------------------------------------------------------------------
# 6. NK-cell identification: canonical markers, positive AND negative
#    evidence required - never a single-marker call.
# ------------------------------------------------------------------------

log_msg("\n--- NK-cell identification ---")

nk_positive_markers <- c("NKG7", "GNLY", "KLRD1", "NCAM1", "KLRB1", "GZMB", "PRF1", "FCGR3A")
t_cell_markers       <- c("CD3D", "CD3E", "CD3G")

present_nk_markers <- intersect(nk_positive_markers, rownames(merged))
present_t_markers  <- intersect(t_cell_markers, rownames(merged))

log_msg("NK marker genes found in data: ", paste(present_nk_markers, collapse = ", "),
        " (", length(present_nk_markers), " / ", length(nk_positive_markers), " expected)")
log_msg("T-cell marker genes found in data: ", paste(present_t_markers, collapse = ", "),
        " (", length(present_t_markers), " / ", length(t_cell_markers), " expected)")

if (length(present_nk_markers) == 0) {
  stop("None of the expected NK marker genes were found in this dataset's gene symbols. ",
       "Cannot proceed with NK identification - check whether features.tsv used gene symbols ",
       "or a different identifier (e.g. Ensembl IDs).")
}

merged <- AddModuleScore(merged, features = list(present_nk_markers), name = "NK_score")
if (length(present_t_markers) > 0) {
  merged <- AddModuleScore(merged, features = list(present_t_markers), name = "Tcell_score")
} else {
  merged$Tcell_score1 <- NA_real_
  log_msg("WARNING: no T-cell markers found in data - NK calls below rely on NK score alone ",
          "and have NOT been checked against T-cell/NKT contamination. Treat with extra caution.")
}

cluster_scores <- merged@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    mean_NK_score = mean(NK_score1, na.rm = TRUE),
    mean_Tcell_score = mean(Tcell_score1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_NK_score))

write_csv(cluster_scores, file.path(qc_dir, "cluster_nk_tcell_scores.csv"))

log_msg("\nPer-cluster NK/T-cell module scores (results/qc/xiong_scrna/cluster_nk_tcell_scores.csv):")
print(as.data.frame(cluster_scores))

# A cluster is called NK if its mean NK score is in the upper tier
# (above the overall median NK score across clusters, a relative - not
# absolute - threshold) AND its mean T-cell score is at or below the
# overall median T-cell score. Both conditions must hold.
median_nk <- median(cluster_scores$mean_NK_score, na.rm = TRUE)
median_t  <- if (length(present_t_markers) > 0) median(cluster_scores$mean_Tcell_score, na.rm = TRUE) else NA_real_

cluster_scores <- cluster_scores %>%
  mutate(
    nk_candidate = mean_NK_score > median_nk &
      (is.na(median_t) | mean_Tcell_score <= median_t)
  )

nk_clusters <- cluster_scores$seurat_clusters[cluster_scores$nk_candidate]

log_msg("\nCluster(s) called as NK-cell candidates (NK score > median AND T-cell score <= median): ",
        if (length(nk_clusters) > 0) paste(nk_clusters, collapse = ", ") else "NONE")

merged$cell_type_call <- ifelse(as.character(Idents(merged)) %in% as.character(nk_clusters),
                                "NK cell (candidate)", "Other / not NK")

nk_feature_plot <- FeaturePlot(merged, features = "NK_score1") +
  labs(title = "Xiong GSE223138 - NK marker module score")
ggsave(file.path(qc_plot_dir, "umap_nk_score.png"), nk_feature_plot, width = 7, height = 6, dpi = 150)

umap_nk_call <- DimPlot(merged, reduction = "umap", group.by = "cell_type_call") +
  labs(title = "Xiong GSE223138 - NK-cell candidate calls")
ggsave(file.path(qc_plot_dir, "umap_nk_candidate_calls.png"), umap_nk_call, width = 7, height = 6, dpi = 150)

log_msg("Wrote plots/umap_nk_score.png and plots/umap_nk_candidate_calls.png.")

# NK cell counts by real disease status (first genuinely real NK-level
# number this project has produced)
nk_by_group <- merged@meta.data %>%
  filter(!is.na(pd_status)) %>%
  group_by(sample_id, pd_status) %>%
  summarise(
    n_cells_total = n(),
    n_nk_candidate = sum(cell_type_call == "NK cell (candidate)"),
    pct_nk = 100 * n_nk_candidate / n_cells_total,
    .groups = "drop"
  )

write_csv(nk_by_group, file.path(qc_dir, "nk_cell_counts_by_sample.csv"))

log_msg("\nNK-candidate cell counts per sample (results/qc/xiong_scrna/nk_cell_counts_by_sample.csv):")
print(as.data.frame(nk_by_group))

# ------------------------------------------------------------------------
# 7. Top marker genes per cluster (useful for manually auditing the NK
#    call, and for later DE work) - descriptive only, not hypothesis tested.
# ------------------------------------------------------------------------

log_msg("\n--- Top marker genes per cluster (for manual audit) ---")
all_markers <- tryCatch(
  FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE),
  error = function(e) {
    log_msg("FindAllMarkers failed: ", conditionMessage(e))
    NULL
  }
)
if (!is.null(all_markers)) {
  write_csv(all_markers, file.path(qc_dir, "cluster_marker_genes.csv"))
  log_msg("Wrote cluster_marker_genes.csv (", nrow(all_markers), " marker gene x cluster rows).")
}

# ------------------------------------------------------------------------
# 8. Save processed object and write report
# ------------------------------------------------------------------------

saveRDS(merged, file.path(validated_dir, "xiong_scrna_processed.rds"))
log_msg("\nWrote results/validated_data/xiong_scrna_processed.rds (QC-filtered, clustered, NK-annotated).")

writeLines(qc_log, file.path(qc_dir, "qc_report.txt"))

# ------------------------------------------------------------------------
# 9. Console summary
# ------------------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("XIONG GSE223138 - scRNA-seq PROCESSING SUMMARY\n")
cat("====================================================\n\n")
cat("Samples loaded:        ", length(seurat_list), " / 6\n")
cat("Cells before QC:       ", n_before, "\n")
cat("Cells after QC:        ", n_after, " (", sprintf("%.1f%%", 100 * n_after / n_before), " retained)\n", sep = "")
cat("Clusters found:        ", n_clusters, " (resolution = ", CLUSTER_RESOLUTION, ")\n", sep = "")
cat("NK-candidate clusters: ", if (length(nk_clusters) > 0) paste(nk_clusters, collapse = ", ") else "NONE", "\n")
cat("Total NK-candidate cells: ", sum(merged$cell_type_call == "NK cell (candidate)"), " / ", ncol(merged), "\n")
cat("\nOutputs:\n")
cat("  results/qc/xiong_scrna/qc_metrics_pre_filter.csv\n")
cat("  results/qc/xiong_scrna/qc_adaptive_thresholds_per_sample.csv\n")
cat("  results/qc/xiong_scrna/cluster_nk_tcell_scores.csv\n")
cat("  results/qc/xiong_scrna/nk_cell_counts_by_sample.csv\n")
cat("  results/qc/xiong_scrna/cluster_marker_genes.csv\n")
cat("  results/qc/xiong_scrna/plots/ (QC violins, UMAPs by cluster/sample/disease/NK score)\n")
cat("  results/qc/xiong_scrna/qc_report.txt\n")
cat("  results/validated_data/xiong_scrna_processed.rds\n")
cat("\nREVIEW BEFORE PROCEEDING: qc_violin_pre_filter.png, umap_by_sample.png (batch check),\n")
cat("and cluster_nk_tcell_scores.csv (confirm NK call makes biological sense).\n")
cat("====================================================\n")

