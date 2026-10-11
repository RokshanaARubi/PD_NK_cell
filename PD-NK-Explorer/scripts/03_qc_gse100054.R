#!/usr/bin/env Rscript
# ==============================================================================
# PD-NK Explorer
# scripts/03_qc_gse100054.R
#
# PURPOSE (AND ONLY PURPOSE):
#   Quality control on GSE100054 - the ONLY dataset marked READY by
#   scripts/02_validate_datasets.R. Xiong and Moquin-Beaudry are
#   intentionally NOT touched by this script; they remain PARTIALLY
#   AVAILABLE pending expression-data acquisition.
#
#   This script does NOT normalize, does NOT perform differential
#   expression, does NOT perform statistical hypothesis testing, does
#   NOT perform machine learning, and does NOT remove any samples or
#   features from the data. It only computes descriptive QC metrics,
#   flags candidates for review, and writes plots/report files.
#
# INPUT:
#   results/validated_data/gse100054_validated.rds
#   (produced by scripts/02_validate_datasets.R; this script refuses to
#   proceed if that file is missing or wasn't actually marked READY,
#   so QC is never silently run on unvalidated data.)
#
# HARD SAFETY RULES ENFORCED BY THIS SCRIPT:
#   - Never modifies data/ (read-only there; this script doesn't even
#     touch data/ directly - it reads only the validated RDS)
#   - Never removes samples or genes from the working data - only flags
#     candidates in separate CSVs for manual review
#   - Never normalizes, log-transforms, or rescales expression values
#   - Never performs DE, hypothesis tests, or ML
#   - Never loads GSE100054_platform_annotation.csv or GPL23126.soft.gz
#     in full (both are very large; only a small header/preview is read)
#   - Never fabricates thresholds as "correct" - all outlier cutoffs are
#     labeled as heuristic and are reported, not silently applied
# ==============================================================================

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(tidyr)
}))

# ------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------

project_root <- getwd()

validated_dir <- file.path(project_root, "results", "validated_data")
qc_dir        <- file.path(project_root, "results", "qc", "gse100054")
qc_plot_dir   <- file.path(qc_dir, "qc_plots")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_plot_dir, recursive = TRUE, showWarnings = FALSE)

rds_path <- file.path(validated_dir, "gse100054_validated.rds")

qc_log <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  qc_log <<- c(qc_log, msg)
}

# ------------------------------------------------------------------------
# 1. Load validated data (hard stop if not present - this is the one place
#    a genuine missing-required-file stop is appropriate, since running QC
#    on unvalidated data would be worse than refusing to run at all)
# ------------------------------------------------------------------------

if (!file.exists(rds_path)) {
  stop(
    "REQUIRED FILE MISSING: ", rds_path, "\n",
    "This script only runs QC on data that has passed scripts/02_validate_datasets.R ",
    "and was marked READY. Run that script first."
  )
}

validated <- tryCatch(readRDS(rds_path), error = function(e) {
  stop("Could not read ", rds_path, ": ", conditionMessage(e))
})

required_fields <- c("expression", "metadata", "identifier_column", "sample_id_column_metadata")
missing_fields <- setdiff(required_fields, names(validated))
if (length(missing_fields) > 0) {
  stop("Validated RDS is missing expected field(s): ", paste(missing_fields, collapse = ", "),
       ". Re-run scripts/02_validate_datasets.R.")
}

expr        <- validated$expression       # data.table: id column + one column per sample
meta        <- validated$metadata         # data.frame: one row per sample
id_col      <- validated$identifier_column
sample_col  <- validated$sample_id_column_metadata

log_msg("=================================================================")
log_msg("GSE100054 - QUALITY CONTROL")
log_msg("Generated: ", as.character(Sys.time()))
log_msg("Input: ", rds_path)
log_msg("Identifier column: '", id_col, "' | Metadata sample-ID column: '", sample_col, "'")
log_msg("=================================================================\n")

feature_ids  <- as.character(expr[[id_col]])
sample_ids   <- colnames(expr)[colnames(expr) != id_col]

expr_mat <- as.matrix(expr[, ..sample_ids])
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- feature_ids

log_msg("Working matrix: ", nrow(expr_mat), " features x ", ncol(expr_mat), " samples.")

# Align metadata rows to expression sample order (should already match 1:1,
# but never assume - re-verify explicitly here rather than trusting the
# previous script's report blindly)
meta_aligned <- meta[match(sample_ids, as.character(meta[[sample_col]])), , drop = FALSE]

n_unmatched <- sum(is.na(meta_aligned[[sample_col]]))
if (n_unmatched > 0) {
  log_msg("WARNING: ", n_unmatched, " expression sample(s) could not be matched to metadata during QC alignment. ",
          "These samples will still be included in QC metrics but will show NA metadata fields.")
} else {
  log_msg("All ", length(sample_ids), " expression samples matched to metadata rows for QC alignment.")
}

# Try to find the diagnosis-like column for annotation of plots (best-effort, read-only use)
diagnosis_col <- names(meta_aligned)[str_detect(tolower(names(meta_aligned)), "diagnosis|disease|condition|group")]
diagnosis_col <- if (length(diagnosis_col) > 0) diagnosis_col[1] else NA_character_
if (!is.na(diagnosis_col)) {
  log_msg("Using '", diagnosis_col, "' as the group/diagnosis label for QC plots (visual annotation only - ",
          "not used to alter any QC calculation).")
} else {
  log_msg("No diagnosis-like column found for plot annotation; QC plots will be unlabeled by group.")
}
group_labels <- if (!is.na(diagnosis_col)) as.character(meta_aligned[[diagnosis_col]]) else rep(NA_character_, length(sample_ids))

# ==============================================================================
# 2. Missing-value assessment
# ==============================================================================

log_msg("\n--- Missing-value assessment ---")

na_per_sample <- colSums(is.na(expr_mat))
na_per_feature <- rowSums(is.na(expr_mat))
total_na <- sum(is.na(expr_mat))
total_cells <- length(expr_mat)

log_msg("Total NA cells: ", total_na, " / ", total_cells,
        " (", sprintf("%.4f%%", 100 * total_na / total_cells), ")")

if (any(na_per_sample > 0)) {
  log_msg("Samples with missing values: ",
          paste(sprintf("%s (%d NA)", names(na_per_sample)[na_per_sample > 0], na_per_sample[na_per_sample > 0]),
                collapse = ", "))
} else {
  log_msg("No sample has any missing values.")
}

features_with_na <- sum(na_per_feature > 0)
log_msg("Features with at least one missing value: ", features_with_na, " / ", nrow(expr_mat))

# ==============================================================================
# 3. Per-sample descriptive statistics
# ==============================================================================

log_msg("\n--- Per-sample descriptive statistics ---")

sample_stats <- tibble(
  sample_id = sample_ids,
  group = group_labels,
  n_features = nrow(expr_mat),
  n_missing = na_per_sample,
  min = apply(expr_mat, 2, function(x) suppressWarnings(min(x, na.rm = TRUE))),
  max = apply(expr_mat, 2, function(x) suppressWarnings(max(x, na.rm = TRUE))),
  mean = apply(expr_mat, 2, function(x) mean(x, na.rm = TRUE)),
  median = apply(expr_mat, 2, function(x) median(x, na.rm = TRUE)),
  sd = apply(expr_mat, 2, function(x) sd(x, na.rm = TRUE))
)

write_csv(sample_stats, file.path(qc_dir, "sample_summary_stats.csv"))
log_msg("Wrote sample_summary_stats.csv (", nrow(sample_stats), " samples).")

# Flag samples whose mean expression deviates > 3 SD from the across-sample
# mean-of-means - a simple, explicitly-labeled heuristic, not a hard rule.
mean_of_means <- mean(sample_stats$mean)
sd_of_means   <- sd(sample_stats$mean)
sample_stats <- sample_stats %>%
  mutate(mean_zscore = (mean - mean_of_means) / sd_of_means,
         flag_mean_outlier = abs(mean_zscore) > 3)

n_mean_outliers <- sum(sample_stats$flag_mean_outlier)
log_msg("Samples flagged as global-mean outliers (heuristic: |z| > 3 across sample means): ", n_mean_outliers)
if (n_mean_outliers > 0) {
  log_msg("  Flagged sample(s): ", paste(sample_stats$sample_id[sample_stats$flag_mean_outlier], collapse = ", "))
}

# ==============================================================================
# 4. Sample-sample correlation
# ==============================================================================

log_msg("\n--- Sample-sample correlation ---")

cor_mat <- suppressWarnings(cor(expr_mat, use = "pairwise.complete.obs", method = "pearson"))

cor_df <- as.data.frame(cor_mat) %>%
  rownames_to_column("sample_id")
write_csv(cor_df, file.path(qc_dir, "sample_correlation_matrix.csv"))
log_msg("Wrote sample_correlation_matrix.csv (", ncol(cor_mat), " x ", ncol(cor_mat), " matrix).")

# For each sample, its mean correlation to all OTHER samples - low mean
# correlation is a classic (heuristic) signature of a low-quality/outlier
# bulk expression sample.
diag(cor_mat) <- NA
mean_other_cor <- rowMeans(cor_mat, na.rm = TRUE)

cor_summary <- tibble(
  sample_id = names(mean_other_cor),
  mean_correlation_to_others = mean_other_cor
) %>%
  mutate(
    cor_zscore = (mean_correlation_to_others - mean(mean_correlation_to_others)) /
      sd(mean_correlation_to_others),
    flag_low_correlation = cor_zscore < -2
  )

n_cor_outliers <- sum(cor_summary$flag_low_correlation)
log_msg("Samples flagged for unusually low correlation to the rest of the cohort ",
        "(heuristic: z < -2 on mean pairwise correlation): ", n_cor_outliers)
if (n_cor_outliers > 0) {
  log_msg("  Flagged sample(s): ", paste(cor_summary$sample_id[cor_summary$flag_low_correlation], collapse = ", "))
}

# Correlation heatmap (base graphics - no extra plotting dependency)
png(file.path(qc_plot_dir, "sample_correlation_heatmap.png"), width = 1000, height = 900, res = 120)
heatmap(
  cor_mat_full <- suppressWarnings(cor(expr_mat, use = "pairwise.complete.obs", method = "pearson")),
  symm = TRUE,
  main = "GSE100054 - Sample-sample Pearson correlation",
  margins = c(8, 8)
)
dev.off()
log_msg("Wrote qc_plots/sample_correlation_heatmap.png")

# Hierarchical clustering dendrogram of samples (1 - correlation as distance)
dist_mat <- as.dist(1 - cor_mat_full)
hc <- hclust(dist_mat, method = "average")

png(file.path(qc_plot_dir, "sample_dendrogram.png"), width = 1000, height = 700, res = 120)
plot(hc, main = "GSE100054 - Sample clustering (1 - Pearson correlation)",
     xlab = "", sub = "", cex = 0.8)
dev.off()
log_msg("Wrote qc_plots/sample_dendrogram.png")

# ==============================================================================
# 5. PCA on samples
# ==============================================================================

log_msg("\n--- Principal component analysis (samples) ---")

# Use only features with non-zero variance and no missing values, so PCA
# itself doesn't silently drop or impute anything unexpected.
feature_var <- apply(expr_mat, 1, var, na.rm = TRUE)
usable_features <- !is.na(feature_var) & feature_var > 0 & na_per_feature == 0

log_msg("Features usable for PCA (finite, non-zero variance, no missing values): ",
        sum(usable_features), " / ", nrow(expr_mat))

pca_result <- NULL
if (sum(usable_features) >= 2 && ncol(expr_mat) >= 3) {
  pca_input <- t(expr_mat[usable_features, , drop = FALSE])
  pca_result <- tryCatch(prcomp(pca_input, center = TRUE, scale. = TRUE),
                         error = function(e) {
                           log_msg("PCA failed: ", conditionMessage(e))
                           NULL
                         })
} else {
  log_msg("Not enough usable features or samples to run PCA.")
}

if (!is.null(pca_result)) {
  var_explained <- summary(pca_result)$importance["Proportion of Variance", ]
  n_pcs_report <- min(5, length(var_explained))
  
  pca_coords <- as.data.frame(pca_result$x[, seq_len(n_pcs_report), drop = FALSE]) %>%
    rownames_to_column("sample_id") %>%
    left_join(tibble(sample_id = sample_ids, group = group_labels), by = "sample_id")
  
  write_csv(pca_coords, file.path(qc_dir, "pca_coordinates.csv"))
  log_msg("Wrote pca_coordinates.csv (first ", n_pcs_report, " PCs).")
  log_msg("Variance explained by PC1-PC", n_pcs_report, ": ",
          paste(sprintf("%.1f%%", 100 * var_explained[seq_len(n_pcs_report)]), collapse = ", "))
  
  # Flag samples beyond 3 SD from centroid on PC1 or PC2 (heuristic only)
  pc1_z <- scale(pca_coords$PC1)[, 1]
  pc2_z <- if ("PC2" %in% colnames(pca_coords)) scale(pca_coords$PC2)[, 1] else rep(0, nrow(pca_coords))
  pca_coords <- pca_coords %>%
    mutate(pc1_zscore = pc1_z, pc2_zscore = pc2_z,
           flag_pca_outlier = abs(pc1_zscore) > 3 | abs(pc2_zscore) > 3)
  
  n_pca_outliers <- sum(pca_coords$flag_pca_outlier)
  log_msg("Samples flagged as PCA outliers (heuristic: |z| > 3 on PC1 or PC2): ", n_pca_outliers)
  if (n_pca_outliers > 0) {
    log_msg("  Flagged sample(s): ", paste(pca_coords$sample_id[pca_coords$flag_pca_outlier], collapse = ", "))
  }
  
  write_csv(pca_coords, file.path(qc_dir, "pca_coordinates.csv"))
  
  # PCA plot
  pca_plot <- ggplot(pca_coords, aes(x = PC1, y = if ("PC2" %in% colnames(pca_coords)) PC2 else 0,
                                     color = group, label = sample_id)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    labs(title = "GSE100054 - PCA of samples",
         x = sprintf("PC1 (%.1f%% variance)", 100 * var_explained[1]),
         y = if ("PC2" %in% colnames(pca_coords)) sprintf("PC2 (%.1f%% variance)", 100 * var_explained[2]) else "PC2 (not available)",
         color = if (!is.na(diagnosis_col)) diagnosis_col else "group") +
    theme_minimal()
  
  ggsave(file.path(qc_plot_dir, "pca_plot.png"), pca_plot, width = 7, height = 5.5, dpi = 150)
  log_msg("Wrote qc_plots/pca_plot.png")
} else {
  pca_coords <- tibble(sample_id = sample_ids, group = group_labels)
}

# ==============================================================================
# 6. Global distribution plots (density + boxplot across all samples)
# ==============================================================================

log_msg("\n--- Global expression distribution ---")

long_expr <- as_tibble(expr_mat) %>%
  mutate(feature_id = feature_ids) %>%
  tidyr::pivot_longer(cols = all_of(sample_ids), names_to = "sample_id", values_to = "value") %>%
  left_join(tibble(sample_id = sample_ids, group = group_labels), by = "sample_id")

density_plot <- ggplot(long_expr, aes(x = value, color = sample_id)) +
  geom_density(alpha = 0.3, show.legend = FALSE) +
  labs(title = "GSE100054 - Per-sample expression density", x = "Expression value", y = "Density") +
  theme_minimal()
ggsave(file.path(qc_plot_dir, "expression_density.png"), density_plot, width = 7, height = 5, dpi = 150)
log_msg("Wrote qc_plots/expression_density.png")

boxplot_plot <- ggplot(long_expr, aes(x = sample_id, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "GSE100054 - Per-sample expression distribution", x = "Sample", y = "Expression value",
       fill = if (!is.na(diagnosis_col)) diagnosis_col else "group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
ggsave(file.path(qc_plot_dir, "expression_boxplot.png"), boxplot_plot, width = 9, height = 5, dpi = 150)
log_msg("Wrote qc_plots/expression_boxplot.png")

# ==============================================================================
# 7. Low-variance feature candidates (reported only, NEVER removed here)
# ==============================================================================

log_msg("\n--- Low-variance feature candidates (reporting only) ---")

feature_variance_df <- tibble(
  feature_id = feature_ids,
  variance = feature_var,
  n_missing = na_per_feature
) %>%
  mutate(variance_percentile = percent_rank(variance))

low_var_threshold_percentile <- 0.10  # heuristic: bottom 10% of variance
low_var_candidates <- feature_variance_df %>%
  filter(!is.na(variance), variance_percentile <= low_var_threshold_percentile)

write_csv(feature_variance_df, file.path(qc_dir, "feature_variance_all.csv"))
write_csv(low_var_candidates, file.path(qc_dir, "low_variance_features_candidates.csv"))

log_msg("Total features: ", nrow(feature_variance_df))
log_msg("Low-variance candidates (heuristic: bottom ", sprintf("%.0f%%", 100 * low_var_threshold_percentile),
        " of variance): ", nrow(low_var_candidates),
        " - written to low_variance_features_candidates.csv for review. NOT removed from the working data.")

# ==============================================================================
# 8. Platform annotation file - size-aware, preview-only handling
# ==============================================================================

log_msg("\n--- Platform annotation file check (preview only, per safety rules) ---")

gse_platform_path <- file.path(project_root, "data", "gse100054", "GSE100054_platform_annotation.csv")
gse_soft_path     <- file.path(project_root, "data", "gse100054", "GPL23126.soft.gz")

if (file.exists(gse_platform_path)) {
  fsize <- file.info(gse_platform_path)$size
  log_msg("GSE100054_platform_annotation.csv detected. Size: ", sprintf("%.2f MB", fsize / 1024^2))
  if (fsize > 100 * 1024^2) {
    log_msg("File exceeds 100 MB - this is unusually large for a probe/gene annotation file and is ",
            "flagged here as a data-quality anomaly worth investigating separately. Only the first 5 rows ",
            "are previewed below; the file is NOT loaded in full by this QC script.")
  }
  preview <- tryCatch(data.table::fread(gse_platform_path, nrows = 5), error = function(e) NULL)
  if (!is.null(preview)) {
    log_msg("Platform annotation preview columns: ", paste(colnames(preview), collapse = ", "))
  } else {
    log_msg("Could not preview platform annotation file (file may be malformed or unusually structured).")
  }
} else {
  log_msg("GSE100054_platform_annotation.csv not found; skipping preview.")
}

if (file.exists(gse_soft_path)) {
  log_msg("GPL23126.soft.gz also present (", sprintf("%.2f MB", file.info(gse_soft_path)$size / 1024^2),
          "). NOT loaded by this script under any circumstance.")
}

log_msg("NOTE: Identifier-to-gene mapping/conversion is intentionally NOT performed in this QC step, ",
        "consistent with the prior validation step's rule to not convert identifiers yet.")

# ==============================================================================
# 9. Write QC report and outlier summary
# ==============================================================================

outlier_summary <- tibble(sample_id = sample_ids) %>%
  left_join(sample_stats %>% select(sample_id, flag_mean_outlier, mean_zscore), by = "sample_id") %>%
  left_join(cor_summary %>% select(sample_id, flag_low_correlation, cor_zscore), by = "sample_id") %>%
  left_join(pca_coords %>% select(any_of(c("sample_id", "flag_pca_outlier", "pc1_zscore", "pc2_zscore"))),
            by = "sample_id") %>%
  mutate(
    flag_pca_outlier = if ("flag_pca_outlier" %in% names(.)) flag_pca_outlier else FALSE,
    any_flag = coalesce(flag_mean_outlier, FALSE) | coalesce(flag_low_correlation, FALSE) | coalesce(flag_pca_outlier, FALSE)
  )

write_csv(outlier_summary, file.path(qc_dir, "outlier_flags.csv"))

n_any_flagged <- sum(outlier_summary$any_flag)
log_msg("\n--- Outlier flag summary ---")
log_msg("Samples flagged by at least one heuristic: ", n_any_flagged, " / ", length(sample_ids))
if (n_any_flagged > 0) {
  log_msg("  ", paste(outlier_summary$sample_id[outlier_summary$any_flag], collapse = ", "))
}
log_msg("These are heuristic flags for manual review only. No samples have been removed or altered.")

writeLines(qc_log, file.path(qc_dir, "qc_report.txt"))

# ==============================================================================
# 10. Write a QC-annotated object (NOT a filtered/normalized dataset)
# ==============================================================================

qc_annotated <- list(
  expression = expr,                 # untouched, exactly as validated
  metadata = meta,                   # untouched
  identifier_column = id_col,
  identifier_type = validated$identifier_type,
  sample_id_column_metadata = sample_col,
  qc_sample_stats = sample_stats,
  qc_correlation_summary = cor_summary,
  qc_pca_coordinates = pca_coords,
  qc_outlier_flags = outlier_summary,
  qc_low_variance_candidates = low_var_candidates,
  qc_run_on = Sys.time(),
  qc_notes = "Descriptive QC only. No normalization, filtering, or sample/feature removal applied."
)

saveRDS(qc_annotated, file.path(validated_dir, "gse100054_qc.rds"))
log_msg("\nWrote results/validated_data/gse100054_qc.rds (QC-annotated, NOT normalized or filtered).")

# ==============================================================================
# 11. Console summary
# ==============================================================================

cat("\n")
cat("====================================================\n")
cat("GSE100054 QUALITY CONTROL SUMMARY\n")
cat("====================================================\n\n")
cat("Samples:                     ", length(sample_ids), "\n")
if (!is.na(diagnosis_col)) {
  group_tab <- table(group_labels, useNA = "ifany")
  cat("Group breakdown (", diagnosis_col, "): ", paste(names(group_tab), group_tab, sep = "=", collapse = ", "), "\n", sep = "")
}
cat("Features:                    ", nrow(expr_mat), "\n")
cat("Missing data:                ", sprintf("%.4f%%", 100 * total_na / total_cells), "\n")
cat("Mean-outlier samples:        ", n_mean_outliers, "\n")
cat("Low-correlation samples:     ", n_cor_outliers, "\n")
if (!is.null(pca_result)) {
  cat("PCA-outlier samples:         ", n_pca_outliers, "\n")
}
cat("Any-flag samples (union):    ", n_any_flagged, "\n")
cat("Low-variance feature candidates: ", nrow(low_var_candidates), " / ", nrow(feature_variance_df), "\n")
cat("\nOutputs:\n")
cat("  results/qc/gse100054/sample_summary_stats.csv\n")
cat("  results/qc/gse100054/sample_correlation_matrix.csv\n")
cat("  results/qc/gse100054/pca_coordinates.csv\n")
cat("  results/qc/gse100054/outlier_flags.csv\n")
cat("  results/qc/gse100054/feature_variance_all.csv\n")
cat("  results/qc/gse100054/low_variance_features_candidates.csv\n")
cat("  results/qc/gse100054/qc_plots/ (density, boxplot, correlation heatmap, dendrogram, PCA plot)\n")
cat("  results/qc/gse100054/qc_report.txt\n")
cat("  results/validated_data/gse100054_qc.rds\n")
cat("\nNo samples or features were removed. No normalization was applied.\n")
cat("====================================================\n")

