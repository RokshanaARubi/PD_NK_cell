

pkgs <- c("dplyr", "readr", "stringr", "tibble", "tidyr")
miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(miss) > 0) stop("Missing packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tibble); library(tidyr)
})

root <- getwd()
P <- function(...) file.path(root, ...)
out_dir <- P("results", "final"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need <- c(P("results", "de", "xiong_nk", "de_primary_PD_vs_Control.csv"),
          P("results", "validated_data", "gse100054_validation.rds"))
for (f in need) if (!file.exists(f)) stop("REQUIRED FILE MISSING: ", f)
rd <- function(f) if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL
rr <- function(f) if (file.exists(f)) readRDS(f) else NULL

de_xiong   <- rd(P("results", "de", "xiong_nk", "de_primary_PD_vs_Control.csv"))
nk_counts  <- rd(P("results", "qc", "xiong_scrna", "nk_cell_counts_by_sample_v2.csv"))
val        <- rr(P("results", "validated_data", "gse100054_validation.rds"))
fun        <- rr(P("results", "validated_data", "functional_analysis.rds"))
roc        <- rr(P("results", "validated_data", "roc_analysis.rds"))

# ------------------------------------------------------------------------
# 1. Final per-gene verdicts
# ------------------------------------------------------------------------
cv <- bind_rows(lapply(names(val$candidate_validation), function(s) val$candidate_validation[[s]]))
cv <- cv %>% distinct(gene, .keep_all = TRUE)    # tier2 first, then wide-set genes
verd <- cv %>%
  mutate(final_status = case_when(
    verdict == "Not measurable on platform" ~ "Not measurable in bulk",
    FDR_xiong < 0.05 & concordant & FDR_candset < 0.05 & robust_to_NK_adjustment %in% TRUE ~
      "REPRODUCIBLE NK BIOMARKER",
    concordant & P_xiong < 0.01 & P_bulk < 0.05 ~ "Suggestive (nominal in both datasets, not FDR-significant)",
    str_detect(verdict, "^CONTRADICTED") ~ "Contradicted by bulk data",
    TRUE ~ "Not replicated"),
    robust_to_NK_adjustment = robust_to_NK_adjustment %in% TRUE) %>%
  select(gene, set, tier, logFC_xiong, P_xiong, FDR_xiong, logFC_bulk, P_bulk, FDR_candset,
         logFC_bulk_NKadj, P_bulk_NKadj, concordant, robust_to_NK_adjustment, verdict, final_status) %>%
  arrange(factor(final_status, levels = c("REPRODUCIBLE NK BIOMARKER",
                                          "Suggestive (nominal in both datasets, not FDR-significant)",
                                          "Not replicated", "Contradicted by bulk data",
                                          "Not measurable in bulk")), P_xiong)
write_csv(verd, file.path(out_dir, "biomarker_verdicts.csv"))
st_tab <- table(verd$final_status)

# ------------------------------------------------------------------------
# 2. Key numbers
# ------------------------------------------------------------------------
n_nk <- if (!is.null(nk_counts)) sum(nk_counts$n_NK) else NA
n_cells <- if (!is.null(nk_counts)) sum(nk_counts$n_cells) else NA
glob <- val$set_level_tests
roc_tbl <- roc$roc_summary
best <- if (!is.null(roc_tbl)) roc_tbl %>% arrange(desc(AUC)) %>% slice(1) else NULL
cam <- roc$camera_robustness
cam_sig <- if (!is.null(cam)) cam %>% filter(padj_xiong < 0.05) else NULL
cam_ok <- if (!is.null(cam_sig)) sum(cam_sig$P_camera_bulk < 0.05 & cam_sig$same_direction) else NA
prox <- roc$composition_proxies

summ <- tribble(
  ~step, ~metric, ~value,
  "2-4", "Cells after QC (Xiong)", format(n_cells, big.mark = ","),
  "4",   "NK cells called (clusters 12, 3, 5)", format(n_nk, big.mark = ","),
  "5",   "Donors (Control / PD)", "2 / 4",
  "5",   "Genes tested (NK pseudobulk)", format(nrow(de_xiong), big.mark = ","),
  "5",   "Genes with FDR < 0.05 (PD vs Control)", as.character(sum(de_xiong$FDR < 0.05)),
  "6",   "Candidates: Tier 1 / Tier 2 (exploratory)",
  paste0(sum(str_detect(verd$tier, "^Tier 1"), na.rm = TRUE), " / ", sum(str_detect(verd$tier, "^Tier 2"), na.rm = TRUE)),
  "7",   "Tier-2 genes measurable in GSE100054", as.character(sum(val$candidate_validation$tier2_primary$verdict != "Not measurable on platform")),
  "7",   "Tier-2 direction concordance", paste0(glob$pct_concordant[glob$set == "tier2_primary"], "%  (binomial P = ", signif(glob$binom_P_vs_50pct[glob$set == "tier2_primary"], 3), ")"),
  "7",   "Genes with final status REPRODUCIBLE", as.character(sum(verd$final_status == "REPRODUCIBLE NK BIOMARKER")),
  "7",   "Genes 'Suggestive' (nominal both)", as.character(sum(str_detect(verd$final_status, "^Suggestive"))),
  "8",   "Xiong pathways padj<0.05 (fgsea, anti-conservative)", as.character(if (!is.null(fun)) sum(fun$gsea_xiong_primary$padj < 0.05, na.rm = TRUE) else NA),
  "8",   "Xiong Hallmark/KEGG sets replicating in bulk (camera, same direction, P<0.05)", as.character(cam_ok),
  "9",   "Best signature AUC", if (!is.null(best)) paste0(best$signature, ": ", best$AUC, " [", best$CI_low, "-", best$CI_high, "], perm P = ", best$perm_P_one_sided) else "NA",
  "9",   "Signatures with permutation P < 0.05", if (!is.null(roc_tbl)) paste0(sum(roc_tbl$perm_P_one_sided < 0.05, na.rm = TRUE), " of ", nrow(roc_tbl)) else "NA",
  "10",  "Overall conclusion", if (sum(verd$final_status == "REPRODUCIBLE NK BIOMARKER") > 0)
    "At least one gene passed all reproducibility criteria (see biomarker_verdicts.csv)" else
      "NOT REPLICATED: no NK candidate passed the reproducibility criteria")
write_csv(summ, file.path(out_dir, "final_summary.csv"))

# ------------------------------------------------------------------------
# 3. Plain-language report
# ------------------------------------------------------------------------
sug <- verd %>% filter(str_detect(final_status, "^Suggestive")) %>% pull(gene)
comp_txt <- if (!is.null(prox)) paste(sprintf("%s (P=%.3g)", prox$cell_type_proxy, prox$wilcox_P), collapse = "; ") else "n/a"
lines <- c(
  "# PD-NK Explorer - Final report",
  "",
  paste0("Generated: ", Sys.time()),
  "",
  "**Research question.** Can NK-cell-specific transcriptomic signatures, identified from single-cell RNA-seq and validated across independent bulk RNA-seq datasets, serve as reproducible peripheral biomarkers for Parkinson's disease?",
  "",
  "## Data actually used",
  "- Discovery: Xiong et al. GSE223138, PBMC scRNA-seq, 6 donors (2 Control, 2 Early PD, 2 Late PD).",
  "- Independent validation: GSE100054, bulk PBMC microarray (GPL23126), 9 Control vs 10 PD.",
  "- NOT used: Moquin-Beaudry PRJNA1145007 (only SRA run metadata available; no processed matrix), PPMI (controlled access).",
  "",
  "## Results",
  sprintf("1. After QC %s cells were retained; %s were called NK (clusters 12, 3, 5).",
          format(n_cells, big.mark = ","), format(n_nk, big.mark = ",")),
  sprintf("2. Donor-level pseudobulk edgeR (2 Control vs 4 PD): %d of %s genes reached FDR < 0.05. All candidates are therefore exploratory (Tier 2: nominal P < 0.01 and |log2FC| >= 1).",
          sum(de_xiong$FDR < 0.05), format(nrow(de_xiong), big.mark = ",")),
  sprintf("3. In GSE100054 the Tier-2 candidates showed %s direction concordance (chance = 50%%), binomial P = %s.",
          paste0(glob$pct_concordant[glob$set == "tier2_primary"], "%"),
          signif(glob$binom_P_vs_50pct[glob$set == "tier2_primary"], 3)),
  sprintf("4. Final gene status: %s.", paste(sprintf("%s = %d", names(st_tab), as.integer(st_tab)), collapse = "; ")),
  if (length(sug)) sprintf("   Suggestive genes (nominal in both datasets, not FDR-significant, to be treated as hypotheses only): %s.", paste(sug, collapse = ", ")) else NULL,
  sprintf("5. Pathway analysis: fgsea reported many significant pathways, but gene-permutation GSEA is anti-conservative for correlated sets. A correlation-aware camera test in the bulk data confirmed %s of the Xiong-significant Hallmark/KEGG sets.", cam_ok),
  if (!is.null(best)) sprintf("6. ROC (direction fixed a priori, signatures defined from Xiong only): best AUC = %s (95%% CI %s-%s), permutation P = %s; %s of %s signatures had P < 0.05. With n = 19 these AUCs are indistinguishable from chance.",
                              best$AUC, best$CI_low, best$CI_high, best$perm_P_one_sided,
                              sum(roc_tbl$perm_P_one_sided < 0.05, na.rm = TRUE), nrow(roc_tbl)) else NULL,
  paste0("7. Bulk PBMC cell-composition proxies, PD vs Control: ", comp_txt, ". Composition differences, not NK-intrinsic expression, are the most visible blood signal."),
  "",
  "## Conclusion",
  if (sum(verd$final_status == "REPRODUCIBLE NK BIOMARKER") == 0)
    "In the available public data, NK-cell transcriptomic signatures from Xiong GSE223138 did NOT replicate in GSE100054 and cannot be called reproducible peripheral biomarkers for Parkinson's disease. This is a negative / inconclusive result, not evidence that no NK signal exists."
  else "Some genes passed all criteria; see biomarker_verdicts.csv and treat them as candidates for prospective validation.",
  "",
  "## Limitations",
  "- Discovery had only 2 Control and 4 PD donors (low power); stage contrasts used 2 donors per group.",
  "- Validation set is a different technology (microarray) and bulk (mixed cell types), so NK-specific changes are diluted.",
  "- Probe-to-gene mapping used gene_assignment symbols; lncRNA-type genes (AC/AL/LINC) are often absent from the array.",
  "- No replication cohort at single-cell level (Moquin-Beaudry raw data not processed). No PPMI.",
  "- Batch/donor effects in the scRNA-seq object were not formally integrated.",
  "",
  "## Next steps that could change the conclusion",
  "- Process Moquin-Beaudry (PRJNA1145007) from FASTQ and repeat steps 3-7 as a true scRNA-seq replication.",
  "- Larger single-cell or sorted-NK cohorts; PPMI blood RNA-seq once access is granted.")
writeLines(lines[!vapply(lines, is.null, logical(1))], file.path(out_dir, "final_report.md"))

# ------------------------------------------------------------------------
# 4. UI bundle
# ------------------------------------------------------------------------
ui_bundle <- list(
  summary = summ, verdicts = verd, nk_counts_by_sample = nk_counts,
  xiong_de = de_xiong,
  stage_de = lapply(setNames(c("Early_vs_Control", "Late_vs_Control", "Late_vs_Early"),
                             c("Early_vs_Control", "Late_vs_Control", "Late_vs_Early")),
                    function(n) rd(P("results", "de", "xiong_nk", paste0("de_stage_", n, ".csv")))),
  bulk_de = val$bulk_de_modelA, set_level_tests = glob,
  gsea_xiong = if (!is.null(fun)) fun$gsea_xiong_primary else NULL,
  pathway_replication = if (!is.null(fun)) fun$pathway_replication else NULL,
  camera_robustness = cam,
  roc_summary = roc_tbl, signature_scores = roc$signature_scores,
  composition_proxies = prox,
  bulk_group = val$group, generated = Sys.time())
saveRDS(ui_bundle, P("results", "validated_data", "ui_results.rds"))

cat("\n", strrep("=", 70), "\nSTEP 10 - FINAL SUMMARY\n", strrep("=", 70), "\n", sep = "")
print(as.data.frame(summ), right = FALSE)
cat("\nFinal gene status counts:\n"); print(st_tab)
cat("\nWrote: results/final/final_summary.csv, biomarker_verdicts.csv, final_report.md\n")
cat("       results/validated_data/ui_results.rds   (used by the Shiny app)\n")
cat(strrep("=", 70), "\n")

