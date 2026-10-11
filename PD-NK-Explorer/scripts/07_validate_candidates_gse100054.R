#!/usr/bin/env Rscript
# ==============================================================================
# PD-NK Explorer
# scripts/07_validate_candidates_gse100054.R
#
# METHODOLOGY STEP 7: Independent bulk dataset validation
#
# QUESTION:
#   Do the NK-cell candidate genes found in Xiong scRNA-seq (script 06) show
#   the SAME direction (and ideally significance) of change in PD vs Control in
#   an independent bulk PBMC dataset (GSE100054)?
#
# WHAT IT DOES
#   1. Loads Xiong candidates (Tier 2 = primary set; nominal P<0.01 = wider set)
#   2. Maps GSE100054 probes -> gene symbols (auto-detects the annotation
#      columns; collapses multi-probe genes by highest mean expression)
#   3. Reports which candidates are even MEASURABLE on this platform
#   4. limma (trend + robust) PD vs Control on the gene-level matrix
#        Model A: ~ group
#        Model B: ~ group + NK_proxy   (adjusts for NK-cell abundance in PBMC)
#   5. Per-candidate verdict: direction concordance + replication P-values
#   6. Set-level replication tests (binomial sign test, limma geneSetTest)
#   7. Saves everything the later scripts (08 functional, 09 ROC) and the
#      Shiny "Cross-Dataset Validation" tab need.
#
# DOES NOT: touch data/, re-run Xiong analysis, do pathway analysis or ROC.
#
# READS
#   results/de/xiong_nk/nk_candidate_biomarkers.csv
#   results/de/xiong_nk/de_primary_PD_vs_Control.csv
#   results/validated_data/gse100054_validated.rds       (from script 02)
#   data/gse100054/GSE100054_platform_annotation.csv     (from script 01)
# WRITES
#   results/validation/gse100054/*.csv, *.txt, plots/*.png
#   results/validated_data/gse100054_validation.rds      (for Shiny + 08/09)
# ==============================================================================

pkgs <- c("data.table", "dplyr", "stringr", "tibble", "readr", "ggplot2", "tidyr", "limma")
miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(miss) > 0) stop("Missing packages: ", paste(miss, collapse = ", "),
                           "\n(limma is Bioconductor: BiocManager::install('limma'))")
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(tibble)
  library(readr); library(ggplot2); library(tidyr); library(limma)
})
set.seed(42)

# ---- USER-ADJUSTABLE SETTINGS ------------------------------------------------
GROUP_COL_OVERRIDE <- NULL   # e.g. "disease state:ch1" if auto-detection picks the wrong column
P_NOMINAL          <- 0.05   # nominal replication threshold in the bulk data
FDR_CANDIDATE      <- 0.05   # BH FDR across the tested candidate set
NK_PROXY_GENES     <- c("NKG7", "GNLY", "KLRD1", "NCAM1", "FCGR3A", "PRF1", "GZMB", "KLRF1", "NCR1")
# ------------------------------------------------------------------------------

root     <- getwd()
de_dir   <- file.path(root, "results", "de", "xiong_nk")
val_rds  <- file.path(root, "results", "validated_data", "gse100054_validated.rds")
annot_f  <- file.path(root, "data", "gse100054", "GSE100054_platform_annotation.csv")
out_dir  <- file.path(root, "results", "validation", "gse100054")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

LOG <- character(0)
say <- function(...) { m <- paste0(...); cat(m, "\n"); LOG <<- c(LOG, m) }

for (f in c(file.path(de_dir, "nk_candidate_biomarkers.csv"),
            file.path(de_dir, "de_primary_PD_vs_Control.csv"), val_rds, annot_f)) {
  if (!file.exists(f)) stop("REQUIRED FILE MISSING: ", f)
}

say(strrep("=", 70))
say("STEP 7 - VALIDATION OF XIONG NK CANDIDATES IN GSE100054 (bulk PBMC)")
say("Generated: ", as.character(Sys.time()))
say(strrep("=", 70))

# ------------------------------------------------------------------------
# 1. Xiong candidates
# ------------------------------------------------------------------------
xi_all  <- read_csv(file.path(de_dir, "de_primary_PD_vs_Control.csv"), show_col_types = FALSE)
xi_tier <- read_csv(file.path(de_dir, "nk_candidate_biomarkers.csv"),  show_col_types = FALSE)

xi_all$gene_uc  <- toupper(xi_all$gene)
sets <- list(
  tier2_primary = xi_tier %>% transmute(gene = toupper(gene), logFC_xiong = logFC,
                                        P_xiong = PValue, FDR_xiong = FDR, tier),
  p01_wide      = xi_all %>% filter(PValue < 0.01) %>%
    transmute(gene = gene_uc, logFC_xiong = logFC, P_xiong = PValue,
              FDR_xiong = FDR, tier = "nominal P<0.01")
)
say("\nCandidate sets from Xiong NK pseudobulk:")
for (s in names(sets)) say("  ", s, ": ", nrow(sets[[s]]), " genes (",
                           sum(sets[[s]]$logFC_xiong > 0), " up / ",
                           sum(sets[[s]]$logFC_xiong < 0), " down in PD)")
say("NOTE: no Xiong gene passed FDR<0.05, so ALL candidates are exploratory. ",
    "This step is the first real test of whether any of them are real.")

# ------------------------------------------------------------------------
# 2. Load GSE100054 + align metadata + define groups
# ------------------------------------------------------------------------
v <- readRDS(val_rds)
expr <- v$expression; meta <- as.data.frame(v$metadata)
id_col <- v$identifier_column; sample_col <- v$sample_id_column_metadata

sample_ids <- setdiff(colnames(expr), id_col)
probe_ids  <- as.character(expr[[id_col]])
X <- as.matrix(expr[, ..sample_ids]); storage.mode(X) <- "numeric"; rownames(X) <- probe_ids
meta <- meta[match(sample_ids, as.character(meta[[sample_col]])), , drop = FALSE]
if (anyNA(meta[[sample_col]])) stop("Some expression samples have no metadata row.")
say("\nGSE100054: ", nrow(X), " probes x ", ncol(X), " samples")

# platform info (printed so you can verify - see note in the chat message)
plat_cols <- grep("platform", names(meta), ignore.case = TRUE, value = TRUE)
for (pc in plat_cols) say("  metadata '", pc, "': ", paste(unique(meta[[pc]]), collapse = " | "))

# auto-detect the group column by VALUES (not just names)
detect_group <- function(meta) {
  best <- NULL; best_score <- 0
  for (cn in names(meta)) {
    vals <- tolower(as.character(meta[[cn]]))
    is_pd  <- str_detect(vals, "parkinson|\\bpd\\b")
    is_ctl <- str_detect(vals, "control|healthy|\\bhc\\b|normal") & !is_pd
    sc <- min(sum(is_pd), sum(is_ctl))
    if (sc > best_score) { best <- cn; best_score <- sc }
  }
  best
}
gcol <- if (!is.null(GROUP_COL_OVERRIDE)) GROUP_COL_OVERRIDE else detect_group(meta)
if (is.null(gcol)) stop("Could not find a PD/Control column. Set GROUP_COL_OVERRIDE. Columns: ",
                        paste(names(meta), collapse = " | "))
vals <- tolower(as.character(meta[[gcol]]))
grp <- case_when(str_detect(vals, "parkinson|\\bpd\\b") ~ "PD",
                 str_detect(vals, "control|healthy|\\bhc\\b|normal") ~ "Control",
                 TRUE ~ NA_character_)
say("Group column used: '", gcol, "'  (values: ", paste(unique(meta[[gcol]]), collapse = " | "), ")")
if (anyNA(grp)) say("WARNING: ", sum(is.na(grp)), " sample(s) not classifiable and will be DROPPED.")
keep <- !is.na(grp)
X <- X[, keep, drop = FALSE]; meta <- meta[keep, , drop = FALSE]
group <- factor(grp[keep], levels = c("Control", "PD"))
say("Samples used: ", paste(names(table(group)), table(group), sep = "=", collapse = ", "))
if (min(table(group)) < 3) stop("Fewer than 3 samples in a group - cannot validate.")

# ------------------------------------------------------------------------
# 3. Scale check (script 03 never normalised / transformed anything)
# ------------------------------------------------------------------------
rng <- range(X, na.rm = TRUE)
say("\nExpression range: ", round(rng[1], 2), " to ", round(rng[2], 2),
    " | median sample median: ", round(median(apply(X, 2, median, na.rm = TRUE)), 2))
if (rng[2] > 50) {
  say("Values look NON-log (max > 50) -> applying log2(x + 1).")
  X <- log2(pmax(X, 0) + 1)
} else {
  say("Values look already log2-scale -> used as is (verify this against the GEO processing description).")
}
sample_med <- apply(X, 2, median, na.rm = TRUE)
say("Between-sample spread of medians (SD): ", round(sd(sample_med), 3),
    "  (large values would suggest normalisation problems)")

# ------------------------------------------------------------------------
# 4. Probe -> gene symbol mapping (auto-detected, defensive)
# ------------------------------------------------------------------------
say("\n--- Probe -> gene mapping ---")
hdr <- names(fread(annot_f, nrows = 0))
say("Annotation columns (", length(hdr), "): ", paste(head(hdr, 25), collapse = ", "),
    if (length(hdr) > 25) ", ..." else "")

id_cands  <- hdr[tolower(hdr) %in% tolower(c("ID", "probeset_id", "transcript_cluster_id", "ID_REF",
                                             "Probe_Id", "ProbeID", "ILMN_Probe_ID", "Probe_ID"))]
sym_cands <- hdr[str_detect(tolower(hdr), "^gene_assignment$|symbol|^gene$|gene_name|ilmn_gene")]
if (length(id_cands) == 0 || length(sym_cands) == 0)
  stop("Could not find ID/symbol columns in annotation. Columns are listed above; ",
       "edit id_cands / sym_cands in this script.")
ann <- fread(annot_f, select = unique(c(id_cands, sym_cands)), showProgress = FALSE)

# choose the ID column that actually matches the expression row names
ov <- sapply(id_cands, function(cn) length(intersect(as.character(ann[[cn]]), rownames(X))))
id_use <- id_cands[which.max(ov)]
say("ID column used: '", id_use, "' (", max(ov), " / ", nrow(X), " expression probes matched)")
if (max(ov) == 0) stop("No annotation ID column overlaps expression probe IDs.")

# choose symbol column (prefer an explicit symbol column, else gene_assignment)
sym_use <- if (any(str_detect(tolower(sym_cands), "symbol"))) sym_cands[str_detect(tolower(sym_cands), "symbol")][1] else sym_cands[1]
say("Symbol column used: '", sym_use, "'")

raw_sym <- as.character(ann[[sym_use]])
parse_syms <- function(x, is_assignment) {
  lapply(strsplit(x, "\\s*///\\s*|;|\\|"), function(ent) {
    if (is_assignment) ent <- vapply(strsplit(ent, "\\s*//\\s*"), function(p) if (length(p) >= 2) p[2] else NA_character_, "")
    ent <- toupper(trimws(ent))
    unique(ent[!is.na(ent) & nzchar(ent) & ent != "---" & ent != "NA"])
  })
}
syms <- parse_syms(raw_sym, is_assignment = (tolower(sym_use) == "gene_assignment"))
map <- data.table(probe = rep(as.character(ann[[id_use]]), lengths(syms)), symbol = unlist(syms))
map <- unique(map[probe %in% rownames(X)])
say("Probes with >=1 gene symbol: ", uniqueN(map$probe), " / ", nrow(X),
    " | distinct symbols: ", uniqueN(map$symbol))

# collapse to gene level: keep the probe with the highest mean (phenotype-blind)
pm <- rowMeans(X, na.rm = TRUE)
map[, mean_expr := pm[probe]]
setorder(map, symbol, -mean_expr)
best <- map[, .SD[1], by = symbol]
G <- X[best$probe, , drop = FALSE]; rownames(G) <- best$symbol
say("Gene-level matrix: ", nrow(G), " genes x ", ncol(G), " samples")

# ------------------------------------------------------------------------
# 5. Candidate coverage (what is even measurable here?)
# ------------------------------------------------------------------------
cov_list <- lapply(names(sets), function(s) {
  d <- sets[[s]]; d$set <- s; d$on_platform <- d$gene %in% rownames(G); d
})
coverage <- bind_rows(cov_list)
write_csv(coverage, file.path(out_dir, "candidate_coverage.csv"))
say("\n--- Candidate coverage on GSE100054 ---")
for (s in names(sets)) {
  d <- coverage %>% filter(set == s)
  say("  ", s, ": ", sum(d$on_platform), " / ", nrow(d), " measurable")
  miss_g <- d$gene[!d$on_platform]
  if (length(miss_g)) say("    not measurable: ", paste(head(miss_g, 40), collapse = ", "),
                          if (length(miss_g) > 40) ", ..." else "")
}
say("Genes like AC#####.#, LINC#####, etc. are often missing on arrays, which shrinks the testable set.")

# ------------------------------------------------------------------------
# 6. NK-abundance proxy (PBMC bulk = mixture of cell types)
# ------------------------------------------------------------------------
nk_present <- intersect(NK_PROXY_GENES, rownames(G))
nk_proxy <- NULL
if (length(nk_present) >= 3) {
  Z <- t(scale(t(G[nk_present, , drop = FALSE])))
  nk_proxy <- colMeans(Z, na.rm = TRUE)
  say("\nNK proxy built from ", length(nk_present), " genes: ", paste(nk_present, collapse = ", "))
  tt <- t.test(nk_proxy ~ group)
  say("NK proxy PD vs Control: mean Control = ", round(tt$estimate[1], 3),
      ", PD = ", round(tt$estimate[2], 3), ", P = ", signif(tt$p.value, 3),
      "  (descriptive; a difference here means NK abundance, not NK-intrinsic expression, may drive signals)")
} else {
  say("\nFewer than 3 NK proxy genes on platform -> Model B (NK-adjusted) will be skipped.")
}

# ------------------------------------------------------------------------
# 7. limma: Model A (group) and Model B (group + NK proxy)
# ------------------------------------------------------------------------
fit_model <- function(design) {
  eBayes(lmFit(G, design), trend = TRUE, robust = TRUE)
}
dA <- model.matrix(~ group)
fitA <- fit_model(dA)
resA <- topTable(fitA, coef = "groupPD", number = Inf, sort.by = "none") %>%
  rownames_to_column("gene") %>% as_tibble()
write_csv(resA %>% arrange(P.Value), file.path(out_dir, "bulk_de_all_genes_modelA.csv"))
say("\nModel A (PD vs Control) genome-wide: FDR<0.05 = ", sum(resA$adj.P.Val < 0.05),
    " | nominal P<0.05 = ", sum(resA$P.Value < 0.05), " of ", nrow(resA), " genes")

resB <- NULL
if (!is.null(nk_proxy)) {
  dB <- model.matrix(~ group + nk_proxy)
  fitB <- fit_model(dB)
  resB <- topTable(fitB, coef = "groupPD", number = Inf, sort.by = "none") %>%
    rownames_to_column("gene") %>% as_tibble()
  write_csv(resB %>% arrange(P.Value), file.path(out_dir, "bulk_de_all_genes_modelB_NKadj.csv"))
  say("Model B (NK-adjusted): FDR<0.05 = ", sum(resB$adj.P.Val < 0.05),
      " | nominal P<0.05 = ", sum(resB$P.Value < 0.05))
}

# ------------------------------------------------------------------------
# 8. Per-candidate verdicts
# ------------------------------------------------------------------------
verdict <- function(d) {
  d <- d %>%
    left_join(resA %>% select(gene, probe_logFC_bulk = logFC, t_bulk = t, P_bulk = P.Value,
                              FDR_genomewide = adj.P.Val), by = "gene")
  if (!is.null(resB)) {
    d <- d %>% left_join(resB %>% select(gene, logFC_bulk_NKadj = logFC, P_bulk_NKadj = P.Value),
                         by = "gene")
  } else { d$logFC_bulk_NKadj <- NA_real_; d$P_bulk_NKadj <- NA_real_ }
  d <- d %>% rename(logFC_bulk = probe_logFC_bulk)
  testable <- !is.na(d$P_bulk)
  d$FDR_candset <- NA_real_
  d$FDR_candset[testable] <- p.adjust(d$P_bulk[testable], "BH")
  d$concordant <- sign(d$logFC_xiong) == sign(d$logFC_bulk)
  d$concordant_NKadj <- sign(d$logFC_xiong) == sign(d$logFC_bulk_NKadj)
  d$verdict <- case_when(
    !testable                                              ~ "Not measurable on platform",
    d$concordant & d$FDR_candset < FDR_CANDIDATE            ~ "VALIDATED (concordant, FDR<0.05 in candidate set)",
    d$concordant & d$P_bulk < P_NOMINAL                     ~ "Replicated (concordant, nominal P<0.05)",
    d$concordant                                           ~ "Concordant only (not significant)",
    !d$concordant & d$P_bulk < P_NOMINAL                    ~ "CONTRADICTED (opposite direction, P<0.05)",
    TRUE                                                   ~ "Discordant (not significant)")
  d$robust_to_NK_adjustment <- d$concordant_NKadj %in% TRUE & d$P_bulk_NKadj < P_NOMINAL
  d
}
cand_val <- lapply(names(sets), function(s) verdict(sets[[s]]) %>% mutate(set = s, .before = 1)) %>%
  setNames(names(sets))
for (s in names(cand_val)) {
  write_csv(cand_val[[s]] %>% arrange(P_bulk), file.path(out_dir, paste0("candidate_validation_", s, ".csv")))
}

say("\n--- Per-candidate verdicts ---")
for (s in names(cand_val)) {
  say(s, ":")
  tb <- table(cand_val[[s]]$verdict)
  for (nm in names(tb)) say("   ", sprintf("%-55s %d", nm, tb[[nm]]))
}
top_show <- cand_val$tier2_primary %>% filter(!is.na(P_bulk)) %>% arrange(P_bulk) %>%
  select(gene, logFC_xiong, logFC_bulk, P_bulk, FDR_candset, verdict) %>% head(15)
say("\nTop tier-2 candidates by bulk P-value:")
print(as.data.frame(top_show %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))

# ------------------------------------------------------------------------
# 9. Set-level replication tests (more powerful than gene-by-gene)
# ------------------------------------------------------------------------
say("\n--- Set-level replication tests ---")
tstat <- setNames(fitA$t[, "groupPD"], rownames(fitA$t))
global_rows <- list()
for (s in names(cand_val)) {
  d <- cand_val[[s]] %>% filter(!is.na(P_bulk))
  n <- nrow(d); k <- sum(d$concordant)
  bt <- if (n >= 5) binom.test(k, n, 0.5, alternative = "greater")$p.value else NA_real_
  up <- d$gene[d$logFC_xiong > 0]; dn <- d$gene[d$logFC_xiong < 0]
  gs_up <- if (length(up) >= 3) tryCatch(geneSetTest(names(tstat) %in% up, tstat, alternative = "up"), error = function(e) NA_real_) else NA_real_
  gs_dn <- if (length(dn) >= 3) tryCatch(geneSetTest(names(tstat) %in% dn, tstat, alternative = "down"), error = function(e) NA_real_) else NA_real_
  sp <- if (n >= 5) suppressWarnings(cor(d$logFC_xiong, d$logFC_bulk, method = "spearman")) else NA_real_
  global_rows[[s]] <- tibble(set = s, n_testable = n, n_concordant = k,
                             pct_concordant = round(100 * k / max(n, 1), 1),
                             binom_P_vs_50pct = bt, spearman_logFC = sp,
                             n_up = length(up), geneSetTest_P_up = gs_up,
                             n_down = length(dn), geneSetTest_P_down = gs_dn)
  say(sprintf("  %-14s testable=%d concordant=%d (%.1f%%) binomial P=%s | Spearman(logFC)=%s | set test: up P=%s, down P=%s",
              s, n, k, 100 * k / max(n, 1), signif(bt, 3), round(sp, 3), signif(gs_up, 3), signif(gs_dn, 3)))
}
global_tests <- bind_rows(global_rows)
write_csv(global_tests, file.path(out_dir, "set_level_replication_tests.csv"))
say("Interpretation guide: ~50% concordance = no replication (chance). ",
    "Clearly >50% with small P = the Xiong signature has real support even if single genes are underpowered.")

# ------------------------------------------------------------------------
# 10. Signature score per sample (feeds script 09 ROC)
#     Built ONLY from Xiong direction; no GSE100054 labels used for gene choice.
# ------------------------------------------------------------------------
sig_scores <- NULL
sig_genes <- cand_val$tier2_primary %>% filter(!is.na(P_bulk))
if (nrow(sig_genes) >= 3) {
  Zg <- t(scale(t(G[sig_genes$gene, , drop = FALSE])))
  w <- sign(sig_genes$logFC_xiong)          # +1 up in PD, -1 down in PD
  sig_scores <- tibble(sample_id = colnames(G), group = as.character(group),
                       signature_score = as.numeric(crossprod(w, Zg) / length(w)),
                       NK_proxy = if (!is.null(nk_proxy)) nk_proxy else NA_real_)
  write_csv(sig_scores, file.path(out_dir, "signature_scores_per_sample.csv"))
  wt <- tryCatch(wilcox.test(signature_score ~ group, data = sig_scores), error = function(e) NULL)
  say("\nXiong tier-2 signature score (direction-weighted mean z): Wilcoxon PD vs Control P = ",
      if (!is.null(wt)) signif(wt$p.value, 3) else "NA",
      "  [descriptive here; formal ROC with cross-validation is script 09]")
}

# ------------------------------------------------------------------------
# 11. Plots
# ------------------------------------------------------------------------
sc <- cand_val$tier2_primary %>% filter(!is.na(P_bulk))
if (nrow(sc) > 0) {
  p1 <- ggplot(sc, aes(logFC_xiong, logFC_bulk, colour = verdict)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text(data = sc %>% filter(P_bulk < P_NOMINAL), aes(label = gene),
              size = 3, vjust = -0.7, check_overlap = TRUE, show.legend = FALSE) +
    labs(title = "Xiong NK (scRNA pseudobulk) vs GSE100054 (bulk PBMC)",
         subtitle = "Points in top-right / bottom-left quadrants = same direction in both datasets",
         x = "log2FC PD vs Control, Xiong NK cells", y = "log2FC PD vs Control, GSE100054 bulk") +
    theme_minimal() + theme(legend.position = "bottom", legend.direction = "vertical")
  ggsave(file.path(plot_dir, "xiong_vs_bulk_logFC.png"), p1, width = 8, height = 7.5, dpi = 150)
  
  top_g <- sc %>% arrange(P_bulk) %>% head(12) %>% pull(gene)
  long <- as.data.frame(t(G[top_g, , drop = FALSE])) %>%
    rownames_to_column("sample_id") %>% mutate(group = group) %>%
    pivot_longer(all_of(top_g), names_to = "gene", values_to = "expr") %>%
    mutate(gene = factor(gene, levels = top_g))
  p2 <- ggplot(long, aes(group, expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_jitter(width = 0.15, size = 1.3) +
    facet_wrap(~ gene, scales = "free_y") +
    scale_fill_manual(values = c(Control = "steelblue", PD = "firebrick")) +
    labs(title = "GSE100054: top Xiong candidates by bulk P-value", x = NULL, y = "log2 expression") +
    theme_minimal() + theme(legend.position = "none")
  ggsave(file.path(plot_dir, "top_candidates_boxplots.png"), p2, width = 10, height = 8, dpi = 150)
}
if (!is.null(sig_scores)) {
  p3 <- ggplot(sig_scores, aes(group, signature_score, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_jitter(width = 0.15, size = 2) +
    scale_fill_manual(values = c(Control = "steelblue", PD = "firebrick")) +
    labs(title = "Xiong NK tier-2 signature score in GSE100054", x = NULL, y = "Signature score") +
    theme_minimal() + theme(legend.position = "none")
  ggsave(file.path(plot_dir, "signature_score_by_group.png"), p3, width = 5, height = 5, dpi = 150)
}
vt <- bind_rows(cand_val) %>% count(set, verdict)
p4 <- ggplot(vt, aes(verdict, n, fill = set)) + geom_col(position = "dodge") + coord_flip() +
  labs(title = "Validation verdicts", x = NULL, y = "Number of candidate genes") + theme_minimal()
ggsave(file.path(plot_dir, "verdict_counts.png"), p4, width = 9, height = 4.5, dpi = 150)
say("\nWrote plots to results/validation/gse100054/plots/")

# ------------------------------------------------------------------------
# 12. Save UI/downstream object + report
# ------------------------------------------------------------------------
saveRDS(list(
  candidate_validation = cand_val,
  coverage = coverage,
  bulk_de_modelA = resA, bulk_de_modelB_NKadj = resB,
  set_level_tests = global_tests,
  signature_scores = sig_scores,
  gene_expression = G, group = group, sample_metadata = meta,
  group_column = gcol, nk_proxy_genes = nk_present, nk_proxy = nk_proxy,
  settings = list(P_NOMINAL = P_NOMINAL, FDR_CANDIDATE = FDR_CANDIDATE),
  run_on = Sys.time()
), file.path(root, "results", "validated_data", "gse100054_validation.rds"))
writeLines(LOG, file.path(out_dir, "validation_report.txt"))

cat("\n", strrep("=", 70), "\n", sep = "")
cat("STEP 7 SUMMARY\n")
cat(strrep("=", 70), "\n")
for (s in names(cand_val)) {
  vv <- cand_val[[s]]
  cat(sprintf("%-14s | candidates %d | measurable %d | validated %d | replicated(nominal) %d | contradicted %d\n",
              s, nrow(vv), sum(!is.na(vv$P_bulk)),
              sum(str_detect(vv$verdict, "^VALIDATED")), sum(str_detect(vv$verdict, "^Replicated")),
              sum(str_detect(vv$verdict, "^CONTRADICTED"))))
}
cat("\nOutputs: results/validation/gse100054/  and  results/validated_data/gse100054_validation.rds\n")
cat("NEXT: script 08 (functional analysis) and 09 (ROC) should use genes from candidate_validation_*.csv\n")
cat(strrep("=", 70), "\n")

