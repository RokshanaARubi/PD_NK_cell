
pkgs <- c("dplyr", "readr", "stringr", "tibble", "tidyr", "ggplot2", "pROC", "limma", "msigdbr")
miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(miss) > 0) stop("Missing packages: ", paste(miss, collapse = ", "),
                           "\n  install.packages('pROC')")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tibble)
  library(tidyr); library(ggplot2); library(pROC); library(limma); library(msigdbr)
})
set.seed(42)

# ---- SETTINGS ----------------------------------------------------------------
PADJ_SET     <- 0.05    # Xiong pathway padj used to DEFINE signature gene sets
MIN_GENES    <- 5       # min measurable genes for a signature
N_PERM       <- 5000
# ------------------------------------------------------------------------------

root    <- getwd()
val_rds <- file.path(root, "results", "validated_data", "gse100054_validation.rds")
fun_rds <- file.path(root, "results", "validated_data", "functional_analysis.rds")
out_dir <- file.path(root, "results", "roc", "gse100054")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
for (f in c(val_rds, fun_rds)) if (!file.exists(f)) stop("REQUIRED FILE MISSING: ", f)

v   <- readRDS(val_rds)
fun <- readRDS(fun_rds)
G     <- v$gene_expression
group <- v$group                       # factor Control / PD
cat("GSE100054 gene-level matrix: ", nrow(G), " genes x ", ncol(G), " samples (",
    paste(names(table(group)), table(group), sep = "=", collapse = ", "), ")\n", sep = "")

zrow <- function(M) { Z <- t(scale(t(M))); Z[is.na(Z)] <- 0; Z }   # label-free scaling

# ------------------------------------------------------------------------
# 1. Signatures defined from Xiong ONLY
# ------------------------------------------------------------------------
sigs <- list()

# (a) Tier-2 gene signature from script 07 (direction-weighted z)
t2 <- v$candidate_validation$tier2_primary %>% filter(!is.na(P_bulk))
if (nrow(t2) >= MIN_GENES) {
  w <- sign(t2$logFC_xiong)
  sigs$Tier2_gene_signature <- as.numeric(crossprod(w, zrow(G[t2$gene, , drop = FALSE])) / length(w))
}

# (b) Pathway-derived signatures: leading-edge genes of Xiong Hallmark/KEGG sets
gx <- fun$gsea_xiong_primary %>% filter(collection %in% c("Hallmark", "KEGG"), padj < PADJ_SET)
split_le <- function(d) unique(toupper(unlist(strsplit(d$leadingEdge, ";"))))
up_genes <- split_le(gx %>% filter(NES > 0))
dn_genes <- split_le(gx %>% filter(NES < 0))
both <- intersect(up_genes, dn_genes); up_genes <- setdiff(up_genes, both); dn_genes <- setdiff(dn_genes, both)
up_m <- intersect(up_genes, rownames(G)); dn_m <- intersect(dn_genes, rownames(G))
cat("Pathway signature genes (Xiong Hallmark/KEGG padj<", PADJ_SET, "): up-in-PD ",
    length(up_m), " | down-in-PD ", length(dn_m), "  [from ", nrow(gx), " sets]\n", sep = "")
if (length(up_m) >= MIN_GENES) sigs$Pathway_UP_score   <- colMeans(zrow(G[up_m, , drop = FALSE]))
if (length(dn_m) >= MIN_GENES) sigs$Pathway_DOWN_score_inverted <- -colMeans(zrow(G[dn_m, , drop = FALSE]))
if (all(c("Pathway_UP_score", "Pathway_DOWN_score_inverted") %in% names(sigs)))
  sigs$Pathway_composite <- (sigs$Pathway_UP_score + sigs$Pathway_DOWN_score_inverted) / 2

# (c) NK proxy (control: is NK abundance itself discriminative?)
if (!is.null(v$nk_proxy)) sigs$NK_abundance_proxy_PDhigh <- as.numeric(v$nk_proxy)

sig_df <- as_tibble(sigs) %>% mutate(sample_id = colnames(G), group = as.character(group), .before = 1)
write_csv(sig_df, file.path(out_dir, "signature_scores_all.csv"))

# ------------------------------------------------------------------------
# 2. ROC (direction fixed: PD > Control), CI, permutation P
# ------------------------------------------------------------------------
auc_fast <- function(score, y) {            # y logical (TRUE = PD)
  r <- rank(score); n1 <- sum(y); n0 <- sum(!y)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
y <- group == "PD"
rows <- list(); rocs <- list()
for (nm in names(sigs)) {
  s <- sigs[[nm]]
  ro <- roc(response = group, predictor = s, levels = c("Control", "PD"),
            direction = "<", quiet = TRUE)
  ci <- suppressWarnings(tryCatch(as.numeric(ci.auc(ro, method = "delong")), error = function(e) rep(NA_real_, 3)))
  obs <- auc_fast(s, y)
  perm <- replicate(N_PERM, auc_fast(s, sample(y)))
  p_perm <- (sum(perm >= obs) + 1) / (N_PERM + 1)            # one-sided: AUC > 0.5
  rows[[nm]] <- tibble(signature = nm, AUC = round(as.numeric(auc(ro)), 3),
                       CI_low = round(ci[1], 3), CI_high = round(ci[3], 3),
                       perm_P_one_sided = signif(p_perm, 3))
  rocs[[nm]] <- ro
}
roc_tbl <- bind_rows(rows)

# ------------------------------------------------------------------------
# 3. Leave-one-out CV for a 2-feature model (UP score + DOWN score)
#    Features are Xiong-defined and label-free, so no selection leakage.
# ------------------------------------------------------------------------
loo_tbl <- NULL
if (all(c("Pathway_UP_score", "Pathway_DOWN_score_inverted") %in% names(sigs))) {
  D <- tibble(y = as.integer(y), up = sigs$Pathway_UP_score, dn = sigs$Pathway_DOWN_score_inverted)
  pred <- numeric(nrow(D))
  for (i in seq_len(nrow(D))) {
    fit <- suppressWarnings(glm(y ~ up + dn, data = D[-i, ], family = binomial()))
    pred[i] <- predict(fit, newdata = D[i, ], type = "response")
  }
  ro_loo <- roc(response = group, predictor = pred, levels = c("Control", "PD"), direction = "<", quiet = TRUE)
  rocs$LOOCV_logistic_UP_DOWN <- ro_loo
  loo_tbl <- tibble(signature = "LOOCV_logistic_UP_DOWN", AUC = round(as.numeric(auc(ro_loo)), 3),
                    CI_low = NA_real_, CI_high = NA_real_,
                    perm_P_one_sided = signif((sum(replicate(N_PERM, auc_fast(sample(pred), y)) >= auc_fast(pred, y)) + 1) / (N_PERM + 1), 3))
  roc_tbl <- bind_rows(roc_tbl, loo_tbl)
}
write_csv(roc_tbl, file.path(out_dir, "roc_summary.csv"))
cat("\n=== ROC (GSE100054; PD = positive; direction fixed a priori) ===\n")
print(as.data.frame(roc_tbl))

# ------------------------------------------------------------------------
# 4. Composition check: is the pathway score just NK abundance?
# ------------------------------------------------------------------------
comp_tbl <- NULL
if (!is.null(v$nk_proxy) && "Pathway_composite" %in% names(sigs)) {
  ct <- suppressWarnings(cor.test(sigs$Pathway_composite, as.numeric(v$nk_proxy), method = "spearman"))
  comp_tbl <- tibble(test = "Spearman: Pathway_composite vs NK proxy", rho = unname(ct$estimate), P = ct$p.value)
  cat(sprintf("\nPathway composite vs NK proxy: rho = %.3f, P = %.3g\n", comp_tbl$rho, comp_tbl$P))
  write_csv(comp_tbl, file.path(out_dir, "composition_check.csv"))
}
# neutrophil / lymphocyte-type proxies (label-free; descriptive)
proxy_sets <- list(Neutrophil = c("S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "MNDA"),
                   T_cell     = c("CD3D", "CD3E", "CD3G", "IL7R", "TRAC"),
                   Monocyte   = c("CD14", "LYZ", "CSF1R", "CST3"),
                   B_cell     = c("MS4A1", "CD79A", "CD79B"))
prox <- lapply(proxy_sets, function(g) { g <- intersect(g, rownames(G)); if (length(g) >= 2) colMeans(zrow(G[g, , drop = FALSE])) else NULL })
prox <- Filter(Negate(is.null), prox)
prox_tbl <- bind_rows(lapply(names(prox), function(nm) {
  wt <- suppressWarnings(wilcox.test(prox[[nm]] ~ group))
  tibble(cell_type_proxy = nm, mean_Control = mean(prox[[nm]][!y]), mean_PD = mean(prox[[nm]][y]), wilcox_P = wt$p.value)
}))
cat("\nBulk cell-composition proxies (PD vs Control):\n"); print(as.data.frame(prox_tbl %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
write_csv(prox_tbl, file.path(out_dir, "composition_proxies.csv"))

# ------------------------------------------------------------------------
# 5. Robustness: limma::camera (accounts for inter-gene correlation)
# ------------------------------------------------------------------------
cat("\n=== Robustness: camera on bulk for top Xiong Hallmark/KEGG sets ===\n")
new_api <- "collection" %in% names(formals(msigdbr))
getset <- function(...) tryCatch(msigdbr(species = "Homo sapiens", ...), error = function(e) NULL)
H <- if (new_api) getset(collection = "H") else getset(category = "H")
K <- if (new_api) getset(collection = "C2", subcollection = "CP:KEGG_LEGACY") else getset(category = "C2", subcategory = "CP:KEGG")
dfs <- Filter(Negate(is.null), list(H, K))
cam_tbl <- NULL
if (length(dfs) > 0) {
  df <- bind_rows(dfs)
  sets <- lapply(split(toupper(df$gene_symbol), df$gs_name), unique)
  idx <- ids2indices(sets, rownames(G)); idx <- idx[lengths(idx) >= 10 & lengths(idx) <= 500]
  design <- model.matrix(~ group)
  cam <- camera(G, idx, design, contrast = 2, inter.gene.cor = NULL, sort = FALSE)  # estimates correlation per set
  cam_tbl <- cam %>% rownames_to_column("pathway") %>% as_tibble() %>%
    rename(P_camera_bulk = PValue, FDR_camera_bulk = FDR, Direction_bulk = Direction) %>%
    inner_join(fun$gsea_xiong_primary %>% filter(collection %in% c("Hallmark", "KEGG")) %>%
                 select(pathway, NES_xiong = NES, P_xiong = pval, padj_xiong = padj),
               by = "pathway") %>%
    mutate(same_direction = (NES_xiong > 0) == (Direction_bulk == "Up")) %>%
    arrange(P_xiong)
  write_csv(cam_tbl, file.path(out_dir, "camera_robustness_bulk.csv"))
  top <- cam_tbl %>% filter(padj_xiong < 0.05)
  cat(sprintf("Xiong-significant Hallmark/KEGG sets: %d | camera P<0.05 in bulk AND same direction: %d | camera FDR<0.05 same direction: %d\n",
              nrow(top), sum(top$P_camera_bulk < 0.05 & top$same_direction),
              sum(top$FDR_camera_bulk < 0.05 & top$same_direction)))
  print(as.data.frame(top %>% select(pathway, NES_xiong, P_xiong, Direction_bulk, P_camera_bulk, FDR_camera_bulk) %>%
                        head(15) %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
} else cat("[skip] msigdbr gene sets not available for camera check\n")

# ------------------------------------------------------------------------
# 6. Plots
# ------------------------------------------------------------------------
roc_df <- bind_rows(lapply(names(rocs), function(nm)
  tibble(signature = nm, FPR = 1 - rocs[[nm]]$specificities, TPR = rocs[[nm]]$sensitivities)))
lab <- roc_tbl %>% mutate(label = sprintf("%s (AUC %.2f)", signature, AUC))
roc_df <- roc_df %>% left_join(lab %>% select(signature, label), by = "signature")
p1 <- ggplot(roc_df, aes(FPR, TPR, colour = label)) +
  geom_abline(linetype = 2, colour = "grey60") + geom_step(linewidth = 0.9) +
  labs(title = "ROC in GSE100054 (PD vs Control)",
       subtitle = "Signatures defined from Xiong NK data only; direction fixed a priori",
       x = "1 - specificity", y = "Sensitivity", colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(plot_dir, "roc_curves.png"), p1, width = 7.5, height = 7.5, dpi = 150)

p2 <- ggplot(roc_tbl, aes(reorder(signature, AUC), AUC)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high)) + coord_flip(ylim = c(0, 1)) +
  labs(title = "AUC with 95% CI (DeLong)", x = NULL, y = "AUC (0.5 = chance)") + theme_minimal()
ggsave(file.path(plot_dir, "auc_summary.png"), p2, width = 8, height = 4.5, dpi = 150)

long <- sig_df %>% pivot_longer(-c(sample_id, group), names_to = "signature", values_to = "score")
p3 <- ggplot(long, aes(group, score, fill = group)) + geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.5) + facet_wrap(~ signature, scales = "free_y") +
  scale_fill_manual(values = c(Control = "steelblue", PD = "firebrick")) +
  labs(title = "Signature scores by group (GSE100054)", x = NULL, y = "Score") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(plot_dir, "signature_scores_boxplots.png"), p3, width = 10, height = 7, dpi = 150)

# ------------------------------------------------------------------------
# 7. Save + summary
# ------------------------------------------------------------------------
saveRDS(list(roc_summary = roc_tbl, signature_scores = sig_df, composition_proxies = prox_tbl,
             composition_check = comp_tbl, camera_robustness = cam_tbl,
             pathway_signature_genes = list(up = up_m, down = dn_m),
             settings = list(PADJ_SET = PADJ_SET, N_PERM = N_PERM), run_on = Sys.time()),
        file.path(root, "results", "validated_data", "roc_analysis.rds"))

best <- roc_tbl %>% arrange(desc(AUC)) %>% slice(1)
cat("\n", strrep("=", 70), "\nSTEP 9 SUMMARY\n", strrep("=", 70), "\n", sep = "")
cat("Samples: ", paste(names(table(group)), table(group), sep = "=", collapse = ", "), "\n")
cat("Best signature: ", best$signature, " AUC = ", best$AUC, " [", best$CI_low, "-", best$CI_high,
    "] perm P = ", best$perm_P_one_sided, "\n", sep = "")
cat("Signatures with perm P < 0.05: ", sum(roc_tbl$perm_P_one_sided < 0.05, na.rm = TRUE), " of ", nrow(roc_tbl), "\n")
cat("Outputs: results/roc/gse100054/  and  results/validated_data/roc_analysis.rds\n")
cat("NEXT: script 10 (final reproducible-biomarker summary) + wire results into the UI\n")
cat(strrep("=", 70), "\n")