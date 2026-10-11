# PD-NK Explorer - Final report

Generated: 2026-10-10 01:16:27.211556

**Research question.** Can NK-cell-specific transcriptomic signatures, identified from single-cell RNA-seq and validated across independent bulk RNA-seq datasets, serve as reproducible peripheral biomarkers for Parkinson's disease?

## Data actually used
- Discovery: Xiong et al. GSE223138, PBMC scRNA-seq, 6 donors (2 Control, 2 Early PD, 2 Late PD).
- Independent validation: GSE100054, bulk PBMC microarray (GPL23126), 9 Control vs 10 PD.
- NOT used: Moquin-Beaudry PRJNA1145007 (only SRA run metadata available; no processed matrix), PPMI (controlled access).

## Results
1. After QC 55,669 cells were retained; 11,117 were called NK (clusters 12, 3, 5).
2. Donor-level pseudobulk edgeR (2 Control vs 4 PD): 0 of 12,296 genes reached FDR < 0.05. All candidates are therefore exploratory (Tier 2: nominal P < 0.01 and |log2FC| >= 1).
3. In GSE100054 the Tier-2 candidates showed 42.6% direction concordance (chance = 50%), binomial P = 0.879.
4. Final gene status: Contradicted by bulk data = 7; Not measurable in bulk = 21; Not replicated = 115; Suggestive (nominal in both datasets, not FDR-significant) = 9.
   Suggestive genes (nominal in both datasets, not FDR-significant, to be treated as hypotheses only): GNAQ, PRKCQ, JADE2, IKZF1, BCL11B, DDHD2, IL7, FGL2, CD109.
5. Pathway analysis: fgsea reported many significant pathways, but gene-permutation GSEA is anti-conservative for correlated sets. A correlation-aware camera test in the bulk data confirmed 0 of the Xiong-significant Hallmark/KEGG sets.
6. ROC (direction fixed a priori, signatures defined from Xiong only): best AUC = 0.633 (95% CI 0.36-0.906), permutation P = 0.175; 0 of 6 signatures had P < 0.05. With n = 19 these AUCs are indistinguishable from chance.
7. Bulk PBMC cell-composition proxies, PD vs Control: Neutrophil (P=0.035); T_cell (P=0.0653); Monocyte (P=0.0172); B_cell (P=0.035). Composition differences, not NK-intrinsic expression, are the most visible blood signal.

## Conclusion
In the available public data, NK-cell transcriptomic signatures from Xiong GSE223138 did NOT replicate in GSE100054 and cannot be called reproducible peripheral biomarkers for Parkinson's disease. This is a negative / inconclusive result, not evidence that no NK signal exists.

## Limitations
- Discovery had only 2 Control and 4 PD donors (low power); stage contrasts used 2 donors per group.
- Validation set is a different technology (microarray) and bulk (mixed cell types), so NK-specific changes are diluted.
- Probe-to-gene mapping used gene_assignment symbols; lncRNA-type genes (AC/AL/LINC) are often absent from the array.
- No replication cohort at single-cell level (Moquin-Beaudry raw data not processed). No PPMI.
- Batch/donor effects in the scRNA-seq object were not formally integrated.

## Next steps that could change the conclusion
- Process Moquin-Beaudry (PRJNA1145007) from FASTQ and repeat steps 3-7 as a true scRNA-seq replication.
- Larger single-cell or sorted-NK cohorts; PPMI blood RNA-seq once access is granted.
