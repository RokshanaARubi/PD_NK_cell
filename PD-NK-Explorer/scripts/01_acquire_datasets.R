
required_packages <- c("GEOquery", "rentrez", "httr2", "curl", "jsonlite")

missing_pkgs <- required_packages[!vapply(required_packages, requireNamespace,
                                          quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install with:\n",
    '  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")\n',
    '  BiocManager::install("GEOquery")\n',
    '  install.packages(c("rentrez", "httr2", "curl", "jsonlite"))'
  )
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(rentrez)
})

# ---- Create directory structure automatically ----
data_dirs <- c(
  "data",
  "data/xiong",
  "data/moquin_beaudry",
  "data/gse100054",
  "data/ppmi"
)
for (d in data_dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message("Created directory: ", d)
  }
}

# ---- Acquisition report accumulator ----
# Each dataset appends one row here; printed as the final report.
acquisition_log <- data.frame(
  Dataset = character(0),
  Repository_Accession = character(0),
  Files_Acquired = character(0),
  Location = character(0),
  Status = character(0),
  stringsAsFactors = FALSE
)

log_entry <- function(dataset, accession, files, location, status) {
  acquisition_log <<- rbind(
    acquisition_log,
    data.frame(Dataset = dataset, Repository_Accession = accession,
               Files_Acquired = files, Location = location, Status = status,
               stringsAsFactors = FALSE)
  )
}


## ---------------------------------------------------------------------------
## 1. XIONG ET AL. 2024 (GSE223138) — PRIMARY DISCOVERY, PBMC scRNA-seq
## ---------------------------------------------------------------------------
# Verified: npj Parkinsons Dis. 10:185 (2024), doi:10.1038/s41531-024-00790-3
# 6 PBMC samples (2 Control, 2 PD-Early, 2 PD-Late), 10x Genomics scRNA-seq.
# Publicly accessible GEO series (no login required as of verification date).

xiong_accession <- "GSE223138"
xiong_dir <- "data/xiong"

xiong_suppfiles_present <- length(list.files(xiong_dir, recursive = TRUE)) > 0

if (xiong_suppfiles_present) {
  message("[Xiong/", xiong_accession, "] Supplementary files already present in ",
          xiong_dir, " — skipping re-download.")
  xiong_status <- "Already present (skipped)"
} else {
  message("[Xiong/", xiong_accession, "] Downloading GEO series metadata and supplementary files...")
  xiong_status <- tryCatch({
    # Series-level metadata: sample titles, characteristics (donor/sample IDs,
    # disease status/stage if encoded in GEO characteristics fields), platform.
    xiong_gse <- GEOquery::getGEO(xiong_accession, GSEMatrix = TRUE, destdir = xiong_dir)
    xiong_pdata <- Biobase::pData(xiong_gse[[1]])
    write.csv(xiong_pdata, file.path(xiong_dir, "GSE223138_sample_metadata.csv"),
              row.names = FALSE)
    
    # Supplementary files (expected: per-sample 10x-style barcode/feature/matrix
    # files or similar scRNA-seq outputs — GEOquery pulls whatever the
    # submitters attached; verify exact contents after download, see
    # "WHAT I MUST VERIFY MANUALLY" below).
    GEOquery::getGEOSuppFiles(xiong_accession, baseDir = "data", makeDirectory = FALSE,
                              filter_regex = NULL)
    # getGEOSuppFiles writes to data/GSE223138 by default when makeDirectory=FALSE
    # with baseDir="data"; normalize into our target folder if needed.
    default_supp_dir <- file.path("data", xiong_accession)
    if (dir.exists(default_supp_dir) && !identical(normalizePath(default_supp_dir),
                                                   normalizePath(xiong_dir))) {
      file.copy(list.files(default_supp_dir, full.names = TRUE), xiong_dir,
                recursive = TRUE, overwrite = FALSE)
      unlink(default_supp_dir, recursive = TRUE)
    }
    "Downloaded"
  }, error = function(e) {
    message("  ERROR downloading ", xiong_accession, ": ", conditionMessage(e))
    "FAILED - see console error"
  })
}

xiong_files <- list.files(xiong_dir, recursive = TRUE)
log_entry(
  dataset = "Xiong et al. 2024 (Discovery scRNA-seq)",
  accession = xiong_accession,
  files = paste(length(xiong_files), "file(s)"),
  location = xiong_dir,
  status = xiong_status
)

if (xiong_status == "FAILED - see console error" || length(xiong_files) == 0) {
  warning("[Xiong/", xiong_accession, "] No files were acquired. Downstream steps ",
          "that depend on this dataset will fail until this is resolved manually.")
}


## ---------------------------------------------------------------------------
## 2. MOQUIN-BEAUDRY ET AL. 2025 (SRP525054 / PRJNA1145007)
##    INDEPENDENT scRNA-seq REPLICATION, PBMC — RAW READS ONLY
## ---------------------------------------------------------------------------
# Verified: Brain 2025;148(8):2847-2860, doi:10.1093/brain/awaf066
# (preprint: bioRxiv 10.1101/2024.07.26.605020)
#
# IMPORTANT: this is NOT a GEO series with a processed matrix. Only raw
# sc-RNAseq + TCR-seq FASTQ files are deposited, under NCBI SRA/BioProject.
# Final published cohort = 14 PD + 10 CTRL + 3 MSA + 3 PSP (78,876 cells).
#
# This script only retrieves the RUN-LEVEL METADATA (a small table listing
# every SRA run accession, so you know exactly what exists and can select
# PD/CTRL runs). It deliberately does NOT bulk-download the raw FASTQs:
# scRNA-seq FASTQs for ~79k cells across ~20 multiplexed donors are
# typically hundreds of GB to multiple TB, and turning them into an
# analyzable matrix requires running Cell Ranger yourself (out of scope
# for "acquisition and organization" and definitely out of scope for QC/
# analysis, which this script must not do).

moquin_bioproject <- "PRJNA1145007"
moquin_sra_study <- "SRP525054"
moquin_dir <- "data/moquin_beaudry"
moquin_meta_file <- file.path(moquin_dir, "moquin_beaudry_SRA_run_metadata.csv")

if (file.exists(moquin_meta_file)) {
  message("[Moquin-Beaudry/", moquin_sra_study, "] Run metadata already present at ",
          moquin_meta_file, " — skipping re-download.")
  moquin_status <- "Metadata already present (skipped)"
} else {
  message("[Moquin-Beaudry/", moquin_sra_study, "] Querying NCBI SRA for run-level metadata...")
  moquin_status <- tryCatch({
    search_res <- rentrez::entrez_search(
      db = "sra",
      term = paste0(moquin_bioproject, "[BioProject]"),
      retmax = 5000
    )
    if (length(search_res$ids) == 0) {
      stop("No SRA records found for BioProject ", moquin_bioproject,
           " — accession may have changed; verify manually at ",
           "https://www.ncbi.nlm.nih.gov/bioproject/", moquin_bioproject)
    }
    # Fetch run/sample summaries in batches (esummary) rather than full XML,
    # to keep this step lightweight (metadata only, no sequence data).
    summaries <- rentrez::entrez_summary(db = "sra", id = search_res$ids)
    if (!is.list(summaries[[1]])) summaries <- list(summaries)  # normalize single-result case
    
    moquin_meta <- do.call(rbind, lapply(summaries, function(s) {
      data.frame(
        sra_uid       = if (!is.null(s$uid)) s$uid else NA,
        experiment    = if (!is.null(s$expxml)) s$expxml else NA,
        run_accession = if (!is.null(s$runs)) s$runs else NA,
        stringsAsFactors = FALSE
      )
    }))
    write.csv(moquin_meta, moquin_meta_file, row.names = FALSE)
    
    writeLines(
      c(
        "This file lists a starting point for the raw SRA metadata for BioProject",
        paste0(moquin_bioproject, " (SRA study ", moquin_sra_study, ")."),
        "",
        "The 'run_accession' / 'experiment' columns are raw NCBI eSummary XML",
        "fragments, not a clean run table -- parse them or, more practically,",
        "download the canonical SraRunTable.txt directly from the SRA Run",
        "Selector for this BioProject:",
        paste0("  https://www.ncbi.nlm.nih.gov/Traces/study/?acc=", moquin_bioproject),
        "",
        "This script intentionally does NOT download the raw FASTQ/BAM files",
        "themselves (potentially hundreds of GB to multiple TB). Once you have",
        "identified the runs you need from the run table, retrieve them with",
        "the SRA Toolkit (not an R package -- install separately), e.g.:",
        "",
        "  prefetch <SRR_accession>",
        "  fasterq-dump <SRR_accession> --split-files -O data/moquin_beaudry/fastq/",
        "",
        "These raw reads will require Cell Ranger (10x Genomics, GEX + VDJ modes)",
        "to become an analyzable count matrix -- this is real downstream analysis",
        "work, not data acquisition, and is intentionally out of scope here."
      ),
      file.path(moquin_dir, "README_moquin_beaudry_ACCESS_NOTES.txt")
    )
    "Metadata downloaded (raw FASTQ NOT downloaded - see README in this folder)"
  }, error = function(e) {
    message("  ERROR querying SRA for ", moquin_bioproject, ": ", conditionMessage(e))
    "FAILED - see console error"
  })
}

log_entry(
  dataset = "Moquin-Beaudry et al. 2025 (Independent scRNA-seq replication)",
  accession = paste0(moquin_sra_study, " / ", moquin_bioproject),
  files = if (file.exists(moquin_meta_file)) "1 metadata file (+ README; raw FASTQ NOT fetched)" else "0 files",
  location = moquin_dir,
  status = moquin_status
)


## ---------------------------------------------------------------------------
## 3. GSE100054 — PBMC EXTERNAL VALIDATION (MICROARRAY, not RNA-seq)
## ---------------------------------------------------------------------------
# Verified: Illumina HumanHT-12 V3.0 Expression BeadChip, 10 PD + 9 HC = 19
# samples. Multiple independent reanalysis papers confirm this platform and
# sample count. Public GEO series, no login required.

gse100054_accession <- "GSE100054"
gse_dir <- "data/gse100054"
gse_expr_file <- file.path(gse_dir, "GSE100054_expression_matrix.csv")
gse_pdata_file <- file.path(gse_dir, "GSE100054_sample_metadata.csv")

if (file.exists(gse_expr_file) && file.exists(gse_pdata_file)) {
  message("[GSE100054] Expression matrix and metadata already present — skipping re-download.")
  gse_status <- "Already present (skipped)"
} else {
  message("[GSE100054] Downloading GEO series matrix (expression + phenotype + platform)...")
  gse_status <- tryCatch({
    gse_obj <- GEOquery::getGEO(gse100054_accession, GSEMatrix = TRUE, destdir = gse_dir)
    eset <- gse_obj[[1]]
    
    # Expression matrix as deposited in the GEO series matrix (this is
    # whatever processing level the submitters chose -- typically normalized
    # signal for BeadChip arrays, NOT raw/background-uncorrected intensities;
    # verify this manually, see notes below).
    expr_mat <- Biobase::exprs(eset)
    write.csv(expr_mat, gse_expr_file)
    
    # Sample-level phenotype metadata: should include disease/control labels.
    pdata <- Biobase::pData(eset)
    write.csv(pdata, gse_pdata_file, row.names = FALSE)
    
    # Platform annotation (GPL) -- probe-to-gene mapping.
    fdata <- Biobase::fData(eset)
    write.csv(fdata, file.path(gse_dir, "GSE100054_platform_annotation.csv"),
              row.names = FALSE)
    
    "Downloaded"
  }, error = function(e) {
    message("  ERROR downloading ", gse100054_accession, ": ", conditionMessage(e))
    "FAILED - see console error"
  })
}

log_entry(
  dataset = "GSE100054 (PBMC external validation, microarray)",
  accession = gse100054_accession,
  files = if (file.exists(gse_expr_file)) "expression matrix + sample metadata + platform annotation" else "0 files",
  location = gse_dir,
  status = gse_status
)


## ---------------------------------------------------------------------------
## 4. PPMI WHOLE-BLOOD RNA-SEQ — LARGE-SCALE EXTERNAL VALIDATION
##    CONTROLLED ACCESS: THIS SCRIPT NEVER ATTEMPTS TO DOWNLOAD PPMI DATA.
## ---------------------------------------------------------------------------
# Verified: PPMI RNA-seq (whole blood, ~1,570 subjects, longitudinal) is
# distributed exclusively through the LONI Image & Data Archive (IDA) after
# PPMI database registration; FASTQ/BAM additionally require submitting a
# signed Genetic Data Request Form to ppmi@loni.usc.edu. There is no public
# API, anonymous download, or Bioconductor package that legitimately
# retrieves this data.
#
# This block ONLY:
#   (a) documents the expected directory structure for files YOU have
#       already obtained through the legitimate PPMI/LONI process, and
#   (b) checks whether those files are already present, reporting status
#       without fabricating or attempting any download.

ppmi_dir <- "data/ppmi"
ppmi_expected_subdirs <- c(
  "clinical_metadata",     # e.g. participant status, UPDRS, visit schedule (from LONI IDA download, not an API)
  "rnaseq_counts",         # gene-level counts / TPM tables, if you were granted these directly
  "rnaseq_fastq",          # raw FASTQ, only if you completed the Genetic Data Request Form process
  "cbc_differential"       # complete blood count with differential, for blood-composition adjustment
)

for (sub in ppmi_expected_subdirs) {
  full <- file.path(ppmi_dir, sub)
  if (!dir.exists(full)) dir.create(full, recursive = TRUE)
}

if (!file.exists(file.path(ppmi_dir, "README_PPMI_ACCESS.txt"))) {
  writeLines(
    c(
      "PPMI whole-blood RNA-seq data is CONTROLLED ACCESS. This script never",
      "downloads it automatically and never will -- doing so would require",
      "credentials this script cannot legitimately possess.",
      "",
      "TO OBTAIN THIS DATA YOURSELF:",
      "  1. Register for a PPMI database account at https://www.ppmi-info.org",
      "  2. Sign the applicable Data Use Agreement.",
      "  3. Clinical/demographic data and CBC-with-differential: available",
      "     for browsing/export directly within the PPMI database once",
      "     registered (LONI IDA).",
      "  4. Processed RNA-seq quantifications (TPM, gene counts) and raw",
      "     FASTQ/BAM: submit a Genetic Data Request Form to",
      "     ppmi@loni.usc.edu (found in the Genetic Data Download section",
      "     of the PPMI database). Requests are filled in the order received.",
      "",
      "EXPECTED FILE PLACEMENT ONCE YOU HAVE LEGITIMATELY OBTAINED FILES:",
      "  data/ppmi/clinical_metadata/   -> participant status, UPDRS, visit info",
      "  data/ppmi/rnaseq_counts/       -> gene-level counts / TPM tables",
      "  data/ppmi/rnaseq_fastq/        -> raw FASTQ (only if you requested these)",
      "  data/ppmi/cbc_differential/    -> complete blood count w/ differential",
      "",
      "Once files are placed in these folders, re-run this script -- it will",
      "detect and report them without attempting to fetch anything itself."
    ),
    file.path(ppmi_dir, "README_PPMI_ACCESS.txt")
  )
}

ppmi_files_found <- list.files(ppmi_dir, recursive = TRUE)
ppmi_files_found <- ppmi_files_found[!grepl("README_PPMI_ACCESS.txt$", ppmi_files_found)]

ppmi_status <- if (length(ppmi_files_found) > 0) {
  "User-supplied files detected"
} else {
  "NOT ACQUIRED - controlled access, manual registration required (see README)"
}

log_entry(
  dataset = "PPMI whole-blood RNA-seq (Large-scale external validation)",
  accession = "No public accession (LONI IDA controlled access)",
  files = paste(length(ppmi_files_found), "user-supplied file(s)"),
  location = ppmi_dir,
  status = ppmi_status
)


## ---------------------------------------------------------------------------
## FINAL ACQUISITION REPORT
## ---------------------------------------------------------------------------

cat("\n")
cat(strrep("=", 78), "\n")
cat("DATA ACQUISITION REPORT\n")
cat(strrep("=", 78), "\n")
print(acquisition_log, row.names = FALSE)
cat(strrep("=", 78), "\n")

# Hard stop if a REQUIRED, publicly-downloadable dataset failed outright.
# (PPMI is expected to be "not acquired" on a fresh checkout -- that's not
# a script failure, it's a controlled-access dataset waiting on the user.)
hard_failures <- acquisition_log$Dataset[
  grepl("FAILED", acquisition_log$Status) &
    !grepl("PPMI", acquisition_log$Dataset)
]
if (length(hard_failures) > 0) {
  stop(
    "Acquisition failed for the following publicly-available dataset(s): ",
    paste(hard_failures, collapse = "; "),
    "\nSee console output above for the underlying error. Fix and re-run ",
    "this script -- it will skip datasets that already succeeded."
  )
}

message("\nAcquisition step complete. See report above for per-dataset status.")

