

required_packages <- c("curl", "stringr")
missing_pkgs <- required_packages[!vapply(required_packages, requireNamespace,
                                          quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install with: install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(curl)
  library(stringr)
})

# ------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------

project_root <- getwd()
xiong_dir    <- file.path(project_root, "data", "xiong")
raw_dir      <- file.path(xiong_dir, "raw")
meta_path    <- file.path(xiong_dir, "GSE223138_sample_metadata.csv")

if (!file.exists(meta_path)) {
  stop(
    "REQUIRED FILE MISSING: ", meta_path, "\n",
    "This script reads the real supplementary_file URLs from Xiong's already-",
    "downloaded sample metadata. Run scripts/01_acquire_datasets.R first."
  )
}

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("geo_accession", "supplementary_file_1", "supplementary_file_2", "supplementary_file_3")
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols) > 0) {
  stop("Expected column(s) not found in ", meta_path, ": ", paste(missing_cols, collapse = ", "),
       ". This script expects the exact structure produced by scripts/01_acquire_datasets.R.")
}

cat("=================================================================\n")
cat("XIONG GSE223138 - SCRNA-SEQ SUPPLEMENTARY FILE ACQUISITION\n")
cat("Generated: ", as.character(Sys.time()), "\n")
cat("=================================================================\n\n")
cat("Samples found in metadata: ", nrow(meta), "\n\n")

# ------------------------------------------------------------------------
# 1. Download helper (per URL): skip if present, verify filename pattern,
#    never overwrite, log clearly on failure without halting other samples.
# ------------------------------------------------------------------------

expected_suffix_pattern <- "barcodes\\.tsv\\.gz$|features\\.tsv\\.gz$|genes\\.tsv\\.gz$|matrix\\.mtx\\.gz$"

# NCBI's FTP server occasionally returns transient errors (e.g. "450"/"451 -
# file temporarily unavailable") when hit with several back-to-back
# connections, especially for larger files like matrix.mtx.gz. This is a
# server-side hiccup, not a bad URL - so each file gets a few retries with
# an increasing delay before being reported as a real failure.
MAX_ATTEMPTS <- 4
RETRY_DELAY_SECONDS <- c(5, 15, 30)  # delay before attempts 2, 3, 4

download_one <- function(url, dest_dir) {
  if (is.na(url) || !nzchar(url)) {
    return(list(status = "SKIPPED - no URL in metadata", path = NA_character_))
  }
  
  fname <- basename(url)
  
  if (!str_detect(tolower(fname), expected_suffix_pattern)) {
    return(list(status = paste0("REFUSED - filename does not match expected 10x pattern: ", fname),
                path = NA_character_))
  }
  
  dest_path <- file.path(dest_dir, fname)
  
  if (file.exists(dest_path) && file.info(dest_path)$size > 0) {
    return(list(status = "Already present (skipped)", path = dest_path))
  }
  
  last_error <- NA_character_
  
  for (attempt in seq_len(MAX_ATTEMPTS)) {
    if (attempt > 1) {
      delay <- RETRY_DELAY_SECONDS[attempt - 1]
      cat("    retrying in ", delay, "s (attempt ", attempt, "/", MAX_ATTEMPTS, ")...\n", sep = "")
      Sys.sleep(delay)
    }
    
    result <- tryCatch({
      curl::curl_download(url, destfile = dest_path, quiet = TRUE, mode = "wb")
      if (file.exists(dest_path) && file.info(dest_path)$size > 0) {
        list(status = "Downloaded", path = dest_path)
      } else {
        list(status = "FAILED - download produced an empty or missing file", path = NA_character_)
      }
    }, error = function(e) {
      list(status = paste0("FAILED - ", conditionMessage(e)), path = NA_character_)
    })
    
    if (result$status == "Downloaded") return(result)
    last_error <- result$status
    # Clean up any partial file before retrying
    if (file.exists(dest_path)) unlink(dest_path)
  }
  
  list(status = paste0(last_error, " (after ", MAX_ATTEMPTS, " attempts)"), path = NA_character_)
}

# ------------------------------------------------------------------------
# 2. Download all three files for each of the 6 GSM samples
# ------------------------------------------------------------------------

acquisition_log <- data.frame(
  gsm = character(0), file_role = character(0), filename = character(0),
  status = character(0), stringsAsFactors = FALSE
)

for (i in seq_len(nrow(meta))) {
  gsm <- meta$geo_accession[i]
  sample_dir <- file.path(raw_dir, gsm)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("--- ", gsm, " ---\n")
  
  files_for_sample <- list(
    barcodes = meta$supplementary_file_1[i],
    features = meta$supplementary_file_2[i],
    matrix   = meta$supplementary_file_3[i]
  )
  
  for (role in names(files_for_sample)) {
    url <- files_for_sample[[role]]
    res <- download_one(url, sample_dir)
    cat("  [", role, "] ", res$status, "\n", sep = "")
    
    acquisition_log <- rbind(acquisition_log, data.frame(
      gsm = gsm, file_role = role,
      filename = if (!is.na(res$path)) basename(res$path) else basename(ifelse(is.na(url), "NA", url)),
      status = res$status,
      stringsAsFactors = FALSE
    ))
  }
  cat("\n")
}

# ------------------------------------------------------------------------
# 3. Per-sample completeness check (all 3 files present = sample ready)
# ------------------------------------------------------------------------

cat("=================================================================\n")
cat("PER-SAMPLE COMPLETENESS\n")
cat("=================================================================\n")

sample_status <- data.frame(gsm = unique(acquisition_log$gsm), complete = NA)
for (i in seq_len(nrow(sample_status))) {
  gsm <- sample_status$gsm[i]
  rows <- acquisition_log[acquisition_log$gsm == gsm, ]
  ok <- all(rows$status %in% c("Downloaded", "Already present (skipped)"))
  sample_status$complete[i] <- ok
  cat(gsm, ": ", if (ok) "COMPLETE (all 3 files present)" else "INCOMPLETE - see log above", "\n", sep = "")
}

n_complete <- sum(sample_status$complete)
cat("\nSamples fully acquired: ", n_complete, " / ", nrow(sample_status), "\n")

# ------------------------------------------------------------------------
# 4. Write machine-readable acquisition log
# ------------------------------------------------------------------------

out_dir <- file.path(project_root, "results", "data_acquisition")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(acquisition_log, file.path(out_dir, "xiong_scrna_supplementary_acquisition_log.csv"), row.names = FALSE)

cat("\nWrote results/data_acquisition/xiong_scrna_supplementary_acquisition_log.csv\n")
cat("Files saved under: data/xiong/raw/<GSM_accession>/\n")

if (n_complete == nrow(sample_status)) {
  cat("\nAll 6 Xiong samples now have real barcodes/features/matrix.mtx files locally.\n")
  cat("NEXT STEP: scRNA-seq QC, clustering, and NK-cell identification (methodology steps 2-4)\n")
  cat("can now run on real data - this is the next script to write.\n")
} else {
  cat("\nSome samples are incomplete. Re-run this script after checking your network/VPN -\n")
  cat("it will only retry the files that failed, not re-download what already succeeded.\n")
}
cat("=================================================================\n")

