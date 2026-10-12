

suppressWarnings(suppressMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tools)
  library(tibble)
}))

# ------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------

project_root <- getwd()

data_dir        <- file.path(project_root, "data")
val_out_dir     <- file.path(project_root, "results", "data_validation")
validated_dir   <- file.path(project_root, "results", "validated_data")

dir.create(val_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(validated_dir, recursive = TRUE, showWarnings = FALSE)

NOT_VERIFIED <- "NOT VERIFIED"
NOT_AVAILABLE <- "NOT AVAILABLE"

# Simple null-coalescing helper (avoids depending on rlang's %||%)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Per-dataset log buffers (written to their own .txt files at the end)
xiong_log   <- character(0)
moquin_log  <- character(0)
gse_log     <- character(0)

log_to <- function(buffer_name, ...) {
  msg <- paste0(...)
  cat(msg, "\n")
  assign(buffer_name, c(get(buffer_name, envir = .GlobalEnv), msg), envir = .GlobalEnv)
}

# ------------------------------------------------------------------------
# Generic safe helpers
# ------------------------------------------------------------------------

safe_read_csv <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readr::read_csv(path, show_col_types = FALSE, guess_max = 5000, ...),
           error = function(e) {
             message("  [WARN] Failed to read ", path, ": ", conditionMessage(e))
             NULL
           })
}

safe_fread <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  tryCatch(data.table::fread(path, ...),
           error = function(e) {
             message("  [WARN] Failed to fread ", path, ": ", conditionMessage(e))
             NULL
           })
}

peek_gz_lines <- function(path, n = 200) {
  if (!file.exists(path)) return(character(0))
  tryCatch({
    con <- gzfile(path, open = "r")
    on.exit(close(con), add = TRUE)
    readLines(con, n = n, warn = FALSE)
  }, error = function(e) character(0))
}

file_exists_any <- function(dir, patterns) {
  # Returns TRUE if any file in `dir` (recursive) matches any regex in `patterns`
  if (!dir.exists(dir)) return(FALSE)
  files <- list.files(dir, recursive = TRUE, full.names = FALSE)
  if (length(files) == 0) return(FALSE)
  any(sapply(patterns, function(p) any(str_detect(tolower(files), p))))
}

list_matching_files <- function(dir, patterns) {
  if (!dir.exists(dir)) return(character(0))
  files <- list.files(dir, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(character(0))
  hits <- files[sapply(basename(tolower(files)), function(f) {
    any(sapply(patterns, function(p) str_detect(f, p)))
  })]
  hits
}

pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), x * 100)

# Identifier-type heuristics (never converts, only classifies)
classify_identifier <- function(ids) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) return("unknown - no identifiers to inspect")
  
  sample_ids <- ids[seq_len(min(200, length(ids)))]
  
  frac_ensembl <- mean(str_detect(sample_ids, "^ENS[A-Z]*G[0-9]+"), na.rm = TRUE)
  frac_entrez  <- mean(str_detect(sample_ids, "^[0-9]+$"), na.rm = TRUE)
  frac_affy    <- mean(str_detect(sample_ids, "_at$|_s_at$|_x_at$"), na.rm = TRUE)
  frac_symbol  <- mean(str_detect(sample_ids, "^[A-Za-z][A-Za-z0-9\\-\\.]{0,14}$") &
                         !str_detect(sample_ids, "^[0-9]+$"), na.rm = TRUE)
  
  scores <- c(
    ensembl = frac_ensembl,
    entrez  = frac_entrez,
    affy_probe = frac_affy,
    gene_symbol = frac_symbol
  )
  
  best <- names(scores)[which.max(scores)]
  best_score <- max(scores)
  
  if (best_score < 0.5) {
    paste0("other/unknown - no single ID pattern dominates (top guess: ", best,
           " at ", pct(best_score), " of sampled IDs)")
  } else {
    paste0(best, " (", pct(best_score), " of sampled IDs match pattern; heuristic only)")
  }
}

# Disease-label content heuristic: inspect VALUES, not just column names
classify_disease_values <- function(values) {
  v <- tolower(as.character(values))
  v <- v[!is.na(v) & nzchar(v)]
  if (length(v) == 0) return(list(status = "NOT VERIFIED", detail = "no values to inspect"))
  
  has_pd      <- any(str_detect(v, "parkinson|^pd$|\\bpd\\b"))
  has_control <- any(str_detect(v, "control|healthy|\\bhc\\b|normal"))
  has_other   <- any(!str_detect(v, "parkinson|^pd$|\\bpd\\b|control|healthy|\\bhc\\b|normal|^na$|^unknown$|^$"))
  has_missing <- any(str_detect(v, "^na$|^unknown$|^$|^n/a$"))
  
  groups <- c()
  if (has_pd) groups <- c(groups, "Parkinson's disease / PD")
  if (has_control) groups <- c(groups, "healthy control / control")
  if (has_other) groups <- c(groups, "other disease group(s) present")
  if (has_missing) groups <- c(groups, "unknown/missing values present")
  
  status <- if (has_pd && has_control) {
    "POSSIBLE"  # column name + value pattern both suggestive, but never auto-VERIFIED without manual confirmation
  } else if (has_pd || has_control) {
    "POSSIBLE"
  } else {
    "NOT VERIFIED"
  }
  
  list(status = status, detail = paste(groups, collapse = "; "))
}

# ==============================================================================
# 1. XIONG GSE223138
# ==============================================================================

log_to("xiong_log", "=================================================================")
log_to("xiong_log", "XIONG GSE223138 - DATASET VALIDATION")
log_to("xiong_log", "Generated: ", as.character(Sys.time()))
log_to("xiong_log", "=================================================================\n")

xiong_dir <- file.path(data_dir, "xiong")

xiong_meta_path   <- file.path(xiong_dir, "GSE223138_sample_metadata.csv")
xiong_series_path <- file.path(xiong_dir, "GSE223138_series_matrix.txt.gz")
xiong_soft_path   <- file.path(xiong_dir, "GPL20301.soft.gz")

xiong_expression_available <- NOT_AVAILABLE
xiong_metadata_available   <- NOT_AVAILABLE
xiong_sample_count_meta    <- NA_integer_
xiong_sample_count_expr    <- NA_integer_
xiong_cell_level_status    <- NOT_VERIFIED
xiong_disease_status       <- NOT_VERIFIED
xiong_disease_detail       <- NOT_VERIFIED
xiong_identifier_type      <- NA_character_
xiong_gene_or_probe_count  <- NA_integer_
xiong_overall_status       <- "NOT READY"
xiong_notes                <- character(0)

# --- 1a. Sample metadata ---
if (file.exists(xiong_meta_path)) {
  xiong_meta <- safe_read_csv(xiong_meta_path)
  if (!is.null(xiong_meta)) {
    xiong_metadata_available <- "AVAILABLE"
    xiong_sample_count_meta <- nrow(xiong_meta)
    log_to("xiong_log", "Loaded sample metadata: ", nrow(xiong_meta), " rows x ", ncol(xiong_meta), " columns.")
    
    target_cols <- c("cell type:ch1", "disease state:ch1", "tissue:ch1")
    for (tc in target_cols) {
      if (tc %in% colnames(xiong_meta)) {
        uniq_vals <- unique(xiong_meta[[tc]])
        log_to("xiong_log", "Column '", tc, "' unique values: ", paste(uniq_vals, collapse = " | "))
      } else {
        # try loose match in case of naming variants
        loose <- colnames(xiong_meta)[str_detect(tolower(colnames(xiong_meta)), fixed(tolower(tc)))]
        if (length(loose) > 0) {
          uniq_vals <- unique(xiong_meta[[loose[1]]])
          log_to("xiong_log", "Column '", loose[1], "' (matched loosely for '", tc, "') unique values: ",
                 paste(uniq_vals, collapse = " | "))
        } else {
          log_to("xiong_log", "Column '", tc, "' NOT FOUND in metadata. Reporting NOT VERIFIED.")
        }
      }
    }
    
    # Disease-label validation using actual VALUES, not just the column name
    disease_col <- if ("disease state:ch1" %in% colnames(xiong_meta)) "disease state:ch1" else NA
    if (!is.na(disease_col)) {
      dv <- classify_disease_values(xiong_meta[[disease_col]])
      xiong_disease_status <- dv$status
      xiong_disease_detail <- dv$detail
      log_to("xiong_log", "Disease label content check ('", disease_col, "'): status=", dv$status,
             " | groups found: ", dv$detail)
    }
    
    # Supplementary-file fields
    supp_cols <- colnames(xiong_meta)[str_detect(tolower(colnames(xiong_meta)), "supplementary")]
    if (length(supp_cols) > 0) {
      log_to("xiong_log", "Supplementary-file field(s) found in metadata: ", paste(supp_cols, collapse = ", "))
      for (sc in supp_cols) {
        vals <- unique(xiong_meta[[sc]])
        log_to("xiong_log", "  '", sc, "' values: ", paste(vals, collapse = " | "))
      }
    } else {
      log_to("xiong_log", "No 'supplementary_file'-type columns found in metadata.")
    }
  } else {
    log_to("xiong_log", "Metadata file exists but could not be parsed.")
  }
} else {
  log_to("xiong_log", "Sample metadata file NOT FOUND: ", xiong_meta_path)
}

# --- 1b. Series matrix: metadata-only vs. actual data table ---
if (file.exists(xiong_series_path)) {
  lines <- peek_gz_lines(xiong_series_path, n = 500)
  if (length(lines) == 0) {
    log_to("xiong_log", "Series matrix file exists but could not be read (peek failed).")
  } else {
    has_table_begin <- any(str_detect(lines, fixed("!series_matrix_table_begin")))
    has_table_end   <- any(str_detect(lines, fixed("!series_matrix_table_end")))
    
    if (has_table_begin) {
      begin_idx <- which(str_detect(lines, fixed("!series_matrix_table_begin")))[1]
      # the line right after table_begin should be the header row (e.g. "ID_REF" + GSM columns)
      header_idx <- begin_idx + 1
      if (header_idx <= length(lines)) {
        header_line <- lines[header_idx]
        n_fields <- length(str_split(header_line, "\t")[[1]])
        # Check if the very next line after header is table_end (i.e., an EMPTY data table)
        next_line_idx <- header_idx + 1
        empty_table <- (next_line_idx <= length(lines)) &&
          str_detect(lines[next_line_idx], fixed("!series_matrix_table_end"))
        
        if (empty_table) {
          log_to("xiong_log", "Series matrix contains a data-table header but NO data rows ",
                 "(table_begin immediately followed by table_end). This confirms the series matrix ",
                 "holds ONLY metadata/summary information, not actual expression values.")
          xiong_expression_available <- NOT_AVAILABLE
        } else {
          log_to("xiong_log", "Series matrix appears to contain a populated data table ",
                 "(header + at least one data row detected within first ", length(lines), " lines previewed). ",
                 "Header field count (approx): ", n_fields, ". This suggests some expression-like values ",
                 "MAY be present in the series matrix itself - requires full-file confirmation beyond this preview.")
          xiong_expression_available <- "POSSIBLY AVAILABLE (unconfirmed - only first lines previewed)"
        }
      }
    } else {
      log_to("xiong_log", "No '!series_matrix_table_begin' marker found within the first ", length(lines),
             " lines previewed. Series matrix likely contains only sample-level metadata/summary fields ",
             "(GEO header block), consistent with scRNA-seq studies where raw/processed data are deposited ",
             "as separate supplementary files rather than embedded in the series matrix.")
      xiong_expression_available <- NOT_AVAILABLE
    }
  }
} else {
  log_to("xiong_log", "Series matrix file NOT FOUND: ", xiong_series_path)
}

# --- 1c. Search for actual cell-level scRNA-seq files (never inferred from names alone) ---
cell_level_patterns <- c(
  "matrix\\.mtx", "barcodes", "features", "genes\\.tsv", "genes\\.txt",
  "\\.h5$", "\\.h5ad$", "\\.hdf5$", "\\.rds$", "\\.tar$", "\\.tar\\.gz$", "\\.zip$"
)
found_cell_level_files <- list_matching_files(xiong_dir, cell_level_patterns)

if (length(found_cell_level_files) > 0) {
  log_to("xiong_log", "Candidate cell-level / archive files found: ",
         paste(basename(found_cell_level_files), collapse = ", "),
         " -- presence of these filenames does NOT by itself confirm valid scRNA-seq content; ",
         "further inspection required before use.")
  xiong_cell_level_status <- "FILES PRESENT BUT CONTENT NOT VERIFIED"
} else {
  log_to("xiong_log", "No matrix.mtx / barcodes / features / h5 / h5ad / rds / tar / zip files found locally. ",
         "Cell-level scRNA-seq expression data are NOT AVAILABLE in the current download.")
  xiong_cell_level_status <- "NOT AVAILABLE LOCALLY"
  xiong_expression_available <- NOT_AVAILABLE
}

# --- 1d. Overall status for Xiong ---
if (xiong_metadata_available == "AVAILABLE" &&
    xiong_cell_level_status == "NOT AVAILABLE LOCALLY" &&
    xiong_expression_available == NOT_AVAILABLE) {
  xiong_overall_status <- "PARTIALLY AVAILABLE"
  xiong_notes <- c(xiong_notes, "Sample-level metadata present; cell-level expression matrix not present locally. Supplementary scRNA-seq files must be obtained before this dataset can be used for NK-cell identification.")
} else if (xiong_metadata_available != "AVAILABLE") {
  xiong_overall_status <- "NOT READY"
  xiong_notes <- c(xiong_notes, "Sample metadata missing or unreadable.")
} else if (str_detect(xiong_expression_available, "POSSIBLY")) {
  xiong_overall_status <- "PARTIALLY AVAILABLE"
  xiong_notes <- c(xiong_notes, "Series matrix may contain some tabular values; requires full confirmation.")
} else {
  xiong_overall_status <- "NOT READY"
}

log_to("xiong_log", "\nFINAL STATUS - Xiong GSE223138: ", xiong_overall_status)
writeLines(xiong_log, file.path(val_out_dir, "xiong_validation.txt"))

# ==============================================================================
# 2. MOQUIN-BEAUDRY
# ==============================================================================

log_to("moquin_log", "=================================================================")
log_to("moquin_log", "MOQUIN-BEAUDRY - DATASET VALIDATION")
log_to("moquin_log", "Generated: ", as.character(Sys.time()))
log_to("moquin_log", "=================================================================\n")

moquin_dir <- file.path(data_dir, "moquin_beaudry")
moquin_sra_path <- file.path(moquin_dir, "moquin_beaudry_SRA_run_metadata.csv")
moquin_readme_path <- file.path(moquin_dir, "README_moquin_beaudry_ACCESS_NOTES.txt")

moquin_expression_available <- NOT_AVAILABLE
moquin_metadata_available   <- NOT_AVAILABLE
moquin_sample_count_meta    <- NA_integer_
moquin_run_accessions       <- character(0)
moquin_overall_status       <- "NOT READY"
moquin_notes                <- character(0)

if (file.exists(moquin_readme_path)) {
  readme_txt <- tryCatch(readLines(moquin_readme_path, warn = FALSE), error = function(e) character(0))
  if (length(readme_txt) > 0) {
    log_to("moquin_log", "README_moquin_beaudry_ACCESS_NOTES.txt contents:")
    for (l in readme_txt) log_to("moquin_log", "  ", l)
  }
} else {
  log_to("moquin_log", "No README/access-notes file found (optional).")
}

if (file.exists(moquin_sra_path)) {
  moquin_sra <- safe_read_csv(moquin_sra_path)
  if (!is.null(moquin_sra)) {
    moquin_metadata_available <- "AVAILABLE (SRA run metadata only)"
    moquin_sample_count_meta <- nrow(moquin_sra)
    log_to("moquin_log", "Loaded SRA run metadata: ", nrow(moquin_sra), " rows x ", ncol(moquin_sra), " columns.")
    log_to("moquin_log", "Columns present: ", paste(colnames(moquin_sra), collapse = ", "))
    
    if ("run_accession" %in% colnames(moquin_sra)) {
      moquin_run_accessions <- unique(as.character(moquin_sra$run_accession))
      log_to("moquin_log", "Number of SRA runs: ", length(moquin_run_accessions))
      log_to("moquin_log", "run_accession values: ", paste(moquin_run_accessions, collapse = ", "))
    } else {
      log_to("moquin_log", "'run_accession' column NOT FOUND. Cannot enumerate run accessions.")
    }
  } else {
    log_to("moquin_log", "SRA metadata file exists but could not be parsed.")
  }
} else {
  log_to("moquin_log", "SRA run metadata file NOT FOUND: ", moquin_sra_path)
}

# Search for a processed expression matrix locally (explicitly: do NOT download FASTQ)
matrix_patterns <- c(
  "matrix\\.mtx", "barcodes", "features", "\\.h5$", "\\.h5ad$", "\\.hdf5$",
  "\\.rds$", "count", "expr", "\\.tar$", "\\.tar\\.gz$", "\\.zip$"
)
moquin_matrix_candidates <- list_matching_files(moquin_dir, matrix_patterns)
# Exclude the SRA metadata / README themselves from being misclassified as a matrix
moquin_matrix_candidates <- moquin_matrix_candidates[
  !basename(moquin_matrix_candidates) %in% c(basename(moquin_sra_path), basename(moquin_readme_path))
]

if (length(moquin_matrix_candidates) > 0) {
  log_to("moquin_log", "Candidate processed-matrix files found: ",
         paste(basename(moquin_matrix_candidates), collapse = ", "),
         " -- content NOT verified by this script; inspect manually before use.")
  moquin_expression_available <- "CANDIDATE FILE(S) PRESENT - CONTENT NOT VERIFIED"
} else {
  log_to("moquin_log", "No processed expression matrix, 10x files, h5/h5ad, or Seurat object found locally. ",
         "Only SRA run-level metadata is currently available.")
  log_to("moquin_log", "Per project safety rules, FASTQ/SRA raw sequencing data will NOT be downloaded by this ",
         "or any automated step without explicit separate authorization.")
  moquin_expression_available <- NOT_AVAILABLE
}

# Overall status
if (moquin_expression_available == NOT_AVAILABLE && moquin_metadata_available != NOT_AVAILABLE) {
  moquin_overall_status <- "PARTIALLY AVAILABLE"
  moquin_notes <- c(moquin_notes, "Only SRA run accessions available; no processed expression matrix. FASTQ processing required before use, not performed here.")
} else if (moquin_metadata_available == NOT_AVAILABLE) {
  moquin_overall_status <- "NOT READY"
  moquin_notes <- c(moquin_notes, "SRA run metadata missing or unreadable.")
} else if (str_detect(moquin_expression_available, "CANDIDATE")) {
  moquin_overall_status <- "PARTIALLY AVAILABLE"
  moquin_notes <- c(moquin_notes, "Candidate matrix file(s) present but unverified.")
}

log_to("moquin_log", "\nFINAL STATUS - Moquin-Beaudry: ", moquin_overall_status)
writeLines(moquin_log, file.path(val_out_dir, "moquin_validation.txt"))

# ==============================================================================
# 3. GSE100054
# ==============================================================================

log_to("gse_log", "=================================================================")
log_to("gse_log", "GSE100054 - DATASET VALIDATION")
log_to("gse_log", "Generated: ", as.character(Sys.time()))
log_to("gse_log", "=================================================================\n")

gse_dir <- file.path(data_dir, "gse100054")

gse_expr_path       <- file.path(gse_dir, "GSE100054_expression_matrix.csv")
gse_meta_path       <- file.path(gse_dir, "GSE100054_sample_metadata.csv")
gse_platform_path   <- file.path(gse_dir, "GSE100054_platform_annotation.csv")
gse_series_path     <- file.path(gse_dir, "GSE100054_series_matrix.txt.gz")
gse_soft_path       <- file.path(gse_dir, "GPL23126.soft.gz")  # NEVER loaded fully - see bug note below

# --- Explicit guard against the known prior bug ---
# The previous inspection script incorrectly treated GPL23126.soft.gz (a platform
# annotation file) as a candidate expression matrix because it was the largest file
# in the directory. This script NEVER uses file size to choose the expression matrix.
# The expression matrix path is hard-coded to the known, explicitly-named file.
if (file.exists(gse_soft_path)) {
  log_to("gse_log", "NOTE: GPL23126.soft.gz detected. This is a PLATFORM ANNOTATION file and will NOT ",
         "be loaded into memory or treated as an expression matrix (known prior-script bug, now guarded against).")
}

gse_expression_available <- NOT_AVAILABLE
gse_metadata_available   <- NOT_AVAILABLE
gse_dims                 <- NA_character_
gse_n_samples_expr       <- NA_integer_
gse_n_features           <- NA_integer_
gse_id_col               <- NA_character_
gse_n_unique_ids         <- NA_integer_
gse_n_duplicate_ids      <- NA_integer_
gse_expr_min             <- NA_real_
gse_expr_max             <- NA_real_
gse_expr_mean            <- NA_real_
gse_expr_median          <- NA_real_
gse_expr_na_count        <- NA_integer_
gse_expr_type_heuristic  <- NOT_VERIFIED
gse_identifier_type      <- NA_character_
gse_disease_status       <- NOT_VERIFIED
gse_disease_detail       <- NOT_VERIFIED
gse_sample_id_col_meta   <- NA_character_
gse_matched_samples      <- NA_integer_
gse_unmatched_expr       <- character(0)
gse_unmatched_meta       <- character(0)
gse_overall_status       <- "NOT READY"
gse_notes                <- character(0)
gse_expr_dt              <- NULL
gse_meta_df              <- NULL

# --- 3a. Expression matrix ---
if (file.exists(gse_expr_path)) {
  gse_expr_dt <- safe_fread(gse_expr_path)
  if (!is.null(gse_expr_dt) && nrow(gse_expr_dt) > 0) {
    gse_expression_available <- "AVAILABLE"
    gse_dims <- paste0(nrow(gse_expr_dt), " rows x ", ncol(gse_expr_dt), " columns")
    log_to("gse_log", "Loaded GSE100054_expression_matrix.csv: ", gse_dims)
    
    # Identifier column: assume first column unless it is clearly numeric-only expression data
    first_col_name <- colnames(gse_expr_dt)[1]
    first_col_vals <- gse_expr_dt[[1]]
    first_col_is_character_like <- is.character(first_col_vals) ||
      (is.factor(first_col_vals)) ||
      mean(is.na(suppressWarnings(as.numeric(as.character(first_col_vals[1:min(20, length(first_col_vals))]))))) > 0.5
    
    if (first_col_is_character_like) {
      gse_id_col <- first_col_name
      ids <- as.character(first_col_vals)
      gse_n_unique_ids <- length(unique(ids))
      gse_n_duplicate_ids <- sum(duplicated(ids))
      gse_n_features <- nrow(gse_expr_dt)
      gse_n_samples_expr <- ncol(gse_expr_dt) - 1
      
      log_to("gse_log", "Identifier column detected: '", gse_id_col, "'")
      log_to("gse_log", "Unique identifiers: ", gse_n_unique_ids, " | Duplicate identifiers: ", gse_n_duplicate_ids)
      log_to("gse_log", "Rows (features): ", gse_n_features, " | Columns (samples, excluding ID col): ", gse_n_samples_expr)
      log_to("gse_log", "Interpretation: rows appear to represent genes/probes, columns appear to represent samples ",
             "(heuristic based on identifier column position and shape; not independently confirmed beyond this file).")
      
      gse_identifier_type <- classify_identifier(ids)
      log_to("gse_log", "Identifier type (heuristic, NOT converted): ", gse_identifier_type)
      
      numeric_mat <- gse_expr_dt[, -1, with = FALSE]
      numeric_mat <- suppressWarnings(as.matrix(numeric_mat))
      storage.mode(numeric_mat) <- "numeric"
      
      gse_expr_na_count <- sum(is.na(numeric_mat))
      finite_vals <- numeric_mat[is.finite(numeric_mat)]
      
      if (length(finite_vals) > 0) {
        gse_expr_min    <- min(finite_vals)
        gse_expr_max    <- max(finite_vals)
        gse_expr_mean   <- mean(finite_vals)
        gse_expr_median <- median(finite_vals)
        
        log_to("gse_log", "Expression value summary -> min: ", round(gse_expr_min, 3),
               " | max: ", round(gse_expr_max, 3),
               " | mean: ", round(gse_expr_mean, 3),
               " | median: ", round(gse_expr_median, 3),
               " | NA count: ", gse_expr_na_count)
        
        has_negative <- any(finite_vals < 0)
        has_non_integer <- any(finite_vals != floor(finite_vals))
        
        gse_expr_type_heuristic <- if (has_negative) {
          "HEURISTIC ONLY: likely log-ratio / normalized microarray intensity (negative values present)"
        } else if (has_non_integer && gse_expr_max < 30) {
          "HEURISTIC ONLY: likely log-transformed or normalized expression (non-integer, small dynamic range)"
        } else if (has_non_integer) {
          "HEURISTIC ONLY: likely normalized/continuous intensity values (non-integer, wide dynamic range)"
        } else if (!has_non_integer && gse_expr_max > 30) {
          "HEURISTIC ONLY: possibly raw or minimally processed counts (non-negative integers, wide range)"
        } else {
          NOT_VERIFIED
        }
        log_to("gse_log", "Expression type assessment (heuristic only, NOT a confirmed fact): ", gse_expr_type_heuristic)
      } else {
        log_to("gse_log", "No finite numeric values found in expression matrix body.")
      }
    } else {
      log_to("gse_log", "First column does not look like an identifier column (appears numeric). ",
             "Identifier column detection: NOT VERIFIED. Manual inspection required.")
    }
  } else {
    log_to("gse_log", "Expression matrix file exists but could not be parsed or is empty.")
  }
} else {
  log_to("gse_log", "Expression matrix file NOT FOUND: ", gse_expr_path)
}

# --- 3b. Sample metadata ---
if (file.exists(gse_meta_path)) {
  gse_meta_df <- safe_read_csv(gse_meta_path)
  if (!is.null(gse_meta_df)) {
    gse_metadata_available <- "AVAILABLE"
    log_to("gse_log", "\nLoaded GSE100054_sample_metadata.csv: ", nrow(gse_meta_df), " rows x ", ncol(gse_meta_df), " columns.")
    log_to("gse_log", "Metadata columns: ", paste(colnames(gse_meta_df), collapse = ", "))
    
    # Identify disease/control-related columns by NAME first, then verify by VALUE
    candidate_disease_cols <- colnames(gse_meta_df)[
      str_detect(tolower(colnames(gse_meta_df)), "disease|diagnosis|condition|group|status|phenotype")
    ]
    
    if (length(candidate_disease_cols) > 0) {
      log_to("gse_log", "Candidate disease/phenotype column(s) by name: ", paste(candidate_disease_cols, collapse = ", "))
      for (cc in candidate_disease_cols) {
        uniq_vals <- unique(gse_meta_df[[cc]])
        log_to("gse_log", "  Column '", cc, "' unique values: ", paste(uniq_vals, collapse = " | "))
        dv <- classify_disease_values(gse_meta_df[[cc]])
        log_to("gse_log", "    Value-content check: status=", dv$status, " | groups found: ", dv$detail)
        if (dv$status == "POSSIBLE" && (gse_disease_status == NOT_VERIFIED)) {
          gse_disease_status <- dv$status
          gse_disease_detail <- paste0(cc, ": ", dv$detail)
        }
      }
    } else {
      log_to("gse_log", "No obviously named disease/phenotype columns found. Disease labels: NOT VERIFIED.")
    }
    
    # Identify a sample-ID column in metadata
    id_candidates <- colnames(gse_meta_df)[
      str_detect(tolower(colnames(gse_meta_df)), "sample|geo_accession|^id$|gsm")
    ]
    if (length(id_candidates) > 0) {
      gse_sample_id_col_meta <- id_candidates[1]
      log_to("gse_log", "Metadata sample-ID column guessed as: '", gse_sample_id_col_meta, "'")
    } else {
      log_to("gse_log", "Could not confidently identify a sample-ID column in metadata by name.")
    }
  } else {
    log_to("gse_log", "Sample metadata file exists but could not be parsed.")
  }
} else {
  log_to("gse_log", "Sample metadata file NOT FOUND: ", gse_meta_path)
}

# --- 3c. Sample-ID consistency check (expression vs. metadata) ---
if (!is.null(gse_expr_dt) && !is.null(gse_meta_df) && !is.na(gse_id_col)) {
  expr_sample_ids <- colnames(gse_expr_dt)[-1]
  
  if (!is.na(gse_sample_id_col_meta)) {
    meta_sample_ids <- as.character(gse_meta_df[[gse_sample_id_col_meta]])
  } else {
    # fall back: try every column and pick the one with the highest overlap with expr column names
    overlaps <- sapply(gse_meta_df, function(col) {
      length(intersect(as.character(col), expr_sample_ids))
    })
    if (length(overlaps) > 0 && max(overlaps, na.rm = TRUE) > 0) {
      best_col <- names(overlaps)[which.max(overlaps)]
      meta_sample_ids <- as.character(gse_meta_df[[best_col]])
      gse_sample_id_col_meta <- best_col
      log_to("gse_log", "Fallback sample-ID match: column '", best_col, "' has the highest overlap with expression column names.")
    } else {
      meta_sample_ids <- character(0)
    }
  }
  
  if (length(meta_sample_ids) > 0) {
    matched <- intersect(expr_sample_ids, meta_sample_ids)
    unmatched_expr <- setdiff(expr_sample_ids, meta_sample_ids)
    unmatched_meta <- setdiff(meta_sample_ids, expr_sample_ids)
    
    gse_matched_samples <- length(matched)
    gse_unmatched_expr <- unmatched_expr
    gse_unmatched_meta <- unmatched_meta
    
    log_to("gse_log", "\nSample-ID consistency check:")
    log_to("gse_log", "  Expression samples: ", length(expr_sample_ids))
    log_to("gse_log", "  Metadata samples: ", length(meta_sample_ids))
    log_to("gse_log", "  Matched samples: ", length(matched))
    log_to("gse_log", "  Unmatched expression samples: ", length(unmatched_expr),
           if (length(unmatched_expr) > 0) paste0(" (", paste(head(unmatched_expr, 10), collapse = ", "),
                                                  if (length(unmatched_expr) > 10) ", ..." else "", ")") else "")
    log_to("gse_log", "  Unmatched metadata samples: ", length(unmatched_meta),
           if (length(unmatched_meta) > 0) paste0(" (", paste(head(unmatched_meta, 10), collapse = ", "),
                                                  if (length(unmatched_meta) > 10) ", ..." else "", ")") else "")
  } else {
    log_to("gse_log", "Could not extract metadata sample IDs for matching. Sample-ID consistency: NOT VERIFIED.")
  }
} else {
  log_to("gse_log", "Sample-ID consistency check skipped: expression matrix and/or metadata unavailable.")
}

# --- 3d. Platform annotation (use CSV already created; never load .soft.gz fully) ---
if (file.exists(gse_platform_path)) {
  log_to("gse_log", "\nPlatform annotation CSV found (GSE100054_platform_annotation.csv). ",
         "Will be used in place of GPL23126.soft.gz if platform info is needed downstream. Not loaded in full here ",
         "since it is not required for structural validation.")
} else {
  log_to("gse_log", "\nPlatform annotation CSV not found. GPL23126.soft.gz exists but will NOT be loaded ",
         "(per safety rules) unless explicitly required in a future step.")
}

# --- 3e. Overall status for GSE100054 ---
if (gse_expression_available == "AVAILABLE" && gse_metadata_available == "AVAILABLE" &&
    !is.na(gse_matched_samples) && gse_matched_samples > 0) {
  gse_overall_status <- if (gse_matched_samples == length(colnames(gse_expr_dt)) - 1) {
    "READY"
  } else {
    "PARTIALLY AVAILABLE"
  }
  if (gse_overall_status == "PARTIALLY AVAILABLE") {
    gse_notes <- c(gse_notes, "Not all expression samples have matching metadata (or vice versa); resolve before use.")
  }
} else if (gse_expression_available == "AVAILABLE" && gse_metadata_available != "AVAILABLE") {
  gse_overall_status <- "PARTIALLY AVAILABLE"
  gse_notes <- c(gse_notes, "Expression matrix present but sample metadata missing/unreadable.")
} else {
  gse_overall_status <- "NOT READY"
  gse_notes <- c(gse_notes, "Expression matrix and/or metadata missing or unreadable.")
}

log_to("gse_log", "\nFINAL STATUS - GSE100054: ", gse_overall_status)
writeLines(gse_log, file.path(val_out_dir, "gse100054_validation.txt"))

# ==============================================================================
# 4. Sample matching report (machine-readable, all datasets)
# ==============================================================================

sample_matching_rows <- list()

if (!is.null(gse_expr_dt) && !is.na(gse_id_col)) {
  expr_ids_gse <- colnames(gse_expr_dt)[-1]
  meta_ids_gse <- if (!is.na(gse_sample_id_col_meta) && !is.null(gse_meta_df)) {
    as.character(gse_meta_df[[gse_sample_id_col_meta]])
  } else {
    character(0)
  }
  all_ids <- union(expr_ids_gse, meta_ids_gse)
  if (length(all_ids) > 0) {
    sample_matching_rows[["gse100054"]] <- tibble(
      dataset = "gse100054",
      sample_id = all_ids,
      in_expression = all_ids %in% expr_ids_gse,
      in_metadata = all_ids %in% meta_ids_gse
    )
  }
}

sample_matching_report <- if (length(sample_matching_rows) > 0) {
  bind_rows(sample_matching_rows)
} else {
  tibble(dataset = character(0), sample_id = character(0),
         in_expression = logical(0), in_metadata = logical(0))
}

write_csv(sample_matching_report, file.path(val_out_dir, "sample_matching_report.csv"))

# ==============================================================================
# 5. Master validation summary CSV
# ==============================================================================

validation_summary <- tibble(
  dataset = c("xiong", "moquin_beaudry", "gse100054"),
  repository_accession = c("GSE223138", "PRJNA1145007 (SRA BioProject, unconfirmed against local files)", "GSE100054"),
  data_type = c("scRNA-seq (PBMC)", "scRNA-seq (PBMC)", "bulk transcriptomics (PBMC)"),
  expression_available = c(xiong_expression_available, moquin_expression_available, gse_expression_available),
  metadata_available = c(xiong_metadata_available, moquin_metadata_available, gse_metadata_available),
  sample_count_expression = c(NA_integer_, NA_integer_, gse_n_samples_expr),
  sample_count_metadata = c(xiong_sample_count_meta, moquin_sample_count_meta, nrow(gse_meta_df %||% tibble())),
  matched_sample_count = c(NA_integer_, NA_integer_, gse_matched_samples),
  gene_or_probe_count = c(NA_integer_, NA_integer_, gse_n_features),
  identifier_type = c(NA_character_, NA_character_, gse_identifier_type),
  disease_label_status = c(xiong_disease_status, NOT_VERIFIED, gse_disease_status),
  cell_level_status = c(xiong_cell_level_status, "NOT AVAILABLE LOCALLY (SRA metadata only)", "NOT APPLICABLE (bulk data)"),
  expression_status = c(xiong_expression_available, moquin_expression_available, gse_expression_available),
  overall_status = c(xiong_overall_status, moquin_overall_status, gse_overall_status),
  notes = c(
    paste(xiong_notes, collapse = " | "),
    paste(moquin_notes, collapse = " | "),
    paste(gse_notes, collapse = " | ")
  )
)

write_csv(validation_summary, file.path(val_out_dir, "dataset_validation_summary.csv"))

# ==============================================================================
# 6. Write validated RDS objects ONLY for datasets that are READY
# ==============================================================================

if (xiong_overall_status == "READY") {
  # Not expected to trigger given current files, but included for reproducibility.
  saveRDS(list(metadata = xiong_meta), file.path(validated_dir, "xiong_validated.rds"))
  log_to("xiong_log", "Wrote results/validated_data/xiong_validated.rds (status=READY).")
} else {
  cat("Xiong: overall_status =", xiong_overall_status, "-> RDS NOT created (only READY datasets get an RDS).\n")
}

if (moquin_overall_status == "READY") {
  log_to("moquin_log", "Wrote results/validated_data/moquin_validated.rds (status=READY).")
} else {
  cat("Moquin-Beaudry: overall_status =", moquin_overall_status, "-> RDS NOT created.\n")
}

if (gse_overall_status == "READY" && !is.null(gse_expr_dt) && !is.null(gse_meta_df)) {
  saveRDS(
    list(
      expression = gse_expr_dt,
      metadata = gse_meta_df,
      identifier_column = gse_id_col,
      identifier_type = gse_identifier_type,
      sample_id_column_metadata = gse_sample_id_col_meta,
      validated_on = Sys.time()
    ),
    file.path(validated_dir, "gse100054_validated.rds")
  )
  cat("GSE100054: overall_status = READY -> Wrote results/validated_data/gse100054_validated.rds\n")
} else {
  cat("GSE100054: overall_status =", gse_overall_status, "-> RDS NOT created (only READY datasets get an RDS).\n")
}

# ==============================================================================
# 7. Final console report
# ==============================================================================

cat("\n")
cat("====================================================\n")
cat("PD-NK EXPLORER DATASET VALIDATION\n")
cat("====================================================\n\n")

cat("Xiong GSE223138:\n")
cat("  Expression:            ", xiong_expression_available, "\n")
cat("  Metadata:              ", xiong_metadata_available, " (", ifelse(is.na(xiong_sample_count_meta), "NA", xiong_sample_count_meta), " samples)\n", sep = "")
cat("  Cell-level scRNA-seq:  ", xiong_cell_level_status, "\n")
cat("  Disease labels:        ", xiong_disease_status, " (", xiong_disease_detail, ")\n", sep = "")
cat("  Overall:               ", xiong_overall_status, "\n\n")

cat("Moquin-Beaudry:\n")
cat("  Expression:            ", moquin_expression_available, "\n")
cat("  Metadata:              ", moquin_metadata_available, "\n")
cat("  SRA runs:              ", length(moquin_run_accessions), "\n")
cat("  Overall:               ", moquin_overall_status, "\n\n")

cat("GSE100054:\n")
cat("  Expression:            ", gse_expression_available, " (", ifelse(is.na(gse_dims), "dims NOT VERIFIED", gse_dims), ")\n", sep = "")
cat("  Metadata:              ", gse_metadata_available, "\n")
cat("  Samples matched:       ", ifelse(is.na(gse_matched_samples), "NOT VERIFIED",
                                        paste0(gse_matched_samples, "/", length(colnames(gse_expr_dt %||% tibble())[-1]))), "\n")
cat("  Identifier type:       ", ifelse(is.na(gse_identifier_type), "NOT VERIFIED", gse_identifier_type), "\n")
cat("  Disease labels:        ", gse_disease_status, " (", gse_disease_detail, ")\n", sep = "")
cat("  Overall:               ", gse_overall_status, "\n\n")

cat("PPMI:\n")
cat("  NOT INCLUDED - controlled access; will be added later.\n\n")

cat("====================================================\n")
cat("Outputs written to:\n")
cat("  results/data_validation/dataset_validation_summary.csv\n")
cat("  results/data_validation/xiong_validation.txt\n")
cat("  results/data_validation/moquin_validation.txt\n")
cat("  results/data_validation/gse100054_validation.txt\n")
cat("  results/data_validation/sample_matching_report.csv\n")
cat("  results/validated_data/*.rds (only for READY datasets)\n")
cat("====================================================\n")

