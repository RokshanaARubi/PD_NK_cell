

suppressWarnings(suppressMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tools)
}))

# ------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------

# Assumes the script is run from the project root (PD-NK-Explorer/).
project_root <- getwd()

data_dir    <- file.path(project_root, "data")
out_dir     <- file.path(project_root, "results", "data_inspection")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

NOT_VERIFIED <- "NOT VERIFIED"

log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  log_lines <<- c(log_lines, msg)
  cat(msg, "\n")
}

safe_read_lines <- function(path, n = 5) {
  # Reads at most n lines from a plain-text file, safely.
  out <- tryCatch({
    con <- file(path, open = "r")
    on.exit(close(con), add = TRUE)
    readLines(con, n = n, warn = FALSE)
  }, error = function(e) NA_character_)
  out
}

safe_read_gz_lines <- function(path, n = 5) {
  # Reads at most n lines from a .gz file WITHOUT extracting it to disk.
  out <- tryCatch({
    con <- gzfile(path, open = "r")
    on.exit(close(con), add = TRUE)
    readLines(con, n = n, warn = FALSE)
  }, error = function(e) NA_character_)
  out
}

str_trunc_safe <- function(x, width = 200) {
  if (is.na(x)) return(NA_character_)
  if (nchar(x) <= width) return(x)
  paste0(substr(x, 1, width), " ... [truncated]")
}

human_bytes <- function(x) {
  if (is.na(x)) return(NA_character_)
  units <- c("B", "KB", "MB", "GB", "TB")
  i <- 1
  val <- x
  while (val >= 1024 && i < length(units)) {
    val <- val / 1024
    i <- i + 1
  }
  sprintf("%.2f %s", val, units[i])
}

# Recognized delimited-text extensions we attempt to preview.
delim_exts <- c("csv", "tsv", "txt")

# ------------------------------------------------------------------------
# 1. File inventory (all three datasets)
# ------------------------------------------------------------------------

log_msg("=================================================================")
log_msg("PD-NK Explorer :: Dataset Inspection Report")
log_msg("Generated: ", as.character(Sys.time()))
log_msg("Project root: ", project_root)
log_msg("=================================================================\n")

dataset_dirs <- list(
  xiong          = file.path(data_dir, "xiong"),
  moquin_beaudry = file.path(data_dir, "moquin_beaudry"),
  gse100054      = file.path(data_dir, "gse100054")
)

inventory_list <- list()

for (ds_name in names(dataset_dirs)) {
  ds_path <- dataset_dirs[[ds_name]]
  
  if (!dir.exists(ds_path)) {
    log_msg("[", ds_name, "] Directory does not exist: ", ds_path)
    next
  }
  
  files <- list.files(ds_path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  
  if (length(files) == 0) {
    log_msg("[", ds_name, "] Directory exists but contains no files: ", ds_path)
    next
  }
  
  info <- file.info(files)
  
  df <- data.frame(
    dataset       = ds_name,
    filename      = basename(files),
    extension     = tolower(tools::file_ext(files)),
    size_bytes    = info$size,
    size_human    = vapply(info$size, human_bytes, character(1)),
    relative_path = sub(paste0("^", project_root, "/?"), "", files),
    full_path     = files,
    stringsAsFactors = FALSE
  )
  
  inventory_list[[ds_name]] <- df
}

dataset_file_inventory <- if (length(inventory_list) > 0) {
  bind_rows(inventory_list)
} else {
  data.frame(
    dataset = character(0), filename = character(0), extension = character(0),
    size_bytes = numeric(0), size_human = character(0),
    relative_path = character(0), full_path = character(0)
  )
}

write_csv(dataset_file_inventory, file.path(out_dir, "dataset_file_inventory.csv"))
log_msg("Wrote: results/data_inspection/dataset_file_inventory.csv (",
        nrow(dataset_file_inventory), " files across all datasets)\n")

# ------------------------------------------------------------------------
# 2. Helper: classify a file by likely role, based on name/extension only
#    (never based on content assumptions)
# ------------------------------------------------------------------------

classify_file <- function(filename) {
  fname_lower <- tolower(filename)
  
  role <- character(0)
  
  if (str_detect(fname_lower, "barcodes")) role <- c(role, "10x_barcodes")
  if (str_detect(fname_lower, "features|genes\\.tsv|genes\\.txt")) role <- c(role, "10x_features_genes")
  if (str_detect(fname_lower, "matrix\\.mtx")) role <- c(role, "10x_sparse_matrix")
  if (str_detect(fname_lower, "\\.h5ad")) role <- c(role, "h5ad_object")
  if (str_detect(fname_lower, "\\.h5$|\\.hdf5$")) role <- c(role, "h5_object")
  if (str_detect(fname_lower, "\\.rds$")) role <- c(role, "r_serialized_object_possibly_seurat")
  if (str_detect(fname_lower, "seurat")) role <- c(role, "seurat_object_named")
  if (str_detect(fname_lower, "\\.loom$")) role <- c(role, "loom_object")
  if (str_detect(fname_lower, "count|counts|expr|expression|matrix")) role <- c(role, "possible_expression_matrix")
  if (str_detect(fname_lower, "meta|metadata|pheno|clinical|sample_?sheet|coldata")) role <- c(role, "possible_metadata")
  if (str_detect(fname_lower, "cell.?info|cell.?meta|obs\\.")) role <- c(role, "possible_cell_metadata")
  if (str_detect(fname_lower, "donor|subject|patient")) role <- c(role, "possible_donor_metadata")
  if (str_detect(fname_lower, "platform|gpl|annot")) role <- c(role, "possible_platform_annotation")
  if (str_detect(fname_lower, "series_matrix|soft")) role <- c(role, "geo_series_level_file")
  if (str_detect(fname_lower, "supp|readme|md5|checksum")) role <- c(role, "supplementary_or_admin_file")
  if (str_detect(fname_lower, "\\.fastq|\\.fq(\\.gz)?$|\\.sra$")) role <- c(role, "raw_sequencing_reads")
  
  if (length(role) == 0) role <- "unclassified"
  paste(role, collapse = "; ")
}

# ------------------------------------------------------------------------
# 3. Helper: preview a delimited text file (plain or gzipped) safely
# ------------------------------------------------------------------------

preview_delim_file <- function(path, n_preview = 5) {
  ext <- tolower(tools::file_ext(path))
  is_gz <- str_detect(tolower(path), "\\.gz$")
  
  raw_lines <- if (is_gz) {
    safe_read_gz_lines(path, n = n_preview)
  } else {
    safe_read_lines(path, n = n_preview)
  }
  
  if (length(raw_lines) == 0 || all(is.na(raw_lines))) {
    return(list(ok = FALSE, header = NA_character_, preview_lines = NA,
                n_cols_guess = NA_integer_, delimiter_guess = NA_character_))
  }
  
  header <- raw_lines[1]
  
  delim_guess <- if (str_detect(header, "\t")) {
    "\t"
  } else if (str_detect(header, ",")) {
    ","
  } else if (str_detect(header, ";")) {
    ";"
  } else {
    " "
  }
  
  n_cols_guess <- length(str_split(header, delim_guess, simplify = FALSE)[[1]])
  
  list(
    ok = TRUE,
    header = header,
    preview_lines = raw_lines,
    n_cols_guess = n_cols_guess,
    delimiter_guess = delim_guess
  )
}

# ------------------------------------------------------------------------
# 4. Helper: attempt to read a metadata-like delimited file into a data.frame
#    (bounded read; will not attempt on huge files)
# ------------------------------------------------------------------------

MAX_METADATA_READ_BYTES <- 200 * 1024 * 1024  # 200 MB safety cap

try_read_metadata_table <- function(path) {
  size <- tryCatch(file.info(path)$size, error = function(e) NA_real_)
  if (is.na(size) || size > MAX_METADATA_READ_BYTES) {
    return(NULL)
  }
  
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path, ignore.case = TRUE)))
  
  df <- tryCatch({
    if (str_detect(tolower(path), "\\.gz$")) {
      con <- gzfile(path)
      if (ext == "csv") {
        read_csv(con, guess_max = 2000, show_col_types = FALSE, progress = FALSE)
      } else {
        read_tsv(con, guess_max = 2000, show_col_types = FALSE, progress = FALSE)
      }
    } else {
      if (ext == "csv") {
        read_csv(path, guess_max = 2000, show_col_types = FALSE, progress = FALSE)
      } else {
        read_tsv(path, guess_max = 2000, show_col_types = FALSE, progress = FALSE)
      }
    }
  }, error = function(e) NULL)
  
  df
}

# ------------------------------------------------------------------------
# 5. Helper: search a set of column names for topic-related keywords
# ------------------------------------------------------------------------

search_terms <- c(
  "parkinson", "\\bpd\\b", "control", "healthy", "early", "late",
  "stage", "donor", "sample", "cell.?type", "cluster", "\\bnk\\b"
)

search_metadata_terms <- function(col_names) {
  if (length(col_names) == 0) return(character(0))
  hits <- lapply(search_terms, function(term) {
    matched <- col_names[str_detect(tolower(col_names), term)]
    if (length(matched) > 0) {
      paste0(term, " -> ", paste(matched, collapse = ", "))
    } else {
      NA_character_
    }
  })
  hits <- unlist(hits)
  hits[!is.na(hits)]
}

# ==============================================================================
# 6. XIONG ET AL. 2024 (PBMC scRNA-seq, PRIMARY DISCOVERY)
# ==============================================================================

log_msg("-----------------------------------------------------------------")
log_msg("XIONG ET AL. 2024 - INSPECTION")
log_msg("-----------------------------------------------------------------")

xiong_files <- inventory_list[["xiong"]]
xiong_status <- NOT_VERIFIED
xiong_expr_status <- NOT_VERIFIED
xiong_cell_meta_status <- NOT_VERIFIED
xiong_donor_meta_status <- NOT_VERIFIED
xiong_disease_labels_status <- NOT_VERIFIED
xiong_stage_labels_status <- NOT_VERIFIED
xiong_nk_annotation_status <- NOT_VERIFIED
xiong_metadata_columns_df <- data.frame(column_name = character(0), source_file = character(0))
xiong_is_cell_level <- NOT_VERIFIED
xiong_has_cell_ids <- NOT_VERIFIED
xiong_has_donor_ids <- NOT_VERIFIED

if (is.null(xiong_files) || nrow(xiong_files) == 0) {
  log_msg("No files found in data/xiong/. Status: INCOMPLETE.")
  xiong_status <- "INCOMPLETE - no files present"
} else {
  
  xiong_files$role_guess <- vapply(xiong_files$filename, classify_file, character(1))
  log_msg("Files found: ", nrow(xiong_files))
  for (i in seq_len(nrow(xiong_files))) {
    log_msg("  - ", xiong_files$relative_path[i],
            " | ext=", xiong_files$extension[i],
            " | size=", xiong_files$size_human[i],
            " | role_guess=", xiong_files$role_guess[i])
  }
  
  # Detect GEO series-level-only files vs. actual cell-level matrices
  has_series_matrix <- any(str_detect(tolower(xiong_files$filename), "series_matrix|\\.soft"))
  has_10x_triplet <- any(str_detect(tolower(xiong_files$filename), "matrix\\.mtx")) &&
    any(str_detect(tolower(xiong_files$filename), "barcodes")) &&
    any(str_detect(tolower(xiong_files$filename), "features|genes"))
  has_h5 <- any(str_detect(tolower(xiong_files$filename), "\\.h5$|\\.h5ad$|\\.hdf5$"))
  has_rds_seurat <- any(str_detect(tolower(xiong_files$filename), "\\.rds$|seurat"))
  has_large_matrix_csv_tsv <- any(
    xiong_files$extension %in% c("csv", "tsv", "txt", "gz") &
      str_detect(tolower(xiong_files$filename), "count|matrix|expr")
  )
  
  if (has_10x_triplet || has_h5 || has_rds_seurat) {
    xiong_is_cell_level <- "LIKELY cell-level scRNA-seq data present (10x/h5/h5ad/Seurat files detected)"
    xiong_expr_status <- "LIKELY PRESENT (cell-level format detected)"
  } else if (has_large_matrix_csv_tsv) {
    xiong_is_cell_level <- "POSSIBLE cell-level or pseudobulk matrix (plain-text matrix-like file detected; format not confirmed) - NOT VERIFIED"
    xiong_expr_status <- "POSSIBLY PRESENT - format not confirmed, see inspection_report.txt"
  } else if (has_series_matrix && !has_10x_triplet && !has_h5 && !has_rds_seurat) {
    xiong_is_cell_level <- "LIKELY GEO SERIES-LEVEL / PROCESSED SUMMARY DATA ONLY - actual cell-level scRNA-seq matrix NOT confirmed present"
    xiong_expr_status <- "NOT CONFIRMED - only series-level/summary files detected"
  } else {
    xiong_is_cell_level <- NOT_VERIFIED
    xiong_expr_status <- NOT_VERIFIED
  }
  
  # Try to read any metadata-like delimited files for column inspection
  candidate_meta_files <- xiong_files %>%
    filter(extension %in% delim_exts | str_detect(tolower(filename), "\\.gz$")) %>%
    filter(str_detect(tolower(role_guess), "metadata|donor|cell"))
  
  if (nrow(candidate_meta_files) == 0) {
    # fall back: any delimited file, since role guess may miss metadata files
    candidate_meta_files <- xiong_files %>%
      filter(extension %in% delim_exts | str_detect(tolower(filename), "\\.gz$"))
  }
  
  all_cols <- character(0)
  col_source <- character(0)
  
  for (i in seq_len(nrow(candidate_meta_files))) {
    fpath <- candidate_meta_files$full_path[i]
    preview <- preview_delim_file(fpath, n_preview = 5)
    if (isTRUE(preview$ok)) {
      log_msg("  Preview of ", candidate_meta_files$relative_path[i], ":")
      log_msg("    Header: ", preview$header)
      log_msg("    Guessed delimiter: '", preview$delimiter_guess, "' | Guessed n_cols: ", preview$n_cols_guess)
    }
    
    df <- try_read_metadata_table(fpath)
    if (!is.null(df)) {
      cols <- colnames(df)
      all_cols <- c(all_cols, cols)
      col_source <- c(col_source, rep(candidate_meta_files$relative_path[i], length(cols)))
      
      if (str_detect(tolower(candidate_meta_files$filename[i]), "cell")) {
        xiong_cell_meta_status <- paste0("PRESENT (", candidate_meta_files$relative_path[i], ", ",
                                         nrow(df), " rows x ", ncol(df), " cols)")
      }
      if (str_detect(tolower(candidate_meta_files$filename[i]), "donor|subject|sample")) {
        xiong_donor_meta_status <- paste0("PRESENT (", candidate_meta_files$relative_path[i], ", ",
                                          nrow(df), " rows x ", ncol(df), " cols)")
      }
    }
  }
  
  if (length(all_cols) > 0) {
    xiong_metadata_columns_df <- data.frame(
      column_name = all_cols,
      source_file = col_source,
      stringsAsFactors = FALSE
    ) %>% distinct()
    
    hits <- search_metadata_terms(unique(all_cols))
    if (length(hits) > 0) {
      log_msg("  Metadata keyword search hits:")
      for (h in hits) log_msg("    ", h)
    } else {
      log_msg("  Metadata keyword search: no matches found in available column names.")
    }
    
    # Disease / stage / NK status based on literal column-name matches only
    cols_lower <- tolower(all_cols)
    xiong_disease_labels_status <- if (any(str_detect(cols_lower, "parkinson|\\bpd\\b|disease|diagnosis|condition|group"))) {
      "POSSIBLE - matching column name(s) found (verify values manually)"
    } else {
      "NOT VERIFIED - no matching column name found"
    }
    xiong_stage_labels_status <- if (any(str_detect(cols_lower, "stage|early|late|hoehn|updrs"))) {
      "POSSIBLE - matching column name(s) found (verify values manually)"
    } else {
      "NOT VERIFIED - no matching column name found"
    }
    xiong_nk_annotation_status <- if (any(str_detect(cols_lower, "\\bnk\\b|cell.?type|celltype|cluster|annotation"))) {
      "POSSIBLE - matching column name(s) found (verify values manually; does not confirm NK label values exist)"
    } else {
      "NOT VERIFIED - no matching column name found"
    }
    xiong_has_cell_ids <- if (any(str_detect(cols_lower, "barcode|cell.?id|cellid"))) "LIKELY PRESENT" else NOT_VERIFIED
    xiong_has_donor_ids <- if (any(str_detect(cols_lower, "donor|subject|patient|sample.?id"))) "LIKELY PRESENT" else NOT_VERIFIED
  } else {
    log_msg("  No metadata-like delimited files could be read for Xiong.")
  }
  
  write_csv(xiong_metadata_columns_df, file.path(out_dir, "xiong_metadata_columns.csv"))
  
  # Overall status
  xiong_status <- if (xiong_expr_status != NOT_VERIFIED &&
                      !str_detect(xiong_expr_status, "NOT CONFIRMED")) {
    "READY FOR FURTHER INSPECTION (expression-like data present; see details)"
  } else {
    "INCOMPLETE / NEEDS VERIFICATION - see inspection_report.txt"
  }
}

log_msg("")

# ==============================================================================
# 7. MOQUIN-BEAUDRY ET AL. 2024/2025 (PBMC scRNA-seq, INDEPENDENT REPLICATION)
# ==============================================================================

log_msg("-----------------------------------------------------------------")
log_msg("MOQUIN-BEAUDRY ET AL. 2024/2025 - INSPECTION")
log_msg("-----------------------------------------------------------------")

moquin_files <- inventory_list[["moquin_beaudry"]]
moquin_status <- NOT_VERIFIED
moquin_expr_status <- NOT_VERIFIED
moquin_cell_meta_status <- NOT_VERIFIED
moquin_donor_meta_status <- NOT_VERIFIED
moquin_disease_labels_status <- NOT_VERIFIED
moquin_stage_labels_status <- NOT_VERIFIED
moquin_nk_annotation_status <- NOT_VERIFIED
moquin_metadata_columns_df <- data.frame(column_name = character(0), source_file = character(0))
moquin_format_guess <- NOT_VERIFIED

if (is.null(moquin_files) || nrow(moquin_files) == 0) {
  log_msg("No files found in data/moquin_beaudry/. Status: INCOMPLETE.")
  moquin_status <- "INCOMPLETE - no files present"
} else {
  
  moquin_files$role_guess <- vapply(moquin_files$filename, classify_file, character(1))
  log_msg("Files found: ", nrow(moquin_files))
  for (i in seq_len(nrow(moquin_files))) {
    log_msg("  - ", moquin_files$relative_path[i],
            " | ext=", moquin_files$extension[i],
            " | size=", moquin_files$size_human[i],
            " | role_guess=", moquin_files$role_guess[i])
  }
  
  fnames_lower <- tolower(moquin_files$filename)
  
  is_10x     <- any(str_detect(fnames_lower, "matrix\\.mtx")) &&
    any(str_detect(fnames_lower, "barcodes")) &&
    any(str_detect(fnames_lower, "features|genes"))
  is_h5ad    <- any(str_detect(fnames_lower, "\\.h5ad$"))
  is_h5      <- any(str_detect(fnames_lower, "\\.h5$|\\.hdf5$")) && !is_h5ad
  is_seurat  <- any(str_detect(fnames_lower, "\\.rds$")) | any(str_detect(fnames_lower, "seurat"))
  is_sra     <- any(str_detect(fnames_lower, "\\.sra$|\\.fastq|\\.fq(\\.gz)?$"))
  is_proc_mx <- any(str_detect(fnames_lower, "count|matrix|expr")) &
    any(moquin_files$extension %in% c("csv", "tsv", "txt"))
  
  moquin_format_guess <- if (is_10x) {
    "10x matrix files (matrix.mtx + barcodes + features/genes)"
  } else if (is_h5ad) {
    "h5ad (AnnData) file"
  } else if (is_h5) {
    "h5/hdf5 file (format/schema not further verified)"
  } else if (is_seurat) {
    "Likely Seurat object (.rds) - content not verified without loading in R with Seurat installed"
  } else if (is_sra) {
    "SRA/FASTQ raw sequencing reads (no processed count matrix confirmed)"
  } else if (is_proc_mx) {
    "Likely processed count matrix in delimited text format (not confirmed raw vs normalized)"
  } else {
    "OTHER / NOT VERIFIED - could not confidently classify based on filenames"
  }
  log_msg("  Format classification: ", moquin_format_guess)
  
  moquin_expr_status <- if (is_10x || is_h5ad || is_h5 || is_seurat || is_proc_mx) {
    paste0("LIKELY PRESENT (", moquin_format_guess, ")")
  } else if (is_sra) {
    "RAW READS ONLY - no processed expression matrix confirmed"
  } else {
    NOT_VERIFIED
  }
  
  candidate_meta_files <- moquin_files %>%
    filter(extension %in% delim_exts | str_detect(tolower(filename), "\\.gz$")) %>%
    filter(str_detect(tolower(role_guess), "metadata|donor|cell"))
  
  if (nrow(candidate_meta_files) == 0) {
    candidate_meta_files <- moquin_files %>%
      filter(extension %in% delim_exts | str_detect(tolower(filename), "\\.gz$"))
  }
  
  all_cols <- character(0)
  col_source <- character(0)
  
  for (i in seq_len(nrow(candidate_meta_files))) {
    fpath <- candidate_meta_files$full_path[i]
    preview <- preview_delim_file(fpath, n_preview = 5)
    if (isTRUE(preview$ok)) {
      log_msg("  Preview of ", candidate_meta_files$relative_path[i], ":")
      log_msg("    Header: ", preview$header)
    }
    
    df <- try_read_metadata_table(fpath)
    if (!is.null(df)) {
      cols <- colnames(df)
      all_cols <- c(all_cols, cols)
      col_source <- c(col_source, rep(candidate_meta_files$relative_path[i], length(cols)))
      
      if (str_detect(tolower(candidate_meta_files$filename[i]), "cell")) {
        moquin_cell_meta_status <- paste0("PRESENT (", candidate_meta_files$relative_path[i], ", ",
                                          nrow(df), " rows x ", ncol(df), " cols)")
      }
      if (str_detect(tolower(candidate_meta_files$filename[i]), "donor|subject|sample")) {
        moquin_donor_meta_status <- paste0("PRESENT (", candidate_meta_files$relative_path[i], ", ",
                                           nrow(df), " rows x ", ncol(df), " cols)")
      }
    }
  }
  
  if (length(all_cols) > 0) {
    moquin_metadata_columns_df <- data.frame(
      column_name = all_cols,
      source_file = col_source,
      stringsAsFactors = FALSE
    ) %>% distinct()
    
    hits <- search_metadata_terms(unique(all_cols))
    if (length(hits) > 0) {
      log_msg("  Metadata keyword search hits:")
      for (h in hits) log_msg("    ", h)
    } else {
      log_msg("  Metadata keyword search: no matches found in available column names.")
    }
    
    cols_lower <- tolower(all_cols)
    moquin_disease_labels_status <- if (any(str_detect(cols_lower, "parkinson|\\bpd\\b|disease|diagnosis|condition|group"))) {
      "POSSIBLE - matching column name(s) found (verify values manually)"
    } else {
      "NOT VERIFIED - no matching column name found"
    }
    moquin_stage_labels_status <- if (any(str_detect(cols_lower, "stage|early|late|hoehn|updrs"))) {
      "POSSIBLE - matching column name(s) found (verify values manually)"
    } else {
      "NOT VERIFIED - no matching column name found"
    }
    moquin_nk_annotation_status <- if (any(str_detect(cols_lower, "\\bnk\\b|cell.?type|celltype|cluster|annotation"))) {
      "POSSIBLE - matching column name(s) found (verify values manually)"
    } else {
      "NOT VERIFIED - no matching column name found"
    }
  } else {
    log_msg("  Note: metadata for h5ad/h5/Seurat objects, if present, cannot be listed here without")
    log_msg("        loading them via appropriate packages (anndata/hdf5r/Seurat). Not attempted in this script.")
  }
  
  write_csv(moquin_metadata_columns_df, file.path(out_dir, "moquin_metadata_columns.csv"))
  
  moquin_status <- if (moquin_expr_status != NOT_VERIFIED && !str_detect(moquin_expr_status, "RAW READS ONLY")) {
    "READY FOR FURTHER INSPECTION (expression-like data present; see details)"
  } else {
    "INCOMPLETE / NEEDS VERIFICATION - see inspection_report.txt"
  }
}

log_msg("")

# ==============================================================================
# 8. GSE100054 (PBMC bulk transcriptomics, INDEPENDENT PBMC VALIDATION)
# ==============================================================================

log_msg("-----------------------------------------------------------------")
log_msg("GSE100054 - INSPECTION")
log_msg("-----------------------------------------------------------------")

gse_files <- inventory_list[["gse100054"]]
gse_status <- NOT_VERIFIED
gse_expr_status <- NOT_VERIFIED
gse_sample_meta_status <- NOT_VERIFIED
gse_disease_labels_status <- NOT_VERIFIED
gse_platform_status <- NOT_VERIFIED
gse_expr_type_status <- NOT_VERIFIED
gse_metadata_columns_df <- data.frame(column_name = character(0), source_file = character(0))
gse_dims <- NOT_VERIFIED
gse_genes_as_rows <- NOT_VERIFIED

if (is.null(gse_files) || nrow(gse_files) == 0) {
  log_msg("No files found in data/gse100054/. Status: INCOMPLETE.")
  gse_status <- "INCOMPLETE - no files present"
} else {
  
  gse_files$role_guess <- vapply(gse_files$filename, classify_file, character(1))
  log_msg("Files found: ", nrow(gse_files))
  for (i in seq_len(nrow(gse_files))) {
    log_msg("  - ", gse_files$relative_path[i],
            " | ext=", gse_files$extension[i],
            " | size=", gse_files$size_human[i],
            " | role_guess=", gse_files$role_guess[i])
  }
  
  fnames_lower <- tolower(gse_files$filename)
  gse_platform_status <- if (any(str_detect(fnames_lower, "gpl"))) {
    matched <- gse_files$filename[str_detect(fnames_lower, "gpl")]
    paste0("Platform-related filename(s) detected: ", paste(matched, collapse = ", "),
           " (platform ID not fabricated - inspect file contents directly to confirm)")
  } else {
    "NOT VERIFIED - no platform/GPL-named file detected"
  }
  
  # Identify the most likely expression matrix file: largest matrix/series file
  matrix_candidates <- gse_files %>%
    filter(str_detect(tolower(role_guess), "expression_matrix|series_level")) 
  if (nrow(matrix_candidates) == 0) {
    matrix_candidates <- gse_files %>% filter(extension %in% c(delim_exts, "gz"))
  }
  
  expr_file_used <- NA_character_
  
  if (nrow(matrix_candidates) > 0) {
    # pick the largest as most likely to be the full matrix
    matrix_candidates <- matrix_candidates %>% arrange(desc(size_bytes))
    top_candidate <- matrix_candidates[1, ]
    expr_file_used <- top_candidate$relative_path
    
    preview <- preview_delim_file(top_candidate$full_path, n_preview = 5)
    if (isTRUE(preview$ok)) {
      log_msg("  Preview of largest candidate matrix file (", top_candidate$relative_path, "):")
      log_msg("    Header: ", str_trunc_safe(preview$header))
      log_msg("    Guessed n_cols (from header): ", preview$n_cols_guess)
      
      df <- try_read_metadata_table(top_candidate$full_path)
      if (!is.null(df)) {
        gse_dims <- paste0(nrow(df), " rows x ", ncol(df), " cols")
        log_msg("    Parsed dimensions: ", gse_dims)
        
        # Heuristic: if first column values look like gene symbols/IDs and are mostly unique,
        # and column count >> row count relative to typical sample counts, likely genes x samples.
        first_col_vals <- df[[1]]
        looks_like_gene_ids <- is.character(first_col_vals) &&
          mean(str_detect(first_col_vals[1:min(20, length(first_col_vals))], "^[A-Za-z0-9_.\\-]+$"), na.rm = TRUE) > 0.7
        
        if (looks_like_gene_ids && ncol(df) < nrow(df)) {
          gse_genes_as_rows <- "LIKELY YES - first column looks like gene identifiers; rows >> columns"
        } else if (looks_like_gene_ids) {
          gse_genes_as_rows <- "POSSIBLE - first column looks like gene identifiers, but row/column count relationship inconclusive"
        } else {
          gse_genes_as_rows <- NOT_VERIFIED
        }
        
        gse_expr_status <- paste0("LIKELY PRESENT (", expr_file_used, ", ", gse_dims, ")")
        
        # crude heuristic for expression type based on value range of a numeric subsample
        numeric_cols <- df %>% select(where(is.numeric))
        if (ncol(numeric_cols) > 0) {
          sample_vals <- unlist(numeric_cols[1:min(50, nrow(numeric_cols)), ], use.names = FALSE)
          sample_vals <- sample_vals[!is.na(sample_vals)]
          if (length(sample_vals) > 0) {
            has_negative <- any(sample_vals < 0)
            has_non_integer <- any(sample_vals != floor(sample_vals))
            max_val <- max(sample_vals)
            
            gse_expr_type_status <- if (has_negative) {
              "LIKELY LOG-RATIO / NORMALIZED MICROARRAY INTENSITY (negative values present) - NOT VERIFIED further"
            } else if (!has_non_integer && max_val > 20) {
              "POSSIBLE RAW COUNTS (non-negative integers, wide range) - NOT VERIFIED further"
            } else if (has_non_integer) {
              "POSSIBLE NORMALIZED EXPRESSION OR MICROARRAY INTENSITY (non-integer values) - NOT VERIFIED further"
            } else {
              NOT_VERIFIED
            }
          }
        }
      } else {
        log_msg("    File too large or unreadable for full parse; dimensions NOT VERIFIED without a bounded/chunked read.")
      }
    }
  }
  
  # Metadata columns (sample-level) search across delimited files
  candidate_meta_files <- gse_files %>%
    filter(extension %in% delim_exts | str_detect(tolower(filename), "\\.gz$")) %>%
    filter(str_detect(tolower(role_guess), "metadata"))
  
  if (nrow(candidate_meta_files) == 0) {
    candidate_meta_files <- gse_files %>%
      filter(extension %in% delim_exts | str_detect(tolower(filename), "\\.gz$"))
  }
  
  all_cols <- character(0)
  col_source <- character(0)
  
  for (i in seq_len(nrow(candidate_meta_files))) {
    fpath <- candidate_meta_files$full_path[i]
    df <- try_read_metadata_table(fpath)
    if (!is.null(df)) {
      cols <- colnames(df)
      all_cols <- c(all_cols, cols)
      col_source <- c(col_source, rep(candidate_meta_files$relative_path[i], length(cols)))
      gse_sample_meta_status <- paste0("PRESENT (", candidate_meta_files$relative_path[i], ", ",
                                       nrow(df), " rows x ", ncol(df), " cols)")
    }
  }
  
  if (length(all_cols) > 0) {
    gse_metadata_columns_df <- data.frame(
      column_name = all_cols,
      source_file = col_source,
      stringsAsFactors = FALSE
    ) %>% distinct()
    
    hits <- search_metadata_terms(unique(all_cols))
    if (length(hits) > 0) {
      log_msg("  Metadata keyword search hits:")
      for (h in hits) log_msg("    ", h)
    } else {
      log_msg("  Metadata keyword search: no matches found in available column names.")
    }
    
    cols_lower <- tolower(all_cols)
    gse_disease_labels_status <- if (any(str_detect(cols_lower, "parkinson|\\bpd\\b|disease|diagnosis|condition|group"))) {
      "POSSIBLE - matching column name(s) found (verify values manually)"
    } else {
      "NOT VERIFIED - no matching column name found"
    }
  } else {
    log_msg("  No separate sample-metadata delimited file identified/readable for GSE100054.")
    log_msg("  Sample metadata may be embedded in a GEO series matrix file header (not parsed by this script).")
  }
  
  write_csv(gse_metadata_columns_df, file.path(out_dir, "gse100054_metadata_columns.csv"))
  
  gse_status <- if (gse_expr_status != NOT_VERIFIED) {
    "READY FOR FURTHER INSPECTION (expression-like data present; see details)"
  } else {
    "INCOMPLETE / NEEDS VERIFICATION - see inspection_report.txt"
  }
}

log_msg("")

# ==============================================================================
# 9. Write full text report
# ==============================================================================

report_path <- file.path(out_dir, "inspection_report.txt")
writeLines(log_lines, report_path)

# ==============================================================================
# 10. Console summary
# ==============================================================================

cat("\n")
cat("=================================================================\n")
cat("DATASET INSPECTION SUMMARY\n")
cat("=================================================================\n\n")

cat("Xiong:\n")
cat("  Status:            ", xiong_status, "\n")
cat("  Expression data:   ", xiong_expr_status, "\n")
cat("  Cell metadata:     ", xiong_cell_meta_status, "\n")
cat("  Donor metadata:    ", xiong_donor_meta_status, "\n")
cat("  Disease labels:    ", xiong_disease_labels_status, "\n")
cat("  Stage labels:      ", xiong_stage_labels_status, "\n")
cat("  NK annotation:     ", xiong_nk_annotation_status, "\n")
cat("  Cell-level data?:  ", xiong_is_cell_level, "\n\n")

cat("Moquin-Beaudry:\n")
cat("  Status:            ", moquin_status, "\n")
cat("  Format guess:      ", moquin_format_guess, "\n")
cat("  Expression data:   ", moquin_expr_status, "\n")
cat("  Cell metadata:     ", moquin_cell_meta_status, "\n")
cat("  Donor metadata:    ", moquin_donor_meta_status, "\n")
cat("  Disease labels:    ", moquin_disease_labels_status, "\n")
cat("  Stage labels:      ", moquin_stage_labels_status, "\n")
cat("  NK annotation:     ", moquin_nk_annotation_status, "\n\n")

cat("GSE100054:\n")
cat("  Status:            ", gse_status, "\n")
cat("  Expression data:   ", gse_expr_status, "\n")
cat("  Dimensions:        ", gse_dims, "\n")
cat("  Genes as rows?:    ", gse_genes_as_rows, "\n")
cat("  Sample metadata:   ", gse_sample_meta_status, "\n")
cat("  Disease labels:    ", gse_disease_labels_status, "\n")
cat("  Platform:          ", gse_platform_status, "\n")
cat("  Expression type:   ", gse_expr_type_status, "\n\n")

cat("Full report written to: results/data_inspection/inspection_report.txt\n")
cat("File inventory written to: results/data_inspection/dataset_file_inventory.csv\n")
cat("Metadata column listings written to:\n")
cat("  - results/data_inspection/xiong_metadata_columns.csv\n")
cat("  - results/data_inspection/moquin_metadata_columns.csv\n")
cat("  - results/data_inspection/gse100054_metadata_columns.csv\n")
cat("=================================================================\n")

