#' Import AlphaGenome Results
#'
#' Reads a CSV file and creates an AlphaVarSet object.
#' It standardizes column names based on a provided mapping.
#'
#' @param file Path to the CSV file.
#' @param col_map A named list mapping internal names to CSV headers. 
#'        Defaults to: list(score = "quantile_score", modality = "output_type", 
#'        variant_id = "variant_id", gene = "gene_name").
#' @param sep Separator for the file (default ",").
#'
#' @return An \code{AlphaVarSet} object.
#' @importFrom readr read_delim cols
#' @importFrom dplyr rename select mutate any_of
#' @importFrom rlang inform abort warn
#' @importFrom utils modifyList
#' @importFrom stats setNames
#' @export
av_import_csv <- function(file, col_map = list(), sep = ",") {
  
  # 1. Define Default Mapping
  default_map <- list(
    score      = "quantile_score",
    modality   = "output_type",
    variant_id = "variant_id",
    gene       = "gene_name",
    coords     = "scored_interval" 
  )
  
  # Merge user map with defaults
  final_map <- utils::modifyList(default_map, col_map)
  
  # 2. Check File
  if (!file.exists(file)) {
    rlang::abort(paste("File not found:", file))
  }
  
  rlang::inform(paste("Reading", basename(file), "..."))
  
  # 3. Read Data
  raw_data <- readr::read_delim(file, delim = sep, show_col_types = FALSE)
  
  # 4. Validate Mandatory Columns in CSV
  csv_cols <- names(raw_data)
  
  # We need to flip the map for checking: we look for the CSV names
  # map: score = "quantile_score" -> we look for "quantile_score"
  mandatory_internal <- c("score", "modality", "variant_id")
  mandatory_csv <- unlist(final_map[mandatory_internal])
  
  missing_csv_cols <- setdiff(mandatory_csv, csv_cols)
  
  if (length(missing_csv_cols) > 0) {
    rlang::abort(paste0(
      "The CSV is missing required columns based on your mapping.\n",
      "Missing: ", paste(missing_csv_cols, collapse = ", "), "\n",
      "Found in CSV: ", paste(head(csv_cols, 5), collapse = ", ")
    ))
  }
  
  # 5. Perform Renaming
  # dplyr::rename takes new = old.
  # Our final_map is list(new = old). We convert to named vector.
  rename_vec <- unlist(final_map)
  
  # Robust Renaming: Only rename columns that actually exist in the CSV
  # This prevents errors if optional columns (like coords) are missing
  rename_vec <- rename_vec[rename_vec %in% csv_cols]
  
  clean_data <- raw_data %>%
    dplyr::rename(any_of(rename_vec))
  
  # 6. Type Conversion (SAFE)
  # Now we are sure 'score' exists because we checked mandatory cols above
  if ("score" %in% names(clean_data)) {
    clean_data <- clean_data %>%
      dplyr::mutate(score = as.numeric(score))
  } else {
    rlang::abort("Renaming failed internally. 'score' column missing after rename.")
  }
    
  if (any(is.na(clean_data$score))) {
    n_na <- sum(is.na(clean_data$score))
    rlang::warn(paste(n_na, "scores became NA. Check if score column is numeric."))
  }
  
  # 7. Create Object
  obj <- AlphaVarSet(
    data = clean_data,
    meta = list(
      filename = basename(file),
      imported_at = Sys.time(),
      col_map_used = final_map
    )
  )
  
  rlang::inform("Import successful.")
  return(obj)
}