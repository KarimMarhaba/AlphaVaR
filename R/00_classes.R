#' Constructor for AlphaVarSet objects
#'
#' Internal function to create the object. Use av_import_csv() instead.
#'
#' @param data A tibble containing the variant data.
#' @param meta A list containing metadata (filename, date).
#'
#' @return An object of class AlphaVarSet.
#' @keywords internal
new_AlphaVarSet <- function(data, meta = list()) {
  structure(
    list(
      data = data,
      meta = meta,
      stats = list() # Placeholder for future results
    ),
    class = "AlphaVarSet"
  )
}

#' Validator for AlphaVarSet objects
#'
#' Checks if the object is structurally sound.
#'
#' @param x An AlphaVarSet object.
#' @return The object if valid, error otherwise.
#' @keywords internal
validate_AlphaVarSet <- function(x) {
  # 1. Check class
  if (!inherits(x, "AlphaVarSet")) {
    stop("Object is not of class 'AlphaVarSet'", call. = FALSE)
  }
  
  # 2. Check essential internal columns
  # Whatever the input CSV looked like, INSIDE the package we expect these names:
  required_cols <- c("variant_id", "score", "modality")
  missing_cols <- setdiff(required_cols, names(x$data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Invalid AlphaVarSet: Missing internal columns: ", 
      paste(missing_cols, collapse = ", "),
      ". Did the import mapping fail?",
      call. = FALSE
    )
  }
  
  # 3. Check data types
  if (!is.numeric(x$data$score)) {
    stop("Invalid AlphaVarSet: 'score' column must be numeric.", call. = FALSE)
  }
  
  return(x)
}

#' Helper to create the object safely
#' @keywords internal
AlphaVarSet <- function(data, meta = list()) {
  validate_AlphaVarSet(new_AlphaVarSet(data, meta))
}

#' Print method for AlphaVarSet
#' 
#' Gives a nice summary in the console when typing the object name.
#' @param x An AlphaVarSet object.
#' @param ... Additional arguments.
#' @export
print.AlphaVarSet <- function(x, ...) {
  n_vars <- nrow(x$data)
  n_mods <- length(unique(x$data$modality))
  cat("=== AlphaVarSet Object ===\n")
  cat("Variants:  ", format(n_vars, big.mark = ","), "\n")
  cat("Modalities:", n_mods, "\n")
  cat("Source:    ", x$meta$filename %||% "Unknown", "\n")
  if (length(x$stats) > 0) {
    cat("Computed:  ", paste(names(x$stats), collapse = ", "), "\n")
  }
  invisible(x)
}