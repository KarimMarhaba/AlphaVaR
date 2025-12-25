#' Create Analysis Report
#'
#' Generates a comprehensive HTML report from an \code{AlphaVarSet} object.
#' Automatically runs core statistics if they haven't been computed yet.
#'
#' @param obj An \code{AlphaVarSet} object.
#' @param output_file Path where the HTML file should be saved.
#' @param open Logical. Open the file in browser after creating?
#' @param custom_template Optional path to a specific Rmd template.
#'        Useful for testing or custom branding.
#'
#' @return The path to the created file (invisibly).
#' @importFrom rmarkdown render
#' @importFrom tools file_path_as_absolute
#' @export
av_create_report <- function(obj, output_file = "AlphaVariant_Report.html", open = TRUE, custom_template = NULL) {
  validate_AlphaVarSet(obj)
  
  # 1. Check if Stats are present, if not, compute them on the fly
  if (is.null(obj$stats$enrichment)) {
    message("Statistics not found. Computing now (defaults)...")
    obj <- av_calc_enrichment(obj)
  }
  
  # 2. Locate Template
  if (!is.null(custom_template)) {
    template_path <- custom_template
  } else {
    # Try system install location first
    template_path <- system.file("rmd", "report_master.Rmd", package = "AlphaVaR")
    
    # Fallback for development (if package not installed yet)
    if (template_path == "") {
      # Check typical relative locations during dev
      candidates <- c(
        "inst/rmd/report_master.Rmd",       # Working dir = Package Root
        "../../inst/rmd/report_master.Rmd"  # Working dir = tests/testthat
      )
      
      found <- Filter(file.exists, candidates)
      if (length(found) > 0) {
        template_path <- found[1]
      }
    }
  }
  
  if (template_path == "" || !file.exists(template_path)) {
    stop("Report template not found. Try reinstalling the package or providing 'custom_template'.")
  }
  
  message(paste("Rendering report to:", output_file))
  
  # 3. Render
  rmarkdown::render(
    input = template_path,
    output_file = basename(output_file), # Render in local dir first
    output_dir = dirname(output_file),
    params = list(obj = obj),
    quiet = TRUE
  )
  
  # 4. Open
  if (open) {
    if (interactive()) {
      utils::browseURL(output_file)
    } else {
      message("Report created. Open manually.")
    }
  }
  
  invisible(output_file)
}