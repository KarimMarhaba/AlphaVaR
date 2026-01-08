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

#' Export HTML Report as PDF
#'
#' Converts the interactive HTML report into a static PDF file using headless Chrome.
#' This ensures that the layout matches the HTML view exactly.
#'
#' @param obj An object of class \code{AlphaVarSet}.
#' @param output_file Character string. Path to the output PDF file.
#' @param open Logical. If \code{TRUE}, opens the PDF after creation. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{av_create_report} or internal rendering.
#'
#' @details
#' This function requires the \code{pagedown} package and a strictly local installation 
#' of Chrome or Chromium. It generates a temporary HTML report first, then prints it to PDF.
#'
#' @return The path to the generated PDF file (invisibly).
#' @export
#'
#' @examples
#' \dontrun{
#'   my_data <- av_import_csv("variants.csv")
#'   av_export_pdf(my_data, "report.pdf")
#' }
av_export_pdf <- function(obj, output_file, open = TRUE, ...) {
  # 1. Check Class
  if (!inherits(obj, "AlphaVarSet")) {
    stop("Input must be of class 'AlphaVarSet'")
  }
  
  # 2. Check Dependency
  if (!requireNamespace("pagedown", quietly = TRUE)) {
    stop("Package 'pagedown' is required for PDF export. Please install it using install.packages('pagedown').")
  }
  
  # 3. Validation of output path
  if (!grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
    warning("output_file does not end with .pdf. Appending extension.")
    output_file <- paste0(output_file, ".pdf")
  }
  
  message("Generating intermediate HTML report...")
  
  # 4. Create temporary HTML using the existing logic
  # We use a temp file for the HTML so we don't clutter the user's directory
  temp_html <- tempfile(fileext = ".html")
  
  tryCatch({
    av_create_report(obj, output_file = temp_html, open = FALSE, ...)
  }, error = function(e) {
    stop("Failed to generate intermediate HTML report: ", e$message)
  })
  
  message("Converting HTML to PDF via pagedown (this requires Chrome/Chromium)...")
  
  # 5. Convert to PDF
  tryCatch({
    pagedown::chrome_print(input = temp_html, output = output_file)
  }, error = function(e) {
    stop("PDF conversion failed. Ensure Chrome is installed. Error: ", e$message)
  })
  
  # 6. Cleanup
  if (file.exists(temp_html)) {
    unlink(temp_html)
  }
  
  message("PDF Report successfully created: ", output_file)
  
  # 7. Open if requested
  if (open) {
    utils::browseURL(output_file)
  }
  
  invisible(output_file)
}