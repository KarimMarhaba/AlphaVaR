#' Launch AlphaVaR Interactive Explorer
#'
#' Launches a local Shiny application to visualize and explore variant scores interactively.
#' Users can upload CSV files, view volcano plots, and filter results.
#'
#' @return This function does not return; it interrupts R execution to run the app.
#' @export
#' @importFrom shiny runApp
#' @importFrom tools file_path_sans_ext
av_run_shiny <- function() {
  app_dir <- system.file("shiny", "alpha_app", package = "AlphaVaR")

  if (app_dir == "") {
    stop("Could not find the Shiny app directory. Try re-installing the package.")
  }

  shiny::runApp(app_dir, display.mode = "normal")
}