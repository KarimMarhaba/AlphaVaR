test_that("Report generation works", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("DT")
  skip_if_not_installed("plotly")
  
  if (!rmarkdown::pandoc_available()) {
    skip("Pandoc not found - skipping report generation test.")
  }

  # Mock Data
  dat <- dplyr::tibble(
    variant_id = 1:20, score = runif(20), modality = "ATAC", gene = paste0("G", 1:20)
  )
  obj <- AlphaVarSet(dat)
  
  tmp_file <- tempfile(fileext = ".html")
  
  # --- ROBUST TEMPLATE FINDER ---
  # Versuche verschiedene Pfade, je nachdem wo der Test läuft
  candidates <- c(
    "../../inst/rmd/report_master.Rmd",  # testthat interactive
    "inst/rmd/report_master.Rmd",        # package root
    system.file("rmd", "report_master.Rmd", package = "AlphaVariantR") # installed
  )
  
  # Nimm den ersten Pfad, der existiert
  valid_paths <- Filter(file.exists, candidates)
  
  if (length(valid_paths) == 0) {
    skip("Template file not found in test environment.")
  }
  
  local_template <- valid_paths[1]
  
  # Run function
  out <- av_create_report(obj, output_file = tmp_file, open = FALSE, custom_template = local_template)
  
  expect_true(file.exists(tmp_file))
  
  unlink(tmp_file)
})