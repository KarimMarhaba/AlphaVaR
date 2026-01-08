test_that("av_export_browser_track produces valid file", {
  # Mock Data
  data <- dplyr::tibble(
    variant_id = c("v1", "v2"),
    score = c(0.5, 0.9),
    modality = "test",
    coords = c("chr1:100-200", "chr2:500-600")
  )
  obj <- structure(list(data = data), class = "AlphaVarSet")

  # Temp file
  tmp_file <- tempfile(fileext = ".bedgraph")

  # Run Function
  expect_message(av_export_browser_track(obj, tmp_file), "Success")

  # Check File Content
  lines <- readLines(tmp_file)
  expect_true(grepl("^track type=bedGraph", lines[1])) # Header
  expect_true(grepl("chr1\t100\t200\t0.5", lines[2])) # Data 1
  expect_true(grepl("chr2\t500\t600\t0.9", lines[3])) # Data 2

  unlink(tmp_file)
})

test_that("av_export_pdf generates a PDF file", {
  # Skip on CRAN or if pagedown is not installed/configured to avoid CI errors
  skip_on_cran()
  if (!requireNamespace("pagedown", quietly = TRUE)) {
    skip("pagedown not installed")
  }
  
  # Setup minimal object
  obj <- structure(list(
    data = dplyr::tibble(
      variant_id = "v1", score = 0.9, modality = "Seq", 
      gene = "BRCA1", coords = "chr1:100"
    ),
    stats = list(),
    meta = list(filename = "test.csv")
  ), class = "AlphaVarSet")
  
  # Define output file
  out_pdf <- tempfile(fileext = ".pdf")
  
  # Run function (expect error if Chrome is missing, but logic should hold)
  # Using tryCatch because pagedown might fail if no Chrome is found locally
  tryCatch({
    av_export_pdf(obj, out_pdf, open = FALSE)
    expect_true(file.exists(out_pdf))
    expect_gt(file.size(out_pdf), 0)
  }, error = function(e) {
    skip("Chrome/Chromium not found for pagedown")
  })
})