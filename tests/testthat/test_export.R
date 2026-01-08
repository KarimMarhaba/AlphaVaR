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