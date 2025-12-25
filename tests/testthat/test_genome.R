test_that("Coordinate parsing works", {
  # Test standard format
  vec <- c("chr1:100-200", "X:500-600")
  parsed <- parse_coords(vec)
  
  expect_equal(parsed$chrom[1], "chr1")
  expect_equal(parsed$start[1], 100)
  expect_equal(parsed$end[1], 200)
  
  expect_equal(parsed$chrom[2], "X")
})

test_that("Locus plot generation works", {
  # Mock data with coordinates
  dat <- dplyr::tibble(
    variant_id = c("v1", "v2"),
    score = c(0.8, 0.5),
    modality = c("ATAC", "ATAC"),
    gene = c("TargetGene", "TargetGene"),
    coords = c("chr1:1000-1010", "chr1:1050-1060") # This column is vital now
  )
  
  obj <- AlphaVarSet(dat)
  
  p <- av_plot_locus(obj, gene_name = "TargetGene")
  expect_true(inherits(p, "ggplot"))
  
  # Test error if gene missing
  expect_error(av_plot_locus(obj, "WrongGene"), "No variants found")
})

test_that("Locus plot catches missing coord column", {
  # Mock data WITHOUT coordinates
  dat <- dplyr::tibble(
    variant_id = "v1", score = 0.5, modality = "A", gene = "G"
  )
  obj <- AlphaVarSet(dat)
  
  expect_error(av_plot_locus(obj, "G"), "No 'coords' column found")
})