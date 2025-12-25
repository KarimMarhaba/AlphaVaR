test_that("Visualization functions return ggplots", {
  # Mock Data
  dat <- dplyr::tibble(
    variant_id = 1:20,
    score = runif(20),
    modality = rep(c("A", "B"), 10)
  )
  obj <- AlphaVarSet(dat)
  
  # Test Ridge (without stats)
  p1 <- av_plot_ridges(obj)
  expect_true(inherits(p1, "ggplot"))
  
  # Test ECDF
  p2 <- av_plot_ecdf(obj)
  expect_true(inherits(p2, "ggplot"))
  
  # Test Ranking
  p3 <- av_plot_ranking(obj)
  expect_true(inherits(p3, "ggplot"))
})

test_that("Volcano plot handles missing stats", {
  dat <- dplyr::tibble(
    variant_id = 1:10, score = 0.5, modality = "A"
  )
  obj <- AlphaVarSet(dat)
  
  # Must fail because av_calc_enrichment wasn't run
  expect_error(av_plot_volcano(obj), "Please run av_calc_enrichment")
  
  # Run stats and try again
  obj <- av_calc_enrichment(obj)
  p <- av_plot_volcano(obj)
  expect_true(inherits(p, "ggplot"))
})