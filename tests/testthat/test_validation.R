# tests/testthat/test_validation.R

test_that("av_compare_scores calculates correlation and joins data correctly", {
  # 1. Setup Dummy Data
  dummy_data <- tibble::tibble(
    variant_id = c("var1", "var2", "var3", "var4", "var5"),
    score = c(0.9, 0.8, 0.1, 0.2, 0.5), 
    modality = "All",
    gene = "GeneA",
    coords = "chr1:100"
  )
  
  # We create the object manually if no exported constructor exists,
  # or use the constructor function. 
  # Since AlphaVarSet is a function according to the error message, we call it:
  obj <- AlphaVarSet(data = dummy_data) 
  
  # If AlphaVarSet() is not exported or is named differently (e.g., new_AlphaVarSet),
  # we can manually force the structure for the test (fallback):
  if (!exists("obj")) {
    obj <- structure(
      list(data = dummy_data, stats = list(), meta = list()),
      class = "AlphaVarSet"
    )
  }
  
  # 2. External Data Setup
  external_scores <- tibble::tibble(
    var_id = c("var1", "var2", "var3", "var4", "var6"), 
    cadd_phred = c(20, 18, 2, 5, 10)
  )
  
  # 3. Test Execution
  result <- av_compare_scores(
    obj, 
    external_df = external_scores, 
    join_by = c("variant_id" = "var_id"), 
    external_score_col = "cadd_phred",
    method = "spearman"
  )
  
  # 4. Checks
  expect_type(result, "list")
  expect_true(all(c("joined_data", "stats", "correlation") %in% names(result)))
  expect_equal(nrow(result$joined_data), 4) # var5 (internal) and var6 (external) are dropped
  expect_equal(result$stats$estimate[["rho"]], 1)
})

test_that("av_plot_compare_scatter returns a ggplot object", {
  # Minimal setup for the plot test
  joined_df <- tibble::tibble(
    variant_id = c("v1", "v2"),
    score = c(1, 2),
    external_score = c(1, 2)
  )
  
  comp_result <- list(
    joined_data = joined_df,
    correlation = list(method = "Spearman", estimate = 1, p.value = 0.05),
    labels = list(x = "AlphaGenome", y = "External")
  )
  
  p <- av_plot_compare_scatter(comp_result)
  
  # Check if it is a ggplot object
  expect_s3_class(p, "ggplot")
})

test_that("av_plot_compare_venn returns a ggplot object", {
  # Setup: Mock output from av_analyze_overlap
  overlap_result <- list(
    cutoff_percentile = 0.05,
    counts = c(
      intersection = 10,
      alpha_only = 40,
      external_only = 40,
      neither = 100
    ),
    cutoffs = c(internal = 0.8, external = 15)
  )
  
  # Test Execution
  p <- av_plot_compare_venn(overlap_result)
  
  # Checks
  expect_s3_class(p, "ggplot")
  
  # Check if the labels were correctly transferred to the plot data
  # We look for numbers in the plot layer (text)
  # Note: ggplot builds are complex to test; we only check the class and absence of errors.
})