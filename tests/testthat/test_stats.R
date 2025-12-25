test_that("QC summary works", {
  # Mock Data
  dat <- dplyr::tibble(
    variant_id = 1:10,
    score = c(rep(0.1, 5), rep(0.9, 5)), # 5 Low, 5 High
    modality = c(rep("A", 5), rep("B", 5)),
    gene = "G1"
  )
  obj <- AlphaVarSet(dat)
  
  qc <- av_qc_summary(obj)
  expect_equal(qc$n_variants, 10)
  expect_equal(qc$mean_score, 0.5)
})

test_that("Enrichment calculation works", {
  # Create a scenario where Modality B is clearly enriched (all 0.9)
  # Modality A is depleted (all 0.1)
  dat <- dplyr::tibble(
    variant_id = 1:20,
    score = c(rep(0.1, 10), rep(0.9, 10)),
    modality = c(rep("ModA", 10), rep("ModB", 10))
  )
  obj <- AlphaVarSet(dat)
  
  # Run Stats
  obj <- av_calc_enrichment(obj, threshold = 0.5)
  
  # Check structure
  expect_true("enrichment" %in% names(obj$stats))
  expect_true("wilcoxon" %in% names(obj$stats))
  
  # Check Logic: ModB should have high Log2OR (positive)
  res <- obj$stats$enrichment
  mod_b <- res[res$group == "ModB", ]
  mod_a <- res[res$group == "ModA", ]
  
  expect_gt(mod_b$log2OR, 0) # Should be positive
  expect_lt(mod_a$log2OR, 0) # Should be negative
  
  # Wilcoxon check
  w_res <- obj$stats$wilcoxon
  diff_b <- w_res[w_res$group == "ModB", ]$median_diff
  expect_gt(diff_b, 0)
})