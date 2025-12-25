test_that("Gene aggregation works correctly", {
  # Mock Data: Gene A has scores 0.1 and 0.9. Gene B has 0.2.
  dat <- dplyr::tibble(
    variant_id = 1:3,
    score = c(0.1, 0.9, 0.2),
    modality = "ATAC",
    gene = c("GeneA", "GeneA", "GeneB")
  )
  obj <- AlphaVarSet(dat)
  
  # Test MAX method
  res_max <- av_aggregate_genes(obj, method = "max")
  expect_equal(res_max$gene[1], "GeneA")
  expect_equal(res_max$gene_score[1], 0.9)
  
  # Test MEAN method
  res_mean <- av_aggregate_genes(obj, method = "mean")
  expect_equal(res_mean$gene[1], "GeneA")
  expect_equal(res_mean$gene_score[1], 0.5) # (0.1+0.9)/2
  
  # Test COUNT method (Threshold 0.5)
  res_count <- av_aggregate_genes(obj, method = "count_hits", threshold = 0.5)
  expect_equal(res_count$gene_score[1], 1) # Only 0.9 is a hit
})

test_that("Drill down returns correct variants", {
  dat <- dplyr::tibble(
    variant_id = c("v1", "v2"),
    score = c(0.5, 0.6),
    modality = "ATAC",
    gene = c("Target", "Other")
  )
  obj <- AlphaVarSet(dat)
  
  res <- av_get_variants_in_gene(obj, "Target")
  expect_equal(nrow(res), 1)
  expect_equal(res$variant_id, "v1")
  
  # Test non-existent gene warning
  expect_warning(av_get_variants_in_gene(obj, "MissingGene"), "No variants found")
})