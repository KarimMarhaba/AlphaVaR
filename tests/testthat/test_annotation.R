test_that("av_annotate_genes adds correct columns and identifies regions", {
  # 1. Mock AlphaVarSet
  # Create a minimal set with 2 variants:
  # Var1: Should fall into a "Promoter" region (simulating a hit near a TSS)
  # Var2: Should be "Intergenic" (far away from any defined gene)
  
  mock_data <- tibble::tibble(
    variant_id = c("v1", "v2"),
    score = c(0.9, 0.1),
    modality = c("H3K27ac", "H3K27ac"),
    gene = c(NA, NA),
    # Simulated coordinates.
    # We will use TP53 (chr17) as a known target if the hg38 package is loaded.
    coords = c("chr17:7661779-7661779", "chr1:1-1") 
  )
  
  obj <- structure(list(data = mock_data, stats = list(), meta = list()), class = "AlphaVarSet")
  
  # 2. Dependency Check
  # We skip this test on systems lacking the genomic infrastructure to avoid build failures.
  skip_if_not_installed("GenomicFeatures")
  skip_if_not_installed("GenomicRanges")
  
  # 3. Conditional Execution
  # We only run the logic if the specific HG38 annotation package is available.
  if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    
    # Run the function
    res <- av_annotate_genes(obj, txdb_package = "TxDb.Hsapiens.UCSC.hg38.knownGene")
    
    # 4. Assertions
    expect_s3_class(res, "AlphaVarSet")
    
    # Check for structural integrity (new columns added)
    expected_cols <- c("nearest_gene", "dist_to_tss", "region_type")
    expect_true(all(expected_cols %in% names(res$data)))
    
    # Check logic validation for TP53 hit (Var1)
    # Note: Exact classification depends on the specific build of TxDb, 
    # but it should definitely not be intergenic if it sits on the TSS.
    tp53_hit <- res$data[1, ]
    expect_true(tp53_hit$region_type %in% c("promoter", "gene_body"))
    
    # Check logic for the random coordinate (Var2)
    random_hit <- res$data[2, ]
    expect_equal(random_hit$region_type, "intergenic")
    
  } else {
    # If the package is missing, ensure the function handles the error gracefully
    expect_error(av_annotate_genes(obj, "NonExistentPackage"), "is required but not installed")
  }
})