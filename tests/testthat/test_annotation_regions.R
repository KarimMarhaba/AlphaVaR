test_that("av_annotate_regions flags variants in custom regions", {
  # 1. Mock AlphaVarSet
  mock_data <- tibble::tibble(
    variant_id = c("v1_in_HAR", "v2_outside"),
    score = c(0.9, 0.1),
    modality = c("test", "test"),
    gene = c("A", "B"),
    coords = c("chr1:100-100", "chr1:500-500") # v1 at 100, v2 at 500
  )
  obj <- structure(list(data = mock_data, stats = list(), meta = list()), class = "AlphaVarSet")
  
  # 2. Mock Region Data (simulating a BED file import)
  # Region is chr1:50-150 (covers v1)
  # We use rtracklayer::import logic, but here we can pass a GRanges directly to test logic
  # or mock the file import if we want to be strict. 
  # For simplicity/speed in unit tests, we allow passing GRanges directly too.
  
  har_regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(50, 150))
  
  # 3. Run Function
  # We test the capability to accept a GRanges object to avoid file I/O in tests
  res <- av_annotate_regions(obj, regions = har_regions, label = "is_HAR")
  
  # 4. Assertions
  expect_true("is_HAR" %in% names(res$data))
  expect_true(res$data$is_HAR[1])  # v1 should be TRUE
  expect_false(res$data$is_HAR[2]) # v2 should be FALSE
})