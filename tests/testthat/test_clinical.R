test_that("av_annotate_clinvar works with mock VCF", {
  # 1. Setup Mock Data
  # Create a minimal AlphaVarSet
  mock_data <- tibble::tibble(
    variant_id = c("v1", "v2", "v3"),
    score = c(0.9, 0.1, 0.5),
    modality = "ATAC",
    coords = c("chr1:100-100", "chr1:200-200", "chr1:500-500") # v1 hits, v2 hits benign, v3 misses
  )
  obj <- structure(list(data = mock_data, stats = list(), meta = list()), class = "AlphaVarSet")

  # 2. Create a temporary Mock VCF file
  # ClinVar format lines
  vcf_lines <- c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance\">",
    "##INFO=<ID=CLNDN,Number=.,Type=String,Description=\"Disease name\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\trsv1\tA\tT\t.\t.\tCLNSIG=Pathogenic;CLNDN=Severe_Disease",
    "chr1\t200\trsv2\tC\tG\t.\t.\tCLNSIG=Benign;CLNDN=Trait_A"
  )

  tmp_vcf <- tempfile(fileext = ".vcf")
  writeLines(vcf_lines, tmp_vcf)

  # 3. Run Function
  # Skip if VariantAnnotation not installed (standard CRAN check politeness)
  skip_if_not_installed("VariantAnnotation")

  res <- av_annotate_clinvar(obj, vcf_file = tmp_vcf)

  # 4. Assertions
  expect_true("is_pathogenic" %in% colnames(res$data))
  expect_true("clinvar_trait" %in% colnames(res$data))

  # Check v1 (Pathogenic match)
  v1 <- dplyr::filter(res$data, variant_id == "v1")
  expect_true(v1$is_pathogenic)
  expect_equal(v1$clinvar_id, "rsv1")
  expect_match(v1$clinvar_trait, "Severe_Disease")

  # Check v2 (Benign match)
  v2 <- dplyr::filter(res$data, variant_id == "v2")
  expect_false(v2$is_pathogenic) # Should be false as it is Benign
  expect_match(v2$clinical_significance, "Benign")

  # Check v3 (No match)
  v3 <- dplyr::filter(res$data, variant_id == "v3")
  expect_true(is.na(v3$clinvar_id))
  expect_false(v3$is_pathogenic)

  # Clean up
  unlink(tmp_vcf)
})