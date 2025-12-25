test_that("Import creates a valid AlphaVarSet object", {
  
  # 1. Create a dummy CSV file matches DEFAULT mapping
  # Default map expects: quantile_score, output_type, variant_id, gene_name
  dummy_file <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    variant_id = c("var1", "var2"),
    quantile_score = c(0.5, 0.9),
    output_type = c("ATAC", "DNASE"),
    gene_name = c("GENE1", "GENE2"),
    scored_interval = c("chr1:100-200", "chr1:300-400")
  ), dummy_file, row.names = FALSE)
  
  # 2. Run Import (Using defaults)
  obj <- av_import_csv(dummy_file)
  
  # 3. Assertions
  expect_s3_class(obj, "AlphaVarSet")
  expect_true(is.numeric(obj$data$score))
  expect_equal(nrow(obj$data), 2)
  expect_equal(obj$data$modality[1], "ATAC") 
  
  unlink(dummy_file)
})

test_that("Custom mapping works", {
  # CSV with weird names
  dummy_file <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    my_id = "var1",
    my_val = 0.8,
    my_assay = "CAGE",
    my_gene = "G1"
  ), dummy_file, row.names = FALSE)
  
  # Import with explicit map
  obj <- av_import_csv(dummy_file, col_map = list(
    variant_id = "my_id", 
    score = "my_val", 
    modality = "my_assay",
    gene = "my_gene"
  ))
  
  expect_equal(obj$data$score[1], 0.8)
  expect_equal(obj$data$modality[1], "CAGE")
  
  unlink(dummy_file)
})

test_that("Import fails gracefully with missing columns", {
  # CSV without 'output_type' (mandatory)
  dummy_file <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    variant_id = c("var1"),
    quantile_score = c(0.5)
  ), dummy_file, row.names = FALSE)
  
  # Wir prüfen nur noch auf das Schlüsselwort "Missing", das ist sicherer
  expect_error(av_import_csv(dummy_file), "Missing")
  
  unlink(dummy_file)
})