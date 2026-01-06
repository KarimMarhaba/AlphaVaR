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
  
  # KORREKTUR: S3 Konstruktor statt R6 $new()
  # Wir erstellen das Objekt manuell, falls kein exportierter Konstruktor existiert,
  # oder nutzen die Konstruktor-Funktion. 
  # Da AlphaVarSet laut Fehlermeldung eine Funktion ist, rufen wir sie auf:
  obj <- AlphaVarSet(data = dummy_data) 
  
  # Falls AlphaVarSet() nicht exportiert ist oder anders heißt (z.B. new_AlphaVarSet),
  # können wir für den Test auch die Struktur manuell erzwingen (Fallback):
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
  expect_equal(nrow(result$joined_data), 4) # var5 (intern) und var6 (extern) fallen raus
  expect_equal(result$stats$estimate[["rho"]], 1)
})

test_that("av_plot_compare_scatter returns a ggplot object", {
  # Minimales Setup für den Plot-Test
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
  
  # Prüfen ob es ein ggplot Objekt ist
  expect_s3_class(p, "ggplot")
})