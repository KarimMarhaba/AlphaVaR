#' Aggregate Scores by Gene
#'
#' Summarizes variant scores at the gene level to identify the most affected genes.
#'
#' @param obj An \code{AlphaVarSet} object.
#' @param method Aggregation method: "max" (highest single score), "mean" (average burden), 
#'        or "count_hits" (number of variants above threshold).
#' @param threshold Used only if method is "count_hits".
#'
#' @return A tibble with one row per gene, sorted by the aggregated score.
#' @importFrom dplyr group_by summarise arrange desc filter
#' @importFrom rlang arg_match
#' @export
av_aggregate_genes <- function(obj, method = c("max", "mean", "sum", "count_hits"), threshold = 0.5) {
  validate_AlphaVarSet(obj)
  dat <- obj$data
  
  # Check if gene info exists
  if (!"gene" %in% names(dat)) {
    stop("No 'gene' column found in data. Did you map it correctly during import?")
  }
  
  method <- rlang::arg_match(method)
  
  # Filter out NAs in gene names
  clean_dat <- dat %>%
    dplyr::filter(!is.na(gene))
  
  # Aggregation Logic
  res <- clean_dat %>%
    dplyr::group_by(gene)
    
  if (method == "max") {
    res <- res %>% 
      dplyr::summarise(
        gene_score = max(abs(score), na.rm = TRUE),
        n_variants = dplyr::n()
      )
  } else if (method == "mean") {
    res <- res %>% 
      dplyr::summarise(
        gene_score = mean(abs(score), na.rm = TRUE),
        n_variants = dplyr::n()
      )
  } else if (method == "sum") {
    res <- res %>% 
      dplyr::summarise(
        gene_score = sum(abs(score), na.rm = TRUE),
        n_variants = dplyr::n()
      )
  } else if (method == "count_hits") {
    res <- res %>% 
      dplyr::summarise(
        gene_score = sum(abs(score) >= threshold, na.rm = TRUE),
        n_variants = dplyr::n()
      )
  }
  
  # Final Formatting
  res <- res %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(gene_score))
    
  return(res)
}

#' Get Variants for a Specific Gene
#'
#' Extracts all variants associated with a specific gene name.
#'
#' @param obj An \code{AlphaVarSet} object.
#' @param gene_name Name of the gene to inspect (case sensitive).
#'
#' @return A tibble containing only variants for that gene.
#' @importFrom dplyr filter arrange desc
#' @export
av_get_variants_in_gene <- function(obj, gene_name) {
  validate_AlphaVarSet(obj)
  dat <- obj$data
  
  if (!"gene" %in% names(dat)) stop("No 'gene' column found.")
  
  res <- dat %>%
    dplyr::filter(gene == gene_name) %>%
    dplyr::arrange(dplyr::desc(abs(score)))
    
  if (nrow(res) == 0) {
    warning(paste("No variants found for gene:", gene_name))
  }
  
  return(res)
}