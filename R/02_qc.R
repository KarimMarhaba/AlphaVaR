#' Generate QC Summary
#'
#' Provides a quick overview of the dataset dimensions and missing values.
#'
#' @param obj An object of class \code{AlphaVarSet}.
#' @return A tibble with summary statistics.
#' @importFrom dplyr summarise n_distinct
#' @importFrom tibble tibble
#' @export
av_qc_summary <- function(obj) {
  validate_AlphaVarSet(obj)
  
  dat <- obj$data
  
  res <- tibble::tibble(
    n_variants = nrow(dat),
    n_modalities = dplyr::n_distinct(dat$modality),
    n_genes = if ("gene" %in% names(dat)) dplyr::n_distinct(dat$gene, na.rm = TRUE) else 0,
    mean_score = mean(dat$score, na.rm = TRUE),
    missing_scores = sum(is.na(dat$score))
  )
  
  return(res)
}

#' Plot Score Distribution (QC)
#'
#' Visualizes the distribution of raw scores to check for artifacts or skew.
#'
#' @param obj An object of class \code{AlphaVarSet}.
#' @param type Plot type: "histogram" or "density".
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density labs theme_minimal
#' @export
av_plot_qc_dist <- function(obj, type = "histogram") {
  validate_AlphaVarSet(obj)
  
  p <- ggplot2::ggplot(obj$data, ggplot2::aes(x = score)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "QC: Score Distribution", x = "Score", y = "Count")
  
  if (type == "histogram") {
    p <- p + ggplot2::geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7)
  } else {
    p <- p + ggplot2::geom_density(fill = "steelblue", alpha = 0.4)
  }
  
  return(p)
}