#' Plot Enrichment Volcano
#'
#' Visualizes Fisher's Exact Test results (Log2 Odds Ratio vs -Log10 FDR).
#' Requires \code{av_calc_enrichment()} to be run first.
#'
#' @param obj An \code{AlphaVarSet} object with computed stats.
#' @param label_top_n Number of top hits to label (default 10).
#' @return A ggplot object.
#' @importFrom ggrepel geom_text_repel
#' @export
av_plot_volcano <- function(obj, label_top_n = 10) {
  validate_AlphaVarSet(obj)
  
  # 1. Check if stats exist
  if (is.null(obj$stats$enrichment)) {
    stop("No enrichment statistics found. Please run av_calc_enrichment() first.")
  }
  
  res <- obj$stats$enrichment
  params <- obj$stats$params
  
  # 2. Classify for coloring
  plot_dat <- res %>%
    dplyr::mutate(
      class = dplyr::case_when(
        p_adj < 0.05 & log2OR > 1 ~ "Enriched",
        p_adj < 0.05 & log2OR < -1 ~ "Depleted",
        TRUE ~ "NS"
      )
    )
  
  # 3. Plot
  ggplot2::ggplot(plot_dat, ggplot2::aes(x = log2OR, y = neg_log10_fdr)) +
    ggplot2::geom_point(ggplot2::aes(color = class, size = n_obs), alpha = 0.7) +
    ggrepel::geom_text_repel(
      data = head(plot_dat, label_top_n), 
      ggplot2::aes(label = group), 
      size = 3,
      max.overlaps = 20
    ) +
    ggplot2::scale_color_manual(values = c("Enriched" = "#B22222", "Depleted" = "#4682B4", "NS" = "grey80")) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.4) +
    ggplot2::labs(
      title = paste0("Enrichment Volcano (Threshold: ", params$threshold, ")"),
      x = "Log2 Odds Ratio",
      y = "-Log10 FDR",
      color = "Status",
      size = "N Variants"
    ) +
    ggplot2::theme_minimal()
}