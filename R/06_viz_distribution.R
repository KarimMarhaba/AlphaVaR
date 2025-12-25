#' Plot Score Distributions (Ridge Plot)
#'
#' Creates a ridge plot (joyplot) of variant scores per modality.
#' Automatically sorts modalities by median shift if \code{av_calc_enrichment} was run.
#'
#' @param obj An \code{AlphaVarSet} object.
#' @param abs_score Logical. If TRUE (default), uses absolute scores.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes labs theme_minimal geom_vline scale_fill_viridis_c
#' @importFrom ggridges geom_density_ridges_gradient
#' @importFrom dplyr mutate arrange
#' @importFrom forcats fct_reorder
#' @export
av_plot_ridges <- function(obj, abs_score = TRUE) {
  validate_AlphaVarSet(obj)
  
  dat <- obj$data
  
  # Prepare Data
  plot_dat <- dat %>%
    dplyr::mutate(val = if(abs_score) abs(score) else score)
  
  # Sorting Logic: Use Wilcoxon results if available
  subtitle_text <- "Sorted alphabetically"
  
  if (!is.null(obj$stats$wilcoxon)) {
    # Extract order from stats
    order_levels <- obj$stats$wilcoxon %>%
      dplyr::arrange(median_diff) %>%
      dplyr::pull(group)
    
    plot_dat$modality <- factor(plot_dat$modality, levels = order_levels)
    subtitle_text <- "Sorted by median shift (Wilcoxon)"
  }
  
  # Plot
  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = val, y = modality, fill = ggplot2::after_stat(x))) +
    ggridges::geom_density_ridges_gradient(scale = 2.0, rel_min_height = 0.01) +
    ggplot2::scale_fill_viridis_c(name = "Score", option = "magma") +
    ggplot2::labs(
      title = "Score Distributions by Modality",
      subtitle = subtitle_text,
      x = if(abs_score) "Absolute Score" else "Raw Score",
      y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
    
  return(p)
}

#' Plot Ranked Boxplots
#'
#' Displays boxplots of scores, sorted by median value or Wilcoxon shift.
#'
#' @param obj An \code{AlphaVarSet} object.
#' @param abs_score Logical. Use absolute scores?
#' @return A ggplot object.
#' @export
av_plot_ranking <- function(obj, abs_score = TRUE) {
  validate_AlphaVarSet(obj)
  dat <- obj$data
  
  plot_dat <- dat %>%
    dplyr::mutate(val = if(abs_score) abs(score) else score)
  
  # Default sorting by median if no stats present
  if (!is.null(obj$stats$wilcoxon)) {
    order_levels <- obj$stats$wilcoxon %>% dplyr::arrange(median_diff) %>% dplyr::pull(group)
    plot_dat$modality <- factor(plot_dat$modality, levels = order_levels)
  } else {
    plot_dat$modality <- forcats::fct_reorder(plot_dat$modality, plot_dat$val, .fun = median)
  }
  
  global_med <- median(plot_dat$val, na.rm = TRUE)
  
  ggplot2::ggplot(plot_dat, ggplot2::aes(x = modality, y = val, fill = modality)) +
    ggplot2::geom_boxplot(outlier.size = 0.3, alpha = 0.8, size = 0.4) +
    ggplot2::geom_hline(yintercept = global_med, linetype = "dashed", color = "firebrick") +
    ggplot2::scale_fill_viridis_d(option = "cividis") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Effect Size Ranking",
      subtitle = "Red line = Global Median",
      x = NULL,
      y = if(abs_score) "Absolute Score" else "Raw Score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

#' Plot Cumulative Distribution (ECDF)
#'
#' @param obj An \code{AlphaVarSet} object.
#' @param abs_score Logical.
#' @return A ggplot object.
#' @export
av_plot_ecdf <- function(obj, abs_score = TRUE) {
  validate_AlphaVarSet(obj)
  dat <- obj$data
  
  plot_dat <- dat %>%
    dplyr::mutate(val = if(abs_score) abs(score) else score)
  
  ggplot2::ggplot(plot_dat, ggplot2::aes(x = val, color = modality)) +
    ggplot2::stat_ecdf(geom = "step", linewidth = 0.8, alpha = 0.8) +
    ggplot2::scale_color_viridis_d() +
    ggplot2::labs(
      title = "Cumulative Score Distribution",
      x = if(abs_score) "Absolute Score" else "Raw Score",
      y = "Cumulative Probability",
      color = "Assay"
    ) +
    ggplot2::theme_minimal()
}