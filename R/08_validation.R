# R/08_validation.R

#' Compare Internal Scores with External Data
#'
#' Joins the internal variant data with an external dataframe (e.g., CADD scores)
#' and calculates correlation statistics.
#'
#' @param obj An object of class \code{AlphaVarSet}.
#' @param external_df A data.frame/tibble containing external scores.
#' @param join_by A named vector indicating the join keys, e.g., \code{c("variant_id" = "ID_column_in_external")}.
#' @param external_score_col Character. The name of the column in \code{external_df} containing the numeric scores.
#' @param method Character. Correlation method: "pearson", "kendall", or "spearman" (default).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{joined_data}: The merged tibble (inner join).
#'   \item \code{stats}: The raw result of \code{cor.test}.
#'   \item \code{correlation}: A simplified list with method, rho/r, and p-value.
#'   \item \code{labels}: Metadata for plotting.
#' }
#' @export
#' @importFrom dplyr inner_join select rename
#' @importFrom stats cor.test
#' @importFrom rlang sym
av_compare_scores <- function(obj, external_df, join_by = "variant_id", external_score_col, method = "spearman") {
  
  # Input Validation
  stopifnot(inherits(obj, "AlphaVarSet"))
  stopifnot(is.data.frame(external_df))
  
  if (!external_score_col %in% colnames(external_df)) {
    stop(paste("Column", external_score_col, "not found in external_df."))
  }
  
  # Prepare Internal Data
  internal_data <- obj$data
  
  # Perform Inner Join
  # Note: This filters out variants not present in both sets
  joined_df <- dplyr::inner_join(internal_data, external_df, by = join_by)
  
  if (nrow(joined_df) < 3) {
    warning("Less than 3 overlapping variants found. Correlation might be meaningless.")
  }
  
  # Extract vectors for correlation
  # Internal score is always 'score' in AlphaVarSet
  vec_internal <- joined_df$score 
  vec_external <- joined_df[[external_score_col]]
  
  # Calculate Correlation
  cor_res <- stats::cor.test(vec_internal, vec_external, method = method)
  
  # Standardize Output for Plotting
  # We rename the external score column to 'external_score' in the output object 
  # to make plotting functions generic, but keep the original joined data clean.
  plot_data <- joined_df
  plot_data$external_score <- vec_external
  
  result <- list(
    joined_data = plot_data,
    stats = cor_res,
    correlation = list(
      method = method,
      estimate = cor_res$estimate[[1]],
      p.value = cor_res$p.value
    ),
    labels = list(
      x = "AlphaGenome Score",
      y = external_score_col
    )
  )
  
  return(result)
}

#' Plot Score Comparison (Scatterplot)
#'
#' Visualizes the correlation between AlphaGenome scores and an external metric.
#'
#' @param comparison_result The list output from \code{av_compare_scores}.
#' @param alpha Numeric. Transparency of points (default 0.6).
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_minimal annotate
av_plot_compare_scatter <- function(comparison_result, alpha = 0.6) {
  
  data <- comparison_result$joined_data
  cor_stats <- comparison_result$correlation
  lbls <- comparison_result$labels
  
  # Create Label Text
  cor_text <- paste0(
    tools::toTitleCase(cor_stats$method), ": ", 
    round(cor_stats$estimate, 3), 
    "\n(p < ", format.pval(cor_stats$p.value, digits = 3), ")"
  )
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$score, y = .data$external_score)) +
    ggplot2::geom_point(alpha = alpha, color = "#2c3e50") +
    ggplot2::geom_smooth(method = "lm", color = "#e74c3c", se = TRUE, linetype = "dashed") +
    ggplot2::labs(
      title = paste("Validation vs.", lbls$y),
      subtitle = paste("N =", nrow(data), "overlapping variants"),
      x = lbls$x,
      y = lbls$y
    ) +
    ggplot2::theme_minimal() +
    # Add correlation text inside the plot
    ggplot2::annotate(
      "text", 
      x = -Inf, y = Inf, 
      label = cor_text, 
      hjust = -0.1, vjust = 1.5, 
      size = 5, fontface = "italic"
    )
    
  return(p)
}

#' Analyze Top-Hit Overlap
#'
#' Calculates the overlap of the top X% variants between AlphaGenome and external scores.
#' Useful for checking if the model prioritizes similar variants ("Venn Diagram logic").
#'
#' @param comparison_result The list output from \code{av_compare_scores}.
#' @param top_percentile Numeric. Cutoff for top hits (e.g., 0.01 for top 1%).
#' @param high_is_good Logical. If TRUE, high scores are top hits. If FALSE (e.g. p-values), low scores are top.
#'
#' @return A list with the confusion matrix and overlap counts.
#' @export
#' @importFrom stats quantile
av_analyze_overlap <- function(comparison_result, top_percentile = 0.05, high_is_good = TRUE) {
  
  df <- comparison_result$joined_data
  
  # Determine cutoffs
  probs <- if(high_is_good) (1 - top_percentile) else top_percentile
  
  cut_internal <- stats::quantile(df$score, probs = probs, names = FALSE)
  cut_external <- stats::quantile(df$external_score, probs = probs, names = FALSE)
  
  # Label Logic
  if (high_is_good) {
    df$top_internal <- df$score >= cut_internal
    df$top_external <- df$external_score >= cut_external
  } else {
    df$top_internal <- df$score <= cut_internal
    df$top_external <- df$external_score <= cut_external
  }
  
  # Venn / Confusion logic
  # Both Top
  n_intersect <- sum(df$top_internal & df$top_external)
  # Only Alpha
  n_alpha_only <- sum(df$top_internal & !df$top_external)
  # Only External
  n_ext_only <- sum(!df$top_internal & df$top_external)
  # Neither
  n_neither <- sum(!df$top_internal & !df$top_external)
  
  return(list(
    cutoff_percentile = top_percentile,
    cutoffs = c(internal = cut_internal, external = cut_external),
    counts = c(
      intersection = n_intersect,
      alpha_only = n_alpha_only,
      external_only = n_ext_only,
      neither = n_neither
    ),
    overlap_fraction = n_intersect / (n_intersect + n_alpha_only + n_ext_only)
  ))
}

#' Plot Venn Diagram of Top Hits
#'
#' Visualizes the overlap of top-ranked variants between AlphaGenome and external scores
#' using a 2-circle Venn diagram.
#'
#' @param overlap_result The list output from \code{av_analyze_overlap}.
#' @param fill_colors Character vector of length 2. Colors for the circles (default: AlphaVaR blue and External grey/red).
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot geom_polygon annotate theme_void labs coord_fixed scale_fill_identity aes
#' @importFrom tibble tibble
av_plot_compare_venn <- function(overlap_result, fill_colors = c("#3498db", "#e74c3c")) {
  
  counts <- overlap_result$counts
  
  # Helper to generate circle coordinates
  # We use simple trigonometry to draw circles without needing 'ggforce'
  get_circle <- function(center_x, center_y, radius = 1.5) {
    theta <- seq(0, 2 * pi, length.out = 200)
    tibble::tibble(
      x = center_x + radius * cos(theta),
      y = center_y + radius * sin(theta)
    )
  }
  
  # Create Data for Plotting
  # Left Circle: AlphaVaR (Center -1, 0)
  c1 <- get_circle(-0.8, 0)
  c1$group <- "AlphaVaR"
  c1$fill <- fill_colors[1]
  
  # Right Circle: External (Center 1, 0)
  c2 <- get_circle(0.8, 0)
  c2$group <- "External"
  c2$fill <- fill_colors[2]
  
  plot_df <- rbind(c1, c2)
  
  # Labels positions
  # Left (Alpha Only), Center (Intersection), Right (External Only)
  label_df <- tibble::tibble(
    x = c(-1.5, 0, 1.5),
    y = c(0, 0, 0),
    label = c(
      paste0("AlphaVaR\nOnly\n", counts["alpha_only"]),
      paste0("Overlap\n", counts["intersection"]),
      paste0("External\nOnly\n", counts["external_only"])
    ),
    color = c("white", "white", "white") # Text color
  )
  
  # Title info
  pct <- overlap_result$cutoff_percentile * 100
  
  # Plotting
  p <- ggplot2::ggplot() +
    # Draw Circles
    ggplot2::geom_polygon(
      data = plot_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$fill),
      alpha = 0.6,
      color = "white"
    ) +
    # Add Text Labels
    ggplot2::annotate(
      "text", x = label_df$x, y = label_df$y, label = label_df$label,
      size = 5, fontface = "bold", color = "white" # Contrast text
    ) +
    # Settings
    ggplot2::scale_fill_identity() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() + # Remove axes/grid
    ggplot2::labs(
      title = paste0("Top ", pct, "% Hits Overlap"),
      subtitle = paste0("Intersection Count: ", counts["intersection"])
    )
  
  return(p)
}