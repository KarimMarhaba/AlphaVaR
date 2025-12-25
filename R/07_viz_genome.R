#' Parse Coordinate Strings
#'
#' Converts strings like "chr1:100-200" into a tibble with chrom, start, end.
#' @keywords internal
parse_coords <- function(coord_vec) {
  # Regex to capture: (chr...) : (start) - (end)
  # Handles cases with or without "chr" prefix
  pattern <- "^([\\w\\d]+):(\\d+)-(\\d+)$"
  
  parts <- stringr::str_match(coord_vec, pattern)
  
  if (any(is.na(parts[, 1]))) {
    warning("Some coordinates could not be parsed. Expected format 'chr:start-end'.")
  }
  
  dplyr::tibble(
    chrom = parts[, 2],
    start = as.numeric(parts[, 3]),
    end   = as.numeric(parts[, 4])
  )
}

#' Plot Genomic Locus (Manhattan-style)
#'
#' Visualizes variants along genomic coordinates for a specific gene or region.
#'
#' @param obj An \code{AlphaVarSet} object.
#' @param gene_name The gene to plot (filters variants for this gene).
#' @param flank_bp Number of base pairs to extend the plotting window (default 500).
#' @param color_by Column to map to color (default "modality").
#'
#' @return A ggplot object.
#' @importFrom dplyr filter bind_cols
#' @importFrom ggplot2 ggplot aes geom_segment geom_point facet_wrap labs theme_minimal theme element_blank
#' @importFrom stringr str_match
#' @export
av_plot_locus <- function(obj, gene_name, flank_bp = 0, color_by = "modality") {
  validate_AlphaVarSet(obj)
  dat <- obj$data
  
  # 1. Check for required column (mapped as 'coords' during import)
  if (!"coords" %in% names(dat)) {
    stop("No 'coords' column found in data. Please map 'coords' to your interval column in av_import_csv().")
  }
  
  # 2. Filter Data for Gene
  sub_dat <- dat %>% 
    dplyr::filter(gene == gene_name)
  
  if (nrow(sub_dat) == 0) stop(paste("No variants found for gene:", gene_name))
  
  # 3. Parse Coordinates
  parsed <- parse_coords(sub_dat$coords)
  
  # Combine with data
  plot_dat <- dplyr::bind_cols(sub_dat, parsed)
  
  # 4. Determine Window
  min_pos <- min(plot_dat$start) - flank_bp
  max_pos <- max(plot_dat$end) + flank_bp
  chrom   <- plot_dat$chrom[1]
  
# 5. Plot
  # Using .data[[var]] for tidy evaluation instead of deprecated aes_string
  ggplot2::ggplot(plot_dat, ggplot2::aes(x = start, xend = start, y = 0, yend = abs(score))) +
    ggplot2::geom_segment(ggplot2::aes(color = .data[[color_by]]), linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = abs(score), color = .data[[color_by]]), size = 2) +
    # Facet by modality to avoid overplotting
    ggplot2::facet_wrap(~ .data[[color_by]], ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = paste("Genomic Context:", gene_name),
      subtitle = paste0("Region: ", chrom, ":", min_pos, "-", max_pos),
      x = paste("Genomic Position (", chrom, ")"),
      y = "Absolute Score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank()
    )
}