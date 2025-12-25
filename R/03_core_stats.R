#' Calculate Enrichment Statistics
#'
#' Performs Fisher's Exact Test and Wilcoxon Rank Sum Test on the variants.
#' Results are stored in the \code{stats} slot of the object.
#'
#' @param obj An object of class \code{AlphaVarSet}.
#' @param group_by The column to group by (default "modality").
#' @param threshold The score threshold for Fisher's test (default 0.5).
#' @param run_wilcoxon Logical, whether to run the Wilcoxon test (time consuming).
#'
#' @return The updated \code{AlphaVarSet} object with results in \code{obj$stats}.
#' @importFrom dplyr select filter mutate group_by summarise n rowwise ungroup arrange desc bind_rows
#' @importFrom stats fisher.test wilcox.test p.adjust median
#' @importFrom tibble tibble
#' @importFrom rlang sym
#' @export
av_calc_enrichment <- function(obj, group_by = "modality", threshold = 0.5, run_wilcoxon = TRUE) {
  validate_AlphaVarSet(obj)
  
  # Check if grouping column exists
  if (!group_by %in% names(obj$data)) {
    stop(paste("Grouping column", group_by, "not found in data."))
  }
  
  # Prepare data (using rlang::sym for dynamic grouping)
  dat <- obj$data
  grp_sym <- rlang::sym(group_by)
  
  # Use absolute scores for stats as defined in logic
  df_proc <- dat %>%
    dplyr::select(group = !!grp_sym, score) %>%
    dplyr::filter(!is.na(group), !is.na(score)) %>%
    dplyr::mutate(
      val = abs(score),
      is_hit = val >= threshold
    )
  
  # --- Fisher Test ---
  n_total <- nrow(df_proc)
  n_hits_global <- sum(df_proc$is_hit)
  
  enrichment <- df_proc %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(n_obs = dplyr::n(), n_hits = sum(is_hit), .groups = "drop") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      a = n_hits, b = n_obs - n_hits,
      c = n_hits_global - n_hits, d = (n_total - n_hits_global) - b,
      # Fisher Test
      ft = list(tryCatch(
        stats::fisher.test(matrix(c(a, c, b, d), nrow = 2)),
        error = function(e) list(estimate = NA, p.value = NA)
      )),
      or = ft$estimate,
      p_raw = ft$p.value,
      # Log2 Odds Ratio with Haldane correction
      log2OR = log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_raw, method = "BH"),
      neg_log10_fdr = -log10(p_adj)
    ) %>%
    dplyr::select(group, n_obs, n_hits, or, log2OR, p_adj, neg_log10_fdr) %>%
    dplyr::arrange(p_adj, dplyr::desc(abs(log2OR)))
  
  # --- Wilcoxon Test ---
  shift <- NULL
  if (run_wilcoxon) {
    # Extract vectors for speed
    vals <- df_proc$val
    grps <- df_proc$group
    u_grps <- unique(grps)
    
    wilcox_list <- lapply(u_grps, function(g) {
      is_g <- grps == g
      x <- vals[is_g]
      y <- vals[!is_g]
      # Robustness check: need at least 5 obs to run test
      if (length(x) < 5 || length(y) < 5) return(NULL)
      
      wt <- stats::wilcox.test(x, y, alternative = "two.sided", exact = FALSE)
      tibble::tibble(
        group = g,
        median_diff = stats::median(x) - stats::median(y),
        p_val = wt$p.value
      )
    })
    
    # --- BUG FIX START ---
    # Bind rows first
    shift_res <- dplyr::bind_rows(wilcox_list)
    
    # Only calculate p_adj if we actually have results (rows > 0)
    if (nrow(shift_res) > 0) {
      shift <- shift_res %>%
        dplyr::mutate(p_adj = stats::p.adjust(p_val, method = "BH")) %>%
        dplyr::arrange(dplyr::desc(median_diff))
    } else {
      # Return empty structure to avoid downstream errors
      shift <- tibble::tibble(
        group = character(), 
        median_diff = numeric(), 
        p_val = numeric(), 
        p_adj = numeric()
      )
    }
    # --- BUG FIX END ---
  }
  
  # Store results in object
  obj$stats$enrichment <- enrichment
  obj$stats$wilcoxon <- shift
  obj$stats$params <- list(threshold = threshold, group_by = group_by)
  
  return(obj)
}