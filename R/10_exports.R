#' Export Variant Scores as Browser Track
#'
#' Exports the variant scores to a BedGraph format file, suitable for viewing
#' in the UCSC Genome Browser or IGV.
#'
#' @param obj An object of class \code{AlphaVarSet}.
#' @param output_file Character string. The path to the output file (e.g., "scores.bedgraph").
#' @param track_name Character string. The name of the track displayed in the browser.
#' @param description Character string. A description for the browser track line.
#' @param color Character string. RGB color for the track (default: "0,0,255" for blue).
#'
#' @return Invisible. The name of the output file.
#' @export
#'
#' @examples
#' \dontrun{
#' av_export_browser_track(my_obj, "output/variant_scores.bedgraph")
#' }
#' @importFrom readr write_lines
#' @importFrom stringr str_match
#' @importFrom dplyr mutate select arrange
#' @importFrom stats na.omit
av_export_browser_track <- function(obj,
                                    output_file,
                                    track_name = "AlphaVaR_Scores",
                                    description = "Variant Effect Scores",
                                    color = "0,0,255") {

  # 1. Validate Input
  if (!inherits(obj, "AlphaVarSet")) {
    stop("Error: Input 'obj' must be of class 'AlphaVarSet'.")
  }

  # 2. Parse Coordinates (Assumes "chr:start-end" format in 'coords' column)
  # We extract chr, start, end for BedGraph format
  # BedGraph format: chrom chromStart chromEnd dataValue
  df <- obj$data

  # Check if coords exist
  if (!"coords" %in% names(df)) {
    stop("Error: 'coords' column missing in data.")
  }

  # Regex to parse "chr1:100-200"
  # This adds temporary columns for export
  parsed_df <- df %>%
    dplyr::mutate(
      chrom = stringr::str_match(coords, "^([^:]+)")[, 2],
      range_part = stringr::str_match(coords, ":([0-9]+-[0-9]+)")[, 2]
    ) %>%
    dplyr::mutate(
      start = as.numeric(stringr::str_split_fixed(range_part, "-", 2)[, 1]),
      end   = as.numeric(stringr::str_split_fixed(range_part, "-", 2)[, 2])
    ) %>%
    dplyr::select(chrom, start, end, score) %>%
    dplyr::arrange(chrom, start)

  # Check for parsing errors
  if (any(is.na(parsed_df$chrom)) || any(is.na(parsed_df$start))) {
    warning("Some coordinates could not be parsed. Ensure 'coords' format is 'chr:start-end'.")
    parsed_df <- na.omit(parsed_df)
  }

  # 3. Create Header (UCSC Format)
  header_line <- sprintf(
    'track type=bedGraph name="%s" description="%s" visibility=full color=%s autoScale=on alwaysZero=on',
    track_name, description, color
  )

  # 4. Write to file
  # We write the header first, then append the data
  readr::write_lines(header_line, output_file)

  # Format data as tab-separated
  readr::write_delim(
    parsed_df,
    output_file,
    delim = "\t",
    col_names = FALSE,
    append = TRUE
  )

  message(sprintf("Success: Browser track exported to %s (%d variants)", output_file, nrow(parsed_df)))
  invisible(output_file)
}