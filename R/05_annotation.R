#' Annotate Variants with Nearest Gene and Genomic Context
#'
#' Assigns biological context to variants by identifying the nearest gene (TSS),
#' calculating the distance to the Transcription Start Site (TSS), and classifying
#' the genomic region.
#'
#' @description
#' This function maps genomic coordinates from the `AlphaVarSet` to a reference
#' transcript database (TxDb). It appends three columns to the `$data` slot:
#' \itemize{
#'   \item \code{nearest_gene}: The Entrez ID (or identifier) of the closest gene.
#'   \item \code{dist_to_tss}: The distance in base pairs to the nearest TSS.
#'   \item \code{region_type}: Classification into "promoter", "gene_body", or "intergenic".
#' }
#'
#' @param obj An object of class \code{AlphaVarSet}.
#' @param txdb_package Character string. Name of the TxDb package to use 
#' (default: "TxDb.Hsapiens.UCSC.hg38.knownGene"). 
#' Can also be a \code{TxDb} object directly.
#' @param promoter_upstream Numeric. Distance upstream of the TSS to define the promoter region (default: 2000).
#' @param promoter_downstream Numeric. Distance downstream of the TSS to define the promoter region (default: 200).
#'
#' @return An \code{AlphaVarSet} object with updated metadata columns.
#' @export
#' @importFrom GenomicRanges GRanges mcols distanceToNearest promoters findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom GenomicFeatures genes
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom methods is
#' @importFrom rlang .data
av_annotate_genes <- function(obj, 
                              txdb_package = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                              promoter_upstream = 2000,
                              promoter_downstream = 200) {
  
  # 1. Validation
  if (!inherits(obj, "AlphaVarSet")) stop("Input must be an AlphaVarSet object.")
  
  # Load TxDb database
  txdb <- .load_txdb(txdb_package)
  
  # 2. Parse User Coordinates to GRanges
  # Converts the string representation "chr:start-end" into a GRanges object
  variants_gr <- .parse_coords_to_gr(obj$data$coords)
  
  # 3. Extract Genomic Features
  # We extract gene definitions to calculate proximity and location
  message("Extracting gene information from annotation database...")
  genes_gr <- GenomicFeatures::genes(txdb)
  
  # Define promoters based on the TSS of the extracted genes
  promoters_gr <- GenomicRanges::promoters(genes_gr, 
                                           upstream = promoter_upstream, 
                                           downstream = promoter_downstream)
  
  # 4. Calculate Nearest Gene (Distance to TSS)
  message("Calculating distances to nearest TSS...")
  
  # Calculate distance to the nearest promoter/TSS (ignoring strand orientation for simple distance)
  hits <- GenomicRanges::distanceToNearest(variants_gr, promoters_gr, ignore.strand = TRUE)
  
  # Initialize result vectors
  n_vars <- nrow(obj$data)
  res_gene_id <- rep(NA_character_, n_vars)
  res_dist <- rep(NA_real_, n_vars)
  res_type <- rep("intergenic", n_vars) # Default state
  
  # Map results back to the original data structure
  query_idx <- S4Vectors::queryHits(hits)
  subject_idx <- S4Vectors::subjectHits(hits)
  
  # Retrieve Gene IDs (typically Entrez IDs in TxDb objects)
  gene_ids <- names(genes_gr) 
  
  res_gene_id[query_idx] <- gene_ids[subject_idx]
  res_dist[query_idx] <- S4Vectors::mcols(hits)$distance
  
  # 5. Determine Region Type 
  # Hierarchy: Promoter > Gene Body > Intergenic
  
  # A. Check Promoter Overlaps
  promoter_ov <- GenomicRanges::findOverlaps(variants_gr, promoters_gr)
  res_type[S4Vectors::queryHits(promoter_ov)] <- "promoter"
  
  # B. Check Gene Body Overlaps 
  # Only check variants that are NOT already classified as promoters
  gene_ov <- GenomicRanges::findOverlaps(variants_gr, genes_gr)
  gene_hits_idx <- S4Vectors::queryHits(gene_ov)
  
  is_promoter <- res_type == "promoter"
  is_in_gene <- seq_len(n_vars) %in% gene_hits_idx
  
  # Assign "gene_body" only if it overlaps a gene but is not in the promoter region
  res_type[!is_promoter & is_in_gene] <- "gene_body"
  
  # 6. Update and Return Object
  obj$data$nearest_gene <- res_gene_id
  obj$data$dist_to_tss <- res_dist
  obj$data$region_type <- res_type
  
  message("Annotation complete.")
  return(obj)
}

#' Annotate Variants with Custom Regions (e.g., HARs)
#'
#' Checks for overlap between variants and a set of genomic regions provided
#' either as a BED file or a GRanges object. Adds a boolean column to the data.
#'
#' @param obj An object of class `AlphaVarSet`.
#' @param regions Path to a BED file (character) or a `GRanges` object.
#' @param label Character string. Name of the new column (default: "in_region").
#'
#' @return An `AlphaVarSet` object with a new boolean column indicating overlap.
#' @export
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom S4Vectors queryHits
#' @importFrom methods is
av_annotate_regions <- function(obj, regions, label = "in_region") {
  
  # 1. Validation
  if (!inherits(obj, "AlphaVarSet")) stop("Input must be an AlphaVarSet object.")
  if (label %in% names(obj$data)) warning(sprintf("Column '%s' already exists and will be overwritten.", label))
  
  # 2. Parse Regions
  region_gr <- NULL
  
  if (methods::is(regions, "GRanges")) {
    region_gr <- regions
  } else if (is.character(regions) && file.exists(regions)) {
    # We suggest rtracklayer for file import
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
      stop("Package 'rtracklayer' is required to read BED files. Please install it.")
    }
    region_gr <- rtracklayer::import(regions)
  } else {
    stop("Argument 'regions' must be a GRanges object or a valid file path.")
  }
  
  # 3. Parse Variant Coords
  # Reuse internal helper from previous task
  variants_gr <- .parse_coords_to_gr(obj$data$coords)
  
  # 4. Find Overlaps
  # minoverlap=1 is standard. strict overlap.
  ov <- GenomicRanges::findOverlaps(variants_gr, region_gr)
  
  # 5. Create Result Column
  # Initialize FALSE
  res_vec <- rep(FALSE, nrow(obj$data))
  
  # Set TRUE for hits
  # queryHits returns indices of variants that overlap
  res_vec[unique(S4Vectors::queryHits(ov))] <- TRUE
  
  # 6. Update Object
  obj$data[[label]] <- res_vec
  
  return(obj)
}

# --- Internal Helper Functions ---

#' Load TxDb object safely
#' @noRd
.load_txdb <- function(pkg_or_obj) {
  if (methods::is(pkg_or_obj, "TxDb")) {
    return(pkg_or_obj)
  } else if (is.character(pkg_or_obj)) {
    if (!requireNamespace(pkg_or_obj, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", pkg_or_obj))
    }
    # Dynamically load the package namespace
    return(eval(parse(text = paste0(pkg_or_obj, "::", pkg_or_obj))))
  } else {
    stop("Argument 'txdb_package' must be a character string or a TxDb object.")
  }
}

#' Parse "chr:start-end" strings into GRanges
#' @noRd
.parse_coords_to_gr <- function(coord_strings) {
  # Expected format: "chr1:100-200" or "chr1:100"
  
  # Split string by delimiters ':' and '-'
  parts <- stringr::str_split_fixed(coord_strings, "[:-]", 3)
  
  seqnames <- parts[, 1]
  start_pos <- as.numeric(parts[, 2])
  
  # Handle potential SNVs where end coordinate might be omitted implies end = start
  end_pos <- as.numeric(parts[, 3])
  end_pos[parts[, 3] == ""] <- start_pos[parts[, 3] == ""]
  
  if (any(is.na(start_pos))) {
    stop("Failed to parse coordinates. Ensure format is 'chr:start-end'.")
  }
  
  GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges = IRanges::IRanges(start = start_pos, end = end_pos)
  )
}