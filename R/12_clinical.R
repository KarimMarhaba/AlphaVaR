#' Annotate Variants with ClinVar Data
#'
#' Enriches the variant data with clinical significance information from a ClinVar VCF file.
#' It maps variants based on genomic coordinates and flags pathogenic entries.
#'
#' @param obj An object of class \code{AlphaVarSet}.
#' @param vcf_file Character string. Path to a local ClinVar VCF file (can be .vcf or .vcf.gz).
#'   Ideally bgzipped and tabix-indexed (.tbi) for performance, but standard VCFs work too.
#' @param genome Character string. The genome build (default "hg38"). Used to ensure coordinate safety.
#'
#' @details
#' This function uses the \code{VariantAnnotation} package to query ClinVar.
#' It adds the following columns to \code{obj$data}:
#' \itemize{
#'   \item \code{clinvar_id}: The ClinVar variation ID.
#'   \item \code{clinical_significance}: The raw clinical significance string (e.g., "Pathogenic", "Benign").
#'   \item \code{clinvar_trait}: The associated disease or trait name.
#'   \item \code{is_pathogenic}: Boolean flag. TRUE if significance contains "Pathogenic" (and not "Conflicting").
#' }
#'
#' @return A modified \code{AlphaVarSet} object with clinical annotations appended to \code{obj$data}.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps mcols
#' @importFrom IRanges subsetByOverlaps
#' @importFrom dplyr mutate left_join select distinct group_by summarise
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom tidyr separate replace_na
#' @importFrom stringr str_detect
#' @export
av_annotate_clinvar <- function(obj, vcf_file, genome = "hg38") {

  # 1. Validation and Checks
  if (!inherits(obj, "AlphaVarSet")) stop("Input must be an AlphaVarSet object.")
  if (!file.exists(vcf_file)) stop("ClinVar VCF file not found at: ", vcf_file)

  # Check for both suggested packages
  if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    stop("Package 'VariantAnnotation' is required. Please install it via BiocManager::install('VariantAnnotation').")
  }
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Package 'Rsamtools' is required for reading tabix indexes. Please install it via BiocManager::install('Rsamtools').")
  }

  message(">>> Phase 4: Annotating variants with ClinVar data...")

  # 2. Parse Object Coordinates to GRanges
  if (!"coords" %in% colnames(obj$data)) stop("Column 'coords' missing in obj$data.")

  # Helper to parse coord string into BED-like columns for GRanges
  coords_df <- obj$data %>%
    dplyr::select(.data$variant_id, .data$coords) %>%
    tidyr::separate(.data$coords, into = c("chr", "range"), sep = ":") %>%
    tidyr::separate(.data$range, into = c("start", "end"), sep = "-", fill = "right") %>%
    dplyr::mutate(end = ifelse(is.na(.data$end), .data$start, .data$end))

  gr_query <- GenomicRanges::makeGRangesFromDataFrame(
    coords_df,
    keep.extra.columns = TRUE,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end"
  )

  # 3. Read ClinVar VCF (Smart Read)
  tbi_file <- paste0(vcf_file, ".tbi")
  has_index <- file.exists(tbi_file)

  # Note: ScanVcfParam is in VariantAnnotation, but TabixFile is in Rsamtools
  param <- VariantAnnotation::ScanVcfParam(
    info = c("CLNSIG", "CLNDN"),
    fixed = "ALT"
  )

  message(if(has_index) "    -> Index found. Performing targeted query (fast)." else "    -> No index found. Scanning full VCF (this may take a moment)...")

  if (has_index) {
    # CORRECTED: Use Rsamtools::TabixFile
    tfile <- Rsamtools::TabixFile(vcf_file)
    vcf <- VariantAnnotation::readVcf(tfile, genome = genome, param = param)
    vcf <- IRanges::subsetByOverlaps(vcf, gr_query)
  } else {
    vcf <- VariantAnnotation::readVcf(vcf_file, genome = genome, param = param)
    vcf <- IRanges::subsetByOverlaps(vcf, gr_query)
  }

  if (length(vcf) == 0) {
    message("    ! No ClinVar overlaps found.")
    obj$data$is_pathogenic <- FALSE
    obj$data$clinical_significance <- NA_character_
    obj$data$clinvar_trait <- NA_character_
    obj$data$clinvar_id <- NA_character_
    return(obj)
  }

  # 4. Extract Data from VCF Object
  info_data <- VariantAnnotation::info(vcf)
  gr_clinvar <- GenomicRanges::granges(vcf)

  # Extract Significance (Flatten list to string)
  clnsig_list <- info_data$CLNSIG
  if (inherits(clnsig_list, "CharacterList") || is.list(clnsig_list)) {
    clnsig_vec <- vapply(clnsig_list, function(x) paste(unique(x), collapse = "|"), character(1))
  } else {
    clnsig_vec <- as.character(clnsig_list)
  }

  # Extract Disease Name
  clndn_list <- info_data$CLNDN
  if (inherits(clndn_list, "CharacterList") || is.list(clndn_list)) {
    clndn_vec <- vapply(clndn_list, function(x) paste(unique(x), collapse = "|"), character(1))
  } else {
    clndn_vec <- as.character(clndn_list)
  }

  GenomicRanges::mcols(gr_clinvar)$clinvar_sig <- clnsig_vec
  GenomicRanges::mcols(gr_clinvar)$clinvar_trait <- clndn_vec
  GenomicRanges::mcols(gr_clinvar)$clinvar_id <- names(vcf)

  # 5. Overlap and Merge
  hits <- GenomicRanges::findOverlaps(gr_query, gr_clinvar)

  matches <- tibble::tibble(
    variant_id = gr_query$variant_id[S4Vectors::queryHits(hits)],
    clinvar_id = gr_clinvar$clinvar_id[S4Vectors::subjectHits(hits)],
    clinical_significance = gr_clinvar$clinvar_sig[S4Vectors::subjectHits(hits)],
    clinvar_trait = gr_clinvar$clinvar_trait[S4Vectors::subjectHits(hits)]
  )

  # Aggregate duplicates
  matches_agg <- matches %>%
    dplyr::group_by(.data$variant_id) %>%
    dplyr::summarise(
      clinvar_id = paste(unique(.data$clinvar_id), collapse = ";"),
      clinical_significance = paste(unique(.data$clinical_significance), collapse = ";"),
      clinvar_trait = paste(unique(.data$clinvar_trait), collapse = ";")
    )

  # 6. Update Object Data
  new_data <- obj$data %>%
    dplyr::left_join(matches_agg, by = "variant_id") %>%
    dplyr::mutate(
      is_pathogenic = stringr::str_detect(.data$clinical_significance, "(?i)Pathogenic") &
                      !stringr::str_detect(.data$clinical_significance, "(?i)Conflicting|Benign")
    ) %>%
    dplyr::mutate(is_pathogenic = tidyr::replace_na(.data$is_pathogenic, FALSE))

  obj$data <- new_data

  message(paste0("    -> Annotation complete. Found ", nrow(matches_agg), " ClinVar hits."))
  message(paste0("    -> ", sum(new_data$is_pathogenic), " variants flagged as Pathogenic."))

  return(obj)
}