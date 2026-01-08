# PROJECT CONTEXT: AlphaVaR (R Package)

## 1. Project Overview
**AlphaVaR** is an R package designed to analyze, visualize, and report variant effect scores derived from AlphaGenome (deep learning models).
**Goal:** Transform raw CSV outputs (variant scores) into biological insights (gene burden, tissue enrichment, genomic context).
**Status:** Functional. Core Stats & Visualization stable. **Phases 1 (Annotation), 2 (Validation), 3 (Accessibility/Viz), and 4 (Clinical Integration) are implemented.**

## 2. Core Data Structure (The "Contract")
The package relies on a central S3 class object: `AlphaVarSet`.
All analysis functions take this object as input and return it (potentially with updated slots).
*Note: Validation functions (Phase 2) return list objects to avoid bloating the main object with external comparison data.*

### Class: `AlphaVarSet`
Constructed via `R/00_classes.R` and `R/01_import.R`.

* **`$data` (tibble):** The main data storage. **Crucial:** Column names are standardized upon import.
    * **Core Columns:**
        * `variant_id` (chr): Unique ID.
        * `score` (num): The numeric effect score (standardized from input).
        * `modality` (chr): The assay or tissue type (standardized from input).
        * `gene` (chr): Target gene name (from input, optional).
        * `coords` (chr): Genomic interval string (e.g., "chr1:100-200").
    * **Annotation Columns (created by `av_annotate_*`):**
        * `nearest_gene` (chr): Entrez ID or Symbol of the nearest gene.
        * `dist_to_tss` (num): Distance to the nearest Transcription Start Site.
        * `region_type` (chr): Classification ('promoter', 'gene_body', 'intergenic').
        * *Dynamic columns:* Boolean flags for custom regions (e.g., `is_HAR`, `in_peak`).
    * **Clinical Columns (created by `av_annotate_clinvar`):**
        * `clinvar_id` (chr): The ClinVar variation ID.
        * `clinical_significance` (chr): Raw significance string (e.g., "Pathogenic|Likely pathogenic").
        * `clinvar_trait` (chr): Associated disease/trait name.
        * `is_pathogenic` (lgl): Boolean flag for pathogenic variants (excluding conflicting/benign).
* **`$stats` (list):** Storage for computed statistics to enable caching.
    * `$enrichment`: Result of Fisher's Exact Test.
    * `$wilcoxon`: Result of Wilcoxon Rank Sum Test.
    * `$params`: Parameters used for calculation.
* **`$meta` (list):** Metadata like filename, import timestamp, mapping used.

## 3. File Architecture & Key Functions
The package follows a strict modular structure in the `R/` directory.

### Data Ingestion
* **`R/01_import.R`**:
    * `av_import_csv(file, col_map)`: Reads CSV, standardizes columns, returns `AlphaVarSet`.

### Statistics (The Engine)
* **`R/03_core_stats.R`**:
    * `av_calc_enrichment(obj)`: Runs Fisher's/Wilcoxon tests. Stores results in `obj$stats`.

### Gene Analysis
* **`R/04_gene_analysis.R`**:
    * `av_aggregate_genes(obj)`: Aggregates variant scores per gene.
    * `av_get_variants_in_gene(obj)`: Drill-down accessor.

### Biological Annotation (Phase 1)
* **`R/05_annotation.R`**:
    * `av_annotate_genes(obj, txdb_package, ...)`: Maps variants to nearest gene/TSS using `GenomicFeatures`. Adds `nearest_gene`, `dist_to_tss`, `region_type`.
    * `av_annotate_regions(obj, regions, label)`: Flags variants overlapping with external regions (BED/GRanges). Adds boolean columns.

### Scientific Validation (Phase 2)
* **`R/08_validation.R`**:
    * `av_compare_scores(obj, external_df, ...)`: Joins internal variants with external scores (e.g., CADD) and calculates correlation. Returns a comparison list (does not modify object).
    * `av_plot_compare_scatter(comparison_result)`: Creates a scatterplot of Internal vs. External scores with correlation stats.
    * `av_analyze_overlap(comparison_result)`: Calculates "Top Hit" overlaps (intersection counts) for Venn-like analysis.
    * `av_plot_compare_venn(overlap_result)`: Visualizes overlap as a Venn diagram (ggplot2).

### Visualization & Exports (Phase 3)
* **`R/06_viz_distributions.R`**: Ridge plots, Ranking boxplots.
* **`R/06_viz_enrichment.R`**: Volcano plots.
* **`R/07_viz_genome.R`**: Locus plots (Manhattan-style).
* **`R/10_exports.R`**:
    * `av_export_browser_track(obj, output_file, ...)`: Exports BedGraph files for UCSC Genome Browser/IGV.
* **`R/11_shiny.R`**:
    * `av_run_shiny()`: Launches the interactive Shiny app located in `inst/shiny/alpha_app`.

### Clinical Integration (Phase 4)
* **`R/12_clinical.R`**:
    * `av_annotate_clinvar(obj, vcf_file, genome)`: Enriches variants with ClinVar data using `VariantAnnotation` (supports Tabix). Adds clinical columns and pathogenicity flags.

### Reporting
* **`R/09_report.R`**:
    * `av_create_report(obj)`: Renders interactive HTML report (RMarkdown).

### Utilities
* **`R/globals.R`**: Global variable declarations.
* **`R/99_utils.R`**: Central utils and pipe imports.

## 4. Dependencies
* **Core:** `dplyr`, `ggplot2`, `readr`, `rlang`, `stringr`, `tibble`, `tidyr`
* **Viz & App:** `ggridges`, `ggrepel`, `shiny`
* **Bioinformatics:**
    * *Imports:* `GenomicRanges`, `GenomicFeatures`, `IRanges`, `S4Vectors`, `methods`
    * *Suggests:* `TxDb.Hsapiens.UCSC.hg38.knownGene`, `rtracklayer` (for BED import), `VariantAnnotation` (for ClinVar), `Rsamtools` (for Tabix indexing)
* **Reporting:** `rmarkdown`, `DT`, `plotly`
* **Dev:** `testthat`, `devtools`, `roxygen2`

## 5. Development Guidelines
1.  **Strict TDD:** Every new feature must have a corresponding test in `tests/testthat/`.
2.  **No Namespace Pollution:** Never use `library()` in functions. Use `package::function()` or `@importFrom`.
3.  **Stability:** The `AlphaVarSet` structure is the source of truth.
4.  **R CMD Check:** Must pass with 0 Errors, 0 Warnings, 0 Notes.