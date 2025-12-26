# PROJECT CONTEXT: AlphaVaR (R Package)

## 1. Project Overview
**AlphaVaR** is an R package designed to analyze, visualize, and report variant effect scores derived from AlphaGenome (deep learning models).
**Goal:** Transform raw CSV outputs (variant scores) into biological insights (gene burden, tissue enrichment, genomic context).
**Status:** Functional, fully tested (testthat), documentation ready (roxygen2), 0 R CMD CHECK errors.

## 2. Core Data Structure (The "Contract")
The package relies on a central S3 class object: `AlphaVarSet`.
All analysis functions take this object as input and return it (potentially with updated slots).

### Class: `AlphaVarSet`
Constructed via `R/00_classes.R` and `R/01_import.R`.

* **`$data` (tibble):** The main data storage. **Crucial:** Column names are standardized upon import.
    * `variant_id` (chr): Unique ID.
    * `score` (num): The numeric effect score (standardized from input, e.g., 'quantile_score').
    * `modality` (chr): The assay or tissue type (standardized from input, e.g., 'output_type').
    * `gene` (chr): Target gene name (optional but recommended).
    * `coords` (chr): Genomic interval string (e.g., "chr1:100-200").
* **`$stats` (list):** Storage for computed statistics to enable caching.
    * `$enrichment`: Result of Fisher's Exact Test.
    * `$wilcoxon`: Result of Wilcoxon Rank Sum Test.
    * `$params`: Parameters used for calculation (thresholds, etc.).
* **`$meta` (list):** Metadata like filename, import timestamp, mapping used.

## 3. File Architecture & Key Functions
The package follows a strict modular structure in the `R/` directory.

### Data Ingestion
* **`R/01_import.R`**:
    * `av_import_csv(file, col_map = list(...))`: Reads CSV. **Key Logic:** Uses `col_map` to rename user columns to internal standards (`score`, `modality`, `gene`). Validates mandatory columns. Returns `AlphaVarSet`.

### Statistics (The Engine)
* **`R/03_core_stats.R`**:
    * `av_calc_enrichment(obj, group_by="modality", threshold=0.5)`: Runs Fisher's Exact Test (Enrichment/Depletion) and Wilcoxon (Shift). Stores results in `obj$stats`. Handles empty result sets gracefully.

### Gene Analysis
* **`R/04_gene_analysis.R`**:
    * `av_aggregate_genes(obj, method="max")`: Aggregates variant scores per gene. Methods: "max", "mean", "sum", "count_hits".
    * `av_get_variants_in_gene(obj, gene_name)`: Drill-down accessor.

### Visualization (ggplot2 based)
* **`R/06_viz_distributions.R`**:
    * `av_plot_ridges(obj)`: Density ridges of scores per modality.
    * `av_plot_ranking(obj)`: Boxplots sorted by median/effect.
* **`R/06_viz_enrichment.R`**:
    * `av_plot_volcano(obj)`: Plots Log2OR vs -Log10FDR. Requires `av_calc_enrichment` to be run first.
* **`R/07_viz_genome.R`**:
    * `av_plot_locus(obj, gene_name)`: "Manhattan-style" locus plot. Parses `$data$coords` string into chrom/start/end on the fly.

### Reporting
* **`R/09_report.R`**:
    * `av_create_report(obj, output_file, custom_template)`: Renders the RMarkdown template.
    * **Template:** `inst/rmd/report_master.Rmd`. Generates an interactive HTML (Plotly, DT).

### Utilities
* **`R/globals.R`**: Handles global variable declarations for `R CMD check` (dplyr/ggplot piping).
* **`R/99_utils.R`**: Central import for pipe `%>%`.

## 4. Dependencies
* **Core:** `dplyr`, `ggplot2`, `readr`, `rlang`, `stringr`
* **Viz:** `ggridges`, `ggrepel`
* **Reporting:** `rmarkdown`, `DT`, `plotly` (suggested/imported)
* **Dev:** `testthat`, `devtools`, `roxygen2`

## 5. Development Guidelines
1.  **Strict TDD:** Every new feature must have a corresponding test in `tests/testthat/`.
2.  **No Namespace Pollution:** Never use `library()` in functions. Use `package::function()` or `@importFrom`.
3.  **Stability:** The `AlphaVarSet` structure is the source of truth. Do not change internal column names (`score`, `modality`) without updating the entire chain.
4.  **R CMD Check:** Must pass with 0 Errors, 0 Warnings, 0 Notes.