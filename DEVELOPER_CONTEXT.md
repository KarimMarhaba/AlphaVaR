# PROJECT CONTEXT: AlphaVaR (R Package)

## 1. Project Overview
**AlphaVaR** is an R package designed to analyze, visualize, and report variant effect scores derived from AlphaGenome (deep learning models).
**Goal:** Transform raw CSV outputs (variant scores) into biological insights (gene burden, tissue enrichment, genomic context).
**Status:** Functional. **Phases 1 (Annotation), 2 (Validation), 3 (Accessibility/Viz), and 4 (Clinical Integration) are fully implemented.** PDF Export capability added.

### Workflow Pipeline
1.  **Ingest:** `av_import_csv()` (01) reads raw data -> creates `AlphaVarSet`.
2.  **Process:** `av_calc_enrichment()` (03) adds stats -> `av_annotate_*` (05/12) adds biology.
3.  **Validate:** `av_qc_summary()` (02) & `av_compare_scores()` (08) check quality.
4.  **Output:** User generates insights via `av_create_report()` (09), `av_run_shiny()` (11), or plots (06/07).

## 2. Core Data Structure (The "Contract")
The package relies on a central S3 class object: `AlphaVarSet`.
All analysis functions take this object as input and return it.

### Class: `AlphaVarSet`
Constructed via `R/00_classes.R` (structure/validation) and `R/01_import.R` (population).
It includes a custom S3 `print` method for console summaries.

* **`$data` (tibble):** The main data storage.
    * **Default Mapping (CSV Header -> Internal Name):**
        * `"quantile_score"` -> `score` (num)
        * `"output_type"` -> `modality` (chr)
        * `"variant_id"` -> `variant_id` (chr)
        * `"gene_name"` -> `gene` (chr, optional)
        * `"scored_interval"` -> `coords` (chr, optional)
    * **Strictly Validated Columns (by `validate_AlphaVarSet`):**
        * `variant_id` (chr): Unique ID.
        * `score` (num): The numeric effect score.
        * `modality` (chr): The assay or tissue type.
    * **Standard Columns (Populated by Import):**
        * `gene` (chr): Target gene name.
        * `coords` (chr): Genomic interval string.
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

### Infrastructure & Classes
* **`R/00_classes.R`**:
    * `AlphaVarSet(data, meta)`: Constructor helper that ensures valid structure.
    * `validate_AlphaVarSet(x)`: Enforces schema integrity (checks for required columns and numeric scores).
    * `print.AlphaVarSet(x)`: S3 method to display variant counts, modalities, and source metadata in the console.

### Data Ingestion
* **`R/01_import.R`**:
    * `av_import_csv(file, col_map = list(), sep = ",")`: 
        * Merges user `col_map` with defaults using `utils::modifyList`.
        * Checks for mandatory CSV columns before renaming.
        * Safely converts `score` to numeric (warns on NAs).
        * Stores metadata: filename, timestamp (`imported_at`), and mapping used (`col_map_used`).

### Quality Control
* **`R/02_qc.R`**:
    * `av_qc_summary(obj)`: Returns basic metrics (variant count, modalities, genes, mean score, missing values) as a tibble to verify import health.
    * `av_plot_qc_dist(obj, type = "histogram")`: Visualizes raw score distribution (histogram or density) to spot artifacts or skew.

### Statistics (The Engine)
* **`R/03_core_stats.R`**:
    * `av_calc_enrichment(obj, group_by = "modality", threshold = 0.5, run_wilcoxon = TRUE)`:
        * **Fisher's Exact Test:** Computes enrichment of high scores (`abs(score) >= threshold`) vs background. Includes Haldane correction for Log2 Odds Ratio to handle zeros.
        * **Wilcoxon Rank Sum:** Compares distribution of one group vs. all others (Group vs Rest). Includes safety checks (requires >= 5 obs).
        * Stores results in `$stats$enrichment` and `$stats$wilcoxon` (tibbles with FDR adjustments).

### Gene Analysis
* **`R/04_gene_analysis.R`**:
    * `av_aggregate_genes(obj, method = "max", threshold = 0.5)`: 
        * Aggregates variant scores per gene using absolute values.
        * Supports methods: `"max"` (single highest impact), `"mean"` (average burden), `"sum"` (total burden), and `"count_hits"` (count of variants >= threshold).
    * `av_get_variants_in_gene(obj, gene_name)`: Drill-down accessor that returns all variants for a specific gene, sorted by absolute score (descending).

### Biological Annotation (Phase 1)
* **`R/05_annotation.R`**:
    * `av_annotate_genes(obj, txdb_package, promoter_upstream = 2000, promoter_downstream = 200)`: 
        * Maps variants to nearest gene using `GenomicFeatures` (TxDb).
        * **Hierarchy:** Promoter > Gene Body > Intergenic.
        * **Outputs:** Adds `nearest_gene` (ID), `dist_to_tss` (numeric), and `region_type`.
        * Internally parses "chr:start-end" coordinate strings.
    * `av_annotate_regions(obj, regions, label)`: 
        * Flags variants overlapping with external regions.
        * Inputs: `GRanges` object or path to BED file (requires `rtracklayer`).
        * **Output:** Adds a boolean column named via `label` (e.g., `is_HAR`).

### Scientific Validation (Phase 2)
* **`R/08_validation.R`**:
    * `av_compare_scores(obj, external_df, ...)`: 
        * Performs an **inner join** with external data (e.g., CADD, PhyloP).
        * Calculates correlation (Spearman/Pearson) and returns a standalone list result (does not modify the main `AlphaVarSet`).
    * `av_plot_compare_scatter(comparison_result)`: 
        * Scatterplot with linear regression line (`geom_smooth`). 
        * Automatically annotates the plot with correlation stats (rho, p-value).
    * `av_analyze_overlap(comparison_result, top_percentile, high_is_good)`: 
        * Calculates intersection of top hits to check for model agreement.
        * **Smart feature:** `high_is_good` parameter handles metrics where low values indicate high impact (e.g., p-values vs. raw scores).
    * `av_plot_compare_venn(overlap_result)`: 
        * Draws a 2-circle Venn diagram using `ggplot2::geom_polygon` (custom implementation to avoid heavy external dependencies).

### Visualization & Exports (Phase 3)
* **`R/06_viz_distributions.R`**:
    * `av_plot_ridges(obj, abs_score = TRUE)`: Ridge plot (joyplot). **Smart feature:** Automatically sorts modalities by biological impact (median shift) if Wilcoxon stats exist; otherwise sorts alphabetically.
    * `av_plot_ranking(obj, abs_score = TRUE)`: Boxplots ranked by effect size. Includes a global median reference line.
    * `av_plot_ecdf(obj, abs_score = TRUE)`: Cumulative distribution function (ECDF) to compare score saturation across assays.
* **`R/06_viz_enrichment.R`**:
    * `av_plot_volcano(obj, label_top_n = 10)`: Visualizes Fisher test results (Log2OR vs -Log10 FDR).
        * **Features:** Color-codes groups as "Enriched"/"Depleted" (FDR < 0.05, |Log2OR| > 1), maps point size to variant count (`n_obs`), and labels top hits via `ggrepel`.
        * **Pre-requisite:** Requires `av_calc_enrichment` results in `$stats`.
* **`R/07_viz_genome.R`**:
    * `av_plot_locus(obj, gene_name, flank_bp = 0, color_by = "modality")`:
        * Visualizes spatial distribution of variants along the genome ("Lollipop" or Manhattan-style track).
        * **Mechanism:** Internally parses `coords` string (e.g., "chr1:100-200") via regex, filters by gene, and facets plots by modality.
        * **Requirement:** The `coords` column must exist (mapped during import).
* **`R/10_exports.R`**:
    * `av_export_browser_track(obj, output_file, track_name = "AlphaVaR_Scores", description, color)`: 
        * Exports variant scores to **BedGraph format** for UCSC Genome Browser or IGV.
        * **Process:** Parses `coords` string to extract chrom/start/end, generates a custom track header (visibility, color), and writes tab-separated values.
* **`R/11_shiny.R`**:
    * `av_run_shiny()`: Launches the interactive Shiny dashboard ("AlphaVaR Explorer").
        * **Features:** Enables interactive exploration (upload CSV, filter results, view Volcano plots) without coding.
        * **Location:** Runs the app bundle stored in `inst/shiny/alpha_app`.

### Clinical Integration (Phase 4)
* **`R/12_clinical.R`**:
    * `av_annotate_clinvar(obj, vcf_file, genome = "hg38")`:
        * Annotates variants using a local ClinVar VCF.
        * **Performance:** Detects `.tbi` index to use `Rsamtools::TabixFile` for fast, targeted queries; falls back to full scan if missing.
        * **Workflow:** Parses `coords` to `GRanges`, overlaps with VCF, and aggregates multiple hits per variant.
        * **Outputs:** Adds `clinvar_id`, `clinical_significance`, `clinvar_trait`, and a calculated `is_pathogenic` boolean flag (strict filter: requires "Pathogenic", excludes "Conflicting").

### Reporting
* **`R/09_report.R`**:
    * `av_create_report(obj, output_file, open = TRUE, custom_template = NULL)`:
        * Renders a comprehensive HTML report via `rmarkdown`.
        * **Smart feature:** Automatically calculates statistics (`av_calc_enrichment`) if missing.
        * **Template Logic:** Auto-detects `report_master.Rmd` in package installation or dev directories.
    * `av_export_pdf(obj, output_file, open = TRUE)`:
        * Converts the HTML report to a static PDF using headless Chrome.
        * **Dependency:** Requires `pagedown` and a local Chrome/Chromium installation.
        * **Workflow:** Generates a temporary HTML, prints to PDF, and cleans up.

### Utilities
* **`R/globals.R`**: 
    * Registers variable names (e.g., `variant_id`, `score`, `gene`) used in non-standard evaluation (NSE) to prevent false positives ("no visible binding") during R CMD Check.
* **`R/99_utils.R`**: 
    * Central import hub (especially for the magrittr pipe `%>%` via `dplyr`) to ensure operators are available package-wide without namespace collisions.

## 4. Dependencies
* **Core:** `dplyr`, `ggplot2`, `readr`, `rlang`, `stringr`, `tibble`, `tidyr`
* **Viz & App:** `ggridges`, `ggrepel`, `shiny`
* **Bioinformatics:**
    * *Imports:* `GenomicRanges`, `GenomicFeatures`, `IRanges`, `S4Vectors`, `methods`
    * *Suggests:* `TxDb.Hsapiens.UCSC.hg38.knownGene`, `rtracklayer` (for BED import), `VariantAnnotation` (for ClinVar), `Rsamtools` (for Tabix indexing)
* **Reporting:** `rmarkdown`, `DT`, `plotly`, `pagedown` (for PDF export)
* **Dev:** `testthat`, `devtools`, `roxygen2`

## 5. Development Guidelines
1.  **Strict TDD:** Every new feature must have a corresponding test in `tests/testthat/`.
2.  **No Namespace Pollution:** Never use `library()` in functions. Use `package::function()` or `@importFrom`.
3.  **Stability:** The `AlphaVarSet` structure is the source of truth.
4.  **R CMD Check:** Must pass with 0 Errors, 0 Warnings, 0 Notes.