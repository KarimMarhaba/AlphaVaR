# ALPHA_VAR ROADMAP: Development Strategy

## Objective
Enhance the R package `AlphaVaR` to become a comprehensive toolkit for analyzing deep-learning derived variant effects.
**Primary Goal:** Enable a medical doctoral thesis focusing on Human Accelerated Regions (HARs) and neuropsychiatric disorders.
**Target Audience:** Geneticists and clinicians with limited coding experience.

## Current State
* **Core:** Import, QC, Statistics (Fisher/Wilcoxon), Basic Plotting, HTML Reporting.
* **Status:** Stable, Tested, Documented.

---

## DEVELOPMENT PHASES

Please implement the following features in the strict order provided. Each phase builds upon the previous one.

### PHASE 1: Biological Annotation (Context is King)
**Goal:** Transform raw genomic coordinates into biological meaning.

1.  **Function `av_annotate_genes(obj, txdb_package = "TxDb.Hsapiens.UCSC.hg38.knownGene")`**
    * **Logic:** Use `GenomicRanges` to find the nearest gene (TSS) for each variant.
    * **Output:** Add columns `nearest_gene`, `dist_to_tss`, `region_type` (promoter, intron, intergenic) to `obj$data`.
    * **Dependency:** `GenomicRanges`, `AnnotationDbi`.
    -> DONE
    -> TODO: im Report anzeigen (grafisch/Tabelle)

2.  **Function `av_annotate_regions(obj, region_bed_file, label = "HAR")`**
    * **Logic:** Check overlaps between variants and a set of regions (e.g., HARs provided as BED).
    * **Output:** Add boolean/string column (e.g., `is_HAR`) to `obj$data`.
    * **Use Case:** User uploads a HAR list; package flags variants inside HARs.

### PHASE 2: Scientific Validation (The "Reviewer Pleaser")
**Goal:** Demonstrate that AlphaGenome adds value over established scores.

3.  **Function `av_compare_scores(obj, comparison_scores, method = "spearman")`**
    * **Input:** User provides a dataframe with `variant_id` and external scores (e.g., CADD, PhyloP).
    * **Logic:** Join with internal data. Calculate correlation.
    * **Viz:** Create Scatterplots (AlphaGenome vs. CADD) and Venn diagrams of top hits (e.g., "Top 1% in AlphaGenome but low in CADD").

### PHASE 3: Accessibility & Visualization (The "Wow" Factor)
**Goal:** Make the tool usable for non-coders and exportable to standard tools.

4.  **Function `av_export_browser_track(obj, output_file, format = "bedgraph")`**
    * **Logic:** Export scores so they can be loaded into UCSC Genome Browser or IGV.
    * **Detail:** Color-code tracks by score intensity.

5.  **Function `av_run_shiny()`**
    * **Logic:** Launch a local Shiny App included in the package.
    * **Features:**
        * Upload CSV.
        * Interactive Volcano Plot (click point -> show gene info).
        * Download filtered tables.

### PHASE 4: Clinical Integration (The "Medical Doctor" Features)
**Goal:** Direct clinical relevance.

6.  **Function `av_annotate_clinvar(obj)`**
    * **Logic:** Optional enrichment using ClinVar data (via API or offline dump).
    * **Output:** Flag variants known to be "Pathogenic" or "Likely Pathogenic".

---

## INSTRUCTIONS FOR AI DEVELOPER
* **Modular Approach:** When implementing a step, create a new file (e.g., `R/05_annotation.R`) unless it fits logically into existing modules.
* **TDD:** Always write the test case (`tests/testthat/`) *before* or immediately after the function implementation.
* **Object Integrity:** Ensure the `AlphaVarSet` object structure remains consistent. New annotations should be added as new columns in `obj$data`.
* **Documentation:** Update `roxygen2` headers immediately.

## NEXT ACTION
Start with **PHASE 1, Task 1 (av_annotate_genes)**.