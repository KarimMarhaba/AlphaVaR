
# AlphaVaR

**AlphaVaR** provides a robust statistical framework to analyze variant
effect scores derived from AlphaGenome. It bridges the gap between raw
deep learning scores and biological insight.

The package covers the full workflow from data ingestion to clinical
reporting: **Annotation** (Phase 1), **Validation** (Phase 2),
**Visualization** (Phase 3), and **Clinical Integration** (Phase 4).

## Workflow Overview

How `AlphaVaR` processes your data:

``` mermaid
graph TD
    %% Nodes and Styles
    classDef input fill:#e1f5fe,stroke:#01579b,stroke-width:2px;
    classDef object fill:#fff9c4,stroke:#fbc02d,stroke-width:2px,stroke-dasharray: 5 5;
    classDef calc fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px;
    classDef output fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px;

    %% 1. INPUT
    UserCSV[(" Input CSV<br/>(Scores, VariantID, Modality)")]:::input
    
    %% 2. INGESTION
    subgraph Ingestion ["1. Ingestion"]
        Import["av_import_csv()"]:::calc
    end

    %% 3. THE OBJECT
    AVSet{{" AlphaVarSet Object<br/>(The Central Hub)"}}:::object

    %% 4. PROCESSING & ANNOTATION
    subgraph Processing ["2. Processing & Annotation"]
        QC["QC & Validation<br/>av_qc_summary()"]:::calc
        Stats["Calc Enrichment<br/>av_calc_enrichment()"]:::calc
        Bio["Annotate Biology<br/>av_annotate_genes()<br/>av_annotate_clinvar()"]:::calc
    end

    %% 5. OUTPUTS
    subgraph Outputs ["3. Outputs & Insights"]
        Report[(" HTML/PDF Report<br/>av_create_report()")]:::output
        App[(" Shiny App Explorer<br/>av_run_shiny()")]:::output
        Tracks[("genome.bedgraph<br/>Browser Export")]:::output
        Plots[("📊 Plots<br/>Volcano, Ridges, Locus")]:::output
    end

    %% Connections
    UserCSV --> Import
    Import --> AVSet
    
    AVSet --> QC
    AVSet --> Stats
    AVSet --> Bio
    
    QC --> AVSet
    Stats --> AVSet
    Bio --> AVSet

    AVSet --> Report
    AVSet --> App
    AVSet --> Tracks
    AVSet --> Plots
```

Understanding the Workflow: The diagram illustrates the central role of
the AlphaVarSet object. Think of it as the “single source of truth.”

    Input: You start by loading raw data into this object (Blue).

    Cycle: You pass the object to analysis functions (Green), which add layers of information (stats, biology) and return the updated object.

    Output: At any point, you can generate reports or visualizations (Purple) based on the current state of the object.

## Installation

You can install the development version directly from GitHub. Because
`AlphaVaR` depends on Bioconductor packages, dependencies will be
handled automatically via the `Remotes` field.

``` r
# install.packages("devtools")
devtools::install_github("KarimMarhaba/AlphaVaR", dependencies = TRUE)
```

**Note:** Ensure you have configured your GitHub Personal Access Token
(PAT) if this is a private repository.

## Quick Start

``` r
library(AlphaVaR)

# 1. Import Data
av_set <- av_import_csv("path/to/your_scores.csv")

# 2. Run Analysis Pipeline
av_set <- av_set %>%
  av_calc_enrichment() %>%
  av_annotate_genes(txdb_package = "TxDb.Hsapiens.UCSC.hg38.knownGene")

# 3. Create Report
av_create_report(av_set, output_file = "Analysis_Report.html")
```

For more details, see
`vignette("workflow_tutorial", package = "AlphaVaR")`.
