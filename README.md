# 🧬 SMARCE1-hESC-RNAseq

**Reproducible bulk RNA-seq analysis pipeline investigating the transcriptomic consequences of SMARCE1 knockout and long 3′UTR deletion in human embryonic stem cells (hESCs).**

---

## Overview

This repository contains the computational workflow used to analyze bulk RNA sequencing data from **wild-type (WT)**, **SMARCE1 knockout (KO)**, and **long 3′UTR deletion (DEL)** human embryonic stem cells.

The pipeline performs:

* Differential gene expression analysis using **DESeq2**
* Gene Ontology (GO) enrichment analysis
* Gene Set Enrichment Analysis (GSEA)
* Semantic similarity reduction using **rrvgo**
* Generation of publication-quality figures and a reproducible **Quarto** report

The analysis is modular, allowing computationally intensive steps to be run once while enabling rapid report rendering during manuscript preparation.

---

## Repository Structure

```text
SMARCE1-hESC-RNAseq/
├── README.md
├── hesc_KO_DEL.qmd              # Main Quarto report
├── _quarto.yml
├── fast_render.Rproj
│
├── scripts/
│   ├── run_all.R
│   ├── 01_run_analysis.R
│   └── 02_run_enrichment.R
│
├── R/
│   └── plotting_functions.R
│
├── data/
│   ├── metadata.csv
│   ├── tx2gene.tsv
│   ├── quant_trimmed/           # Not tracked by Git
│   └── 5wk_del/                 # Not tracked by Git
│
├── objects/                     # Generated cache (ignored)
└── results/                     # Generated output (ignored)
```

---

## Analysis Workflow

```text
Salmon Quantifications
          │
          ▼
      tximport
          │
          ▼
       DESeq2
          │
          ├──────────────┐
          ▼              ▼
   Differential      Shrunken
    Expression      Log2 Fold Change
          │
          ▼
 GO Enrichment + GSEA
          │
          ▼
      rrvgo Summary
          │
          ▼
  Quarto Report & Figures
```

---

## Running the Pipeline

Run the complete analysis:

```r
source("scripts/run_all.R")
```

This performs the complete computational workflow and generates cached analysis objects in the `objects/` directory.

---

## Rendering the Report

After the analysis completes:

```bash
quarto render hesc_KO_DEL.qmd
```

or simply click **Render** in RStudio.

Because the computational analyses are cached, subsequent report renders are substantially faster and do not rerun DESeq2, GO, or GSEA unless the analysis scripts are executed again.

---

## Data

The repository includes:

* `metadata.csv` — sample metadata
* `tx2gene.tsv` — transcript-to-gene mapping used by `tximport`

Large intermediate analysis files, Salmon quantifications, cached objects, and generated results are intentionally excluded from version control.

---

## Software

Developed using R with the following major packages:

* DESeq2
* tximport
* clusterProfiler
* fgsea
* rrvgo
* tidyverse
* Quarto

---

## Reproducibility

The analysis follows a two-stage workflow:

1. **Analysis scripts** perform all computationally intensive calculations and save reusable analysis objects.
2. **Quarto** loads these saved objects to generate figures and tables without repeating expensive computations.

This design improves reproducibility while greatly reducing report rendering time during manuscript preparation.

---

## Author

**Gladys Hazitha Samuel, Ph.D.**

Computational & Molecular Biologist
University of Connecticut

---

## License

This repository is intended for academic research. A formal open-source license will be added prior to public release.
