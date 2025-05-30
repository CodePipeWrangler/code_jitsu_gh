# Epigenetic GWAS with Random Forests

This project builds upon the methods and findings from my dissertation research, which investigated the use of DNA methylation as a basis for genome-wide association studies (GWAS). Specifically, it explores how **5-methylcytosine (5-mC)** variation contributes to phenotypic traits and whether this epigenetic information can be harnessed to accelerate **genetic gains** in breeding programs.

---

## Overview

While traditional GWAS rely on genetic markers (e.g., SNPs), this study focuses on **epigenetic markers** derived from cytosine methylation (5-mC).  
These markers often remain undercharacterized and their contribution to phenotypic variation is less well understood.  
By applying both statistical and machine learning models, this project seeks to:

- Detect significant associations between methylation patterns and phenotypic traits.
- Compare model performance and efficiency.
- Improve interpretability and predictive utility of epigenetic variation.

---

## Key Contributions

- Developed and preprocessed custom input datasets through extensive **data wrangling and feature engineering**.
- **Projected genotypes** for recombinant inbred lines (RILs) by combining parental genotype data with recombination breakpoint maps — enabling downstream methylation analysis at inferred loci.
- Implemented both **joint-linkage mapping** (R) and **random forest models** (Python/R) for GWAS analysis.

*A detailed walkthrough of the dataset creation and RIL genotype projection process will be added soon.*

---

## Methods

- Input data includes:
  - **5-mC methylation profiles** from bisulfite sequencing
  - **Phenotypic trait data** from experimental populations
  - **Inferred RIL genotypes** via recombination-informed parental projection
- Tools and frameworks:
  - R (`lme4`, `data.table`) for statistical modeling
  - Python (`scikit-learn`, `pandas`, `numpy`) for machine learning-based GWAS

---

## Project Structure

```text
epiGWAS/
├── data/           # Methylation profiles, phenotype data, projected genotypes
├── scripts/        # GWAS, ML modeling, and data engineering code
├── results/        # Output files, figures, and summary stats
└── README.md       # Project overview

