# XYomics - R package for the analysis of condition-specific sex differences in omics data

## Table of contents
- [Introduction](#introduction)
- [Content](#content)
- [Data](#data)
- [Requirements](#requirements)
- [License](#license)
- [Instructions](#instructions)

## Introduction

This R package provides tools for analyzing condition-specific molecular sex differences in omics data. It enables identification of sex-specific, sex-dimorphic, and sex-modulated genes through differential expression and interaction analyses, accompanied by pathway enrichment and gene regulatory network analysis capabilities.

## Downloads

### Package File

You can download the latest version of the XYomics package here:

- [XYomics_0.1.1.tar.gz](https://gitlab.com/uniluxembourg/lcsb/bds/xyomics/-/raw/main/XYomics_0.1.1.tar.gz)

### Checksum

MD5 Checksum: `2c1bba2aa8291e4edfbb0e50c86b52a4`

## Content

Five main analysis modules (+ dedicated implementations for single-cell data analysis):

1.  **Differential Expression Analysis**: Identifies sex-specific and sex-dimorphic genes across conditions in bulk and single-cell data.

2.  **Sex Interaction Analysis**: Detects sex-modulated genes through interaction analysis in bulk and single-cell data.

3.  **Pathway Enrichment Analysis**: Performs pathway analysis on sex-biased gene sets for bulk and single-cell data.

4.  **Gene Regulatory Network (GRN) Analysis**: Constructs and analyzes condition-specific gene regulatory networks for single-cell data.

5.  **Plotting and Reporting Functions**: Generates visualizations and comprehensive analysis reports.

## Data

Includes example datasets and supports various omics data types (RNA-seq, microarray, proteomics). Input data should be formatted as described in the vignettes.

## Example output
 Please find below the vignettes for the XYomics R-packages:

- Vignette for sex-dependent analysis of bulk RNA-seq data: [XYomics_bulk_example](https://gitlab.com/uniluxembourg/lcsb/bds/xyomics/-/blob/main/xyomics/vignettes/XYomics_bulk_example.Rmd)
- Vignette for sex-dependent analysis of single-cell RNA-seq data: [XYomics_sc_example](https://gitlab.com/uniluxembourg/lcsb/bds/xyomics/-/blob/main/xyomics/vignettes/XYomics_sc_example.Rmd)

## Requirements

Required R version: ≥ 4.2.0

Key dependencies:

- **Bioconductor:** DESeq2, limma, edgeR, clusterProfiler, org.Hs.eg.db, ReactomePA  
- **CRAN:** ggplot2, dplyr, igraph  

## License

MIT License

## Instructions

### Install from GitLab

- For Ubuntu users, please install the following dependencies before running the installation script:
```
apt-get install libcurl-dev libcurl4-openssl-dev libudunits2-dev libgdal-dev libharfbuzz-dev libfribidi-dev libharfbuzz-dev libfribidi-dev
```

#### 1. Install the dependencies:

```{r}
install_XYpackages <- function() {
  # Helper to check which packages are not installed
  not_installed <- function(pkgs) {
    setdiff(pkgs, rownames(installed.packages()))
  }
  # Ensure BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  # Install Bioconductor packages
  bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "topGO", "edgeR", "DESeq2")
  to_install <- not_installed(bioc_pkgs)
  if (length(to_install) > 0) {
    BiocManager::install(to_install, update = FALSE, ask = FALSE)
    message(paste(to_install, collapse = ", "), " packages added...")
  }
  # Install CRAN helper packages
  cran_pkgs <- c("devtools", "remotes", "DT", "kableExtra")
  to_install_cran <- not_installed(cran_pkgs)
  if (length(to_install_cran) > 0) {
    install.packages(to_install_cran)
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    os_type <- Sys.info()[["sysname"]]
    message("Installing 'sf' package...")

    if (os_type == "Darwin") {
      # macOS
      install.packages("sf", type = "binary")
    } else if (os_type == "Linux" || os_type == "Windows") {
      # Linux or Windows
      install.packages("sf", type = "source")
    } else {
      warning("Unknown OS. Please install 'sf' manually.")
    }
  }
  # Install GitHub packages if missing
  github_pkgs <- list(
    PCSF = "IOR-Bioinformatics/PCSF",
    multienrichjam = "jmw86069/multienrichjam"
  )
  for (pkg in names(github_pkgs)) {
    if (pkg %in% not_installed(pkg)) {
      remotes::install_github(github_pkgs[[pkg]], upgrade = "never")
    }
  }
  # Message if nothing was added
  if (length(to_install) == 0 && length(to_install_cran) == 0 &&
      all(!names(github_pkgs) %in% not_installed(names(github_pkgs)))) {
    message("No new packages added...")
  }
}

install_XYpackages()
```
Install from GitLab:

    devtools::install_url("https://gitlab.com/uniluxembourg/lcsb/bds/xyomics/-/raw/main/XYomics_0.1.1.tar.gz")

Install from Local Download:
1. Download the package file: [XYomics_0.1.1.tar.gz](https://gitlab.com/uniluxembourg/lcsb/bds/xyomics/-/raw/main/XYomics_0.1.1.tar.gz)
2. Install using R:
   
   install.packages("path/to/XYomics_0.1.1.tar.gz", repos = NULL, type = "source")

The package installation is tested on the following operating systems 
- Windows 11
- MAC OS 
- Ubuntu 20.4