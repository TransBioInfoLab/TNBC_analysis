# Introduction 

Repository with supplementary code for article:

> Proteogenomic analysis of triple-negative breast cancer identifies subtype-specific therapeutic vulnerabilities and epigenetic immune suppression

## Authors
 - Antonio Colaprico
 - Brian D. Lehmann 
 - Tiago C. Silva
 - Xi S.Chen

## Folder and files structure

- `01_data_download.R`: Download majority of the data used in the analysis. This is the first file to be executed.
- `data:` contains data used by the analysis code that are not retrieved with
`01_data_download.R`
- `analysis/TCGA:` contains code used to perform majority of TCGA analysis
- `analysis/CPTAC:` contains code used to perform majority of CPTAC analysis
- `analysis/PDTX:` contains code used to perform majority of PDTX analysis
- `analysis/depmap`: contains code used to perform majority of 
DepMap(The Cancer Dependency Map Project) analysis
- `analysis/GDSC`: contains code used to perform majority of GDSC (
Genomics of Drug Sensitivity in Cancer) analysis

# Instructions

- To run the analysis you will need a machine with R 4.0 and at least 32Gb. 

- We expect the working directory to be the upper level where you can see folder 
`analysis` and `data`. 

- Install all required packages, then run `01_data_download.R` first.
The data should be save under `data` directory.

- Then go to the `analysis` folder and select  the appropriate analysis.

## If you are having issues

For any reason, if you are not able to run the code, please, create an issue in GitHub.

# Code to install required packages

```{R, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.13",ask = FALSE) # Install last version of Bioconductor

list.of.packages <- c(
  "readr",
  "readxl",
  "plyr",
  "dplyr",
  "tidyr",
  "GenomicRanges",
  "SummarizedExperiment",
  "myGene",
  "illuminaHumanv4.db",
  "EnsDb.Hsapiens.v75",
  "FDb.InfiniumMethylation.hg19",
  "limma",
  "edgeR",
  "DESeq2",
  "sva",
  "stats",
  "GSVA",
  "GSEABase",
  "GSVAdata",
  "genefilter",
  "maftools",
  "circlize",
  "ComplexHeatmap",
  "TCGAbiolinks",
  "ELMER",
  "flowCore",
  "flowStats",
  "cowplot",
  "gridGraphics",
  "survminer",
  "survminer",
  "ChIPseeker",
  "qvalue",
  "EDASeq",
  "mygene",
  "ggVennDiagram",
  "ggcyto",
  "tidyverse",
  "CancerSubtypes",
  "gdata"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

devtools::install_github("igordot/msigdbr",ref = "60a89f3aa43fe0745905ab58f121d9d28601d1ad")
devtools::install_github("dviraran/xCell")
```

# Session information

```
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7
```
