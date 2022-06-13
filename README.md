
# Introduction 

Repository with supplementary code for article:

> Proteogenomic analysis of triple-negative breast cancer identifies subtype-specific therapeutic vulnerabilities and epigenetic immune suppression

## Authors
 - Antonio Colaprico
 - Brian D. Lehmann 
 - Tiago C. Silva
 - Xi S.Chen

## How to cite this code:

[![DOI](https://zenodo.org/badge/347098639.svg)](https://zenodo.org/badge/latestdoi/347098639)


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

- To run the analysis you will need a machine with R 4.2 and at least 32Gb. 

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
BiocManager::install(version = "3.15",ask = FALSE) # Install last version of Bioconductor

list.of.packages <- c(
  "readr",
  "readxl",
  "plyr",
  "dplyr",
  "tidyr",
  "GenomicRanges",
  "SummarizedExperiment",
  "mygene",
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
**R version 4.1.0 (2021-05-18)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_US.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_US.UTF-8_, _LC_COLLATE=en_US.UTF-8_, _LC_MONETARY=en_US.UTF-8_, _LC_MESSAGES=en_US.UTF-8_, _LC_PAPER=en_US.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_US.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 

* grid 
* parallel 
* stats4 
* stats 
* graphics 
* grDevices 
* utils 
* datasets 
* methods 
* base 


**other attached packages:** 

* msigdbr(v.6.2.1) 
* pander(v.0.6.4) 
* gdata(v.2.18.0) 
* CancerSubtypes(v.1.18.0) 
* NMF(v.0.23.0) 
* cluster(v.2.1.2) 
* rngtools(v.1.5) 
* pkgmaker(v.0.32.2) 
* registry(v.0.5-1) 
* sigclust(v.1.1.0) 
* forcats(v.0.5.1) 
* stringr(v.1.4.0) 
* purrr(v.0.3.4) 
* tibble(v.3.1.4) 
* tidyverse(v.1.3.1) 
* ggcyto(v.1.20.0) 
* flowWorkspace(v.4.4.0) 
* ncdfFlow(v.2.38.0) 
* BH(v.1.75.0-0) 
* RcppArmadillo(v.0.10.6.0.0) 
* ggVennDiagram(v.1.1.4) 
* EDASeq(v.2.26.1) 
* ShortRead(v.1.50.0) 
* GenomicAlignments(v.1.28.0) 
* Rsamtools(v.2.8.0) 
* Biostrings(v.2.60.2) 
* XVector(v.0.32.0) 
* qvalue(v.2.24.0) 
* ChIPseeker(v.1.28.3) 
* survminer(v.0.4.9) 
* ggpubr(v.0.4.0) 
* ggplot2(v.3.3.5) 
* gridGraphics(v.0.5-1) 
* cowplot(v.1.1.1) 
* flowStats(v.4.4.0) 
* flowCore(v.2.4.0) 
* ELMER(v.2.16.0) 
* ELMER.data(v.2.16.0) 
* TCGAbiolinks(v.2.21.5) 
* ComplexHeatmap(v.2.8.0) 
* circlize(v.0.4.13) 
* maftools(v.2.8.0) 
* GSVAdata(v.1.28.0) 
* hgu95a.db(v.3.13.0) 
* GSEABase(v.1.54.0) 
* graph(v.1.70.0) 
* annotate(v.1.70.0) 
* XML(v.3.99-0.7) 
* GSVA(v.1.40.1) 
* sva(v.3.40.0) 
* BiocParallel(v.1.26.2) 
* genefilter(v.1.74.0) 
* mgcv(v.1.8-36) 
* nlme(v.3.1-152) 
* DESeq2(v.1.32.0) 
* edgeR(v.3.34.1) 
* limma(v.3.48.3) 
* FDb.InfiniumMethylation.hg19(v.2.2.0) 
* TxDb.Hsapiens.UCSC.hg19.knownGene(v.3.2.2) 
* EnsDb.Hsapiens.v75(v.2.99.0) 
* ensembldb(v.2.16.4) 
* AnnotationFilter(v.1.16.0) 
* illuminaHumanv4.db(v.1.26.0) 
* org.Hs.eg.db(v.3.13.0) 
* mygene(v.1.28.0) 
* GenomicFeatures(v.1.44.2) 
* AnnotationDbi(v.1.54.1) 
* SummarizedExperiment(v.1.23.4) 
* Biobase(v.2.52.0) 
* MatrixGenerics(v.1.4.3) 
* matrixStats(v.0.60.1) 
* GenomicRanges(v.1.44.0) 
* GenomeInfoDb(v.1.28.4) 
* IRanges(v.2.26.0) 
* S4Vectors(v.0.30.0) 
* BiocGenerics(v.0.38.0) 
* tidyr(v.1.1.3) 
* dplyr(v.1.0.7) 
* plyr(v.1.8.6) 
* readxl(v.1.3.1) 
* readr(v.2.0.1) 
* xCell(v.1.1.0) 


**loaded via a namespace (and not attached):** 

* graphlayouts(v.0.7.1) 
* lattice(v.0.20-44) 
* haven(v.2.4.3) 
* vctrs(v.0.3.8) 
* blob(v.1.2.2) 
* survival(v.3.2-13) 
* DBI(v.1.1.1) 
* R.utils(v.2.10.1) 
* SingleCellExperiment(v.1.14.1) 
* rappdirs(v.0.3.3) 
* gsubfn(v.0.7) 
* jpeg(v.0.1-9) 
* zlibbioc(v.1.38.0) 
* htmlwidgets(v.1.5.3) 
* mvtnorm(v.1.1-2) 
* GlobalOptions(v.0.1.2) 
* irlba(v.2.3.3) 
* DEoptimR(v.1.0-9) 
* tidygraph(v.1.2.0) 
* Rcpp(v.1.0.7) 
* KernSmooth(v.2.23-20) 
* DelayedArray(v.0.18.0) 
* RcppParallel(v.5.1.4) 
* Hmisc(v.4.5-0) 
* fs(v.1.5.0) 
* fastmatch(v.1.1-3) 
* digest(v.0.6.27) 
* png(v.0.1-7) 
* scatterpie(v.0.1.7) 
* DOSE(v.3.18.2) 
* ggraph(v.2.0.5) 
* pkgconfig(v.2.0.3) 
* GO.db(v.3.13.0) 
* gridBase(v.0.4-7) 
* DelayedMatrixStats(v.1.14.3) 
* iterators(v.1.0.13) 
* GetoptLong(v.1.0.5) 
* xfun(v.0.25) 
* zoo(v.1.8-9) 
* tidyselect(v.1.1.1) 
* reshape2(v.1.4.4) 
* pcaPP(v.1.9-74) 
* viridisLite(v.0.4.0) 
* rtracklayer(v.1.52.1) 
* rlang(v.0.4.11) 
* hexbin(v.1.28.2) 
* RVenn(v.1.1.0) 
* glue(v.1.4.2) 
* RColorBrewer(v.1.1-2) 
* modelr(v.0.1.8) 
* RProtoBufLib(v.2.4.0) 
* ggsignif(v.0.6.2) 
* DO.db(v.2.9) 
* jsonlite(v.1.7.2) 
* bit(v.4.0.4) 
* IDPmisc(v.1.1.20) 
* gridExtra(v.2.3) 
* gplots(v.3.1.1) 
* ConsensusClusterPlus(v.1.56.0) 
* stringi(v.1.7.4) 
* TCGAbiolinksGUI.data(v.1.12.0) 
* yulab.utils(v.0.0.2) 
* bitops(v.1.0-7) 
* cli(v.3.0.1) 
* rhdf5filters(v.1.4.0) 
* sqldf(v.0.4-11) 
* RSQLite(v.2.2.8) 
* rrcov(v.1.5-5) 
* data.table(v.1.14.0) 
* flowViz(v.1.56.0) 
* rstudioapi(v.0.13) 
* locfit(v.1.5-9.4) 
* ks(v.1.13.2) 
* VariantAnnotation(v.1.38.0) 
* survMisc(v.0.5.5) 
* R.oo(v.1.24.0) 
* dbplyr(v.2.1.1) 
* lifecycle(v.1.0.0) 
* munsell(v.0.5.0) 
* cellranger(v.1.1.0) 
* R.methodsS3(v.1.8.1) 
* hwriter(v.1.3.2) 
* caTools(v.1.18.2) 
* codetools(v.0.2-18) 
* MultiAssayExperiment(v.1.18.0) 
* htmlTable(v.2.2.1) 
* proto(v.1.0.0) 
* xtable(v.1.8-4) 
* abind(v.1.4-5) 
* farver(v.2.1.0) 
* km.ci(v.0.5-2) 
* aplot(v.0.1.0) 
* biovizBase(v.1.40.0) 
* ggtree(v.3.0.4) 
* BiocIO(v.1.2.0) 
* patchwork(v.1.1.1) 
* dichromat(v.2.0-0) 
* fda(v.5.1.9) 
* Matrix(v.1.3-4) 
* tidytree(v.0.3.4) 
* ellipsis(v.0.3.2) 
* prettyunits(v.1.1.1) 
* lubridate(v.1.7.10) 
* reprex(v.2.0.1) 
* mclust(v.5.4.7) 
* igraph(v.1.2.6) 
* fgsea(v.1.18.0) 
* htmltools(v.0.5.2) 
* BiocFileCache(v.2.0.0) 
* yaml(v.2.2.1) 
* utf8(v.1.2.2) 
* plotly(v.4.9.4.1) 
* fds(v.1.8) 
* hdrcde(v.3.4) 
* aws.s3(v.0.3.21) 
* foreign(v.0.8-81) 
* withr(v.2.4.2) 
* aroma.light(v.3.22.0) 
* bit64(v.4.0.5) 
* foreach(v.1.5.1) 
* robustbase(v.0.93-8) 
* ProtGenerics(v.1.24.0) 
* cytolib(v.2.4.0) 
* GOSemSim(v.2.18.1) 
* rsvd(v.1.0.5) 
* ScaledMatrix(v.1.0.0) 
* memoise(v.2.0.0) 
* evaluate(v.0.14) 
* rio(v.0.5.27) 
* geneplotter(v.1.70.0) 
* tzdb(v.0.1.2) 
* curl(v.4.3.2) 
* fansi(v.0.5.0) 
* checkmate(v.2.0.0) 
* cachem(v.1.0.6) 
* Gviz(v.1.36.2) 
* babelgene(v.21.4) 
* impute(v.1.66.0) 
* rjson(v.0.2.20) 
* openxlsx(v.4.2.4) 
* rstatix(v.0.7.0) 
* ggrepel(v.0.9.1) 
* clue(v.0.3-59) 
* tools(v.4.1.0) 
* rainbow(v.3.6) 
* magrittr(v.2.0.1) 
* RCurl(v.1.98-1.4) 
* car(v.3.0-11) 
* aws.signature(v.0.6.0) 
* ape(v.5.5) 
* ggplotify(v.0.1.0) 
* xml2(v.1.3.2) 
* httr(v.1.4.2) 
* assertthat(v.0.2.1) 
* rmarkdown(v.2.10) 
* boot(v.1.3-28) 
* R6(v.2.5.1) 
* Rhdf5lib(v.1.14.2) 
* nnet(v.7.3-16) 
* progress(v.1.2.2) 
* KEGGREST(v.1.32.0) 
* treeio(v.1.16.2) 
* gtools(v.3.9.2) 
* shape(v.1.4.6) 
* beachmat(v.2.8.1) 
* HDF5Array(v.1.20.0) 
* BiocSingular(v.1.8.1) 
* rhdf5(v.2.36.0) 
* splines(v.4.1.0) 
* carData(v.3.0-4) 
* ggfun(v.0.0.3) 
* colorspace(v.2.0-2) 
* generics(v.0.1.0) 
* base64enc(v.0.1-3) 
* pracma(v.2.3.3) 
* chron(v.2.3-56) 
* pillar(v.1.6.2) 
* Rgraphviz(v.2.36.0) 
* tweenr(v.1.0.2) 
* iCluster(v.2.1.0) 
* GenomeInfoDbData(v.1.2.6) 
* gtable(v.0.3.0) 
* rvest(v.1.0.1) 
* zip(v.2.2.0) 
* restfulr(v.0.0.13) 
* knitr(v.1.33) 
* latticeExtra(v.0.6-29) 
* shadowtext(v.0.0.8) 
* biomaRt(v.2.48.3) 
* fastmap(v.1.1.0) 
* Cairo(v.1.5-12.2) 
* doParallel(v.1.0.16) 
* broom(v.0.7.9) 
* BSgenome(v.1.60.0) 
* scales(v.1.1.1) 
* filelock(v.1.0.2) 
* backports(v.1.2.1) 
* plotrix(v.3.8-1) 
* enrichplot(v.1.12.2) 
* hms(v.1.1.0) 
* ggforce(v.0.3.3) 
* KMsurv(v.0.1-5) 
* polyclip(v.1.10-0) 
* lazyeval(v.0.2.2) 
* Formula(v.1.2-4) 
* crayon(v.1.4.1) 
* MASS(v.7.3-54) 
* downloader(v.0.4) 
* sparseMatrixStats(v.1.4.2) 
* viridis(v.0.6.1) 
* reshape(v.0.8.8) 
* rpart(v.4.1-15) 
* compiler(v.4.1.0) 
