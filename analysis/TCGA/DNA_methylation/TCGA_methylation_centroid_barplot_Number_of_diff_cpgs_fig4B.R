#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Article:
# Proteogenomic analysis of triple-negative breast cancer 
# identifies subtype-specific therapeutic vulnerabilities and 
# epigenetic immune suppression 
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Authors: 
# - Antonio Colaprico
# - Brian D. Lehmann 
# - Tiago C. silva
# - Xi S.Chen
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Date: 29 June 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# BEFORE YOU START
# To run this analysis you need to run previously the following code:
# - 01_download_data.R: to download DNAm data
# - TCGA_methylation_centroid.R: to perform the analysis used in this file
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(dplyr)
library(ggpubr)
library(readr)
library(FDb.InfiniumMethylation.hg19)
library(SummarizedExperiment)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(tidyverse)
library(dplyr)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/TCGA/"
dir.code <- "analysis/TCGA/DNA_methylation//"
dir.plots <- "plots/TCGA/DNA_methylation"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
diffmean.cut <- 1
fdr.cut <- 0.05
cMax <- 12
tssDistance <- 3000
chrToExclude <- c("chrX","chrY","chrM")

dir.output <-  "analysis_results/TCGA/dnam_results"
files <- dir(dir.output,pattern = "centroid_methy_testing",full.names = TRUE)
results <- plyr::ldply(
  .data = files,
  .fun = function(f){
    data <- read.csv(f) 
    data$status <- "Not significant"
    data$status <- NA
    data$status[data$logFC > 1 & data$adj.P.Val < 0.05] <- "Hypermethylated"
    data$status[data$logFC < -1 & data$adj.P.Val < 0.05] <- "Hypomethylated"
    df <- plyr::count( data$status) 
    df
    df$subtype <- stringr::str_extract(pattern = "BL1|BL2|M|LAR",basename(f))
    df %>% dplyr::filter(!is.na(x))
  })

results$subtype <- factor(
  results$subtype,levels = 
    results %>% dplyr::filter(x == "Hypomethylated") %>% dplyr::arrange(freq) %>% pull(subtype)
)
results$x <- factor(results$x,levels = c("Hypomethylated","Hypermethylated"))


p <- ggpubr::ggbarplot(
  results,
  x = "subtype"
  ,y = "freq",
  fill = "x", 
  color = "white",
  position = position_dodge(0.8)
) + xlab("") + ylab("Differently methylated CpGs (n)") + labs(fill = "") +
  scale_fill_manual(values = c("#3366FF","#FFCC66")) +
  rotate_y_text()

ggsave(
  filename = file.path(dir.plots,"Diff_methylated_probes_Figure4B.pdf"),
  plot = p,
  width = 5,
  height = 5
)
