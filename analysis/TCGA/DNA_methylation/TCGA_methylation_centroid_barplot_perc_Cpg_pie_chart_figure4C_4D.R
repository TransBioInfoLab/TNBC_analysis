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
library(ggplot2)
library(readr)
library(FDb.InfiniumMethylation.hg19)
library(SummarizedExperiment)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

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
load(file.path(dir.data,"TCGA_TNBC_192samples_Methylation.Rdata"))

diffmean.cut <- 1
fdr.cut <- 0.05
cMax <- 12
tssDistance <- 3000
chrToExclude <- c("chrX","chrY","chrM")


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

source(file.path(dir.code,'DMRbarplot.R'))

dir.output <-  "analysis_results/TCGA/dnam_results"
files <- dir(dir.output,pattern = "centroid_methy_testing",full.names = TRUE)

plots <- plyr::llply(files,.fun = function(f){
  subtype <- stringr::str_extract(pattern = "BL1|BL2|M|LAR",basename(f))
  
  DMR_results <- f %>% read.csv() %>% 
    as.data.frame 
  
  infiniumMethylation <- FDb.InfiniumMethylation.hg19::get450k() %>% as.data.frame()
  tabRanges <-  infiniumMethylation[DMR_results$X,]
  
  DMR_results <- data.frame(
    Probe = DMR_results$X,
    chr =  tabRanges$seqnames,
    start = tabRanges$start,
    end =  tabRanges$end,
    strand = "*",
    pvalue = DMR_results$P.Value,
    qvalue = DMR_results$adj.P.Val,
    meth.diff = DMR_results$logFC,
    stringsAsFactors = FALSE,
    row.names = DMR_results$X
  ) %>% dplyr::filter(chr %in% paste0("chr",c(1:22)))
  
  DMR_results.filtered <-  DMR_results %>% dplyr::filter(abs(meth.diff) > diffmean.cut & qvalue < fdr.cut) 
  
  DMR_results.filtered.up <- DMR_results.filtered %>% dplyr::filter(meth.diff > 0)
  DMR_results.filtered.down <- DMR_results.filtered %>% dplyr::filter(meth.diff < 0)
  
  probes <- BRCA.met %>% 
    subset(subset = as.character(rownames(BRCA.met)) %in% c(DMR_results.filtered.up$Probe)) %>%
    rowRanges()
  peakAnno <- annotatePeak(probes, tssRegion=c(-tssDistance, tssDistance), TxDb=txdb, annoDb="org.Hs.eg.db")
  toWrite_merged <- merge(
    x = DMR_results.filtered.up,
    y = as.data.frame(peakAnno),
    by.x = "Probe",
    by.y = "probeID"
  )
  
  write.csv(toWrite_merged, file = paste0(dir.plots,"/TNBC_",subtype,"_DMRs_annotationUp.csv" ))
  pdf(paste0(dir.plots,"/DMRs_Pies_TNBC_",subtype,"_DMRs_annotationUp.pdf"),width = 6, height = 6)
  plotAnnoPie(peakAnno) 
  dev.off()
  
  met.sel.Down <- subset(BRCA.met,subset = as.character(rownames(BRCA.met)) %in% c(DMR_results.filtered.down$Probe))
  probes <- rowRanges(met.sel.Down)
  peakAnno <- annotatePeak(probes, tssRegion=c(-tssDistance, tssDistance), TxDb=txdb, annoDb="org.Hs.eg.db")
  
  toWrite_merged <- merge(
    x = DMR_results.filtered.down,
    y = as.data.frame(peakAnno),
    by.x = "Probe",
    by.y = "probeID"
  )
  
  write.csv(toWrite_merged, file = paste0(dir.plots,"/TNBC_",subtype,"_DMRs_annotationDown.csv" ))
  pdf(paste0(dir.plots,"/DMRs_Pies_TNBC_",subtype,"_DMRs_annotationDown.pdf"),width = 6, height = 6)
  plotAnnoPie(peakAnno) 
  dev.off()
  
  p <- DMRbarplot(
    methdiff = DMR_results, 
    logFC = TRUE,
    plot = TRUE, 
    qvalue.cutoff = fdr.cut, 
    meth.cutoff = diffmean.cut,
    exclude = chrToExclude,
    subtitle = paste0("CpGs (", "TNBC others vs TNBC ",subtype,")"),
    cMax
  )
  list("barplot" = p)
},.progress = "time")
names(plots) <- stringr::str_extract(pattern = "BL1|BL2|M|LAR",basename(files))

barplot.percent.cpgs <- ggarrange(
  plots$BL1$barplot,
  plots$BL2$barplot,
  plots$LAR$barplot,
  plots$M$barplot,
  labels = LETTERS[1:4],
  ncol = 4, nrow = 1,
  common.legend = TRUE
)
ggsave(
  barplot.percent.cpgs , 
  filename = file.path(dir.plots,"TNBC_Fig4C_DMR_chr_barplot.pdf"),
  width = 25,
  height = 6
)
