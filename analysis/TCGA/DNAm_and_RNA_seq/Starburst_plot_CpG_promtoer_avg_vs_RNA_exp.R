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
# - TCGA_RNASeq_limma: to perform the analysis used in this file
# - TCGA_methylation_centroid_barplot_perc_Cpg_pie_chart_figure4C_4D.R: 
# to annotate the diff. methylated CpGs
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
library(plyr)
library(ggrepel)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/TCGA/"
dir.code <- "analysis/TCGA/DNA_methylation//"
dir.plots <- "plots/TCGA/DNA_methylation"
dir.output.rnaseq <-  "analysis_results/TCGA/rnaseq_results"
dir.output.dnam <-  "analysis_results/TCGA/dnam_results"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

files.results.dnam <- dir(dir.output.dnam,pattern = "centroid_methy_testing",full.names = TRUE)
file.resuls.rna.seq <- dir(dir.output.rnaseq,full.names = TRUE)

load(file.path(dir.data,"TCGA_TNBC_192samples_Methylation.Rdata"))
probes <- BRCA.met %>% rowRanges()

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tssDistance <- 3000
peakAnno <- annotatePeak(probes, tssRegion=c(-tssDistance, tssDistance), TxDb = txdb, annoDb = "org.Hs.eg.db")

logFC.rna.cut <- 1
logFC.dnam.cut <- 1
fdr.cut <- 0.05

plots <- plyr::llply(file.resuls.rna.seq,.fun = function(f.rna){
  
  subtype <- stringr::str_extract(pattern = "BL1|BL2|M|LAR",string = f.rna)
  f.dnam <- files.results.dnam[grep(subtype,basename(files.results.dnam))]
  
  # For each gene take the mean of all columns
  dnam.results <- read.csv(f.dnam) %>% 
    dplyr::rename(probeID = X, logFC.dnam = logFC, adj.P.Val.dnam = adj.P.Val)
  
  rnaseq.results <- read.csv(f.rna) %>% 
    dplyr::rename(SYMBOL = ID, logFC.rna = logFC,adj.P.Val.rna = adj.P.Val) %>% 
    dplyr::select(logFC.rna, SYMBOL, adj.P.Val.rna)
  
  merged <- dplyr::left_join(
    x = dnam.results,
    y = as.data.frame(peakAnno)
  ) 
  
  promoter.avg <- merged %>% 
    dplyr::filter(grepl("Promoter", annotation)) %>% 
    plyr::ddply("SYMBOL", numcolwise(mean))
  
  dnam.rna.results <- promoter.avg %>% dplyr::left_join(rnaseq.results)
  
  dnam.rna.results$status <- "Not significant"
  dnam.rna.results$status[dnam.rna.results$logFC.dnam > logFC.dnam.cut & dnam.rna.results$logFC.rna < (- 1 * logFC.rna.cut) & dnam.rna.results$adj.P.Val.rna < fdr.cut] <- "HyperMethylated & Gene Down-regulated"
  dnam.rna.results$status[dnam.rna.results$logFC.dnam < (-1 * logFC.dnam.cut) & dnam.rna.results$logFC.rna > logFC.rna.cut & dnam.rna.results$adj.P.Val.rna < fdr.cut] <- "HypoMethylated & Gene Up-regulated"
  
  
  p <- ggpubr::ggscatter(
    data = dnam.rna.results ,
    x = "logFC.dnam",
    y = "logFC.rna",
    fill = "status",
    palette = c("#ffcc00","#3366FF","black"),
    color = "status",
    xlab = "DNA Methylation (log2 FC)",
    ylab = "RNA (log2 FC)",
    title = paste0(subtype," methylation")
  ) + geom_text_repel(
    data = dnam.rna.results %>% dplyr::filter(grepl("Methylated",status)), 
    aes(label = SYMBOL), cex = 2, max.overlaps = 200
  )
  p
})
names(plots) <- subtype <- stringr::str_extract(pattern = "BL1|BL2|M|LAR",string = file.resuls.rna.seq)
merged <- ggpubr::ggarrange(plotlist = plots,ncol = 2,nrow = 2)
ggplot2::ggsave(plot = merged,filename = file.path(dir.plots,"Starburst_plot.pdf"),width = 10,height = 10)
