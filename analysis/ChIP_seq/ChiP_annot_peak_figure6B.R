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
# Date: 28 October 2020
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(ggrepel)
library(plyr)
library(ggpubr)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Paths
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.data <- "data/ChIP_seq"
dir.plots <- "plots/ChIP_seq"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Annotate ChIP peaks with chipseeker
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samplefiles <- dir(dir.data,recursive = TRUE,full.names = TRUE,pattern = "bed|txt")
names <- paste0(basename(dirname(samplefiles)),"_",stringr::str_extract(c("H3K27_DMSOvTAZ|H3K27_TAZvDMSO"),string = samplefiles))
samplefiles <- as.list(samplefiles)
names(samplefiles) <- names

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnnoList <- plyr::llply(
  samplefiles,
  annotatePeak,
  TxDb = txdb,
  tssRegion = c(-1000, 1000),
  verbose = FALSE,
  .progress = "time"
)

peakAnnoList
# Plot distribution of peaks
plotAnnoBar(peakAnnoList)

#################annotate tazemetostat treated sample
plyr::l_ply(names(peakAnnoList),.fun = function(name){
  peakAnnot <- data.frame(peakAnnoList[[name]]@anno)
  
  # Get the entrez IDs
  entrez <- peakAnnot$geneId
  
  # Return the gene symbol for the set of Entrez IDs
  annotations_edb <- AnnotationDbi::select(
    EnsDb.Hsapiens.v75,
    keys = entrez,
    columns = c("GENENAME"),
    keytype = "ENTREZID"
  )
  
  # Change IDs to character type to merge
  annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)
  
  ###fix colum names
  begin <- grep("strand",colnames(peakAnnot)) + 1
  end <- grep("annotation",colnames(peakAnnot)) - 1
  colnames(peakAnnot)[begin:end] = c(
    name,
    "score",
    ".",
    "enrichment",
    "p-value_-log10",
    "q-value_-log10",
    "q-value"
  )[1:(end - begin + 1)]
  
  # Write to file
  peakAnnot %>%
    left_join(annotations_edb, by = c("geneId" = "ENTREZID")) %>%
    readr::write_tsv(
      file = paste0(dir.data,"/",colnames(peakAnnot)[6],"_annotation.tsv"),
    )
},.progress = "time")

# Select promoter peaks and average by gene

files <- dir(dir.data,pattern = "DMSOvTAZ_annotation.tsv",full.names = TRUE)

plots <- plyr::llply(.data = files,.fun = function(f){
  
  DMSO_promoter <- readr::read_tsv(f) %>% 
    dplyr::mutate(geneId = as.character(geneId)) %>% dplyr::filter(annotation == "Promoter")
  
  DMSO_promoter_ave <- ddply(DMSO_promoter,"GENENAME",numcolwise(mean),drop=FALSE)
  
  ###change to negative
  DMSO_promoter_ave[,6] <- -DMSO_promoter_ave[,6]
  
  # Read in RNA testing
  RNA <- read.csv(file.path(dir.data,"../cell_lines/DE_anno_BT549.csv"), row.names=1, header=TRUE) %>% 
    dplyr::rename(GENENAME = gene_id) %>% dplyr::filter(adj.P.Val < 0.05)
  data <- merge(DMSO_promoter_ave, RNA, by.x = "GENENAME")
  
  MHC <- c(
    "NLRC5",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "HLA-H",
    "CIITA",
    "HLA-DPA1",
    "HLA-DQB1",
    "HLA-DMA",
    "HLA-DMB",
    "HLA-DRB1"
  )
  HOX <- read.csv(file.path(dir.data,"HOX_genes.csv"),  header=TRUE)
  ggplot(data, aes(x = logFC, y = enrichment))+
    geom_point(
      color = ifelse(
        data$GENENAME %in% c(HOX$GENES,MHC), 
        ifelse(data$GENENAME %in% MHC,"red","orange"),
        "grey50"
      ),
      size = 0.6
    ) + 
    xlab("RNAseq (logFC) Tazemetosat:DMSO") +
    ylab("H3K27me3 peak Tazemetosat:DMSO") +
    ylim(-7,-1)+
    geom_text_repel(
      data = data[data$GENENAME %in% MHC,], 
      aes(label = GENENAME),
      cex = 4,  
      max.overlaps = Inf,
      box.padding = 1.5,
      color = "black",
      segment.colour  = "black"
    ) +
    geom_text_repel(
      data = data[data$GENENAME %in% HOX$GENES,], 
      aes(label = GENENAME),
      cex = 4,  
      max.overlaps = Inf,
      box.padding = 1.5,
      color = "black",
      segment.colour  = "black"
    ) +
    theme(
      panel.grid = element_line(colour = "white"),
      panel.background = element_rect(fill = "white", colour = "white"),
      axis.line = element_line(size = 0.5, colour = "black"),
      axis.title = element_text(size = 20),
      aspect.ratio=1,
      axis.text = element_text(size = 20),
      panel.border = element_rect(linetype = "solid", fill = NA, size=0.5)
    )
},.inform = TRUE,.progress = "time")


pmerged <- ggarrange(
  plotlist = plots,
  ncol = 3, 
  nrow = 1,
  common.legend = FALSE
)

ggsave(
  pmerged , 
  filename = file.path(dir.plots,paste0("Scatter.pdf")),
  width = 20,
  height = 6
)
