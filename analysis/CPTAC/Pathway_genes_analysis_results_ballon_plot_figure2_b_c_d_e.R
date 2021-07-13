#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# This script creates the ballon plots (Fig2) with RNA-seq, Protein and phospo data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Article:
# Proteogenomic analysis of triple-negative breast cancer 
# identifies TNBC_subtype-specific therapeutic vulnerabilities and 
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Libraries                                                    |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
library(GSEABase)
library(ggpubr)
library(readr)
library(data.table)
library(scales)
library(ggpubr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Directories                                                    |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
dir.cptac <- "data/CPTAC/"
dir.plot.out <- file.path("plots/CPTAC/Ballon_plots/")

# create folder if it does not exist
for(p in grep("dir",ls(),value = T)) {
  dir.create(get(p),recursive = TRUE,showWarnings = FALSE)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Data                                                    |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

CPTAC_phospho_pathway_genes <-  readr::read_csv(
  file.path(dir.cptac, "CPTAC_genelist_prospective.csv")
) %>% as.data.frame() %>% dplyr::rename(protein = GENE, pathway = PATHWAY)

# For is AKT1xS477s, we just need the second part after x
CPTAC_phospho_pathway_genes <- CPTAC_phospho_pathway_genes %>% 
  tidyr::separate(col = "SITE",into = c("Gene","phosphosite"),sep = "x") 

# Change HLA-B to HLA.B
CPTAC_phospho_pathway_genes$protein <- gsub("-",".",CPTAC_phospho_pathway_genes$protein)

Genelist <- CPTAC_phospho_pathway_genes$protein
Genelist <- c(Genelist, c("CCND2", "AURKA", "AURKB"))
print(Genelist)

files.rnaseq <- dir(dir.cptac,full.names = T,pattern = "centroid_RNAseq")
dataGene_merged <- plyr::adply(files.rnaseq,.margins = 1,.fun = function(f){
  # Read with RNAseq results
  data <- readr::read_csv(
    file = f,
    col_types = readr::cols()
  ) %>% as.data.frame %>% dplyr::rename(GeneSymbol = X1,FDR = adj.P.Val)
  data$GeneSymbol <- gsub("-", ".", data$GeneSymbol)
  data <- data[order(abs(data$logFC), decreasing = TRUE), ]
  data <- data %>% 
    dplyr::select(c("GeneSymbol", "logFC", "FDR")) %>% 
    dplyr::filter(GeneSymbol %in% Genelist)
  data$Subtype <- stringr::str_extract(basename(f),"[[:alnum:]]*")
  data
},.id = NULL)

files.proteome <- dir(dir.cptac,full.names = T,pattern = "centroid_proteome")
dataProt_merged <- plyr::adply(files.proteome,.margins = 1,.fun = function(f){
  # Read with RNAseq results
  data <- readr::read_csv(
    file = f,
    col_types = readr::cols()
  ) %>% as.data.frame %>% dplyr::rename(GeneSymbol = X1,FDR = adj.P.Val)
  data$GeneSymbol <- gsub("-", ".", data$GeneSymbol)
  data <- data[order(abs(data$logFC), decreasing = TRUE), ]
  data <- data %>% 
    dplyr::select(c("GeneSymbol", "logFC", "FDR")) %>% 
    dplyr::filter(GeneSymbol %in% Genelist)
  data$Subtype <- stringr::str_extract(basename(f),"[[:alnum:]]*")
  data
},.id = NULL)


files.phospo <- dir(dir.cptac,full.names = T,pattern = "centroid_phosphosite")
dataPhospho_merged <- plyr::adply(files.phospo,.margins = 1,.fun = function(f){
  # Read with RNAseq results
  data <- readr::read_csv(
    file = f,
    col_types = readr::cols()
  ) %>% as.data.frame %>% dplyr::rename(FDR = adj.P.Val)
  data <- data %>% tidyr::separate(col = X1,into = c("GeneSymbol","Gene_site"),sep = "x")
  data$GeneSymbol <- gsub("-", ".", data$GeneSymbol)
  data <- data[order(abs(data$logFC), decreasing = TRUE), ]  
  data <- data %>% 
    dplyr::select(c("GeneSymbol", "Gene_site","logFC", "FDR")) %>% 
    dplyr::filter(GeneSymbol %in% Genelist)
  data$Subtype <- stringr::str_extract(basename(f),"[[:alnum:]]*")
  data
},.id = NULL)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plots
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plots per genes, genes in the x-axis, subtypes in the y-axis
plyr::a_ply(
  .data = unique(CPTAC_phospho_pathway_genes$protein),
  .margins = 1,
  .fun = function(curProtein){
    
    rna <- dataGene_merged[dataGene_merged$GeneSymbol %in% curProtein, ]
    if(nrow(rna) > 0){
      rna$DataType <- "R"
      colnames(rna) <- c("ID","logFC", "FDR", "Subtype","DataType")
    } else {
      rna <- NULL
    }
    
    prot <- dataProt_merged[dataProt_merged$GeneSymbol %in% curProtein, ]
    if(nrow(prot) > 0){
      prot$DataType <- "P"
      colnames(prot) <- c("ID","logFC", "FDR", "Subtype","DataType")
    } else{
      prot <- NULL
    }
    
    phospo <- dataPhospho_merged[dataPhospho_merged$GeneSymbol %in% curProtein, ]
    selected.sites <- CPTAC_phospho_pathway_genes[CPTAC_phospho_pathway_genes$protein == curProtein,]$phosphosite
    idx <- grep(pattern = paste(selected.sites,collapse = "|"),x = phospo$Gene_site,ignore.case = TRUE)
    if(length(idx) == 0){
      phospo <- NULL
    } else {
      phospo <- phospo[idx,]
      phospo$GeneSymbol <- NULL
      colnames(phospo) <- c("ID","logFC", "FDR", "Subtype")
      phospo$DataType <- paste0("pP\n",gsub(paste0(curProtein,"-"),"",phospo$ID))
    }
    
    dataCPTAC_merged_all <- rbind(rna, prot, phospo)
    
    if(is.null(dataCPTAC_merged_all)) return()
    
    dataCPTAC_merged_all$DataType <- factor(dataCPTAC_merged_all$DataType , levels = rev(unique(dataCPTAC_merged_all$DataType)))
    dataCPTAC_merged_all$subtype <- factor(dataCPTAC_merged_all$Subtype, levels = c("BL1", "BL2", "M", "LAR"))
    dataCPTAC_merged_all$logFC <- rescale(dataCPTAC_merged_all$logFC, to = c(-1, 1))
    
    p_gene_box <- ggballoonplot(
      data = dataCPTAC_merged_all,
      x = "Subtype",
      y = "DataType",
      size = 20,
      fill = "logFC",
      rotate.x.text = FALSE,
      show.label = FALSE
    ) + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(
      axis.text.y = element_text(
        angle = 0,
        hjust = 1,
        size = 20
      ),
      axis.text.x = element_text(angle = 90, size = 20)
    ) + ggthemes::theme_base() + xlab("") + ylab("") + 
      theme(
        legend.position = "none",plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()
      ) + ggtitle(curProtein)
    ggsave(
      p_gene_box,
      filename = paste0(dir.plot.out, curProtein, ".pdf"),
      width = 4,
      height = 4
    )
  },.progress = "time")

