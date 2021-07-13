#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# The following script will create heatmap with GSVA results
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
# Date: 21 Jun 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(TCGAbiolinks)
library(GenomicRanges)
library(ggpubr)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/GISTIC/"
dir.plots <- "plots/GISTIC/"
dir.analysis <- "analysis/TCGA/CNV"
dir.output <- "analysis_results/GISTIC/"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
source(file.path(dir.analysis,'VisualizeCNplot.R'))

# get gene information
genes <- TCGAbiolinks:::get.GRCh.bioMart("hg19") %>% 
  dplyr::select(external_gene_name, chromosome_name,start_position,end_position) %>% 
  dplyr::rename(GeneSymbol = external_gene_name, Chr = chromosome_name, Start = start_position, End = end_position) %>% 
  unique


GeneToVisualize <- c(
  "PIK3CA","KDR","KIT","PDGFRA","FGFR1","BCL11A",
  "MLH1","APC","PTEN",
  "CASP8","PIK3CA","EGFR","JAK2",
  "ERRB3","KRAS","CTNNB1","MLH1","APC",
  "CSF1R","PTCH1","PTEN","RB1",
  "FGFR3","KDR","KIT","PDGFRA",
  "PTEN","FGFR2","BRCA2","RB1",
  "MAP2K1","BRCA1","PIK3CA",
  "FGFR2","MLH1","FGFR2","PTEN",
  "BRCA2",
  "CCNE1","AKT3","ATR","MLH1",
  "PIK3R1","ETV6","CREBBP","PALB2","NOTCH2",
  "NOTCH3","BRD4","MYC","B2M","HLA-A","HLA-B","HLA-C",
  "CCND2","FOXM1","MYB","PIM1","BAP1","MYCN","MCL1",
  "CHEK1","NFIB","CDKN2B","STK11","EZH2","SUZ12","SMARCD1",
  "PBRM1","DPF3","SMARCA4","ARID1B"
)

GeneToVisualize <- sort(unique(GeneToVisualize))

FolderList <- list.dirs(dir.data)

cn.plot.list <- plyr::alply(
  .data = c("TNBCall", "BL1", "BL2", "LAR", "M"),
  .margins = 1,
  .fun = function(curSubtype){
    folderCur <- FolderList[grep(curSubtype, FolderList)]
    FolderListCur <- list.files(folderCur)
    
    fileCurScore <- FolderListCur[grep("scores.gistic", FolderListCur)]
    fileCurAmp <- FolderListCur[grep("amp_genes.conf_99", FolderListCur)]
    fileCurDel <- FolderListCur[grep("del_genes.conf_99", FolderListCur)]
    
    
    tabAnno_conf_Amp <- readr::read_delim(
      paste0(folderCur, "/", fileCurAmp),
      "\t",
      escape_double = FALSE,
      trim_ws = TRUE
    )
    
    tabAnno_conf_Del <- readr::read_delim(
      paste0(folderCur, "/", fileCurDel),
      "\t",
      escape_double = FALSE,
      trim_ws = TRUE
    )
    
    scores_filt_call <- readr::read_delim(
      paste0(folderCur, "/", fileCurScore),
      "\t",
      escape_double = FALSE,
      trim_ws = TRUE
    ) %>% as.data.frame %>% dplyr::mutate( Aberration.Kind = 0) %>% 
      dplyr::mutate(Aberration.Kind = replace(Type, Type == "Amp",1)) %>% 
      dplyr::mutate(Aberration.Kind = replace(Type, Type == "Del",-1)) %>%
      dplyr::rename("Region Start [bp]" = Start, "Region End [bp]" = End,"score" = "G-score") %>%
      dplyr::mutate(Ratio = score/`-log10(q-value)`)
    
    FDR.thresh <- 1 / 10 ^ 0.25
    
    scores_filt_call_Amp_full <- scores_filt_call %>% dplyr::filter(Type == "Amp")
    scores_filt_call_Del_full <- scores_filt_call %>% dplyr::filter(Type == "Del")
    scores_filt_FDR_Amp <- scores_filt_call_Amp_full %>% dplyr::filter(`-log10(q-value)` > -log10(FDR.thresh))
    scores_filt_FDR_Del <- scores_filt_call_Del_full %>% dplyr::filter(`-log10(q-value)` > -log10(FDR.thresh))
    
    
    write.csv(
      scores_filt_FDR_Amp,
      file = paste0(dir.output,paste0("TNBC_",curSubtype,"_scores_filt_FDR_Amp.csv"))
    )
    write.csv(
      scores_filt_FDR_Del,
      file = paste0(dir.output,paste0("TNBC_",curSubtype,"_scores_filt_FDR_Del.csv"))
    )
    
    genes_GR <-  makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)
    
    # Adding gene names to the gistic result
    df_GR <- scores_filt_call_Amp_full
    colnames(df_GR)[colnames(df_GR) == "Region Start [bp]"] <- "start"
    colnames(df_GR)[colnames(df_GR) == "Region End [bp]"] <- "end"
    df_GR_new <-  makeGRangesFromDataFrame(df_GR, keep.extra.columns = TRUE)
    hits_Amp <- findOverlaps(genes_GR, df_GR_new, type = "any")
    
    df_GR_ann_Amp <- cbind(df_GR[subjectHits(hits_Amp), ], genes[queryHits(hits_Amp), ])
    df_GR_ann_Amp <- df_GR_ann_Amp[order(df_GR_ann_Amp$score, decreasing = TRUE), ]
    df_GR_ann_Amp <- df_GR_ann_Amp[!duplicated(df_GR_ann_Amp$GeneSymbol), ]
    
    GenesinWidePeak_Amp <- NULL
    if (ncol(tabAnno_conf_Amp) > 1) {
      for (curColumn in grep("[[:alnum:]]*p|[[:alnum:]]q",colnames(tabAnno_conf_Amp),value = T)) {
        curGenesWidePeak <- intersect(
          as.matrix(tabAnno_conf_Amp[, curColumn]),
          df_GR_ann_Amp$GeneSymbol
        )
        print(sort(curGenesWidePeak))
        GenesinWidePeak_Amp <- c(curGenesWidePeak, GenesinWidePeak_Amp)
      }
      
      df_GR_ann_Amp_GeneWidePeak <- df_GR_ann_Amp %>% 
        dplyr::filter(GeneSymbol %in% GenesinWidePeak_Amp)
    } else if (ncol(tabAnno_conf_Amp) == 1) {
      df_GR_ann_Amp_GeneWidePeak <- tabAnno_conf_Amp
    }
    
    
    df_GR <- scores_filt_call_Del_full
    colnames(df_GR)[colnames(df_GR) == "Region Start [bp]"] <- "start"
    colnames(df_GR)[colnames(df_GR) == "Region End [bp]"] <- "end"
    df_GR$Chromosome <- gsub("chr0", "", df_GR$Chromosome)
    df_GR$Chromosome <- gsub("chr", "", df_GR$Chromosome)
    df_GR$Chromosome <- as.numeric(df_GR$Chromosome)
    df_GR <- df_GR[df_GR$Type %in% "Del", ]
    df_GR_new <- makeGRangesFromDataFrame(df_GR, keep.extra.columns = TRUE)
    hits_Del <- findOverlaps(genes_GR, df_GR_new, type = "any")
    df_GR_ann_Del <- cbind(df_GR[subjectHits(hits_Del), ], genes[queryHits(hits_Del), ])
    df_GR_ann_Del <- df_GR_ann_Del[order(df_GR_ann_Del$score, decreasing = TRUE), ]
    df_GR_ann_Del <- df_GR_ann_Del[!duplicated(df_GR_ann_Del$GeneSymbol), ]
    
    GenesinWidePeak_Del <- NULL
    if (ncol(tabAnno_conf_Del) > 1) {
      for (curColumn in grep("[[:alnum:]]*p|[[:alnum:]]q",colnames(tabAnno_conf_Del),value = T)) {
        
        curGenesWidePeak <-intersect(
          as.matrix(
            tabAnno_conf_Del[, curColumn]),
          df_GR_ann_Del$GeneSymbol
        )
        print(sort(curGenesWidePeak))
        GenesinWidePeak_Del <- c(curGenesWidePeak, GenesinWidePeak_Del)
      }
      
      df_GR_ann_Del_GeneWidePeak <- df_GR_ann_Del %>% 
        dplyr::filter(GeneSymbol %in% GenesinWidePeak_Del)
    } else if (ncol(tabAnno_conf_Del) == 1) {
      df_GR_ann_Del_GeneWidePeak <- tabAnno_conf_Del
    }
    
    
    write.csv(
      df_GR_ann_Amp_GeneWidePeak,
      file = paste0(dir.output,paste0("TNBC_",curSubtype,"_scores_GeneWidePeak_Amp.csv"))
    )
    write.csv(
      df_GR_ann_Del_GeneWidePeak,
      file = paste0(dir.output,paste0("TNBC_",curSubtype,"_scores_GeneWidePeak_Del.csv"))
    )
    
    transformationFormulaAxis2 <- rbind(scores_filt_FDR_Amp, scores_filt_FDR_Del)
    transformationFormulaAxis2_thresh <- mean(
      transformationFormulaAxis2$score / transformationFormulaAxis2$`-log10(q-value)`
    )

    p <- VisualizeCNplot(
      scores_filt_call = scores_filt_call,
      GeneToVisualize = GeneToVisualize,
      chrArms = TRUE,
      tabAnno_conf_Amp = tabAnno_conf_Amp,
      tabAnno_conf_Del = tabAnno_conf_Del,
      FDR.thresh = FDR.thresh,
      anno_Amp_GeneWidePeak = df_GR_ann_Amp_GeneWidePeak,
      anno_Del_GeneWidePeak = df_GR_ann_Del_GeneWidePeak,
      show.names = "significant",
      titleplot = paste0("TNBC ", curSubtype, " subtype")
    )
    
    p <- p + ggthemes::theme_base() + 
      labs(title = "",x = "") + 
      theme( axis.text.x = element_text(angle = 90, vjust = 0.5))   
    return(p)
  },.inform = TRUE
)
names(cn.plot.list) <-  c("TNBCall", "BL1", "BL2", "LAR", "M")


pmerged <- ggarrange(
  plotlist = cn.plot.list,
  labels = names(cn.plot.list),
  common.legend = TRUE,
  legend = "none",
  hjust = -1,
  font.label = list(size = 14, face = "bold"),
  ncol = 1,
  nrow = 5
)

ggsave(
  pmerged,
  filename = file.path(dir.plots,"S5h_CNplot_Gistic2_v5_full_bdl_rotate.pdf"),
  width = 12,
  height = 14
)
