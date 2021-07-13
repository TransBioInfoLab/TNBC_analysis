#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# The following script will create Fig1C. TCGA heatmap with RNA-seq, CNA, mut,
# Clinical information and others
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
library(maftools)
library(dplyr)
library(SummarizedExperiment)
library(readr)
library(ComplexHeatmap)
library(circlize)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Folders where we was saved the data
# See script: 01_data_download.R
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data,"TCGA")
dir.Gistic2 <-  file.path(dir.data,"GISTIC/")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o colors                                                        |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# RPPA Heatmap
pal_protein <- colorRampPalette(
  c('#3361A5', '#248AF3', '#14B3FF', 
    '#88CEEF', '#C1D5DC', '#EAD397', 
    '#FDB31A','#E42A2A', '#A31D1D')
)(100)

pal_rna <- colorRampPalette(
  c("#352A86","#343DAE","#0262E0","#1389D2",
    "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D")
)(100)

colors.tnbc <- c(
  "BL1" = "red",
  "BL2" = "orange",
  "LAR" = "green",
  "M" = "blue"
)

colors.pam50 <- c(
  "Basal" = "red",
  "Her2" = "lightpink1",
  "LumA" = "blue",
  "LumB" = "deepskyblue",
  "Normal" = "limegreen",
  "Normal-like" = "limegreen"
)

colors.pathways <- c(
  "Ag presentation" = "#165713",
  "Ag\npresentation" = "#165713",
  "cell cycle" = "#d27400",
  "Cell cycle" = "#d27400",
  "EMT" = "#5ca8ee",
  "EMT/wound repair" = "#5ca8ee",
  "AR" = "#1eff93",
  "AR signaling" = "#1eff93",
  "Immune" = "#0AF9F6",
  "Immune checkpoint" = "blue",
  "Immune\ncheckpoint" = "blue",
  "xCell sig." = "#4363d8",
  "RTK/PI3K/mTOR" = "pink",
  "chromatin" = "gold",
  "adhesion/motility growth factor" = "red",
  "Growth factor" = "red",
  "development" = "grey",
  "Development" = "grey",
  "notch" = "brown",
  "other" = "yellow",
  "growth factor" = "#34b3b0",
  "MAPK" = "black",
  "DNA repair" = "orange",
  "adhesion/motility" = 'blue',
  "TF" = "purple",
  "PI3K" = "pink",
  "Miscellaneous" = "#038074",
  "Misc." = "#038074"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Data                                                                    |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# Load supplemental table 2
metadata <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) 

# Load RNA-SEQ 
load(file.path(dir.data,"TCGA/TCGA_TNBC_192samples_RNASeq.Rdata"))
colnames(dataFilt) <- substr(colnames(dataFilt), 1, 12)
dataFilt_TNBC <- log2(dataFilt[, S2$patient] + 1)

# Gene list for CNA/mutation
gene.cna.mut <- read_csv(file = "data/gene_list_cna_mutation_fig1B.csv") %>% as.data.frame %>% 
  dplyr::rename(Hugo_Symbol = Genes)
rownames(gene.cna.mut) <- gene.cna.mut$Hugo_Symbol
GeneListCN <- unique(gene.cna.mut$Hugo_Symbol)


load(file.path(dir.tcga,"TCGA_TNBC_192samples_mc3_maf.Rdata")) 

mc3.maf.tnbc <- mc3.maf.tnbc %>%
  dplyr::select(c("Hugo_Symbol", "patient", "Variant_Classification")) %>% 
  dplyr::filter(!Variant_Classification %in% c("Silent","3'Flank","3'UTR", "5'Flank","5'UTR","Intron")) %>%
  dplyr::left_join(y = S2) %>% 
  dplyr::filter(Hugo_Symbol %in% GeneListCN) %>% 
  dplyr::left_join(
    y = gene.cna.mut
  )
mc3.maf.tnbc <- mc3.maf.tnbc %>% dplyr::select(c("Hugo_Symbol", "patient", "Variant_Classification"))
# Keep only samples with information in MC3 MAF
S2 <- S2 %>% dplyr::filter(patient %in% mc3.maf.tnbc$patient)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Sort by correlation within each subtype 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
OrderSubtypes <- c("BL1", "BL2", "M", "LAR")
S2 <- plyr::adply(OrderSubtypes,.margins = 1,.fun = function(s) {
  df_sel <- S2 %>% dplyr::filter(subtype %in% s)
  df_sel[order(df_sel[, s], decreasing = TRUE), ]
})
rownames(S2) <-  S2$patient

# Final Heatmap order
sampleOrderMut <- S2$patient

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Adding  copy number results from GISTIC2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
Gistic2.files <- dir(dir.Gistic2,pattern = "all_data_by_genes",full.names = TRUE,recursive = T)

Gistic2.data <- lapply(Gistic2.files,FUN = function(f){readr::read_tsv(f)}) %>% 
  purrr::reduce(inner_join) %>% as.data.frame()

rownames(Gistic2.data) <- Gistic2.data$`Gene Symbol`
Gistic2.data <- Gistic2.data[,-c(1:3)]
colnames(Gistic2.data) <- substr(colnames(Gistic2.data), 1, 12)

missing.samples <- sampleOrderMut[!sampleOrderMut %in% colnames(Gistic2.data)]

# Create an empty matrix for the missing samples
missing.gistic2 <- matrix(0, nrow(Gistic2.data), length(missing.samples)) %>% as.data.frame()
rownames(missing.gistic2) <- rownames(Gistic2.data)
colnames(missing.gistic2) <- sampleOrderMut[!sampleOrderMut %in% colnames(Gistic2.data)]
Gistic2.data <- cbind(Gistic2.data,missing.gistic2)


dataCN_TNBC_melt_v1 <- cbind(Hugo_Symbol = rownames(Gistic2.data),Gistic2.data)
dataCN_TNBC_melt <- reshape2::melt(dataCN_TNBC_melt_v1)
colnames(dataCN_TNBC_melt) <- colnames(mc3.maf.tnbc)
dataCN_TNBC_melt$AMP_DEL <- rep("-", nrow(dataCN_TNBC_melt))

# thresholds
copyNumber_del_thresh <- -.7
copyNumber_amp_thresh <- 1

dataCN_TNBC_melt[dataCN_TNBC_melt$Variant_Classification > copyNumber_amp_thresh, "AMP_DEL"] <- "AMP"
dataCN_TNBC_melt[dataCN_TNBC_melt$Variant_Classification < copyNumber_del_thresh, "AMP_DEL"] <- "DEL"
dataCN_TNBC_melt$Variant_Classification <- dataCN_TNBC_melt$AMP_DEL
dataCN_TNBC_melt <- dataCN_TNBC_melt[dataCN_TNBC_melt$Variant_Classification %in% c("AMP", "DEL"), ]

dataCN_TNBC_melt <-  subset(dataCN_TNBC_melt, select = colnames(mc3.maf.tnbc))

curMAF_merged <- rbind(mc3.maf.tnbc, dataCN_TNBC_melt)
curMAF_merged$patient <- substr(curMAF_merged$patient,1,12)
curMAF <- curMAF_merged

MAF_ch <- matrix(0, length(unique(curMAF$Hugo_Symbol)), length(unique(curMAF$patient)))
MAF_ch <- as.data.frame(MAF_ch)
rownames(MAF_ch) <- unique(curMAF$Hugo_Symbol)
colnames(MAF_ch) <- unique(curMAF$patient)

# order genes by pathways
gene.cna.mut <- gene.cna.mut[gene.cna.mut$Hugo_Symbol %in% rownames(MAF_ch), ]
MAF_ch <- MAF_ch[gene.cna.mut$Hugo_Symbol, ]

# Fill matrix 
for (gene in rownames(MAF_ch)) {
  # For each gene get the samples
  curMAF_sel <- curMAF[curMAF$Hugo_Symbol %in% gene, ]
  sampleList_cur <- unique(curMAF_sel$patient)
  
  # For each patient fill the matrix  
  for (sample in sampleList_cur) {
    curMAF_sel_sub_sample <- curMAF_sel[curMAF_sel$patient %in% sample, ]
    MAF_ch[gene, sample] <- paste(curMAF_sel_sub_sample$Variant_Classification, collapse = ";")
  }
}
MAF_ch[MAF_ch == 0] <- ""
MAF_ch <- MAF_ch[, sampleOrderMut]

# creating alter fun for complexHeatmap OncoPrint
col = c(
  "Missense_Mutation" = "#CC3333",
  "Frame_Shift_Del" = "gold",
  "Splice_Site" = "green",
  "Nonsense_Mutation" = "blue",
  "In_Frame_Del" = "orange",
  "Frame_Shift_Ins" = "#333399",
  "RNA" = "#CC3333",
  "MUT" = "black",
  "DEL" = "cornflowerblue",
  "AMP" = "#FF6666"
)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = "#FFFFFF", col = NA))
  },
  
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = col["DEL"], col = NA))
  },
  
  MUT = function(x, y, w, h) {
    grid.rect(x,  y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = col["MUT"], col = NA))
  },
  
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  RNA = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = col["RNA"], col = NA))
  },
  
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = col["AMP"], col = NA))
  }
  
)



# reduce barcode names
colnames(MAF_ch) <- substr(colnames(MAF_ch), nchar("TCGA-A2-A"), 12)
sampleOrderMut_short <- substr(sampleOrderMut, 9, 12)

ha_op2 = rowAnnotation(foo = anno_text(
  unique(gene.cna.mut$Annotation[1]),
  location = 2,
  rot = -90,
  just = "right",
  gp = gpar(fontsize = 20)
))
MAF_ch_sel <- MAF_ch[gene.cna.mut$Genes[1:5], ]
rownames(MAF_ch_sel)

oncoPrint_class <- function(
  dataMAF,
  classGenes,
  AnnotationPathway,
  legendON = FALSE,
  sampleNamesON = FALSE
) {
  common_classGenes <- intersect(rownames(dataMAF), classGenes)
  infoAnno <- AnnotationPathway[common_classGenes, "Annotation"][1]
  infoColor <- AnnotationPathway[common_classGenes, "Color"][1]
  datalist <- list(infoAnno = c(infoAnno =  infoColor, "B" = "orange"))
  names(datalist) <- infoAnno
  names(datalist[[1]])[1] <- infoAnno
  df <- data.frame(
    infoAnno = c(
      rep(infoAnno, length(common_classGenes)),
      rep("B", 0)
    )
  )
  colnames(df) <- infoAnno
  
  ha_row = HeatmapAnnotation(
    df = df,
    col = datalist,
    show_legend = FALSE,
    show_annotation_name = FALSE,
    which = "row",
    width = unit(1, "cm")
  )
  
  ht_op_class_curr <- oncoPrint(
    dataMAF[common_classGenes, ],
    alter_fun = alter_fun,
    height = unit(length(common_classGenes) / 3,"cm"),
    col = col,
    top_annotation = NULL,
    right_annotation = NULL,
    alter_fun_is_vectorized = FALSE,
    left_annotation = ha_row,
    column_title = "",
    heatmap_legend_param = list(
      title = "Alterations",  
      at = names(alter_fun), 
      labels = names(alter_fun) %>% {gsub("_"," ",.)} %>% {gsub("AMP","Amplification",.)}  %>% {gsub("DEL","Deletion",.)}
    ),
    remove_empty_columns = FALSE,
    remove_empty_rows = FALSE,
    show_pct = FALSE,
    border = TRUE,
    show_column_names = sampleNamesON,
    row_order = rownames(dataMAF[common_classGenes, ]),
    column_order = sampleOrderMut_short,
    show_heatmap_legend = legendON
  )
  return(ht_op_class_curr)
}

ht_op_class_cell_cycle <- oncoPrint_class(
  dataMAF = MAF_ch,
  legendON = TRUE,
  classGenes = gene.cna.mut$Hugo_Symbol[gene.cna.mut$Annotation == unique(gene.cna.mut$Annotation)[1]],
  AnnotationPathway = gene.cna.mut
)
ht_op_class_dna_repair <- oncoPrint_class(
  dataMAF = MAF_ch,
  classGenes = gene.cna.mut$Hugo_Symbol[gene.cna.mut$Annotation == unique(gene.cna.mut$Annotation)[2]],
  AnnotationPathway = gene.cna.mut
)
ht_op_class_pi3k <- oncoPrint_class(
  dataMAF = MAF_ch,
  classGenes = gene.cna.mut$Hugo_Symbol[gene.cna.mut$Annotation == unique(gene.cna.mut$Annotation)[3]],
  AnnotationPathway = gene.cna.mut
)
ht_op_class_notch <- oncoPrint_class(
  dataMAF = MAF_ch,
  classGenes = gene.cna.mut$Hugo_Symbol[gene.cna.mut$Annotation == unique(gene.cna.mut$Annotation)[4]],
  AnnotationPathway = gene.cna.mut
)
ht_op_class_mapk <-oncoPrint_class(
  dataMAF = MAF_ch,
  classGenes = gene.cna.mut$Hugo_Symbol[gene.cna.mut$Annotation == unique(gene.cna.mut$Annotation)[5]],
  AnnotationPathway = gene.cna.mut
)
ht_op_class_epigenetic <- oncoPrint_class(
  dataMAF = MAF_ch,
  classGenes = gene.cna.mut$Hugo_Symbol[gene.cna.mut$Annotation == unique(gene.cna.mut$Annotation)[6]],
  AnnotationPathway = gene.cna.mut
)

ht_op_class_immune <- oncoPrint_class(
  dataMAF = MAF_ch,
  classGenes = gene.cna.mut$Hugo_Symbol[gene.cna.mut$Annotation == unique(gene.cna.mut$Annotation)[7]],
  legendON = FALSE,
  AnnotationPathway = gene.cna.mut,
  sampleNamesON = FALSE
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Adding clinical information
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

ha1 <- HeatmapAnnotation(
  "TNBC Subtype Consensus" = S2$subtype,
  "BL1 TNBCtype Corr." = S2$BL1,
  "BL2 TNBCtype Corr." = S2$BL2,
  "M TNBCtype Corr." = S2$M,
  "LAR TNBCtype Corr." = S2$LAR,
  "BRCA Subtype PAM50" = S2$BRCA_Subtype_PAM50,
  "Age groups" = S2$agegroups,
  "Lymphnodes positive" = S2$number_of_lymphnodes_positive_by_he,
  "TIME class" = S2$TIME_CLASS,
  which = "column",
  gap = unit(1, "mm"),
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    "Lymphnodes positive" = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 10),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(3, "cm"),
      title = "Lymphnodes positive"
    ),
    "BRCA Subtype PAM50" = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 10),
      nrow = 2,
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(3, "cm"),
      title = "BRCA Subtype\nPAM50"
    ),
    "Age groups" = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 10),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(3, "cm"),
      title = "Age groups"
    ),
    "TNBC Subtype Consensus" = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 10),
      nrow = 2,
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(3, "cm"),
      title = "TNBC Subtype\nConsensus"
    ),
    "TIME class" = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 10),
      nrow = 2,
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(3, "cm"),
      title = "TIME class"
    ),
    "M TNBCtype Corr." =  list(
      direction = "horizontal",
      legend_width  = unit(3, "cm"),
      labels_gp = gpar(fontsize = 10),
      title_gp = gpar(fontsize = 10, fontface = "bold")
    ),
    "LAR TNBCtype Corr." =  list(
      direction = "horizontal",
      legend_width  = unit(3, "cm"),
      labels_gp = gpar(fontsize = 10),
      title_gp = gpar(fontsize = 10, fontface = "bold")
    ),
    "BL1 TNBCtype Corr." =  list(
      direction = "horizontal",
      legend_width  = unit(3, "cm"),
      labels_gp = gpar(fontsize = 10),
      title_gp = gpar(fontsize = 10, fontface = "bold")
    ),
    "BL2 TNBCtype Corr." =  list(
      direction = "horizontal",
      legend_width  = unit(3, "cm"),
      labels_gp = gpar(fontsize = 10),
      title_gp = gpar(fontsize = 10, fontface = "bold")
    )
  ), 
  simple_anno_size = unit(0.3,"cm"),
  col = list(
    "Age groups" = c(
      "61-100" = "#221c47",
      "41-60" = "#5f4fc4",
      "0-40" = "#59ddf7"
    ),
    "TNBC Subtype Consensus" = colors.tnbc,
    "BRCA Subtype PAM50" = colors.pam50,
    "TIME class" = c(
      "FI" = "#FCFFA4FF",
      "ID" = "#420A68FF",
      "MR" = "#AE305CFF",
      "SR" = "#F8850FFF"
    ),
    "LAR TNBCtype Corr." = colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "green")),
    "M TNBCtype Corr." = colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "blue")),
    "BL1 TNBCtype Corr." = colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "red")),
    "BL2 TNBCtype Corr." = colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "orange")),
    "Lymphnodes positive" = colorRamp2(c(0, 20), c("white", "black"))
  )
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Adding immune ESTIMATE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

TCGA_immune_ESTIMATE <- rbind(as.numeric(S2$ESTIMATE_StromalScore),as.numeric(S2$ESTIMATE_ImmuneScore))
colnames(TCGA_immune_ESTIMATE) <- S2$patient
rownames(TCGA_immune_ESTIMATE) <- c("ESTIMATE StromalScore","ESTIMATE ImmuneScore")
TCGA_immune_ESTIMATE <- TCGA_immune_ESTIMATE[,sampleOrderMut]

ht_estimate <- Heatmap(
  matrix = TCGA_immune_ESTIMATE,
  name = "ESTIMATE",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width  = unit(4, "cm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  height = unit(0.8,"cm"),
  width  = unit(30,"cm"),
  cluster_column_slices = FALSE,
  column_split = factor(
    S2$subtype,
    levels = unique(S2$subtype)
  ),
  row_names_gp =  gpar(fontsize = 10),
  cluster_columns = FALSE,
  show_column_names = FALSE,
  cluster_rows  = FALSE,
  top_annotation = ha1
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Adding MUSIC Azizi Immune Cells
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
Music_Azizi_tab_sel_Scale <- subset(
  S2,
  select = c(
    "Exhausted T cell, T-regulatory Cell",
    "Naive T cells",
    "Monocytic Lineage",
    "B cells"
  )
) %>% t 
Music_Azizi_tab_sel_Scale <- pheatmap:::scale_rows(x = Music_Azizi_tab_sel_Scale)

orderAzizi <- c(
  "Monocytic Lineage",
  "Exhausted T cell, T-regulatory Cell",
  "T−regulatory Cell"
)
orderAzizi <- rownames(Music_Azizi_tab_sel_Scale)
ht_Azizi <- Heatmap(
  matrix = Music_Azizi_tab_sel_Scale[orderAzizi, sampleOrderMut],
  name = "Azizi scRNA",
  left_annotation = HeatmapAnnotation(
    "Azizi scRNA" = rep("Azizi scRNA", length(orderAzizi)),
    which = "row",
    col = list("Azizi scRNA" = c("Azizi scRNA" = "blue")),
    show_legend = FALSE,
    show_annotation_name = FALSE
  ),
  col =  colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  use_raster = TRUE,
  raster_device = c("png"),
  raster_quality = 2,
  cluster_columns = FALSE,
  cluster_rows  = FALSE
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Adding MUSIC Karaayvaz Immune Cells
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

Karaayvaz_sel <- c(
  "Kaaraayvaz_Epithelial",
  "Kaaraayvaz_Myoepithelial",
  "Kaaraayvaz_Stroma",
  "Kaaraayvaz_Monocyte",
  "Kaaraayvaz_Lymphocyte"
)

Music_tab_sel <- subset(S2, select = Karaayvaz_sel) %>% t
Music_Karaayvaz_tab_sel_Scale <- pheatmap:::scale_rows(Music_tab_sel)

col_music = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
orderKaraayvaz <- Karaayvaz_sel
orderKaraayvaz <- rownames(Music_Karaayvaz_tab_sel_Scale)
ht_Karaayvaz <- Heatmap(
  matrix = Music_Karaayvaz_tab_sel_Scale[orderKaraayvaz, sampleOrderMut],
  name = "Karaayvaz scRNA",
  height = unit(2.0,"cm"),
  col = col_music,
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width  = unit(3, "cm"),
    labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  left_annotation = HeatmapAnnotation(
    "Karaayvaz scRNA" = rep("Karaayvaz scRNA", length(orderKaraayvaz)),
    which = "row",
    col = list("Karaayvaz scRNA" = c("Karaayvaz scRNA" = "blue")),
    show_legend = FALSE,
    show_annotation_name = FALSE
  ),
  use_raster = TRUE,
  row_names_gp =  gpar(fontsize = 10),
  raster_device = c("png"),
  raster_quality = 2,
  cluster_columns = FALSE,
  cluster_rows  = FALSE
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Adding MUSIC Nguygen Immune Cells
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
Nguygen_sel <- c(
  "Nguyen_Basal",
  "Nguyen_Basal_Myoepithelial",
  "Nguyen_Luminal_1_1",
  "Nguyen_Luminal_1_2",
  "Nguyen_Luminal_2"
)

Music_tab_sel <- subset(S2,select = Nguygen_sel) %>% t
Music_Nguygen_tab_sel_Scale <- pheatmap:::scale_rows(Music_tab_sel)
orderNguygen <- rownames(Music_Nguygen_tab_sel_Scale)[c(2,4,5)]
ht_Nguygen <- Heatmap(
  matrix = Music_Nguygen_tab_sel_Scale[orderNguygen, sampleOrderMut],
  name = "Nguygen scRNA",
  height = unit(0.8,"cm"),
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width  = unit(3, "cm"),
    labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  left_annotation = HeatmapAnnotation(
    "Nguygen scRNA" = rep("Nguygen scRNA", length(orderNguygen)),
    which = "row",
    col = list("Nguygen scRNA" = c("Nguygen scRNA" = "blue")),
    show_legend = FALSE,
    show_annotation_name = FALSE
  ),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  use_raster = TRUE,
  row_names_gp =  gpar(fontsize = 10),
  raster_device = c("png"),
  raster_quality = 2,
  cluster_columns = FALSE,
  cluster_rows  = FALSE
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Xcell                                                        |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
orderXCell <- c(
  'xCell_CD4+ memory T-cells',
  'xCell_CD4+ naive T-cells',
  'xCell_CD8+ naive T-cells',
  'xCell_CD8+ T-cells',
  'xCell_CD8+ Tcm',
  'xCell_Tregs',
  'xCell_B-cells',
  'xCell_Plasma cells',
  'xCell_Memory B-cells',
  'xCell_Monocytes',
  'xCell_DC',
  'xCell_Macrophages'
)

RNA_xCell <- subset(S2,select = orderXCell) %>% t

Xcell_tab_sel_Scale <-  pheatmap:::scale_rows(RNA_xCell)

ht_music <- Heatmap(
  matrix = Xcell_tab_sel_Scale[orderXCell, sampleOrderMut],
  name = "xCell",
  height = unit(4,"cm"),
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width  = unit(3, "cm"),
    labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  left_annotation = HeatmapAnnotation(
    "xCell" = rep("xCell", length(orderXCell)),
    which = "row",
    col = list("xCell" = c("xCell" = "#4363d8")),
    show_legend = FALSE,
    show_annotation_name = FALSE
  ),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp =  gpar(fontsize = 10),
  cluster_columns = FALSE,
  cluster_rows  = FALSE
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Gene expression                                                             |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o CPTAC   
CPTAC_GENES <- read_csv("data/CPTAC_GENES_v5.2.csv")
CPTAC_GENES$pathway[which(CPTAC_GENES$pathway == "Immune checkpoint")] <- "Immune\ncheckpoint"
CPTAC_GENES$pathway[which(CPTAC_GENES$pathway == "Ag presentation" )] <- "Ag\npresentation" 
colnames(CPTAC_GENES)[1] <- "GENE_SYMBOL"
CPTAC_GENES <- CPTAC_GENES[!duplicated(CPTAC_GENES$GENE_SYMBOL),]
CPTAC_GENES <- rbind(CPTAC_GENES,data.frame("GENE_SYMBOL" = "PDCD1", "pathway" = "Immune\ncheckpoint"))
CPTAC_GENES <- rbind(CPTAC_GENES,data.frame("GENE_SYMBOL" = "CD14", "pathway" = "Immune"))
CPTAC_GENES <- rbind(CPTAC_GENES,data.frame("GENE_SYMBOL" = "FCGR3A", "pathway" = "Immune"))

GeneImmuneList <- c(
  "CD8A",
  "CD4",
  "PDCD1",
  "CD14",
  "FCGR3A",
  "HLA-A",
  "HLA-B",
  "HLA-C",
  "TAP1",
  "TAP2",
  "B2M",
  "LAG3",
  "IDO1",
  "TIGIT",
  "HAVCR2",
  "CD274",
  "PDCD1LG2",
  "CTLA4",
  "C10orf54",
  "CD276",
  "VTCN1"
)


commonGenes_rna_iMS <- intersect(GeneImmuneList, rownames(dataFilt_TNBC))
setdiff(GeneImmuneList, rownames(dataFilt_TNBC))

dataExpr_TNBC <- dataFilt_TNBC[commonGenes_rna_iMS, ]
rownames(dataFilt_TNBC)[grep("MIC", rownames(dataFilt_TNBC))]

dataExpr_TNBC_scale <- pheatmap:::scale_rows(dataExpr_TNBC)

row.break <- CPTAC_GENES[match(GeneImmuneList,CPTAC_GENES$GENE_SYMBOL),]$pathway
row.annot <- HeatmapAnnotation(
  df = data.frame("Pathways" = row.break),
  which = "row",
  show_legend = FALSE,
  show_annotation_name = FALSE, 
  col = list("Pathways" = colors.pathways)
)

ht_GeneExpr <-  Heatmap(
  matrix = dataExpr_TNBC_scale[, sampleOrderMut],
  name = "RNA-Seq\n(row z-score)",
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_width  = unit(3, "cm"),
    labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  height = unit(10,"cm"),
  col = colorRamp2(seq(-2, 2, by = 4/99), pal_rna),
  row_split = row.break,
  left_annotation = row.annot,
  row_names_gp =  gpar(fontsize = 10),
  cluster_columns = FALSE,
  cluster_rows  = FALSE
)

#===============================================================================
# Drawiwng heatmap
#===============================================================================
pdf(
  "plots/TNBC_Figure1C_from_S2.pdf",
  width = 17,
  height = 20
)
ht_list <- 
  ht_estimate  %v% # estimate
  ht_Karaayvaz %v% # Music ht_Karaayvaz
  ht_Nguygen %v% # Music ht_Nguygen
  ht_music  %v% # xCell
  ht_GeneExpr %v% # RNA
  #-------------------# mut/CNA pathways
  ht_op_class_cell_cycle  %v%  # cell cycle
  ht_op_class_dna_repair  %v%  # DNA repair
  ht_op_class_pi3k  %v%  # PI3K
  ht_op_class_notch %v%   # NOTCH
  ht_op_class_mapk %v%   # MAPK
  ht_op_class_epigenetic  %v% 
  ht_op_class_immune # Epigenetic
draw(
  ht_list,
  column_title = "TCGA", 
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  ht_gap = unit(1, "mm"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
decorate_annotation("Cell cycle pathway", {
  grid.text(
    "Cell\ncycle",
    gp = gpar(fontsize = 12),
    x = unit(-0.2, "cm"),
    rot = 90,
    just = "bottom"
  )
})


decorate_annotation("Immune", {
  grid.text(
    "Immune",
    gp = gpar(fontsize = 12),
    x = unit(-0.2, "cm"),
    rot = 90,
    just = "bottom"
  )
})
decorate_annotation("DNA repair pathway", {
  grid.text(
    "DNA\nrepair",
    gp = gpar(fontsize = 12),
    x = unit(-0.2, "cm"),
    rot = 90,
    just = "bottom"
  )
})
decorate_annotation("Growth factor/PI3K pathway", {
  grid.text("PI3K",
            gp = gpar(fontsize = 12),
            
            x = unit(-0.2, "cm"),
            rot = 90,
            just = "bottom")
})
decorate_annotation("Notch pathway", {
  grid.text("NOTCH",
            gp = gpar(fontsize = 12),
            
            x = unit(-0.2, "cm"),
            rot = 90,
            just = "bottom")
})
decorate_annotation("MAPK pathway", {
  grid.text("MAPK",
            gp = gpar(fontsize = 12),
            
            x = unit(-0.2, "cm"),
            rot = 90,
            just = "bottom")
})
decorate_annotation("PRC2/BAF pathway", {
  grid.text(
    "Epigenetic",
    gp = gpar(fontsize = 12),
    
    x = unit(-0.2, "cm"),
    rot = 90,
    just = "bottom"
  )
})

decorate_annotation("xCell", {
  grid.text(
    "xCell\nsignature",
    gp = gpar(fontsize = 12),
    x = unit(-0.2, "cm"),
    rot = 90,
    just = "bottom"
  )
})

decorate_annotation("Karaayvaz scRNA", {
  grid.text(
    "scRNA\ndeconvolution",
    gp = gpar(fontsize = 12),
    x = unit(-0.2, "cm"),
    y = unit(0.7, "cm"),
    rot = 90,
    just = "bottom"
  )
})


dev.off()
