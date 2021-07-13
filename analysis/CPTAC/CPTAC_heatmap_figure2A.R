#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# This script creates a heatmap to show
# RNA-seq, protein and phospo data for CPTAC samples
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
library(readr)
library(ComplexHeatmap)
library(rlang)
library(circlize)
library(dplyr)
library(TCGAbiolinks)
library(DelayedMatrixStats)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Folders where we was saved the data
# See script: 01_data_download.R
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots/CPTAC"
dir.tcga <- file.path(dir.data,"TCGA")
dir.cptac <- file.path(dir.data,"CPTAC")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o colors                                                        |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# 3) For ATAC Heatmap
pal_atac <- colorRampPalette(
  c('#3361A5', '#248AF3', '#14B3FF', 
    '#88CEEF', '#C1D5DC', '#EAD397', 
    '#FDB31A','#E42A2A', '#A31D1D')
)(100)

pal_rna <- colorRampPalette(
  c("#352A86","#343DAE","#0262E0","#1389D2",
    "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D")
)(100)

pal_phospo <- colorRampPalette(c("blue","white","#82118f"))(100)

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
  "EGFR" = "#429abd",
  "EMT/wound repair" = "#5ca8ee",
  "AR" = "#1eff93",
  "AR signaling" = "#1eff93",
  "Immune" = "#0AF9F6",
  "Immune checkpoint" ="blue",
  "Immune\ncheckpoint" ="blue",
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
  "adhesion/motility" ='blue',
  "TF" = "purple",
  "PI3K" = "pink",
  "xCell sig." = "#4363d8",
  "PRC2\nAUXILLARY" = "#756767",
  "PRC2\nCORE" = "#756767",
  "SWI-SNF" = "green",
  "UTX" = "black"
)

color.bl1.cor <- colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "red"))
color.bl2.cor <- colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "orange"))
color.m.cor <- colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "blue"))
color.lar.cor <- colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "green"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Data                                                                       |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o CPTAC   
CPTAC_GENES <- readxl::read_xlsx(file.path(dir.cptac,"CPTAC_protein_SIG_LFC1.xlsx"))
colnames(CPTAC_GENES)[1] <- "GENE_SYMBOL"

load(file.path(dir.cptac,"cptac_TNBC_linkedomics_27samples.Rdata"))
cptac_meta <-  TNBC_samples_27

# Get TNBC_subtype data
cptac_meta <- cptac_meta %>% group_by(TNBC_subtype) %>% plyr::arrange(factor(TNBC_subtype,levels = c("BL1","BL2","M","LAR")))
cptac_meta$Sample.ID <- cptac_meta$patient

#cptac_meta <- plyr::adply(c("BL1","BL2","M","LAR"), 1 , function(TNBC_subtype) {
#  cptac_meta %>% 
#    dplyr::filter(cptac_meta$TNBC_subtype == !!TNBC_subtype)  %>% 
#    dplyr::arrange(desc(!!(parse_expr(TNBC_subtype))))
#})

# Set samples in the same order
protein <- cptac_TNBC_Prot[,match(cptac_meta$Sample.ID,colnames(cptac_TNBC_Prot))]
rna <- cptac_TNBC_RNA[,match(cptac_meta$Sample.ID,colnames(cptac_TNBC_RNA))]
phospo <- cptac_TNBC_Phospho_median[,match(cptac_meta$Sample.ID,colnames(cptac_TNBC_Phospho_median))]

# Checking 
all(colnames(phospo) == colnames(rna))
all(colnames(phospo) == colnames(protein))

# Set genes in the same order
protein <- protein[rownames(protein) %in% CPTAC_GENES$GENE_SYMBOL,]
protein <- protein[match(CPTAC_GENES$GENE_SYMBOL,rownames(protein)),]

rna <- rna[rownames(rna) %in% CPTAC_GENES$GENE_SYMBOL,]
rna <- rna[match(CPTAC_GENES$GENE_SYMBOL,rownames(rna)),]

phospo <- phospo[rownames(phospo) %in% CPTAC_GENES$GENE_SYMBOL,]
phospo <- phospo[match(CPTAC_GENES$GENE_SYMBOL,rownames(phospo)),]

# data with no NA
common.genes <- CPTAC_GENES$GENE_SYMBOL %>% 
  intersect(rownames(protein)) %>%
  intersect(rownames(rna)) %>%
  intersect(rownames(phospo))

phospo.no.na <- phospo[common.genes,]
rna.no.na <- rna[common.genes,]
protein.no.na <- protein[common.genes,]

all(rownames(protein.no.na) == rownames(rna.no.na))
all(rownames(protein.no.na) == rownames(phospo.no.na))

rna.and.xcell <- rna.no.na
protein.and.xcell <- protein.no.na
phospo.and.xcell <- phospo.no.na

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Heatmap plot                                                  |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
library(ComplexHeatmap)
ha <- columnAnnotation(
  "TNBC subtype Consensus" = cptac_meta$TNBC_subtype[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  "BL1 TNBCtype Corr." = cptac_meta$BL1[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  "BL2 TNBCtype Corr." = cptac_meta$BL2[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  "M TNBCtype Corr." = cptac_meta$M[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  "LAR TNBCtype Corr." = cptac_meta$LAR[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  "PAM50" = cptac_meta$PAM50[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  col = list( 
    "TNBC subtype Consensus"= colors.tnbc,
    "BL1 TNBCtype Corr." = color.bl1.cor,
    "BL2 TNBCtype Corr." = color.bl2.cor,
    "LAR TNBCtype Corr." = color.lar.cor,
    "M TNBCtype Corr." = color.m.cor,
    "PAM50" = colors.pam50
  ),
  simple_anno_size = unit(0.3, "cm"),
  annotation_legend_param = list(
    "BL1 TNBCtype Corr." = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(2, "cm")),
    "LAR TNBCtype Corr." = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(2, "cm")),
    "BL2 TNBCtype Corr." = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(2, "cm")),
    "M TNBCtype Corr." = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(2, "cm")),
    "TNBC subtype Consensus" = list(legend_gp = gpar(fontsize = 10),nrow = 2)
  ),
  show_annotation_name = TRUE,
  show_legend = FALSE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8)
)

ha.no.names <- columnAnnotation(
  "TNBC subtype Consensus" = cptac_meta$TNBC_subtype[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  "BL1 TNBCtype Corr." = cptac_meta$BL1[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  "BL2 TNBCtype Corr." = cptac_meta$BL2[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  "M TNBCtype Corr." = cptac_meta$M[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  "LAR TNBCtype Corr." = cptac_meta$LAR[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  "PAM50" = cptac_meta$PAM50[match(colnames(protein.and.xcell),cptac_meta$Sample.ID)],
  col = list(
    "TNBC subtype Consensus"= colors.tnbc,
    "BL1 TNBCtype Corr." = color.bl1.cor,
    "BL2 TNBCtype Corr." = color.bl2.cor,
    "LAR TNBCtype Corr." = color.lar.cor,
    "M TNBCtype Corr." = color.m.cor,
    "PAM50" = colors.pam50
  ),
  simple_anno_size = unit(0.3, "cm"),
  annotation_legend_param = list(
    "BL1 TNBCtype Corr." = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      legend_width  = unit(2, "cm")),
    "LAR TNBCtype Corr." = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      legend_width  = unit(2, "cm")),
    "BL2 TNBCtype Corr." = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      legend_width  = unit(2, "cm")),
    "M TNBCtype Corr." = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      legend_width  = unit(2, "cm")),
    "TNBC subtype Consensus" = list(legend_gp = gpar(fontsize = 8),nrow = 2)
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 8)
)


col.rna <- colorRamp2(seq(-2, 2, by = 4/99), pal_rna)
col.protein <- colorRamp2(seq(-2, 2, by = 4/99), pal_atac)
col.phospo <- colorRamp2(seq(-2, 2, by = 4/99), pal_phospo)

CPTAC_GENES <- CPTAC_GENES %>% filter(GENE_SYMBOL %in% rownames(protein.and.xcell))
CPTAC_GENES <- CPTAC_GENES[!duplicated(CPTAC_GENES$GENE_SYMBOL),]
CPTAC_GENES <- CPTAC_GENES[match(rownames(protein.and.xcell),CPTAC_GENES$GENE_SYMBOL),]
#----------------------------------------------
# Protein
#----------------------------------------------
row.break <- c(CPTAC_GENES$Subtype)
row.annot <- HeatmapAnnotation(
  df = data.frame("Subtype" = row.break),
  which = "row",
  show_legend = FALSE,
  show_annotation_name = FALSE,
  col = list("Subtype" = colors.tnbc)
)

ht.protein <- Heatmap(
  protein.and.xcell %>% data.matrix() %>% t() %>% scale() %>% t(),
  show_column_names = FALSE,
  show_row_names = FALSE,
  cluster_rows = FALSE,
  col = col.protein, 
  column_title = "Protein",
  cluster_columns = FALSE, 
  row_title_gp =  gpar(fontface = "bold",fontsize = 10),
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  raster_device = c("png"),
  cluster_row_slices = FALSE,
  #column_order = order,
  column_split = factor(cptac_meta$TNBC_subtype,levels = unique(cptac_meta$TNBC_subtype)),
  use_raster = TRUE,
  width = unit(3, "cm"),
  raster_quality = 2,
  top_annotation = ha.no.names,
  row_names_gp =  gpar(fontsize = 8),
  heatmap_legend_param =  list(
    direction = "horizontal",
    legend_width  = unit(3, "cm"),
    labels_gp = gpar(fontsize = 8),
    title_gp = gpar(fontface = "bold",fontsize = 10)
  ),
  name = "Protein\n(row z-score)"
)


#----------------------------------------------
# RNA
#----------------------------------------------

ht.rna <- Heatmap(
  rna.and.xcell %>% data.matrix() %>% t() %>% scale() %>% t(),
  show_column_names = FALSE,
  cluster_rows = FALSE,
  col = col.rna, 
  use_raster = TRUE,
  column_title = "RNA-Seq",
  left_annotation = row.annot,
  row_split = factor(row.break,levels = unique(row.break)),
  cluster_columns = FALSE, 
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  raster_device = c("png"),
  column_split = factor(cptac_meta$TNBC_subtype,levels = unique(cptac_meta$TNBC_subtype)),
  cluster_row_slices = FALSE,
  #column_order = order,
  width = unit(3, "cm"),
  raster_quality = 2,
  top_annotation = ha,
  show_row_names = FALSE,
  row_names_gp =  gpar(fontsize = 8),
  heatmap_legend_param =  list(
    direction = "horizontal",
    legend_width  = unit(3, "cm"),
    labels_gp = gpar(fontsize = 8),
    title_gp = gpar(fontface = "bold",fontsize = 10)),
  name = "RNA-Seq\n(row z-score)"
)

pdf(
  file.path(dir.plots,"CPTAC_TNBC_GENES_prospective_2A.pdf"),
  width = 8, 
  height = 7
)
draw(
  ht.rna + ht.protein,
  newpage = TRUE, 
  column_title = "CPTAC", 
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()