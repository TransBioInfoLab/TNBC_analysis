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
# Date: 23 June 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(readr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Path
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.data <- "data"
dir.plots <- "plots/MET500_metabric"
dir.met500 <- file.path(dir.data,"MET500")
dir.metabric <- file.path(dir.data,"metabric")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Colors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pal_rna <- colorRampPalette(
  c(
    "#352A86","#343DAE","#0262E0","#1389D2",
    "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"
  )
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

colors.biopsy.tissue <- c(
  "Bone Marrow" = "red",
  "Brain" = "#1d964b", # green
  "Liver"= "grey",
  "Lung"= "blue",
  "Lymph Node"= "purple",
  "Skin"= "#b8f000",
  "Soft Tissue" = "orange",
  "Primary breast" = "#ff8fad"
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
  "Immune checkpoint" ="blue",
  "Immune\ncheckpoint" ="blue",
  "xCell sig." = "#4363d8",
  "RTK/PI3K/mTOR" = "pink",
  "chromatin" = "gold",
  "adhesion/motility growth factor" = "red",
  "Growth factor" = "red",
  "growth factor" = "#34b3b0",
  "GF" = "red",
  "development" = "grey",
  "Development" = "grey",
  "notch" = "brown",
  "other" = "yellow",
  "MAPK" = "black",
  "DNA repair" = "orange",
  "adhesion/motility" ='blue',
  "TF" = "purple",
  "PI3K" = "pink",
  "xCell sig." = "#4363d8"
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# xcell heatmap
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
orderXCell <-
  c(
    'CD4+ memory T-cells',
    'CD4+ naive T-cells',
    'CD8+ naive T-cells',
    'CD8+ T-cells',
    'CD8+ Tcm',
    'Tregs',
    'B-cells',
    'Plasma cells',
    'Memory B-cells',
    'Monocytes',
    'DC',
    'Macrophages'
  )

orderHeatmap <- c(
  "xCell sig.",
  "Immune",
  "Ag presentation",
  "Immune checkpoint",
  "Cell cycle",
  "DNA repair",
  "MAPK",
  "Development",
  "EMT/wound repair",
  "TF",
  "Growth factor",
  "AR signaling",
  "PI3K"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(dir.met500,"MET500_TNBC.rda"))

# load(file.path(dir.metabric,"Metabric_TNBC_data_withxCell.Rdata"))
load( file.path(dir.metabric,"metabric_TNBC.rda"))

all(colnames(MET500.tnbc) == colnames(MET500.tnbc.xcell))
all(colnames(MET500.tnbc) == MET500.tnbc.meta$sample_detailed)

MET500.tnbc.with.xcell <- rbind(MET500.tnbc,MET500.tnbc.xcell)
MET500.tnbc.meta$ImmuneScore <- MET500.tnbc.with.xcell["ImmuneScore",] %>% as.numeric()

Metabric.tnbc.meta <- Metabric.tnbc.meta %>% as.data.frame
rownames(Metabric.tnbc.meta) <- Metabric.tnbc.meta$sample
all(colnames(Metabric.tnbc) == colnames(Metabric.tnbc.xcell))
Metabric.tnbc.meta <- Metabric.tnbc.meta[colnames(Metabric.tnbc),,drop = FALSE]
all(colnames(Metabric.tnbc) == rownames(Metabric.tnbc.meta))

Metabric.tnbc.meta_with_xcell <- rbind(Metabric.tnbc, Metabric.tnbc.xcell)
Metabric.tnbc.meta$ImmuneScore  <- Metabric.tnbc.meta_with_xcell["ImmuneScore",] %>% as.numeric()

# what to plot
CPTAC_GENES <- read_csv(file.path(dir.data,"CPTAC_GENES_v5.2.csv"))


colnames(CPTAC_GENES)[1] <- "GENE_SYMBOL"
CPTAC_GENES <- CPTAC_GENES[!duplicated(CPTAC_GENES$GENE_SYMBOL),]

CPTAC_GENES <- rbind(CPTAC_GENES,data.frame(GENE_SYMBOL = orderXCell,pathway = "xCell sig."))

common.genes <- CPTAC_GENES$GENE_SYMBOL %>% 
  intersect(rownames(MET500.tnbc.with.xcell)) %>% 
  intersect(rownames(Metabric.tnbc.meta_with_xcell))
CPTAC_GENES <- CPTAC_GENES %>% dplyr::filter(GENE_SYMBOL %in% common.genes)

CPTAC_GENES <- CPTAC_GENES %>% dplyr::arrange(factor(pathway,levels = orderHeatmap))

CPTAC_GENES$pathway[which(CPTAC_GENES$pathway == "Immune checkpoint")] <- "Immune\ncheckpoint"
CPTAC_GENES$pathway[which(CPTAC_GENES$pathway == "Ag presentation" )] <- "Ag\npresentation" 
CPTAC_GENES$pathway[which(CPTAC_GENES$pathway == "Growth factor" )] <- "GF" 
CPTAC_GENES$pathway[which(CPTAC_GENES$pathway == "AR signaling" )] <- "AR" 


# Improving data
MET500.tnbc.meta$biopsy_tissue <- gsub("_"," ",MET500.tnbc.meta$biopsy_tissue) %>% stringr::str_to_title()
MET500.tnbc.meta <- MET500.tnbc.meta %>% dplyr::filter(subtype != "UNC")
MET500.tnbc.with.xcell <- MET500.tnbc.with.xcell[,colnames(MET500.tnbc.with.xcell) %in% MET500.tnbc.meta$sample_detailed]
MET500.tnbc.with.xcell <-  MET500.tnbc.with.xcell[,colnames(MET500.tnbc.with.xcell) %in% MET500.tnbc.meta$sample_detailed]
MET500.tnbc.with.xcell <- MET500.tnbc.with.xcell[,match(MET500.tnbc.meta$sample_detailed, colnames(MET500.tnbc.with.xcell))]
MET500.tnbc.with.xcell <- MET500.tnbc.with.xcell[,match(MET500.tnbc.meta$sample_detailed, colnames(MET500.tnbc.with.xcell))]
all(colnames(MET500.tnbc.with.xcell) == MET500.tnbc.meta$sample_detailed)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MET500
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MET500.tnbc.with.xcell <- MET500.tnbc.with.xcell[common.genes %in% CPTAC_GENES$GENE_SYMBOL,]
MET500.tnbc.with.xcell <- MET500.tnbc.with.xcell[na.omit(match(common.genes,rownames(MET500.tnbc.with.xcell))),]

all(colnames(MET500.tnbc.with.xcell) == colnames(MET500.tnbc.with.xcell))

# Heatmap plot
library(ComplexHeatmap)
ha <- columnAnnotation(
  "TNBC subtype" = MET500.tnbc.meta$subtype[match(colnames(MET500.tnbc.with.xcell),MET500.tnbc.meta$sample_detailed)],
  "Biopsy tissue" = MET500.tnbc.meta$biopsy_tissue[match(colnames(MET500.tnbc.with.xcell),MET500.tnbc.meta$sample_detailed)],
  "Immune Score" = MET500.tnbc.meta$ImmuneScore[match(colnames(MET500.tnbc.with.xcell),MET500.tnbc.meta$sample_detailed)],
  col = list(
    "TNBC subtype" = colors.tnbc,
    "Biopsy tissue" = colors.biopsy.tissue,
    "Immune Score" = colorRamp2(breaks = c(0,0.6),colors = c("white","black"))
  ),
  annotation_legend_param = list(
    "TNBC subtype" = list(legend_gp = gpar(fontsize = 10),ncol = 1),
    "Biopsy tissue" = list(legend_gp = gpar(fontsize = 10),ncol = 1),
    "Immune Score" = list(title = "Immune\nScore")
  ),
  show_annotation_name = TRUE,
  show_legend = TRUE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 10)
)

order <- unlist(plyr::alply(c("BL1","BL2","M","LAR"), 1 , function(x) {
  idx <- which(MET500.tnbc.meta$subtype[match(colnames(MET500.tnbc.with.xcell),MET500.tnbc.meta$sample_detailed)] == x)
  if(length(idx) == 1) return(idx)
  aux <- na.omit(MET500.tnbc.with.xcell[,idx])
  order <- t(aux) %>% dist %>% hclust(method = "average")
  as.numeric(idx[order$order])
}))


row.annot <- HeatmapAnnotation(
  df = data.frame("Pathway" = CPTAC_GENES$pathway),
  which = "row",
  show_legend = FALSE,
  show_annotation_name = FALSE,
  col = list("Pathway" = colors.pathways)
)

col <- colorRamp2(seq(-2, 2, by = 4/99), pal_rna)

ht.met500 <- Heatmap(
  MET500.tnbc.with.xcell %>% data.matrix()  %>% t() %>% scale() %>% t(),
  show_column_names = FALSE,
  cluster_rows = FALSE,
  show_heatmap_legend = FALSE,
  col = col, 
  #left_annotation = row.annot,
  cluster_columns = FALSE, 
  column_split = factor( 
    MET500.tnbc.meta$subtype[match(colnames(MET500.tnbc.with.xcell),MET500.tnbc.meta$sample_detailed)],
    levels = c("BL1","BL2","M","LAR")),
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  raster_device = c("png"),
  cluster_row_slices = FALSE,
  column_title = paste0("MET500 (n = ", ncol(MET500.tnbc.with.xcell),")"),
  column_order = order,
  row_split =  factor(CPTAC_GENES$pathway, levels = unique(CPTAC_GENES$pathway)),
  width = unit(5, "cm"),
  raster_quality = 2,
  top_annotation = ha,
  row_names_gp =  gpar(fontsize = 10),
  heatmap_legend_param =  list(
    direction = "vertical",
    title = "RNA\n(row z-score)",
    title_gp = gpar(fontface = "bold",fontsize = 10)
  ),
  name = "RNA (row z-score)"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Metabric
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Metabric.tnbc.meta_with_xcell <- Metabric.tnbc.meta_with_xcell[common.genes,]
Metabric.tnbc.meta_with_xcell <- Metabric.tnbc.meta_with_xcell[na.omit(match(common.genes,rownames(Metabric.tnbc.meta_with_xcell))),]

Metabric.tnbc.meta$TNBCtype <- gsub("[A-Z]-","",Metabric.tnbc.meta$TNBCtype)

order <- unlist(plyr::alply(c("BL1","BL2","M","LAR"), 1 , function(x) {
  idx <- which(Metabric.tnbc.meta$TNBCtype[match(colnames(Metabric.tnbc.meta_with_xcell),rownames(Metabric.tnbc.meta))] == x)
  if(length(idx) == 1) return(idx)
  aux <- na.omit(Metabric.tnbc.meta_with_xcell[,idx])
  order <- t(aux) %>% dist %>% hclust(method = "average")
  as.numeric(idx[order$order])
}))

library(ComplexHeatmap)
ha <- columnAnnotation(
  "TNBC subtype" = Metabric.tnbc.meta$TNBCtype[match(colnames(Metabric.tnbc.meta_with_xcell),rownames(Metabric.tnbc.meta))],
  "Biopsy tissue" = rep("Primary breast",ncol(Metabric.tnbc.meta_with_xcell)),
  "Immune Score" = Metabric.tnbc.meta$ImmuneScore[match(colnames(Metabric.tnbc.meta_with_xcell),rownames(Metabric.tnbc.meta))],
  annotation_legend_param = list(
    "TNBC subtype" = list(legend_gp = gpar(fontsize = 10),ncol = 1),
    "Biopsy tissue" = list(legend_gp = gpar(fontsize = 10),ncol = 1),
    "Immune Score" = list(title = "Immune\nScore")
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_name_side = "left",
  col = list(
    "TNBC subtype" = colors.tnbc,
    "Biopsy tissue" = colors.biopsy.tissue,
    "Immune Score" = colorRamp2(breaks = c(0,0.6),colors = c("white","black"))
  ),
  annotation_name_gp = gpar(fontsize = 10))


col <- colorRamp2(seq(-2, 2, by = 4/99), pal_rna)

ht.metabric <- Heatmap(
  Metabric.tnbc.meta_with_xcell %>% data.matrix()  %>% t() %>% scale() %>% t(),
  show_column_names = FALSE,
  cluster_rows = FALSE,
  col = col, 
  column_split = factor(Metabric.tnbc.meta$TNBCtype[match(colnames(Metabric.tnbc.meta_with_xcell),rownames(Metabric.tnbc.meta))],levels = c("BL1","BL2","M","LAR")),
  left_annotation = row.annot,
  cluster_columns = FALSE, 
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  raster_device = c("png"),
  show_heatmap_legend = FALSE,
  cluster_row_slices = FALSE,
  column_order = order,
  column_title = paste0("Metabric (n = ", ncol(Metabric.tnbc.meta_with_xcell),")"),
  row_split =  factor(CPTAC_GENES$pathway, levels = unique(CPTAC_GENES$pathway)),
  width = unit(10, "cm"),
  raster_quality = 2,
  top_annotation = ha,
  row_names_gp =  gpar(fontsize = 6),
  heatmap_legend_param =  list(
    direction = "vertical",
    title_gp = gpar(fontface = "bold",fontsize = 10)
  ),
  name = "RNA (row z-score)"
)


pdf(file.path(dir.plots,"TNBC_MET500_metabric_RNA_zscore_only_immune.pdf"), width = 10, height = 20)
draw(
  ht.metabric + ht.met500,
  newpage = TRUE, 
  #column_title = "MET500 & metabric RNA", 
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

