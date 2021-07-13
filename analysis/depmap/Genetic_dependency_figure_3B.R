#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# The following script downloads all data used to perform TNBC analysis
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
library(readr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(gplots)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/Genetic_dependancy//"
dir.plots <- "plots/Genetic_dependancy"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Colors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colors.tnbc <- c(
  "BL1" = "red",
  "BL2" = "orange",
  "LAR" = "green",
  "M" = "blue"
)

colors.pathways <- c(
  "cell cycle" = "red",
  "EMT" = "#5ca8ee",
  "AR" = "#56a628",
  "Immune" = "#0AF9F6",
  "xCell sig." = "#4363d8",
  "RTK/PI3K/mTOR" = "blue",
  "chromatin" = "#2ae6f7",
  "adhesion/motility growth factor" = "red",
  "development" = "grey",
  "notch" = "brown",
  "other" = "black",
  "growth factor" = "#34b3b0",
  "MAPK" = "orange",
  "DNA repair" = "yellow",
  "adhesion/motility" = 'pink',
  "TF" = "#ceff94"
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- read.csv(
  file.path(dir.data,"RNAi_sgRNA_depmap_merged_select_seq_break.csv"),
  sep = ",",
  header = TRUE,
  row.names = 1
)

matrix <- as.matrix(data[, 1:8])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Heatmap
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ha <- columnAnnotation(
  "TNBC Subtype" = gsub("_[[:alnum:]]*","",colnames(matrix)),
  "Knockdown" = gsub("[[:alnum:]]*_","",colnames(matrix)),
  col = list(
    "TNBC Subtype"  = colors.tnbc,
    "Knockdown" = c("RNAi" = "lightblue", "sgRNA" = "yellow")
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  height = unit(0.3,"cm"),
  annotation_height = unit(0.3,"cm"),
  simple_anno_size = unit(0.3, "cm")
)


row.break <- factor(
  data$pathway %>% as.character(),
  levels = unique(data$pathway) %>% as.character()
)
ha.row <- HeatmapAnnotation(
  "Pathway" = row.break,
  which = "row",
  width = unit(0.2,"cm"),
  annotation_legend_param = list(ncol = 1),
  annotation_width = unit(0.2,"cm"),
  simple_anno_size = unit(0.3, "cm"),
  show_legend = TRUE,
  show_annotation_name = FALSE,
  col = list("Pathway" = colors.pathways)
)

ht <- Heatmap(
  matrix,  
  top_annotation = ha,
  left_annotation = ha.row,
  raster_device = c("png"),
  row_split = row.break,
  raster_quality = 2,
  #row_title = " ",
  height = unit(25,"cm"),
  width = unit(3,"cm"),
  row_names_side = "right",
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  name = "Effect",
  show_column_names = FALSE,
  cluster_row_slices = FALSE,
  row_names_gp =  gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6),
)

pdf(
  file.path(dir.plots,"figure3C_genetic_dependancy.pdf"),
  width = 6,
  height = 20
)
draw(
  ht, 
  heatmap_legend_side = "right",
  merge_legend = TRUE,
  annotation_legend_side = "right",
  newpage = FALSE
)
dev.off()