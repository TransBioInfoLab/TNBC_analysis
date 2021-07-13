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
# Date: 24 June 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ComplexHeatmap)
library(dplyr)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Paths
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.data <- "data/depmap/"
dir.plots <- "plots/depmap"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# heatmap for all sig RNAi
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- readr::read_csv(file.path(dir.data,"depmap_sig_TNBCtype.csv")) %>% 
  data.frame(row.names = .$Gene) %>% dplyr::mutate(Gene = NULL)

ha <- columnAnnotation(
  "Knockdown" = rep("RNAi",4),
  "TNBC subtype" =  gsub("_tvalue","",colnames(data)),
  col = list(
    "TNBC subtype" = c(
      "BL1" = "red",
      "BL2" = "orange",
      "LAR" = "green",
      "M" = "blue",
      "BL3" = "black",
      "UNC" = "grey"
    ),
    "Knockdown" = c("RNAi" = "lightblue", "sgRNA" = "yellow")
  ),
  annotation_legend_param = list("TNBC subtype" = list(legend_gp = gpar(fontsize = 10),nrow = 4)),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10)
)


ht.RNAi <- Heatmap(
  data,
  top_annotation = ha,
  cluster_columns = FALSE,
  column_title = "depmap\nRNAi",
  show_row_names = FALSE,
  width = unit(2,"cm"),
  height = unit(5,"cm"),
  name = "T-value",
  row_title = paste0(nrow(data), " genes"),
  show_column_names = FALSE,
  show_row_dend = FALSE,
  raster_device = c("png"),
  raster_quality = 2
) 

#######heatmap for all sig sgRNA
data <- readr::read_csv(file.path(dir.data,"depmap_sig_sgRNA_TNBCtype.csv"))  %>% 
  data.frame(row.names = .$Gene) %>% dplyr::mutate(Gene = NULL)


ha <- columnAnnotation(
  "Knockdown" = rep("sgRNA",4),
  "TNBC subtype" =  gsub("_tvalue","",colnames(data)),
  col = list(
    "TNBC subtype" = c(
      "BL1" = "red",
      "BL2" = "orange",
      "LAR" = "green",
      "M" = "blue",
      "BL3" = "black",
      "UNC" = "grey"
    ),
    "Knockdown" = c("RNAi" = "lightblue", "sgRNA" = "yellow")
  ),
  annotation_legend_param = list("TNBC subtype" = list(legend_gp = gpar(fontsize = 10),nrow = 4)),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10))

ht.crispr <- Heatmap(
  data,
  top_annotation = ha,
  cluster_columns = FALSE,
  column_title = "sgRNA",
  show_row_names = FALSE,
  width = unit(2,"cm"),
  height = unit(5,"cm"),
  row_title_side = "right",
  row_title = paste0(nrow(data), " genes"),
  name = "T-value",
  show_column_names = FALSE,
  show_row_dend = FALSE,
  raster_device = c("png"),
  raster_quality = 2
) 
pdf(file.path(dir.plots,"depmap_knockdown_figure_5SC.pdf"),height = 4,width = 3)
draw(ht.crispr, merge_legend = TRUE)
draw(ht.RNAi, merge_legend = TRUE)
dev.off()