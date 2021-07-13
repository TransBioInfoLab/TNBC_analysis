#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# This script creates figure 3E
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
# Date: 23 June 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Libraries                                                    |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Directories                                                  |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
dir.data <- "data/"
dir.plot <- "plots/"

# create folder if it does not exist
for(p in grep("dir",ls(),value = T)) {
  dir.create(get(p),recursive = TRUE,showWarnings = FALSE)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Data                                                         |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
type <- readxl::read_xlsx(file.path(dir.data,"summary_heatmap_all_data_BDL.xlsx"),n_max = 1,col_names = F)
tab <- readxl::read_xlsx(file.path(dir.data,"summary_heatmap_all_data_BDL.xlsx"),skip = 1) %>% as.data.frame()

m <- tab[,-c(1:2)]
rownames(m) <- tab$Gene

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Heatmap                                                       |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

colors.tnbc <- c(
  "BL1" = "red",
  "BL2" = "orange",
  "LAR" = "green",
  "M" = "blue"
)
ra <- HeatmapAnnotation(
  "Subtype" = tab$Subtype,
  which = "row",
  col = list("Subtype" = colors.tnbc),
  annotation_legend_param = list(
    "Subtype" = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 10),
      nrow = 1,
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(3, "cm"),
      title = "TNBC Subtype"
    )
  ),annotation_name_side = "top",
  show_annotation_name = TRUE
)

ta <- HeatmapAnnotation(
  "Data type" = type[,-c(1:2)] %>% unlist %>% {gsub("Genomics","Genomic alteration",.)},
  which = "column",
  show_annotation_name = FALSE,
  col = list("Data type" = c(
    "Genomic alteration" = "yellow",
    "Genetic dependancy" = "#db34eb",
    "Pharmacologic dependency" = "#6d98e3",
    "PDX Dependency" = "darkblue"        
  )
  ),
  annotation_legend_param = list(
    "Data type" = list(
      direction = "horizontal",
      labels_gp = gpar(fontsize = 10),
      nrow = 1,
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      legend_width  = unit(3, "cm"),
      title = "Data type"
    )
  )
)


ht <- Heatmap(
  matrix = m,
  col = list("0" = "black", "1" = "red", "NA" = "grey"),
  name = "Legend",
  row_split = tab$Subtype,
  row_title = " ",
  show_row_names = TRUE,
  column_split = type[,-c(1:2)] %>% unlist %>% factor(levels = c("Genomics","Genetic dependancy","Pharmacologic dependency", "PDX Dependency")),
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  rect_gp = gpar(col= "white"),
  column_names_side = "top",
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_names_side = "left",
  left_annotation = ra,
  top_annotation = ta,
  heatmap_legend_param = list(
    direction = "horizontal",
    labels_gp = gpar(fontsize = 10),
    nrow = 1,
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    legend_width  = unit(3, "cm")
  )
)

pdf(file = file.path(dir.plot,"summary_heatmap_all_data_BDL.pdf"),width = 14,height = 10)
draw(
  ht, 
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom", 
  merge_legends = TRUE
)
dev.off()