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
# Date: 21 June 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(ComplexHeatmap)
library(circlize)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/depmap"
dir.plots <- "plots/depmap/"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)


#============================================================
# Data
#============================================================
data <- readr::read_csv(file.path(dir.data,"TNBC_cell_line_consensus.csv")) %>% as.data.frame() 
data[which(data == "UNS")] <- "UNC"
rownames(data) <- data$Cell_line
data$Cell_line <- NULL


#=====================================
# Consensus heatmap
#=====================================

col <-   list(
  "BL1" = "red",
  "BL2" = "orange",
  "LAR" = "green",
  "M" = "blue",
  "UNC" = "black",
  "MSL" = "purple",
  "IM" = "pink",
  "UNS" = "black"
  
)


ha.source <- HeatmapAnnotation(
  "Cell line source" = colnames(data[,-ncol(data)]),
  show_legend = TRUE,
  which = "column",
  annotation_name_side = "left",
  col  = list(
    "Cell line source" = c(
      "GDS"  = "gold",
      "GRAY" = "gray",
      "CCLE" = "purple",
      "VANDERBILT" = "black",
      "GSE10890" = "pink"
    )
  )
)


ht_list <-
  Heatmap(t(data[,-ncol(data)]), 
          name =  "TNBC subtype", 
          col = col,
          column_names_gp = gpar(fontsize = 12),
          heatmap_legend_param = list(legend_direction = "horizontal",
                                      labels_gp = gpar(fontsize = 12), 
                                      title_gp = gpar(fontsize = 12)),
          show_column_names = TRUE,
          column_names_side = "top",
          show_row_names = TRUE,
          row_names_side = "left",
          use_raster = TRUE,
          rect_gp = gpar(col = "white", lwd = 2),
          clustering_method_rows = "average",
          raster_device = c("png"),
          raster_quality = 2,
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          show_row_dend = FALSE,
          #row_title = paste0(nrow(data), " cell lines"),
          row_names_gp = gpar(fontsize = 12),
          #top_annotation = ha,
          width = unit(14, "cm"),
          #column_title =  paste0("Databases"), 
          column_title_gp = gpar(fontsize = 12), 
          row_title_gp = gpar(fontsize = 12)) %v%
  columnAnnotation(
    "Consensus" = data$CONSENSUS,
    show_legend = FALSE,
    annotation_name_side = "left",
    col  = list(
      "Consensus" = c(
        "BL1" = "red",
        "BL2" = "orange",
        "LAR" = "green",
        "M" = "blue",
        "UNC" = "black",
        "MSL" = "purple",
        "IM" = "pink",
        "UNS" = "black"
      )
    )
  )
ht_list


pdf(file = file.path(dir.plots,"consensus_heatmap_cell_lines_flip_S7b.pdf"),width = 10,height = 2.5)
draw(ht_list)
dev.off()
