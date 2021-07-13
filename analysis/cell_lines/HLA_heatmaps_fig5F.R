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
library(ComplexHeatmap)
library(circlize)
library(dplyr)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/cell_lines/"
dir.plots <- "plots/cell_lines"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#============================================================
# Color Palettes
#============================================================

# 1) For RNA Heatmap
pal_rna <- colorRampPalette(
  c("#352A86","#343DAE","#0262E0","#1389D2",
    "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D")
)(100)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
nlrc5 <- readr::read_csv(file.path(dir.data, "NLRC5_HLA_heat.csv"))
ciita <- readr::read_csv(file.path(dir.data, "CIITA_HLA_heat.csv"))

data <- dplyr::full_join(ciita, nlrc5)
matrix <- data[,grep("SAMPLE|CELL|TREAT|COLOR",colnames(data),invert = TRUE)] %>% t
colnames(matrix) <- data$SAMPLE

col.idx <- c(
  grep("CAL51.*UNT",colnames(matrix),value = TRUE),
  grep("CAL51.*TAZMETOSTAT",colnames(matrix),value = TRUE),
  grep("CAL51.*CP1205",colnames(matrix),value = TRUE),
  grep("CAL51.*MAK683",colnames(matrix),value = TRUE),
  grep("CAL120.*UNT",colnames(matrix),value = TRUE),
  grep("CAL120.*TAZMETOSTAT",colnames(matrix),value = TRUE),
  grep("CAL120.*CPI1205",colnames(matrix),value = TRUE),
  grep("CAL120.*MAK683",colnames(matrix),value = TRUE),
  grep("BT549.*UNT",colnames(matrix),value = TRUE),
  grep("BT549.*TAZMETOSTAT",colnames(matrix),value = TRUE),
  grep("BT549.*CP1205",colnames(matrix),value = TRUE),
  grep("BT549.*MAK683",colnames(matrix),value = TRUE)
)
# Change gene and samples order
matrix <- matrix[c(7:11,1:6),col.idx]

ha = HeatmapAnnotation(
  "Treatment" = data$TREAT[match(col.idx,data$SAMPLE)],
  col = list(
    "Treatment" = c(
      "1-UNT" = "#FDE725FF",
      "3-CP1205" = "#35B779FF",
      "2-TAZMETOSTAT" = "#31688EFF",
      "4-MAK683" = "#440154FF"
    )
  ),
  annotation_legend_param =  list(
      legend_direction = "horizontal",
      nrow = 2, 
      ncol = 2,
      labels_gp = gpar(fontsize = 12), 
      title_gp = gpar(fontsize = 12)
    ),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  which = "column",
  annotation_name_gp = gpar(fontsize = 6)
)

ht_list <- Heatmap(
  matrix %>% t %>% scale %>% t,, 
  top_annotation = ha,
  cluster_column_slices = FALSE,
  column_split = factor(data$CELL[match(col.idx,data$SAMPLE)],levels = c("CAL51","CAL120","BT549")),
  row_split = c(rep("MHC-I",5),rep("MHC-II",6)),
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    legend_width  = unit(3, "cm"),
    labels_gp = gpar(fontsize = 12), 
    title_gp = gpar(fontsize = 12)
  ),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_side = "left",
  cluster_columns = FALSE,
  name = "RNA-seq\n(row z-score)")

pdf(file.path(dir.plots,"fig5F.pdf"), width = 5, height = 3.5)
draw(
  ht_list,
  newpage = TRUE, 
  merge_legend = TRUE,
  column_title = "", 
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()
