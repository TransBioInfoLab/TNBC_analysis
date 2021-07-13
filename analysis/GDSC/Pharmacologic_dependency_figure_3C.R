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
# Date: 21 June 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ComplexHeatmap)
library(dplyr)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/GDSC/"
dir.plots <- "plots/GDSC"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
drug.sensitivity <- readr::read_csv(file.path(dir.data,"GDS_less0.1_T.csv")) %>% as.data.frame()

mat <- drug.sensitivity %>% dplyr::select(
  c("BL1_tvalue",  "BL2_tvalue", "M_tvalue","LAR_tvalue")
)

rownames(mat) <- drug.sensitivity$DRUG_NAME

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Heatmap colors
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
colors.tnbc <- c(
  "BL1" = "red",
  "BL2" = "orange",
  "LAR" = "green",
  "M" = "blue"
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
  "Chromatin" = "gold",
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
  "Misc." = "#038074",
  "apoptosis" = "#00533c",
  "RTK" = "pink",
  "Kinase" = "#ffc173",
  "kinase" = "#ffc173",
  "Other"  = "#038074"
)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Draw Heatmap
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
ha <- HeatmapAnnotation(
  "TNBCsubtype" = c("BL1","BL2","M","LAR"), 
  col = list("TNBCsubtype" = colors.tnbc),
  annotation_label = "TNBC subtype"
)
ra <- rowAnnotation("Pathways" = drug.sensitivity$Pathway, show_annotation_name = FALSE, col = list("Pathways" = colors.pathways))
ht <- Heatmap(
  matrix = mat,
  name = "T-value",
  show_row_names = TRUE,
  top_annotation = ha,
  row_names_gp =  gpar(fontsize = 8),
  left_annotation = ra,
  row_split = factor(drug.sensitivity$Pathway,levels = unique(drug.sensitivity$Pathway)),
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,column_title_gp = gpar(fontsize = 8),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 8)
)
pdf(file = file.path(dir.plots,"GDS_less0.1_T.pdf"),width = 4,height = 10)

draw(
  ht,
  newpage = TRUE, 
  column_title = "", 
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()