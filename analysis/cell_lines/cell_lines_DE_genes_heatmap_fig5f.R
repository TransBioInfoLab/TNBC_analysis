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
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Paths
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.data <- "data/cell_lines"
dir.plots <- "plots/cell_lines"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Colors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
drugs.colors <- c(
  "MAK683" = "#FDE725FF" ,
  "CPI1205"  = "#35B779FF",
  "CP1205"  = "#35B779FF",
  "TAZMETOSTAT" = "#31688EFF",
  "UNT" = "#440154FF" 
)


genes.labels <- read.csv(file.path(dir.data,"genes_to_label_FIG5C.csv"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CAL120
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CAL120 <- read.csv(
  file.path(dir.data,"DE_anno_CAL120.csv"),
  sep = ",",
  na.strings = "NA",
  row.names = 1,
  header = TRUE,
  comment.char = "#",
  stringsAsFactors = F
) %>% dplyr::arrange(order(logFC)) %>% dplyr::filter(abs(logFC) > 1.5)

rownames(CAL120) <- CAL120$gene_id
genes.highlight <- genes.labels %>% dplyr::filter(cell_line == "CAL120")  %>% pull(1) %>% as.character()
highlight.idx <- match(genes.highlight,rownames(CAL120))

drugs.CAL120 <- grep("CAL120",colnames(CAL120),value = T) %>% sapply(function(x) { unlist(strsplit(x,"_"))[2]})
CAL120_mat <- as.matrix(CAL120[, grep("CAL120",colnames(CAL120))])

col.order <- c(
  grep("UNT",colnames(CAL120_mat)),
  grep("TAZ",colnames(CAL120_mat)),
  grep("CPI",colnames(CAL120_mat)),
  grep("MAK",colnames(CAL120_mat))
)
drugs.CAL120 <- drugs.CAL120[col.order]
CAL120_mat <- CAL120_mat[,col.order]

ht.CAL120 <- Heatmap(
  CAL120_mat %>% t %>% scale %>% t,
  show_column_names = FALSE,
  show_row_names = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  column_title = "CAL120",
  width = unit(3,"cm"),
  height = unit(7,"cm"),
  row_title = paste0(nrow(CAL120_mat), " DE genes"),
  name = "RNA-Seq row z-score",
  show_row_dend = FALSE,
  raster_device = c("png"),
  raster_quality = 2,
  top_annotation = columnAnnotation(
    "Drugs" = drugs.CAL120,
    col = list("Drugs" = drugs.colors),
    show_legend = FALSE,
    show_annotation_name = FALSE)
) +   rowAnnotation(
  link = anno_mark(
    at = highlight.idx, 
    labels = rownames(CAL120_mat)[highlight.idx], 
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")
  )
)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CAL51
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CAL51 <- read.csv(
  file.path(dir.data,"DE_anno_CAL51.csv"),
  sep = ",",
  na.strings = "NA",
  row.names = 1,
  header = TRUE,
  comment.char = "#",
  stringsAsFactors = F
) %>% dplyr::arrange(order(logFC)) %>% dplyr::filter(abs(logFC) > 1.5)

rownames(CAL51) <- CAL51$gene_id
genes.highlight <-  genes.labels %>% dplyr::filter(cell_line == "CAL51")  %>% pull(1) %>% as.character()
highlight.idx <- match(genes.highlight,rownames(CAL51))

drugs.CAL51 <- grep("CAL51",colnames(CAL51),value = T) %>% sapply(function(x) { unlist(strsplit(x,"_"))[2]})
CAL51_mat <- as.matrix(CAL51[, grep("CAL51",colnames(CAL51))])

col.order <- c(
  grep("UNT",colnames(CAL51_mat)),
  grep("TAZ",colnames(CAL51_mat)),
  grep("CP",colnames(CAL51_mat)),
  grep("MAK",colnames(CAL51_mat)
  )
)

drugs.CAL51 <- drugs.CAL51[col.order]
CAL51_mat <- CAL51_mat[,col.order]


ht.CAL51 <- Heatmap(
  CAL51_mat %>% t %>% scale %>% t,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_title = "CAL51",
  cluster_columns = FALSE,
  width = unit(3,"cm"),
  height = unit(7,"cm"),
  show_heatmap_legend = FALSE,
  row_title = paste0(nrow(CAL51_mat), " DE genes"),
  name = "RNA-Seq row z-score",
  show_row_dend = FALSE,
  raster_device = c("png"),
  raster_quality = 2,
  top_annotation = columnAnnotation(
    "Drugs" = drugs.CAL51,
    col = list("Drugs" = drugs.colors),
    show_legend = FALSE,
    annotation_name_side = "right",
    show_annotation_name = FALSE
  )
) +   rowAnnotation(
  link = anno_mark(
    at = highlight.idx, 
    labels = rownames(CAL51)[highlight.idx], 
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")
  )
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BT549
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BT549 <- read.csv(
  file.path(dir.data,"DE_anno_BT549.csv"),
  sep = ",",
  na.strings = "NA",
  row.names = 1,
  header = TRUE,
  comment.char = "#",
  stringsAsFactors = F
) %>% dplyr::arrange(order(logFC)) %>% dplyr::filter(abs(logFC) > 1.5)

rownames(BT549) <- BT549$gene_id
genes.highlight <-  genes.labels %>% dplyr::filter(cell_line == "BT549")  %>% pull(1) %>% as.character()
highlight.idx <- match(genes.highlight,rownames(BT549))

drugs.BT549 <- grep("BT549",colnames(BT549),value = T) %>% sapply(function(x) { unlist(strsplit(x,"_"))[2]})
BT549_mat <- as.matrix(BT549[, grep("BT549",colnames(BT549))])

col.order <- c(
  grep("UNT",colnames(BT549_mat)),
  grep("TAZ",colnames(BT549_mat)),
  grep("CP",colnames(BT549_mat)),
  grep("MAK",colnames(BT549_mat))
)

drugs.BT549 <- drugs.BT549[col.order]
BT549_mat <- BT549_mat[,col.order]

ht.BT549 <- Heatmap(
  BT549_mat %>% t %>% scale %>% t,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_title = "BT549",
  use_raster = TRUE,
  width = unit(3,"cm"),
  height = unit(7,"cm"),
  heatmap_legend_param = list(
    title = "RNA-Seq\n(row z-score)",
    direction = "horizontal",
    legend_width  = unit(2, "cm"),
    labels_gp = gpar(fontsize = 8),
    title_gp = gpar(fontsize = 10, fontface = "bold")
    ),
  row_title = paste0(nrow(BT549_mat), " DE genes"),
  cluster_columns = FALSE,
  name = "RNA-Seq row z-score",
  show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  raster_device = c("png"),
  raster_quality = 2,
  top_annotation = columnAnnotation(
    "Drugs" = drugs.BT549,
    col = list("Drugs" = drugs.colors),
    show_legend = TRUE,
    annotation_legend_param = list("Drugs" = list(legend_gp = gpar(fontsize = 10),nrow = 1)),
    annotation_name_side = "right",
    show_annotation_name = TRUE
  )
) + rowAnnotation(
  link = anno_mark(
    at = highlight.idx, 
    labels = rownames(BT549_mat)[highlight.idx], 
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")
  )
)

# Put heatmaps together in the same figure
ht1 <- grid.grabExpr(
  draw(
    ht.CAL51, 
    merge_legend = TRUE, 
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
  )
)

ht2 <- grid.grabExpr(
  draw(
    ht.CAL120, 
    merge_legend = TRUE, 
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
  )
)


ht3 <- grid.grabExpr(
  draw(
    ht.BT549, 
    merge_legend = TRUE, 
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
  )
)

pdf(file.path(dir.plots,"EZH2_cell_lines_figure_5D_all.pdf"), height = 5, width = 12)

pushViewport(viewport(x = 0, y = 0.06, width = 0.3, height = 1, just = c("left", "bottom")))
grid.draw(ht1)
popViewport()

pushViewport(viewport(x = 0.3, y = 0, width = 0.3, height = 1, just = c("left", "bottom")))
grid.draw(ht3)
popViewport()

pushViewport(viewport(x = 0.6, y = 0.06, width = 0.3, height = 1, just = c("left", "bottom")))
grid.draw(ht2)
popViewport()
dev.off()


list <-  list(
  "CAL120" = rownames(CAL120_mat),
  "CAL51" = rownames(CAL51_mat),
  "BT549" = rownames(BT549_mat)
)
gplots::venn(list)
