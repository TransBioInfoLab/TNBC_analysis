#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# The following script will perform a consensus clustering
# CC provides quantitative and visual ‘stability’ evidence derived
# from repeated subsampling and clustering. 
# CC reports a consensus of these repetitions, which is robust relative to sampling variability
# Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2881355/
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
# Date: 01 Feb 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(CancerSubtypes)
library(dplyr)
library(circlize)
library(DESeq2)

dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data,"TCGA")
dir.Gistic2 <-  file.path(dir.data,"GISTIC/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data - Please run 01_data_download.R before to have this data available
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matrix_data <- get(load(file.path(dir.tcga,"mbatch_TNBC_legacy_RNA-Seq-gene-normalized-v2_matrix_data.rda")))
matrix_data <- log2(matrix_data + 1)

dataRNA_var <- FSbyVar(matrix_data, cut.type = "topk",value = 5000)
dataFilt_scale <- pheatmap:::scale_rows(matrix_data[rownames(dataRNA_var),])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get metadata in the same order as the data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load supplemental table 2
S2 <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS") 


dataFilt_scale <- dataFilt_scale[,substr(colnames(dataFilt_scale),1,12) %in% S2$patient]
TNBC_pData <- S2[match(substr(colnames(dataFilt_scale),1,12),S2$patient),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get metadata in the same order as the data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

result <- ExecuteCC(
  clusterNum = 5,
  d = dataFilt_scale,
  maxK = 10,
  reps = 500, 
  pItem = 0.8, 
  pFeature = 1,
  verbose = TRUE,
  plot = "png",
  distance = "euclidean",
  clusterAlg = "km",
  title = "TNBC_CC"
)

sil <- silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
p <- plot(sil)

cc_matrix <- result$originalResult[[5]]$consensusMatrix %>% as.data.frame
colnames(cc_matrix) <- substr(names(result$originalResult[[5]]$consensusClass),1,15)
rownames(cc_matrix) <- substr(names(result$originalResult[[5]]$consensusClass),1,15)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adding cluster membership threshold
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We determined membership as the average cluster for each sample from all members 
# of group and selected samples >0.5

# selection of groups
cc_groups <- result$originalResult[[5]]$consensusClass %>% as.data.frame %>% 
  dplyr::rename(CC_group = ".",)  %>% 
  dplyr::mutate(barcode15 = substr(rownames(.),1,15), membership = 0, membership_tresh = CC_group) %>% 
  dplyr::arrange(CC_group)
rownames(cc_groups) <- cc_groups$barcode15

# For each sample we will get the mean distance of all other member of the 
# cluster to this sample
for(sample in cc_groups$barcode15){
  sample.group <- cc_groups[sample,"CC_group"]
  samples.in.same.group <- cc_groups %>% dplyr::filter(CC_group %in% sample.group) %>% dplyr::pull(barcode15)
  mean.from.all.other.samples <- mean(as.numeric(cc_matrix[sample,setdiff(samples.in.same.group,sample)]))
  cc_groups[sample,"membership"] <- mean.from.all.other.samples
}

cc_groups[cc_groups$membership < 0.5,"membership_tresh"] <- "vlow"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Heatmap
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ComplexHeatmap)
m <- result$distanceMatrix %>% as.dist() %>% as.matrix()
diag(m) <- 1
colnames(m) <- rownames(m) <- names(result$group)
group <- result$group
group.membership.values <- cc_groups$membership[match(substr(names(result$group),1,15), cc_groups$barcode15)]
group.membership <- cc_groups$membership_tresh[match(substr(names(result$group),1,15), cc_groups$barcode15)]
group.membership[group.membership != "vlow"] <- "FALSE"
group.membership[group.membership == "vlow"] <- "TRUE"

col = list(
  "TNBC type" = c(
    "BL1" = "red",
    "BL2" = "orange",
    "LAR" = "green",
    "M" = "blue",
    "UNS" = "black"
  ),
  "TNBC type - 6" = c(
    "BL1" = "red",
    "BL2" = "orange",
    "LAR" = "green",
    "M" = "blue",
    "IM" = "yellow",
    "MSL" = "purple",
    "UNS" = "black"
  ),
  "K-means group" = c(
    "1" = "purple",
    "2" = "darkblue", 
    "3" = "black",
    "4" = "gray",
    "5" = "darkgreen",
    "6" = "yellow"
  ),
  "Sample Group membership" = colorRamp2(c(0,1), c("white", "black")),
  "Is low group membership" = c("TRUE" = "orange","FALSE" = "white"),
  "LAR TNBCtype Corr." = colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "green")),
  "M TNBCtype Corr." = colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "blue")),
  "BL1 TNBCtype Corr." = colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "red")),
  "BL2 TNBCtype Corr." = colorRamp2(c(-0.7, 0, 0.7), c("gray", "white", "orange"))
)
hr <- HeatmapAnnotation(
  "K-means group" = group,
  "Sample Group membership" = group.membership.values, 
  "Is low group membership" = group.membership,
  "BL1 TNBCtype Corr." = TNBC_pData$BL1[match(substr(rownames(m),1,12),TNBC_pData$patient)],
  "BL2 TNBCtype Corr." = TNBC_pData$BL2[match(substr(rownames(m),1,12),TNBC_pData$patient)],
  "M TNBCtype Corr." = TNBC_pData$M[match(substr(rownames(m),1,12),TNBC_pData$patient)],
  "LAR TNBCtype Corr."  = TNBC_pData$LAR[match(substr(rownames(m),1,12),TNBC_pData$patient)],
  "TNBC type" = TNBC_pData$subtype[match(substr(rownames(m),1,12),TNBC_pData$patient)],
  which = "row",
  col = col
)

ht <- HeatmapAnnotation(
  "K-means group" = group,
  "Sample Group membership" = group.membership.values, 
  "Is low group membership" = group.membership,
  "BL1 TNBCtype Corr." = TNBC_pData$BL1[match(substr(rownames(m),1,12),TNBC_pData$patient)],
  "BL2 TNBCtype Corr." = TNBC_pData$BL2[match(substr(rownames(m),1,12),TNBC_pData$patient)],
  "M TNBCtype Corr." = TNBC_pData$M[match(substr(rownames(m),1,12),TNBC_pData$patient)],
  "LAR TNBCtype Corr."  = TNBC_pData$LAR[match(substr(rownames(m),1,12),TNBC_pData$patient)],
  "TNBC type" = TNBC_pData$subtype[match(substr(rownames(m),1,12),TNBC_pData$patient)],
  which = "column",
  annotation_name_side = "left",
  show_legend = TRUE,
  col = col, 
  annotation_legend_param = list(
    "K-means group" = list(
      direction = "horizontal"
    ),
    "Sample Group membership" = list(
      direction = "horizontal"
    ),
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
      legend_width  = unit(2, "cm"))
  )
  
)

aux <- data.frame(group)
aux$patient <- substr(rownames(aux),1,12)
aux <- merge(TNBC_pData,aux) %>% unique()
aux <- aux[order(aux$group,aux$subtype),]
order <- match(aux$patient,substr(rownames(m),1,12))

ht <- Heatmap(
  matrix = m,
  top_annotation = ht,
  #left_annotation = hr,
  row_order = order,
  column_order = order,
  heatmap_legend_param = list( 
    direction = "vertical"
  ),
  cluster_rows = FALSE,
  use_raster = TRUE,
  column_title = paste0("K-means using 5000 most variable genes"),
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  heatmap_height = unit(15,"cm"),
  name = "similarity",
  col = circlize::colorRamp2(c(0,1),c("white","black"))
)
ht

## Output
pdf(file.path(dir.plots,"cc_kmeans_k_5_sort_by_tnbc.pdf"),width = 10,height = 6)
ComplexHeatmap::draw(
  ht,
  newpage = TRUE, 
  merge_legends = TRUE
)
dev.off()