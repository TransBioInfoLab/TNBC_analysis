#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# This script performs ELMER analysis on TCGA data.
# ELMER compares distal probes DNA methylation (HM450) and identifies
# putative target genes affected by the change of methylation
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
library(TCGAbiolinks)
library(readxl)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(MultiAssayExperiment)
library(ELMER.data)
library(ELMER)
library(dplyr)
library(png)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Folders where we was saved the data
# See script: 01_data_download.R
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data,"TCGA")
dir.output <-  "analysis_results/TCGA/ELMER_analysis_hg19"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Get metadata
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=

# Load supplemental table 2
metadata <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS") 


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Get data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# RNA-seq
load(file.path(dir.tcga,"TCGA_TNBC_192samples_RNASeq.Rdata"))
rna <- dataFilt
rownames(rna) <- rowRanges(BRCA.exp)$ensembl_gene_id[match(rownames(rna),rowRanges(BRCA.exp)$gene_id)]
rna <- rna[!is.na(rownames(rna)),]

# DNA methylation
met <- get(load(file.path(dir.tcga,"TCGA_TNBC_192samples_Methylation.Rdata")))

#-----------------------------------------------------------
# Create MAE object Only tumor and distal probes
#-----------------------------------------------------------
# Distal probes
genome <- "hg19"
met.platform <- "450K"
distal.probes <- get.feature.probe(genome = genome, met.platform = met.platform)

mae.distal <- createMAE(
  exp = rna, 
  met = met,
  save = FALSE,
  linearize.exp = TRUE,
  filter.probes = distal.probes,
  met.platform = met.platform,
  genome = genome,
  TCGA = TRUE
)

mae.distal$TNBC_subtype <- gsub("[0-9]-","",metadata$subtype[match(mae.distal$patient,metadata$patient)])
for(i in unique(mae.distal$TNBC_subtype)) {
  colData(mae.distal)[[paste0(i,"_vs_others")]] <- "Others"
  colData(mae.distal)[[paste0(i,"_vs_others")]][mae.distal$TNBC_subtype == i] <- i
}

save(mae.distal,file = "data/mae_TNBC_with_groups_distal_probes_hg19.rda")

#-------------------------------------------------------------------------------
# ELMER analyses have 5 main steps:
# 1. Identify distal probes on HM450K or EPIC arrays.
# 2. Identify distal probes with significantly different DNA methylation level between two groups
# 3. Identify putative target genes for differentially methylated distal probes.
# 4. Identify enriched motifs for the distal probes which are significantly 
#    differentially methylated and linked to putative target gene.
# 5. Identify regulatory TFs whose expression associate with DNA methylation at enriched motifs.
#-------------------------------------------------------------------------------
# We wil be performing step 1 to 3 to each TNBC type vs all other ones
# LAR, BL1, BL2, M. UNS (unstable) were removed
#-----------------------------------------------------------
# Perform DNA methylation analysis
#-----------------------------------------------------------
mae.distal <- mae.distal[,mae.distal$TNBC_subtype != "UNS"]
#--------------------------------------------------------------------------------
# Analysis 1: one vs all others
#--------------------------------------------------------------------------------
cores <- 1
mode <- "supervised"

for(tnbc.subtype in unique(mae.distal$TNBC_subtype)) {
  for(direction in c(c("hypo","hyper"))) {
    group.col <- paste0(tnbc.subtype,"_vs_others")
    group1 <- tnbc.subtype
    group2 <- "Others"
    message(group.col, ":", direction, " ", group1, " vs ", group2)
    
    dir.out <- paste0(
      dir.output,"/",
      mode,"/",
      group.col,"-",
      gsub("[[:punct:]]| ", ".", group1),
      "_vs_",
      gsub("[[:punct:]]| ", ".", group2),
      "/", direction
    )
    
    # if (length(dir(dir.out, recursive = T, pattern = "getTF", full.names = T)) > 0) next
    
    sig.diff <- get.diff.meth(
      data = mae.distal, 
      group.col = group.col,
      group1 =  group1,
      group2 = group2,
      mode = mode,
      sig.dif = 0.2,
      diff.dir = direction,
      cores = cores, 
      dir.out = dir.out,
      pvalue = 0.01
    )
    
    nearGenes <- GetNearGenes(
      data = mae.distal, 
      probes = sig.diff$probe, 
      numFlankingGenes = 20
    ) # 10 upstream and 10 dowstream genes
    
    # evaluate correlation bwt gene expression and DNA methylation
    pair <- get.pair(
      data = mae.distal,
      group.col = group.col,
      group1 =  group1,
      group2 = group2,
      nearGenes = nearGenes,
      mode = mode,
      permu.dir = paste0(dir.out,"/permu"),
      permu.size = 10000, # Please set to 100000 to get significant results
      raw.pvalue = 0.05,   
      Pe = 0.05, # Please set to 0.001 to get significant results
      filter.probes = TRUE, # See preAssociationProbeFiltering function
      filter.percentage = 0.05,
      filter.portion = 0.3,
      dir.out = dir.out,
      diff.dir = direction ,
      cores = cores,
      label = direction
    )
    
    tryCatch({
      heatmapPairs(
        data = mae.distal, 
        group.col = group.col,
        group1 =  group1,
        group2 = group2,
        annotation.col = c("TNBC_subtype"),
        pairs = pair,
        filename =  paste0(dir.out,"/heatmap_pairs.pdf")
      )
    }, error = function(e) print(e))
    
    # 4. Identify enriched motifs for the distal probes which are significantly 
    # differentially methylated and linked to putative target gene.
    if(length(unique(pair$Probe)) < 10) next
    
    enriched.motif <- get.enriched.motif(
      data = mae.distal,
      probes = unique(pair$Probe), 
      dir.out = dir.out,
      label = direction,
      min.incidence = 10,
      lower.OR = 1.1
    )
    
    # 5. Identify regulatory TFs whose expression associate with 
    # DNA methylation at enriched motifs.
    TF <- get.TFs(
      data = mae.distal, 
      group.col = group.col,
      group1 =  group1,
      group2 = group2,
      mode = mode,
      diff.dir = direction,
      enriched.motif = enriched.motif,
      dir.out = dir.out,
      cores = cores, 
      label = direction,
      save.plots = FALSE
    )
    
  }
}

#============================================================
# Heatmap
#============================================================
pair.files  <-  dir(
  dir.output,
  pattern = ".*.pairs.significant.csv", 
  full.names = T,
  recursive = T,
  all.files = T
) 

analysis <- file.path(basename(dirname(dirname(pair.files))), basename(dirname(pair.files)))
all.pairs <- NULL
all.pairs <- plyr::adply(.data = pair.files,1, function(f) readr::read_csv(f,col_types = readr::cols()),.id = "Analysis")
all.pairs$Analysis <- gsub("_", " ",gsub("TNBC_subtype-","",paste0(pair.files %>% dirname  %>% basename," in ",pair.files %>% dirname %>% dirname %>% basename)))[all.pairs$Analysis]

all.pairs$id <- paste0(all.pairs$Probe,"_",all.pairs$GeneID)
all.pairs <- all.pairs[!duplicated(all.pairs$id),]

all.pairs.fdr.0.01 <- all.pairs[all.pairs$FDR < 0.01,]
all.pairs.fdr.0.05 <- all.pairs[all.pairs$FDR < 0.05,]
all.pairs.fdr.0.001 <- all.pairs[all.pairs$FDR < 0.001,]

#============================================================
# Color Palettes
#============================================================

# 1) For RNA Heatmap
pal_rna <- colorRampPalette(
  c("#352A86","#343DAE","#0262E0","#1389D2",
    "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D")
)(100)

# 2) For Methylation Heatmap
pal_methylation <- colorRampPalette(
  c("#000436","#021EA9","#1632FB",
    "#6E34FC","#C732D5","#FD619D",
    "#FF9965","#FFD32B","#FFFC5A")
)(100)

# 3) For ATAC Heatmap
pal_atac <- colorRampPalette(
  c('#3361A5', '#248AF3', '#14B3FF', 
    '#88CEEF', '#C1D5DC', '#EAD397', 
    '#FDB31A','#E42A2A', '#A31D1D')
)(100)

# 4) For ATAC Bigwig Heatmap style
pal_atac_bw_heatmap <- colorRampPalette(
  c("white","#2488F0","#7F3F98",
    "#E22929","#FCB31A")
)(100)

#============================================================
# RNA-seq and DNAm plot
#============================================================
group.col <- "TNBC_subtype"

col.order.met <- plyr::alply(
  .data = c("M","BL1","BL2","LAR"), 
  .margins = 1, 
  .fun = function(x) {
    idx <- which(colData(mae.distal)[, group.col] == x)
    aux <- na.omit(assay(getExp(mae.distal)[all.pairs.fdr.0.001$GeneID,idx]))
    order <- t(aux) %>% dist %>% hclust(method = "average")
    as.numeric(idx[order$order])
  }) %>% unlist


# DNA methylation and gene expression annotation
ha = HeatmapAnnotation(
  "TNBC subtype" = mae.distal$TNBC_subtype, 
  col = list(
    "TNBC subtype" = c(
      "BL1" = "red",
      "BL2" = "orange",
      "LAR" = "green",
      "M" = "blue"
    )
  ),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 6)
)

# DNA methylation matrix
ht_list <- Heatmap(
  assay(getMet(mae.distal)[all.pairs.fdr.0.001$Probe,]), 
  name = "Methylation Beta-value", 
  col = colorRamp2(seq(0, 1, by = 1/99), pal_methylation),
  column_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(
    legend_width = unit(5, "cm"),
    legend_direction = "horizontal",
    labels_gp = gpar(fontsize = 12), 
    title_gp = gpar(fontsize = 12)
  ),
  show_column_names = F,
  show_row_names = F,
  use_raster = TRUE,
  clustering_method_rows = "average",
  clustering_distance_rows = "manhattan",
  raster_device = c("png"),
  raster_quality = 2,
  cluster_columns = FALSE,
  #row_order = order.heatmap,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  row_title = paste0(nrow(all.pairs.fdr.0.001), " linked distal probes"),
  column_order = col.order.met,
  row_names_gp = gpar(fontsize = 12),
  top_annotation = ha,
  width = unit(7, "cm"),
  column_title =  paste0("DNA methylation beta values (n = ",length(col.order.met),")"), 
  column_title_gp = gpar(fontsize = 12), 
  row_title_gp = gpar(fontsize = 16)
) 


# RNA-seq matrix
ht_list <- ht_list + 
  Heatmap(
    t(scale (t(assay(getExp(mae.distal)[all.pairs.fdr.0.001$GeneID,])))), 
    name = "RNA-seq (z-score)", 
    col = colorRamp2(seq(-2, 2, by = 4/99), pal_rna),
    column_names_gp = gpar(fontsize = 8),
    show_column_names = FALSE,
    heatmap_legend_param = list(
      legend_width = unit(5, "cm"),
      legend_direction = "horizontal",
      labels_gp = gpar(fontsize = 12), 
      title_gp = gpar(fontsize = 12)
    ),
    show_row_names = FALSE,
    cluster_columns = FALSE,
    #row_order = order.heatmap,
    use_raster = TRUE,
    raster_device = c("png"),
    raster_quality = 2,
    cluster_rows = FALSE,
    column_order = col.order.met,
    row_names_gp = gpar(fontsize = 12),
    top_annotation = ha,
    width = unit(7, "cm"),
    column_title = paste0("RNA-Seq Row-wise z-score (n = ", length(col.order.met),")"), 
    column_title_gp = gpar(fontsize = 12), 
    row_title_gp = gpar(fontsize = 16)
  ) 


pdf(file.path(dir.plots,"ELMER_heatmap_pairs_fdr_0_001_hg19.pdf"), width = 12, height = 4)
draw(
  object = ht_list,
  newpage = TRUE, 
  column_title = "TCGA-BRCA ELMER anti-correlated paired gene-probes", 
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right"
)
dev.off()