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
# Date: 12 Mar 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(genefilter)
library(readxl)
library(limma) 
library(TCGAbiolinks)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(msigdbr)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Folders where we was saved the data
# See script: 01_data_download.R
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data,"TCGA")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# loading pathways list
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
PCA_pathway_list <- readxl::read_xlsx(file.path(dir.data,"PCA pathway list_new_AC.xlsx")) %>% as.data.frame
PCA_pathway_list$subtype <- gsub("TNBC_", "",PCA_pathway_list$subtype)
PCA_pathway_list <- PCA_pathway_list %>% dplyr::filter(.data$dataset %in% c("REACTOME","KEGG","HALLMARK"))


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Get data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
load(file.path(dir.tcga,"TCGA_TNBC_192samples_RNASeq.Rdata"))
colnames(dataFilt) <- substr(colnames(dataFilt),1,12)

# Load supplemental table 2
metadata <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS") 


metadata <- metadata %>% group_by(subtype) %>% plyr::arrange(factor(subtype,levels = c("BL1","BL2","M","LAR")))

dataFilt <- dataFilt[,metadata$patient]


vx = vector(mode="list")
pvals = vector(mode="list")
tabLogFC = vector(mode="list")
k <- 1

for( collection in unique(PCA_pathway_list$dataset)){
    
    PCA_pathway_list_cur <- PCA_pathway_list[PCA_pathway_list$dataset %in% collection, ]
    
    if(collection == "HALLMARK"){
        geneSets <- msigdbr(species = "Homo sapiens",category = "H")
    }
    
    if(collection == "C2_CGP"){
        geneSets <- msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CGP")
    }
    
    if(collection == "REACTOME"){
        geneSets <- msigdbr(species = "Homo sapiens",category = "C2",subcategory = "REACTOME")
    }
    
    if(collection == "KEGG"){
        geneSets <- msigdbr(species = "Homo sapiens",category = "C2",subcategory = "KEGG")
    }
    
    TNBC_cur <- gsva(
        expr = dataFilt, 
        gset.idx.list = split(x = geneSets$gene_symbol, f = geneSets$gs_name), 
        min.sz = 10, 
        max.sz = 500, 
        verbose = TRUE, 
        parallel.sz = 1
    )
    
    designTNBC <- matrix(0,ncol(dataFilt),4)
    colnames(designTNBC) <- unique(metadata$subtype)
    rownames(designTNBC) <- colnames(dataFilt)
    
    for (subtype in unique(metadata$subtype)){
        metadata_cur <- metadata[metadata$subtype %in% subtype,]
        designTNBC[metadata_cur$patient,subtype] <- 1
    }
    
    fitTNBC <- lmFit(TNBC_cur, designTNBC)
    fitTNBC_eBayes <- eBayes(fitTNBC)
    
    topTabBL1 <- topTable(fitTNBC_eBayes, coef="BL1",number = nrow(TNBC_cur))
    topTabBL2 <- topTable(fitTNBC_eBayes, coef="BL2",number = nrow(TNBC_cur))
    topTabLAR <- topTable(fitTNBC_eBayes, coef="LAR",number = nrow(TNBC_cur))
    topTabM <- topTable(fitTNBC_eBayes, coef = "M",number = nrow(TNBC_cur))
    
    commonPathways <- intersect(PCA_pathway_list_cur$pathway, rownames(TNBC_cur))
    
    if (length(commonPathways) == 0){
        curPathways <- paste0(collection,"_",PCA_pathway_list_cur$pathway)
        commonPathways <- intersect(curPathways, rownames(TNBC_cur))
    }
    
    TNBC_cur_sel <- TNBC_cur[commonPathways,]
    
    topTabBL1_sel <- topTabBL1[commonPathways,]
    topTabBL2_sel <- topTabBL2[commonPathways,]
    topTabLAR_sel <- topTabLAR[commonPathways,]
    topTabM_sel <- topTabM[commonPathways,]
    
    PvalMat_sel <- matrix(0,length(commonPathways),length(unique(metadata$subtype)))
    colnames(PvalMat_sel) <- unique(metadata$subtype)
    rownames(PvalMat_sel) <- commonPathways
    
    PvalMat_sel[commonPathways, "BL1"] <- topTabBL1_sel[commonPathways,"adj.P.Val"]
    PvalMat_sel[commonPathways, "BL2"] <- topTabBL2_sel[commonPathways,"adj.P.Val"]
    PvalMat_sel[commonPathways, "LAR"] <- topTabLAR_sel[commonPathways,"adj.P.Val"]
    PvalMat_sel[commonPathways, "M"] <- topTabM_sel[commonPathways,"adj.P.Val"]
    
    logFCMat_sel <- matrix(0,length(commonPathways),length(unique(metadata$subtype)))
    colnames(logFCMat_sel) <- unique(metadata$subtype)
    rownames(logFCMat_sel) <- commonPathways
    
    logFCMat_sel[commonPathways, "BL1"] <- topTabBL1_sel[commonPathways,"logFC"]
    logFCMat_sel[commonPathways, "BL2"] <- topTabBL2_sel[commonPathways,"logFC"]
    logFCMat_sel[commonPathways, "LAR"] <- topTabLAR_sel[commonPathways,"logFC"]
    logFCMat_sel[commonPathways, "M"] <- topTabM_sel[commonPathways,"logFC"]
    
    logFCMat_sel[logFCMat_sel < 0] <- 0
    
    rownames(TNBC_cur_sel) <- 
        gsub("_", " ",
             gsub(
                 paste0(collection,"_"),"",
                 rownames(TNBC_cur_sel)
             )
        )
    
    vx[[k]] <- TNBC_cur_sel
    pvals[[k]] <- PvalMat_sel
    tabLogFC[[k]] <- logFCMat_sel
    
    names(vx)[k] <- collection
    names(pvals)[k] <- collection
    names(tabLogFC)[k] <- collection
    
    k <- k + 1
}

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Draw heatmap
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
ha1 = HeatmapAnnotation(
    "TNBC_subtype" = metadata$subtype,
    which = "column",
    show_annotation_name = FALSE,
    col = list(
        TNBC_subtype = c(
            "BL1" = "red", 
            "BL2" = "orange",
            "LAR" = "green",
            "M" = "blue")
    ),
    annotation_legend_param = list(
        "TNBC_subtype"  = list(
            nrow = 1,
            #direction = "horizontal",
            labels_gp = gpar(fontsize = 10),
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            #legend_width  = unit(3, "cm"),
            title = "TNBC subtype")
    )
)

col_pathways = colorRamp2(c(0,0.5, 1), c("blue", "yellow","red"))

row_ha_1 = rowAnnotation(
    BL1 = anno_barplot(tabLogFC[[1]][,"BL1"],gp = gpar(fill = "red"),axis = FALSE, ylim = c(0,0.35)),
    BL2 = anno_barplot(tabLogFC[[1]][,"BL2"],gp = gpar(fill = "orange"),axis = FALSE, ylim = c(0,0.35)),
    M = anno_barplot(tabLogFC[[1]][,"M"],gp = gpar(fill = "blue"),axis = FALSE, ylim = c(0,0.35)),
    LAR = anno_barplot(tabLogFC[[1]][,"LAR"],gp = gpar(fill = "green"),axis = FALSE, ylim = c(0,0.35)),
    annotation_name_side = "top"
)

ht1 <- Heatmap(
    matrix = vx[[1]],
    name = names(vx)[1],
    use_raster = TRUE,
    raster_device = c("png"),
    raster_quality = 2,
    col = col_pathways, 
    column_split = factor(metadata$subtype,levels = unique(metadata$subtype)),
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE, 
    show_heatmap_legend = FALSE,
    cluster_rows  = TRUE,
    show_row_dend = FALSE,    
    row_title = names(vx)[1],
    top_annotation = ha1,
    width = unit(7,"cm"),
    row_names_gp = gpar(fontsize = 4),
    show_column_names = FALSE,
    right_annotation = row_ha_1
)

row_ha_2 = rowAnnotation(
    BL1 = anno_barplot(tabLogFC[[2]][,"BL1"],gp = gpar(fill = "red"),axis = FALSE, ylim = c(0,0.35)),
    BL2 = anno_barplot(tabLogFC[[2]][,"BL2"],gp = gpar(fill = "orange"),axis = FALSE, ylim = c(0,0.35)),
    M = anno_barplot(tabLogFC[[2]][,"M"],gp = gpar(fill = "blue"),axis = FALSE, ylim = c(0,0.35)),
    LAR = anno_barplot(tabLogFC[[2]][,"LAR"],gp = gpar(fill = "green"),axis = FALSE, ylim = c(0,0.35)),
    show_annotation_name = FALSE
)


ht2 <- Heatmap(
    matrix = vx[[2]],
    use_raster = TRUE,
    raster_device = c("png"),
    raster_quality = 2,
    name = names(vx)[2],
    col = col_pathways,
    cluster_columns = FALSE,
    show_heatmap_legend = FALSE,
    cluster_rows  = TRUE,
    show_row_dend = FALSE,    
    #top_annotation = ha1,
    row_names_gp = gpar(fontsize = 4),
    show_column_names = FALSE,
    row_title = names(vx)[2],
    right_annotation = row_ha_2)


row_ha_3 = rowAnnotation(
    BL1 = anno_barplot(tabLogFC[[3]][,"BL1"],gp = gpar(fill = "red"), ylim = c(0,0.35)),
    BL2 = anno_barplot(tabLogFC[[3]][,"BL2"],gp = gpar(fill = "orange"), ylim = c(0,0.35)),
    M = anno_barplot(tabLogFC[[3]][,"M"],gp = gpar(fill = "blue"), ylim = c(0,0.35)),
    LAR = anno_barplot(tabLogFC[[3]][,"LAR"],gp = gpar(fill = "green"), ylim = c(0,0.35)),
    show_annotation_name = FALSE
)


ht3 <- Heatmap(
    matrix = vx[[3]],
    name = names(vx)[3],
    col = col_pathways,
    cluster_columns = FALSE,
    show_heatmap_legend = TRUE,
    use_raster = TRUE,
    raster_device = c("png"),
    raster_quality = 2,
    heatmap_legend_param = list(
        direction = "horizontal",
        legend_width  = unit(2, "cm"),
        labels_gp = gpar(fontsize = 10),
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        title = "GSVA ES"),
    cluster_rows  = TRUE,
    show_row_dend = FALSE,    
    row_names_gp = gpar(fontsize = 4),
    show_column_names = FALSE,
    row_title = names(vx)[3],
    right_annotation = row_ha_3
)

pdf(file.path(dir.plots,"TNBC_Pathways_heatmap_v7.pdf"), width = 10,height = 10)
ht_list <- ht1 %v% ht2 %v% ht3
draw(
    ht_list, 
    ht_gap = unit(1, "mm"),
    heatmap_legend_side = "top",
    annotation_legend_side = "top")
dev.off()
