#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# This script performs scRNA analysis on TNBC tumors from GSE118390 (Karaayvaz et al. 2018).
# This analysis identifies epitheilial cells through common markers and creates an new seraut
# object of all epithelial cells or epithelial cells from each tumor
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
library(data.table)
library(R.utils)
library(gdata)
library(Seurat)
library(sctransform)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots/scRNA"
dir.scrnaseq <- file.path(dir.data,"scRNA_data")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-----------------------------------------------------------
# Get data and meta
#-----------------------------------------------------------
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118389/suppl/GSE118389_counts_rsem.txt.gz"

file <-  file.path(dir.scrnaseq,basename(url))
if(!file.exists(file)) downloader::download(url, file)

data <- read.table(file)
meta <- read.csv(file.path(dir.scrnaseq,"kara_meta.csv"), header=TRUE, row.names=1, as.is=TRUE)

#-----------------------------------------------------------
# create seraut object
kara.seurat <- CreateSeuratObject(
    data,
    project = "kara",
    meta.data = meta,
    min.cells = 3,
    min.features = 200
  )

# normalize seraut object
kara.seurat <- NormalizeData(kara.seurat)

kara.seurat <- FindVariableFeatures(
    object = kara.seurat,
    mean.function = ExpMean,
    dispersion.function = LogVMR,
    x.low.cutoff = 0.0125,
    x.high.cutoff = 3,
    y.cutoff = 0.5
  )

# Run PCA, UMAP and find clusters
kara.seurat <-ScaleData(kara.seurat)
kara.seurat <- RunPCA(kara.seurat, verbose = FALSE)
kara.seurat <- RunUMAP(kara.seurat, dims = 1:30, verbose = FALSE)
kara.seurat <- RunTSNE(kara.seurat)
kara.seurat <- FindNeighbors(kara.seurat, dims = 1:30, verbose = FALSE)
kara.seurat <- FindClusters(kara.seurat, verbose = FALSE)

###save processed file
saveRDS(kara.seurat, file = file.path(dir.scrnaseq,"kara_seruat.rds"))

#-----------------------------------------------------------
### read in saved file
#kara.seurat=readRDS("kara.seurat.rds")
#-----------------------------------------------------------


#create plot by cluster
DimPlot(kara.seurat, label = TRUE) + NoLegend()
DimPlot(kara.seurat, label = FALSE) 

pdf("kara cluster.pdf", width=4, height=4)
DimPlot(kara.seurat, label = FALSE) 
dev.off()

#create plot by tnbc subtype
Idents(kara.seurat) <- "subtype"

pdf("kara subtype all.pdf", width=4, height=4)
DimPlot(kara.seurat, label = FALSE) 
dev.off()

#create plot by patient id
Idents(kara.seurat) <- "Patient"

pdf("kara patient all.pdf", width=4, height=4)
DimPlot(kara.seurat, label = FALSE) 
dev.off()

Idents(kara.seurat) <- "CellType"
pdf("kara celltype no legend.pdf", width=4, height=4)
DimPlot(kara.seurat, label = FALSE)+ NoLegend() 
dev.off()



#-----------------------------------------------------------
### Identify cell type with common markers

###Bcell
pdf(file.path(dir.plots,"kara_Cd79A.pdf"), width=4, height=4)
FeaturePlot(kara.seurat,reduction = "umap", features = c("CD79A"), pt.size = 0.2,cols= c("grey","red")) 
dev.off()

###endothelial
pdf(file.path(dir.plots,"kara_CD31_PECAM1.pdf"), width=4, height=4)
FeaturePlot(kara.seurat,reduction = "umap", features = c("PECAM1"), pt.size = 0.2,cols= c("grey","red")) 
dev.off()

###Tcell
pdf(file.path(dir.plots,"kara_CD8A.pdf"), width=4, height=4)
FeaturePlot(kara.seurat,reduction = "umap", features = c("CD8A"), pt.size = 0.2,cols= c("grey","red")) 
dev.off()

###Tcell
pdf(file.path(dir.plots,"kara_CD3E.pdf"), width=4, height=4)
FeaturePlot(kara.seurat,reduction = "umap", features = c("CD3E"), pt.size = 0.2,cols= c("grey","red")) 
dev.off()

###myoepithelial
pdf(file.path(dir.plots,"kara_p63.pdf"), width=4, height=4)
FeaturePlot(kara.seurat,reduction = "umap", features = c("TP63"), pt.size = 0.2,cols= c("grey","red")) 
dev.off()

###monocyte
pdf(file.path(dir.plots,"kara_CD14.pdf"), width=4, height=4)
FeaturePlot(kara.seurat,reduction = "umap", features = c("CD14"), pt.size = 0.2,cols= c("grey","red")) 
dev.off()

###epithelial
pdf(file.path(dir.plots,"kara_epcam.pdf"), width=4, height=4)
FeaturePlot(kara.seurat,reduction = "umap", features = c("EPCAM"), pt.size = 0.2,cols= c("grey","red")) 
dev.off()


#-----------------------------------------------------------
###Select only epithelial cell

Idents(kara.seurat) <- "Epi"
pdf("kara epi all.pdf", width=4, height=4)
DimPlot(kara.seurat, label = FALSE) 
dev.off()



#-----------------------------------------------------------
#### create new seraut object on epithelial cells
Idents(kara.seurat) <- "Epi"
kara_epi = subset(kara.seurat, idents = c("Epithelial"), invert = FALSE)

kara_epi <- RunPCA(kara_epi, verbose = FALSE)
kara_epi <- RunUMAP(kara_epi, dims = 1:30, verbose = FALSE)
kara_epi <- RunTSNE(kara_epi)
kara_epi <- FindNeighbors(kara_epi, dims = 1:30, verbose = FALSE)
kara_epi <- FindClusters(kara_epi, verbose = FALSE)


DimPlot(kara_epi, label = TRUE) + NoLegend()

####epithelial cells by patient
pdf(file.path(dir.plots,"kara epi_patient.pdf"), width=4, height=4)
Idents(kara_epi) <- "Patient"
DimPlot(kara_epi, label = TRUE) + NoLegend()
dev.off()

####epithelial cells by subtype
Idents(kara_epi) <- "subtype"

pdf(file.path(dir.plots,"kara epi_subtype.pdf"), width=4, height=4)
DimPlot(kara_epi, label = FALSE, cols=c("orange","red","blue", "green","black")) 
dev.off()


Idents(kara_epi) <- "seurat_clusters"
DimPlot(kara_epi, label = FALSE) 
DimPlot(kara_epi, label = TRUE) + NoLegend()

#-----------------------------------------------------------
######subset epithelial cells by each patient and color by TNBC subtype
Idents(kara_epi) <- "Patient"

###subset by patient 89 
PT089 = subset(kara_epi, idents = c("PT089"), invert = FALSE)
DimPlot(PT089, label = FALSE) 
Idents(PT089) <- "subtype"
pdf(file.path(dir.plots,"PT89_subtype.pdf"), width=4, height=4)
DimPlot(PT089, label = FALSE, cols=c("orange","red","blue", "green","black")) 
dev.off()

#FeaturePlot(object = PT089, features = "LAR")

###subset by patient 39 
PT039 = subset(kara_epi, idents = c("PT039"), invert = FALSE)
DimPlot(PT039, label = FALSE) 
Idents(PT039) <- "subtype"
pdf(file.path(dir.plots,"PT39_subtype.pdf"), width=4, height=4)
DimPlot(PT039, label = FALSE, cols=c("red","blue","orange", "green")) 
dev.off()

###subset by patient 858
PT058 = subset(kara_epi, idents = c("PT058"), invert = FALSE)
DimPlot(PT058, label = FALSE) 
Idents(PT058) <- "subtype"
pdf(file.path(dir.plots,"PT58_subtype.pdf"), width=4, height=4)
DimPlot(PT058, label = FALSE,cols=c("orange","black","red", "green")) 
dev.off()

###subset by patient 81
PT081=subset(kara_epi, idents = c("PT081"), invert = FALSE)
DimPlot(PT081, label = FALSE) 
Idents(PT081) <- "subtype"
pdf(file.path(dir.plots,"PT81_subtype.pdf"), width=4, height=4)
DimPlot(PT081, label = FALSE,cols=c("red","orange", "blue","black")) 
dev.off()

###subset by patient 84
PT084=subset(kara_epi, idents = c("PT084"), invert = FALSE)
DimPlot(PT084, label = FALSE) 
Idents(PT084) <- "subtype"
pdf(file.path(dir.plots,"PT84_subtype.pdf"), width=4, height=4)
DimPlot(PT084, label = FALSE,cols=c("red","orange", "green","blue")) 
dev.off()

###subset by patient126
PT126=subset(kara_epi, idents = c("PT126"), invert = FALSE)
DimPlot(PT126, label = FALSE) 
Idents(PT126) <- "subtype"
pdf(file.path(dir.plots,"PT126_subtype.pdf"), width=4, height=4)
DimPlot(PT126, label = FALSE,cols=c("red","green", "black","orange","blue")) 
dev.off()



