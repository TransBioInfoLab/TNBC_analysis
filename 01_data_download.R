#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# The following script downloads all data used to perform TNBC analysis
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

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Required packages
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
library(TCGAbiolinks)
library(xCell) # Installed from https://github.com/dviraran/xCell
library(readr)
library(data.table)
library(dplyr)
library(EDASeq)
library(illuminaHumanv4.db)
library(readxl)
library(mygene)
library(SummarizedExperiment)
library(GenomicRanges)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Folders where we will save the data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.tcga <- file.path(dir.data,"TCGA")
dir.atac_seq <- file.path(dir.tcga,"ATAC-seq")
dir.cptac <- file.path(dir.data,"CPTAC")
dir.pdtx <- file.path(dir.data,"PDTX")
dir.met500 <- file.path(dir.data,"MET500")
dir.metabric <- file.path(dir.data,"metabric")
dir.depmap <- file.path(dir.data,"depmap")

for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# TNBC subtype inferred using TNBCtype https://cbc.app.vumc.org/tnbc/
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Retrieve TNBCType results
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# TNBC samples were identified by evaluating distribution of ER, PGR and HER2 
# using RNA, protein (RPPA/mass spec) and DNA copy number annotated with clinical 
# assessment (IHC and FISH) provided by TCGA and CPTAC (Supplementary Fig. S1). 
# Load supplemental table 2
S2 <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS") 

tnbc.patient <- S2$patient

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Retrieve TCGA data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#-----------------------------------
# 1.1 - Clinical / pathologic data
#----------------------------------
query.clin <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Clinical", 
  file.type = "xml",
  barcode = tnbc.patient
)

GDCdownload(query.clin, directory = dir.tcga)

BRCA.clin <- GDCprepare_clinic(
  query = query.clin, 
  clinical.info = "patient",
  directory = dir.tcga
)

save(
  BRCA.clin, 
  file = file.path(dir.tcga,"TCGA_TNBC_192samples_Clinical.Rdata")
)

#-----------------------------------
# 1.2 - RNASeq Gene Expression
#----------------------------------
query.rna <- GDCquery(
  project = "TCGA-BRCA", 
  legacy = TRUE,
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type = "results",
  barcode = tnbc.patient,
  sample.type = c("Primary Tumor"),
  experimental.strategy = "RNA-Seq"
)

GDCdownload(query.rna, directory = dir.tcga)
BRCA.exp <- GDCprepare(query = query.rna, directory = dir.tcga)

dataFilt <- BRCA.exp %>% 
  TCGAanalyze_Preprocessing(cor.cut = 0.6) %>%
  TCGAanalyze_Normalization(
    geneInfo = geneInfo,
    method = "gcContent"
  )      %>%  TCGAanalyze_Filtering(
    method = "quantile", 
    qnt.cut =  0.25
  )   

save(
  BRCA.exp, 
  dataFilt,
  file = file.path(dir.tcga,"TCGA_TNBC_192samples_RNASeq.Rdata")
)

#-----------------------------------
# 1.3 - Copy number variation (CNV)
#----------------------------------
query.cnv <- GDCquery(
  project = "TCGA-BRCA",
  data.category =  "Copy number variation",
  legacy = TRUE,
  platform = "Affymetrix SNP Array 6.0",
  data.type = "Copy number segmentation",
  file.type = "nocnv_hg19.seg",
  sample.type = c("Primary Tumor"),
  barcode = tnbc.patient
)

GDCdownload(query.cnv, directory = dir.tcga)

BRCA.cnv <- GDCprepare(query.cnv, directory = dir.tcga) %>% as.data.frame

save(
  BRCA.cnv, 
  file = file.path(dir.tcga,"TCGA_TNBC_192samples_CNV.Rdata")
)

#-----------------------------------
# 1.4 - DNA methylation
#-----------------------------------
query.met <- GDCquery(
  project = "TCGA-BRCA", 
  legacy = TRUE,
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450",
  sample.type = c("Primary Tumor"),
  barcode = tnbc.patient
)

GDCdownload(query.met, directory = dir.tcga)

BRCA.met <- GDCprepare(query.met, directory = dir.tcga)

save(
  BRCA.met, 
  file = file.path(dir.tcga, "TCGA_TNBC_192samples_Methylation.Rdata")
)

#-----------------------------------
# 1.5 - Mutation MAF
#----------------------------------
query.mut <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple nucleotide variation", 
  data.type = "Simple somatic mutation",
  access = "open", 
  legacy = TRUE,
  file.type = "genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.1.1.0.curated.somatic.maf"
)

GDCdownload(query.mut, directory = dir.tcga)

BRCA.mut <- GDCprepare(query.mut, directory = dir.tcga)
BRCA.mut <- BRCA.mut[substr(BRCA.mut$Tumor_Sample_Barcode,1,12) %in% tnbc.patient,]
save(BRCA.mut, file = file.path(dir.tcga, "TCGA_TNBC_192samples_Mutation.Rdata"))

# Get also mutation MC3 MAF (753MB)
# https://gdc.cancer.gov/about-data/publications/mc3-2017
# Download https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc
options(timeout = 1000) # set 1000 second to download the file, default is 60 seconds
downloader::download(
  "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc",
  file.path(dir.tcga,"mc3.v0.2.8.PUBLIC.maf.gz"),
  mode = ifelse(Sys.info()["sysname"] == "Windows","wb","w")
)
mc3.maf <- readr::read_tsv(
  file = file.path(dir.tcga,"mc3.v0.2.8.PUBLIC.maf.gz"),
  progress = TRUE, 
  col_types = readr::cols()
) %>% dplyr::mutate(patient = substr(Tumor_Sample_Barcode,1,12))

mc3.maf.tnbc <- mc3.maf %>% dplyr::filter(patient %in% tnbc.patient)

save(
  mc3.maf.tnbc, 
  file =  file.path(dir.tcga,"TCGA_TNBC_192samples_mc3_maf.Rdata")
)

#-----------------------------------
# 1.6 - microRNA Expression
#----------------------------------
query.miR <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Gene expression",
  data.type = "miRNA gene quantification",
  platform = "Illumina HiSeq",
  file.type = "hg19.mirbase20.mirna.quantification",
  legacy = TRUE,
  barcode = tnbc.patient
)

GDCdownload(query.miR, directory = dir.tcga)

BRCA.miR <- GDCprepare(query.miR, directory = dir.tcga)

save(
  BRCA.miR, 
  file =  file.path(dir.tcga,"TCGA_TNBC_192samples_microRNA.Rdata")
)

#-----------------------------------
# 1.7 - RPPA Protein Expression
# ----------------------------------
# We  have technical replications
# Approach 1: take the avg
query.RPPA <- GDCquery(
  project = "TCGA-BRCA", 
  legacy = TRUE,
  data.category = "Protein expression",
  data.type = "Protein expression quantification",
  platform = "MDA_RPPA_Core", 
  file.type = "expression"
)

query.RPPA$results[[1]] <- query.RPPA$results[[1]][substr(query.RPPA$results[[1]]$cases,1,12) %in% tnbc.patient,]  
query.RPPA$results[[1]] <- query.RPPA$results[[1]][substr(query.RPPA$results[[1]]$cases,14,15) == "01",]
GDCdownload(query.RPPA, directory = dir.tcga)

BRCA.RPPA <- GDCprepare(query.RPPA, directory = dir.tcga)

rppa.samples <- substr(colnames(BRCA.RPPA),1,12)[-1] %>% unique
BRCA.RPPA.avg <- plyr::adply(rppa.samples,.margins = 1,.fun = function(s){
  rowMeans(BRCA.RPPA[,grep(s, colnames(BRCA.RPPA)),drop = FALSE],na.rm = TRUE)
},.id = NULL) %>% t
rownames(BRCA.RPPA.avg) <- BRCA.RPPA$`Composite Element REF`
colnames(BRCA.RPPA.avg) <- substr(colnames(BRCA.RPPA),1,15)[-1] %>% unique

# extracted from mdanderson.org_BRCA.MDA_RPPA_Core.mage-tab.1.5.0.tar.gz
# https://portal.gdc.cancer.gov/legacy-archive/files/ec52b4e4-fc5a-46ca-ac41-37e1a9f22811
RPPA_Core_antibody_annotation <- readr::read_delim(
  file.path(dir.tcga,"mdanderson.org_BRCA.MDA_RPPA_Core.antibody_annotation.txt"),
  "\t", 
  escape_double = FALSE, 
  trim_ws = TRUE
)

save(
  BRCA.RPPA.avg,
  RPPA_Core_antibody_annotation, 
  file = file.path(dir.tcga, "TCGA_TNBC_192samples_RPPA.Rdata")
)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# linkedOmics aux function for CPTAC data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
get_Linked_omics_data <- function(link, destfile){
  linkedomics.url <- "http://linkedomics.org/data_download/CPTAC-BRCA/"
  ProteinUrl <- file.path(
    linkedomics.url,
    link
  )
  if (!file.exists(destfile)) {
    download.file(url = ProteinUrl, destfile = destfile)
  }
  
  data <- fread(destfile) %>% as.data.frame
  rownames(data) <- data[,1]
  data <- data[,-1]
  return(data)
}


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Retrieve protein data from linkedOmics
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dataProt <- get_Linked_omics_data(
  "HS_CPTAC_BRCA_2018_Proteome_Ratio_Norm_gene_Median.cct",
  file.path(dir.cptac,"HS_CPTAC_BRCA_2018_Proteome_Ratio_Norm_gene_Median.cct")
)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Retrieve Phospho data from linkedOmics
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dataPhosphoPhosphosite <- get_Linked_omics_data(
  "HS_CPTAC_BRCA_2018_Phosphoproteome_Ratio_Norm_Site.cct",
  file.path(dir.cptac,"HS_CPTAC_BRCA_2018_Phosphoproteome_Ratio_Norm_Site.cct")
)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Retrieve Phospho gene mean data from linkedOmics
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dataPhospho_gene <- get_Linked_omics_data(
  "HS_CPTAC_BRCA_2018_Phosphoproteome_Ratio_Norm_Gene_median.cct",
  file.path(dir.cptac,"HS_CPTAC_BRCA_2018_Phosphoproteome_Ratio_Norm_Gene_median.cct")
)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Retrieve TCGA RNA-seq gene mean data from linkedOmics
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dataRNA <- get_Linked_omics_data(
  "HS_CPTAC_BRCA_2018_RNA_GENE.cct",
  file.path(dir.cptac, "HS_CPTAC_BRCA_2018_RNA_GENE.cct")
)

cptac_metadata <- readr::read_tsv("http://www.linkedomics.org/data_download/CPTAC-BRCA/HS_CPTAC_BRCA_2018_CLI.tsi")
cptac_TNBCSubtype <- readxl::read_excel(file.path(dir.cptac,"CPTAC prospective revised.xlsx"))
colnames(cptac_TNBCSubtype)[1] <- "patient"
cptac_TNBCSubtype <- cptac_TNBCSubtype %>% dplyr::filter(!is.na(patient))
cptac_TNBCSubtype$TNBC_subtype <- cptac_TNBCSubtype$`TNBC sutype`

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Adding missing samples as NAs
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
cptac_TNBCSubtype <- dplyr::left_join(
  cptac_TNBCSubtype,
  data.frame(
    patient = setdiff(colnames(dataProt),cptac_TNBCSubtype$patient)
  )
) %>% as.data.frame()
rownames(cptac_TNBCSubtype) <- cptac_TNBCSubtype$patient

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Selection of cptac BRCA samples 
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
cptac_BRCA_samples <- colnames(dataProt)
cptac_TNBCSubtype <- cptac_TNBCSubtype[cptac_BRCA_samples,]
cptac_BRCA_RNA <- dataRNA[,cptac_BRCA_samples] 
cptac_BRCA_Prot <- dataProt[,cptac_BRCA_samples] 
cptac_BRCA_Phospho <- dataPhosphoPhosphosite[,cptac_BRCA_samples]

save(
  cptac_BRCA_RNA,
  cptac_BRCA_Prot,
  cptac_BRCA_Phospho,
  cptac_TNBCSubtype,
  file =  file.path(dir.cptac,"cptac_BRCA_prospective_linkedomics_122_allsamples.Rdata")
)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Phospo median or mean
# Prepare phospo sites: since multiple stes maps to one gene, use either mean or median
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
library(dplyr)
cptac_BRCA_Phospho <- cptac_BRCA_Phospho %>% mutate_if(is.character, ~ as.numeric(.x)) %>% data.matrix()
cptac_BRCA_Phospho <- cptac_BRCA_Phospho[,match(colnames(cptac_BRCA_RNA),colnames(cptac_BRCA_Phospho))]
rownames(cptac_BRCA_Phospho) <- rownames(dataPhosphoPhosphosite)

aux <- sapply(rownames(cptac_BRCA_Phospho), function(x){unlist(strsplit(x,"x"))[[1]]}) %>% unique()
cptac_BRCA_Phospho.mean <- plyr::adply(aux,.margins = 1,.fun = function(x){
  colMeans(cptac_BRCA_Phospho[grep(x,rownames(cptac_BRCA_Phospho)),,drop = FALSE],na.rm = TRUE)
},.progress = "time",.id = NULL,.parallel = FALSE)

cptac_BRCA_Phospho.median <- plyr::adply(aux, .margins = 1, .fun = function(x){
  matrixStats::colMedians(cptac_BRCA_Phospho[grep(x,rownames(cptac_BRCA_Phospho)),,drop = FALSE],na.rm = TRUE)
}, .progress = "time", .id = NULL, .parallel = FALSE)

colnames(cptac_BRCA_Phospho.median) <- colnames(cptac_BRCA_Phospho)
colnames(cptac_BRCA_Phospho.mean) <- colnames(cptac_BRCA_Phospho)
rownames(cptac_BRCA_Phospho.mean) <- aux 
rownames(cptac_BRCA_Phospho.median) <- aux


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# selection of all TNBC samples 27 samples
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
cptac_TNBCSubtype <- readxl::read_excel(file.path(dir.cptac,"CPTAC prospective revised.xlsx")) %>% as.data.frame
colnames(cptac_TNBCSubtype)[1] <- "patient"
cptac_TNBCSubtype$TNBC_subtype <- cptac_TNBCSubtype$`TNBC sutype`
cptac_TNBCSubtype <- cptac_TNBCSubtype %>% dplyr::filter(!is.na(patient))
TNBC_samples_27 <- cptac_TNBCSubtype[cptac_TNBCSubtype$TNBC_subtype %in% c("BL1","BL2","LAR","M"),]
rownames(TNBC_samples_27) <- TNBC_samples_27$patient
TNBC_samples_27$patient[-c(1:2)] <- paste0("X",TNBC_samples_27$patient[-c(1:2)])

cptac_TNBC_RNA <- cptac_BRCA_RNA %>% dplyr::select(TNBC_samples_27$patient)
cptac_TNBC_Prot <- cptac_BRCA_Prot %>% dplyr::select(TNBC_samples_27$patient)
cptac_TNBC_Phospho <- cptac_BRCA_Phospho %>% as.data.frame() %>% dplyr::select(TNBC_samples_27$patient)
cptac_TNBC_Phospho_median <- cptac_BRCA_Phospho.median %>% dplyr::select(TNBC_samples_27$patient)

# Add 0 to phospo if it is missing
missinproteins <- setdiff(rownames(cptac_TNBC_Prot),rownames(cptac_BRCA_Phospho.median))
dataPhosho_missing <- matrix(0,length(missinproteins),ncol(cptac_BRCA_Phospho.median)) %>% as.data.frame
colnames(dataPhosho_missing) <- colnames(cptac_BRCA_Phospho.median)
rownames(dataPhosho_missing) <- missinproteins
cptac_BRCA_Phospho.median_full <- rbind(cptac_BRCA_Phospho.median,dataPhosho_missing)

cptac_TNBC_Phospho_median_full <- cptac_BRCA_Phospho.median_full  %>% dplyr::select(TNBC_samples_27$patient)

library(xCell)
cptac_TNBC_xCell_RNA <- as.data.frame(xCellAnalysis(expr = cptac_TNBC_RNA))
cptac_TNBC_xCell_Prot <- as.data.frame(xCellAnalysis(expr = cptac_TNBC_Prot))
cptac_TNBC_xCell_Phospho <- as.data.frame(xCellAnalysis(expr = cptac_TNBC_Phospho_median_full))

save(
  cptac_TNBC_RNA,
  cptac_TNBC_Prot,
  cptac_TNBC_Phospho,
  cptac_TNBC_Phospho_median,
  cptac_TNBC_xCell_RNA,
  cptac_TNBC_xCell_Prot,
  cptac_TNBC_xCell_Phospho,
  TNBC_samples_27,
  file =  file.path(dir.cptac,"cptac_TNBC_linkedomics_27samples.Rdata")
)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Depmap data
# DepMap Public 18Q4 CCLE_depMap_18Q4_TPM_v2.csv
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
depmap <- readr::read_csv("https://depmap.org/portal/download/api/download?file_name=ccle%2Fdepmap-rnaseq-expression-data-ccd0.9%2FCCLE_depMap_18Q4_TPM_v2.csv&bucket=depmap-external-downloads")
save(depmap, file = file.path(dir.depmap,"depmap.rda"))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# PDTX data
# https://caldaslab.cruk.cam.ac.uk/bcape
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
PDTX <- readr::read_tsv("https://ndownloader.figshare.com/files/5968401")
save(PDTX, file = file.path(dir.pdtx,"PDTX.rda"))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# MET500 RNA-seq data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Source https://xenabrowser.net/datapages/?cohort=MET500%20(expression%20centric)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
sourceFile <- "https://ucsc-public-main-xena-hub.s3.us-east-1.amazonaws.com/download/MET500%2FgeneExpression%2FM.mx.txt.gz"
MET500 <- readr::read_tsv(sourceFile)
save(MET500, file = file.path(dir.met500,"MET500.rda"))

MET500.tnbc.meta <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "D-MET500"
) 

# change to genes.symbols
MET500$sample <- getGenes(gsub("\\.[0-9]*","",MET500$sample),fields = "symbol")$symbol
MET500 <- MET500[!is.na(MET500$sample),]

MET500.tnbc <- MET500[,-1] %>% as.matrix()
rownames(MET500.tnbc) <- MET500$sample
colnames(MET500.tnbc) <- gsub("-",".",colnames(MET500.tnbc))
MET500.tnbc <- MET500.tnbc[,colnames(MET500.tnbc)  %in% MET500.tnbc.meta$sample_detailed]
MET500.tnbc.xcell <- xCell::xCellAnalysis(MET500.tnbc)
MET500.tnbc.meta <- MET500.tnbc.meta[match(colnames(MET500.tnbc),MET500.tnbc.meta$sample_detailed),]

save(
  MET500.tnbc.meta,
  MET500.tnbc,
  MET500.tnbc.xcell,
  file = file.path(dir.met500,"MET500_TNBC.rda")
)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Metabric RNA-seq data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Controlled data access - data not included in the code/Github 
# https://www.synapse.org/#!Synapse:syn1757063
# Save file in data/metabric folder
file <- file.path(dir.metabric,"Complete_normalized_expression_data_METABRIC.txt")
file <- "~/TBL Dropbox/Tiago Silva/TNBC/metabric/data/Complete_normalized_expression_data_METABRIC.txt"

options(timeout = 1000) 
downloader::download(
  "https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz",
  file.path(dir.metabric,"brca_metabric.tar.gz"),
  mode = ifelse(Sys.info()["sysname"] == "Windows","wb","w")
)
utils::untar(file.path(dir.metabric,"brca_metabric.tar.gz"),exdir = dir.metabric)

data <- readr::read_tsv(file.path(dir.metabric,"brca_metabric/data_expression_median.txt"))
colnames(data) <- gsub("-","_",colnames(data))
data <- data[!duplicated(data$Hugo_Symbol),]

Metabric.tnbc.meta <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "B-METABRIC"
) %>% dplyr::filter(sample %in% colnames(data))

Metabric.tnbc <- data[,Metabric.tnbc.meta$sample] %>% as.data.frame()
rownames(Metabric.tnbc) <- data$Hugo_Symbol

Metabric.tnbc.xcell <- xCell::xCellAnalysis(Metabric.tnbc)
Metabric.tnbc.meta <- Metabric.tnbc.meta[match(colnames(Metabric.tnbc),Metabric.tnbc.meta$sample),]

save(
  Metabric.tnbc.meta,
  Metabric.tnbc,
  Metabric.tnbc.xcell,
  file = file.path(dir.metabric,"metabric_TNBC.rda")
)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# TCGA BRCA mbatch data used for consensus clustering
# MBatch Omic Browser 
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
options(timeout = 1000) # set 1000 second to download the file, default is 60 seconds
# downloaded 143.3 MB
downloader::download(
  "https://bioinformatics.mdanderson.org/MQA/dszipdata?id=bev-gdc-84810dfecbb341d9603677e906c528df",
  file.path(dir.tcga,"mbatch_TNBC_legacy_RNA-Seq-gene-normalized-v2_matrix_data.zip"),
  mode = ifelse(Sys.info()["sysname"] == "Windows","wb","w")
)
unzip(
  file.path(dir.tcga,"mbatch_TNBC_legacy_RNA-Seq-gene-normalized-v2_matrix_data.zip"),
  files = "original/matrix_data.tsv",
  exdir = dir.tcga
)
mbatch.hg19 <- readr::read_tsv(
  file = file.path(dir.tcga,"original/matrix_data.tsv"),
  progress = TRUE, 
  col_types = readr::cols()
)  %>% as.data.frame()
mbatch.tnbc.hg19 <- mbatch.hg19 %>% 
  dplyr::select(
    c(
      1,
      which(substr(colnames(mbatch.hg19),1,15) %in% paste0(tnbc.patient,"-01")) # get only TNBC primary solid tumors
    )
  )
mbatch.tnbc.hg19$Gene_symbol <- TCGAbiolinks:::GenesCutID(mbatch.tnbc.hg19$X1)
mbatch.tnbc.hg19 <- mbatch.tnbc.hg19[mbatch.tnbc.hg19$Gene_symbol != "?",]
mbatch.tnbc.hg19 <- mbatch.tnbc.hg19[!duplicated(mbatch.tnbc.hg19$Gene_symbol),]
rownames(mbatch.tnbc.hg19) <- mbatch.tnbc.hg19$Gene_symbol
mbatch.tnbc.hg19 <- mbatch.tnbc.hg19 %>% dplyr::mutate(X1 = NULL,Gene_symbol = NULL)
mbatch.tnbc.hg19 <- mbatch.tnbc.hg19 %>% as.matrix
save(
  mbatch.tnbc.hg19,
  file = file.path(dir.tcga,"mbatch_TNBC_legacy_RNA-Seq-gene-normalized-v2_matrix_data.rda")
)

mbatch.brca.hg19.normalized <- mbatch.tnbc.hg19 %>% 
  TCGAanalyze_Normalization(
    geneInfo = geneInfo,
    method = "gcContent"
  ) %>%   TCGAanalyze_Filtering(
    method = "quantile", 
    qnt.cut =  0.25
  )   
save(
  mbatch.brca.hg19.normalized,
  file = file.path(dir.tcga,"mbatch_TNBC_legacy_RNA-Seq-gene-normalized-v2_matrix_data.rda")
)



#============================================================
# TCGA ATAC-seq data
#============================================================
# 1) Get ATAC-seq big wig
# query <- GDCquery_ATAC_seq(tumor = "BRCA",file.type = "bigWigs")
# GDCdownload(query,method = "client")

# Big wigs
url <- "https://api.gdc.cancer.gov/data/f1c06cd3-cf35-41cc-bc75-6db273c94273"
file <- file.path(dir.atac_seq,"BRCA_bigWigs.tgz")
dir.create(dirname(file),recursive = TRUE)


options(timeout = 10000) # set 10000 second to download the file, default is 60 seconds
# Data size 15GB
downloader::download(url,file, mode = ifelse(Sys.info()["sysname"] == "Windows","wb","w"))
untar(tarfile = file,exdir = dirname(file))

# Data
url <- "https://api.gdc.cancer.gov/data/38b8f311-f3a4-4746-9829-b8e3edb9c157"
file <- file.path(dir.atac_seq,"TCGA-ATAC_Cancer_Type-specific_Count_Matrices_log2norm_counts.zip")
dir.create(dirname(file),recursive = TRUE)
# Data size 0.6GB
downloader::download(url,file, mode = ifelse(Sys.info()["sysname"] == "Windows","wb","w"))
unzip(file, exdir = dirname(file))

# 1) Get ATAC-seq normalized count matrix
atac.seq.brca <- readr::read_tsv(file.path(dirname(file),"BRCA_log2norm.txt"))

gdc.file <- "https://api.gdc.cancer.gov/data/7a3d7067-09d6-4acf-82c8-a1a81febf72c"
samples.ids <- readr::read_tsv(gdc.file, col_types = readr::cols())
samples.ids$sample <- substr(samples.ids$Case_ID,1,16)
head(samples.ids)

colnames(atac.seq.brca)[-c(1:5)] <- samples.ids$Case_ID[match(gsub("_","-",colnames(atac.seq.brca)[-c(1:5)]),samples.ids$bam_prefix)]
atac.seq.brca[1:4,1:8]

idx.tnbc <- which(substr(colnames(atac.seq.brca)[-c(1:5)],1,12) %in% S2$barcode) + 5
atac.seq.tnbc <- atac.seq.brca[,c(1:5,idx.tnbc)]

#-------------------------------
# SummarizedExperiment object
#-------------------------------
samples.info <- TCGAbiolinks:::colDataPrepare(unique(colnames(atac.seq.tnbc)[-c(1:5)]))
samples.info$TNBC_subtype <- S2$subtype[match(samples.info$patient,S2$patient)]

# create SE object  
counts <- atac.seq.tnbc[,-c(1:5)]
rowRanges <- makeGRangesFromDataFrame(atac.seq.tnbc)
rowRanges$score <- atac.seq.tnbc$score
rowRanges$name <- atac.seq.tnbc$name

rowRanges$id_peak <-  paste0(
  "peak_id_",
  as.character(seqnames(rowRanges)),
  ":",
  start(rowRanges) - 1,
  "-",
  end(rowRanges)
)

names(rowRanges) <- paste(
  atac.seq.tnbc$name,atac.seq.tnbc$seqnames,
  atac.seq.tnbc$start,atac.seq.tnbc$end, 
  sep = "_"
)

colData <- DataFrame(unique(left_join(samples.info,samples.ids)))

rse <- SummarizedExperiment(
  assays = SimpleList(log2norm = as.matrix(counts)),
  rowRanges = rowRanges, 
  colData = colData
)
save(rse,file = file.path(dir.atac_seq,"TNBC_ATAC_seq_non_merged.rda"), compress = "xz")

# This function will calculate the Means of the peaks for a given group
# in our case we will calculate the mean of the replicates of each patient.
library(dplyr)
groupMeans <- function(mat, groups = NULL, na.rm = TRUE){
  stopifnot(!is.null(groups))
  gm <- lapply(unique(groups), function(x){
    rowMeans(mat[,which(groups == x),drop = F], na.rm=na.rm)
  }) %>% Reduce("cbind",.)
  colnames(gm) <- unique(groups)
  return(gm)
}
matMerged <- groupMeans(
  mat = atac.seq.tnbc[,-c(1:5)], 
  groups = substr(colnames(atac.seq.tnbc)[-c(1:5)],1,12)
)

rse.merged <- SummarizedExperiment(
  assays = SimpleList(log2norm = as.matrix(matMerged)),
  rowRanges = rowRanges, 
  colData = samples.info
)
save(rse.merged,file = file.path(dir.atac_seq,"TNBC_ATAC_seq_merged.rda"), compress = "xz")
#============================================================

