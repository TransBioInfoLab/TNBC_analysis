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
# Methy450_NonHispanic_V2.R
# https://www.programiz.com/r-programming/switch-function R programming
# https://www.programiz.com/ C,C++,Java,R,Python,Koltlin programming languages
# https://www.programiz.com/dsa data structures and algorithsm for C,C++,Java,R,Python,Koltlin programming languages 
library(Hmisc)
library(qvalue)
library(limma)
library(edgeR)
library(hexbin)
library(sva)
library(dplyr)
require(SummarizedExperiment)
options(scipen = 500)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Folders where we was saved the data
# See script: 01_data_download.R
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data,"TCGA")
dir.output <-  "analysis_results/TCGA/dnam_results"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Get data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
load(file.path(dir.tcga,"TCGA_TNBC_192samples_Methylation.Rdata"))

############ 1.data preparation--setup work path,data preprocessing,M_value
dnam <- assay(BRCA.met)
colnames(dnam) <- substr(colnames(dnam), 1,12)
# select !ch. rows
BMIQ.nonrs <- dnam[!grepl("rs.", rownames(dnam)),] %>% na.omit
####### M-value for compuation
BMIQ <- log2(BMIQ.nonrs/(1 - BMIQ.nonrs))

# Load supplemental table 2
metadata <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS") %>% dplyr::filter(!is.na(Stage))

BMIQ <- BMIQ[,colnames(BMIQ) %in% metadata$patient]
metadata <- metadata[match(colnames(BMIQ),metadata$patient),]

########## linear model: M ~ age + stage + Subtype, SVA feature selection, DMA--SVA--info1
DM_limma <- function(
  Data, 
  sample_info, 
  grp_by,
  mult_factor,
  output_dir = NULL, 
  prefix = NULL
){
  
  stopifnot(all(colnames(Data) == sample_info$patient))
  
  cat(paste0("~ 0 + ",grp_by," + Stage + Age"), "\n")
  
  # preparing the data
  sample_info$Stage = factor(sample_info$Stage)
  mod <- model.matrix(formula(paste0("~ 0 + ", mult_factor)), data = sample_info)
  mod0 <- model.matrix(~ Age + Stage, data = sample_info)
  
  Data <- Data %>% na.omit %>% as.matrix
  ###success r
  n.sv <- num.sv(Data,mod,method = "be", B = 20, seed = set.seed(5000))
  n.sv
  ####
  ## success r
  svobj = sva(Data, mod, mod0, n.sv = n.sv)
  
  df <- svobj$sv %>% as.data.frame
  colnames(df) <- paste0("sv", 1:n.sv)
  modSv = cbind(mod,df)
  mod0Sv = cbind(mod0,df)
  
  data.vlmfit = lmFit(Data, modSv)
  Fit.eb <- eBayes(data.vlmfit)
  res.collect = topTable(Fit.eb, coef = 1, adjust.method = "BH", n = Inf, sort.by = "P")
  Temp_id = match(rownames(res.collect), rownames(Data))
  
  if (!is.null(output_dir)) {
    write.csv(res.collect, paste0(output_dir, "/", prefix, "_", "centroid_methy_testing", ".csv"), row.names = TRUE, quote = FALSE)
  }
  res.collect
}

res.BL1 <- DM_limma(
    Data = BMIQ,
    sample_info = metadata,
    grp_by = "BL1",
    mult_factor = "BL1 + Stage + Age",
    output_dir = dir.output,
    prefix = "BL1"
  )

res.BL2 <- DM_limma(
    Data = BMIQ,
    sample_info = metadata,
    grp_by = "BL2",
    mult_factor = "BL2 + Stage + Age",
    output_dir = dir.output,
    prefix = "BL2"
  )

res.LAR <- DM_limma(
    Data = BMIQ,
    sample_info = metadata,
    grp_by = "LAR",
    mult_factor = "LAR + Stage + Age",
    output_dir = dir.output,
    prefix = "LAR"
  )

res.M <- DM_limma(
    Data = BMIQ,
    sample_info = metadata,
    grp_by = "M",
    mult_factor = "M + Stage + Age",
    output_dir = dir.output,
    prefix = "M"
  )
