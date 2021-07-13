#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# RNASeqv2 DEG Analaysis with 172 samples
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
library(limma)
library(edgeR)
library(DESeq2)
library(plyr)
library(sva)
library(dplyr)
library(stats)
require(SummarizedExperiment)
options(scipen = 500)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Folders where we was saved the data
# See script: 01_data_download.R
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data,"TCGA")
dir.output <-  "analysis_results/TCGA/rnaseq_results"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Get data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
load(file.path(dir.tcga,"TCGA_TNBC_192samples_RNASeq.Rdata"))

exp <- assay(BRCA.exp,"raw_count") %>% round
rownames(exp) <- rowRanges(BRCA.exp)$gene_id
colnames(exp) <- substr(colnames(exp), 1,12)


# Load supplemental table 2
metadata <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS") %>% dplyr::filter(!is.na(Stage))

exp <- exp[, metadata$patient]

cpm <- cpm(exp)
lcpm <- cpm(exp, log = TRUE)
keep.exprs <- rowSums(cpm > 1) >= 27
exp <- exp[keep.exprs,]
dim(exp)

voom.exp = voom(exp)$E

## DE test with limma-voom updated
DE_limmaVoom_mult <- function(
  Data,
  sample_info,
  grp_by,
  mult_factor,
  output_dir = NULL,
  prefix = NULL
) {
  
  stopifnot(all(colnames(Data) == sample_info$patient))

  cat(paste0("~ ",grp_by, " + Stage + Age"), "\n")
  
  # preparing the data
  sample_info$Stage = factor(sample_info$Stage)
  Data <- as.matrix(Data)
  Data <- Data[rowSums(Data) != 0,] 
  
  data.des = model.matrix(formula(paste0("~", mult_factor)), data = sample_info)
  data.des0 = model.matrix( ~Age+Stage,data = sample_info)
  
  batch_unsup_sva = svaseq(Data,data.des,data.des0)$sv
  batch = batch_unsup_sva
  adj = "+batch"
  data.des = model.matrix(formula(paste0("~", mult_factor,adj)), data = sample_info)
  
  data.dge = DGEList(counts= Data)
  data.dge = calcNormFactors(data.dge)
  data.v = voom(data.dge, data.des)
  
  data.vlmfit = lmFit(data.v, data.des)
  
  Fit.eb <- eBayes(data.vlmfit)
  
  res.collect = topTable(Fit.eb, coef = grp_by, adjust.method = "BH", n = Inf, sort.by = "P")
  Temp_id = match(res.collect$ID, rownames(Data))
  res.collect = cbind(res.collect, Data[Temp_id, ])
  
  if (!is.null(output_dir)) {
    write.csv(res.collect, paste0(output_dir, "/", prefix, "_", "centroid_testing", ".csv"),row.names=T, quote=F)
  }
  res.collect
}


res.BL1 = DE_limmaVoom_mult(
  Data = exp,
  sample_info = metadata,
  grp_by = "BL1",
  mult_factor = "BL1+Stage+Age",
  output_dir = dir.output,
  prefix = "BL1"
)

res.BL2 = DE_limmaVoom_mult(
  Data = exp,
  sample_info = metadata,
  grp_by = "BL2",
  mult_factor = "BL2+Stage+Age",
  output_dir = dir.output,
  prefix = "BL2"
)

res.LAR = DE_limmaVoom_mult(
  Data = exp,
  sample_info = metadata,
  grp_by = "LAR",
  mult_factor = "LAR+Stage+Age",
  output_dir = dir.output,
  prefix = "LAR"
)

res.M = DE_limmaVoom_mult( 
  Data = exp,
  sample_info = metadata,
  grp_by = "M",
  mult_factor = "M+Stage+Age",
  output_dir = dir.output,
  prefix = "M"
)
