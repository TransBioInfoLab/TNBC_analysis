#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# RPPA DE Analysis
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
library(plyr)
library(sva)
library(dplyr)
options(scipen = 500)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Folders where we was saved the data
# See script: 01_data_download.R
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data,"TCGA")
dir.output <-  "analysis_results/TCGA/rppa_results"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Get data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
load(file.path(dir.tcga,"TCGA_TNBC_192samples_RPPA.Rdata"))

exp <- BRCA.RPPA.avg
exp <- exp[!is.na(rowSums(exp)),] # remove proteins with any NA's
colnames(exp) <- substr(colnames(exp), 1,12)

metadata <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS") %>% dplyr::filter(!is.na(Stage))

exp <- exp[, colnames(exp) %in% metadata$patient]
metadata <- metadata[match(colnames(exp), metadata$patient),]

DM_limma <- function(
  Data, 
  sample_info, 
  grp_by,
  mult_factor,
  output_dir = NULL,
  prefix = NULL, 
  output_name = NULL
){
  
  Temp_id = na.omit(match(colnames(Data), sample_info$patient))
  sample_info = sample_info[Temp_id, ]
  Temp_id = na.omit(match(sample_info$patient, colnames(Data)))
  Data = Data[,Temp_id]
  cat(paste0("~ 0 + ",grp_by, " + Stage + Age"), "\n")
  sample_info$Stage = factor(sample_info$Stage)
  mod = model.matrix(formula(paste0("~", mult_factor)), data = sample_info)
  mod0 = model.matrix(~Age+Stage,data = sample_info)
  
  ############################# Data filtering
  Data <- Data %>% na.omit %>% as.matrix()
  
  ############################# success r
  n.sv = num.sv(Data, mod, method = "be", B = 20, seed = set.seed(5000))
  n.sv
  #############################
  ###### success r
  svobj = sva(Data,mod,mod0,n.sv=n.sv)
  
  df <- svobj$sv %>% as.data.frame
  colnames(df) <- paste0("sv", 1:n.sv)
  modSv = cbind(mod,df)
  mod0Sv = cbind(mod0,df)

  data.vlmfit = lmFit(Data, modSv)
  Fit.eb <- eBayes( data.vlmfit)
  res.collect = topTable(Fit.eb, coef = grp_by, adjust.method = "BH", n = Inf, sort.by = "P")
  Temp_id = match(rownames(res.collect), rownames(Data))
  res.collect = cbind(res.collect, Data[Temp_id, ])
  
  if (!is.null(output_dir)) {
    write.csv(
      x = res.collect, 
      file = paste0(output_dir, "/", prefix, "_", "centroid_RPPA_testing.csv"),
      row.names = T, quote = F
    )
  }
  res.collect
}

res.BL1 = DM_limma(
  Data = exp,
  sample_info = metadata,
  grp_by = "BL1",
  mult_factor = "BL1+Stage+Age",
  output_dir = dir.output,
  prefix = "BL1"
)
res.BL2 = DM_limma(
  Data = exp,
  sample_info = metadata,
  grp_by = "BL2",
  mult_factor = "BL2+Stage+Age",
  output_dir = dir.output,
  prefix = "BL2"
)
res.LAR = DM_limma(
  Data = exp,
  sample_info = metadata,
  grp_by = "LAR",
  mult_factor = "LAR+Stage+Age",
  output_dir = dir.output,
  prefix = "LAR"
)
res.M = DM_limma(
  Data = exp,
  sample_info = metadata,
  grp_by = "M",
  mult_factor = "M+Stage+Age",
  output_dir = dir.output,
  prefix = "M"
)
