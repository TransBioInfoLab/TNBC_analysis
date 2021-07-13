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
## Model 
library(gtools)
library(doParallel)
library(plyr)
library(dplyr)
library(readr)
doParallel::registerDoParallel(cores = 4)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/PDTX/"
dir.plots <- "plots/PDTX"
dir.output <- "analysis_results/PDTX"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
## Models
model_IC50 <- readr::read_csv(file.path(dir.data,"DrugResponsesAUCModels.csv"))
centroids <- read.csv(file.path(dir.data,"PDTX_models_result.csv"), row.names = 1)

model_IC50 <- model_IC50 %>% 
  dplyr::filter(Model %in% rownames(centroids)) %>% 
  mutate(AUC = replace(AUC, AUC == 0, 1e-05))

tt1 <- unique(model_IC50$Drug)

# Aux function
get_tval_pval_lm <- function(y, x, prefix){
  temp <- summary(lm(logit(y) ~ x))
  df <- data.frame("tvalue" = temp$coef[2,3],"pvalue" = temp$coef[2,4])
  colnames(df) <- paste0(prefix,"_",colnames(df))
  return(df)
}

## RNAi 
Models_AUC <- plyr::adply(.data = tt1,.margins = 1,.fun = function(drug){
  temp1 <- model_IC50[model_IC50$Drug == drug, ]
  temp2 <- centroids[temp1$Model,]
  cbind(
    get_tval_pval_lm(y = temp1$AUC, x = temp2$BL1, prefix = "BL1"),
    get_tval_pval_lm(y = temp1$AUC, x = temp2$BL2, prefix = "BL2"),
    get_tval_pval_lm(y = temp1$AUC, x = temp2$M, prefix = "M"),
    get_tval_pval_lm(y = temp1$AUC, x = temp2$LAR, prefix = "LAR")
  )
},.expand = F,.id = "DRUG_NAME", .progress = "time",.parallel = FALSE)
Models_AUC$DRUG_NAME <- tt1
write.csv(Models_AUC, file = file.path(dir.output,"Models_AUC.csv"))


## sample
key <- read.csv(file.path(dir.data,"drug_sample_key.csv"))
sample_IC50 <- readr::read_csv(file.path(dir.data,"DrugResponsesAUCSamples.csv"))
centroids <- read.csv(file.path(dir.data,"PDTX_samples_result.csv"), row.names = 1)
centroids <- centroids[gsub("-",".",key$EXP_ID),]
rownames(centroids) <- as.character(key$DRUG_ID)

sample_IC50 <- sample_IC50 %>% 
  dplyr::filter(ID %in% rownames(centroids)) %>% 
  mutate(AUC = replace(AUC, AUC == 0, 1e-05))

tt1 <- sample_IC50$Drug %>% unique

samples_AUC <- plyr::adply(.data = tt1,.margins = 1,.fun = function(drug){
  temp1 <- sample_IC50 %>% dplyr::filter(.data$Drug == drug)
  temp2 <- centroids[as.character(temp1$ID),]
  cbind(
    get_tval_pval_lm(y = temp1$AUC, x = temp2$BL1, prefix = "BL1"),
    get_tval_pval_lm(y = temp1$AUC, x = temp2$BL2, prefix = "BL2"),
    get_tval_pval_lm(y = temp1$AUC, x = temp2$M, prefix = "M"),
    get_tval_pval_lm(y = temp1$AUC, x = temp2$LAR, prefix = "LAR")
  )
},.expand = F,.id = "DRUG_NAME", .progress = "time",.parallel = FALSE)
samples_AUC$DRUG_NAME <- tt1
write.csv(samples_AUC, file = file.path(dir.output,"Samples_AUC.csv"))