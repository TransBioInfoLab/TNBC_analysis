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
library(dplyr)
library(readr)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/depmap//"
dir.plots <- "plots/depmap/"
dir.analysis <- "analysis/depmap/"
dir.output <- "analysis_results/depmap/"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
RNAi <- read_csv(file.path(dir.data,"depmap_RNAi_TNBC.csv")) %>% as.data.frame()
sgRNA <- read.csv(file.path(dir.data,"sgRNA_gene_effect_TNBC.csv"), row.names = 1)
centroids <- read.csv(file.path(dir.data,"TNBCtype_screen_V2.csv"), row.names = 1)
rownames(RNAi) <- gsub("\\([[:alnum:]]*\\)","",RNAi$Gene) %>% stringr:::str_trim()
RNAi$GENE <- RNAi$Gene <- NULL

all(colnames(sgRNA) %in% colnames(RNAi))
all(colnames(RNAi) %in% rownames(centroids))
all(colnames(sgRNA) %in% rownames(centroids))

centroids <- centroids[rownames(centroids) %in% colnames(RNAi),]
RNAi <- RNAi[,colnames(RNAi) %in% rownames(centroids),]
centroids <- centroids[match(colnames(RNAi),rownames(centroids)),]

# Aux function
get_tval_pval_lm <- function(y,x, prefix){
  temp3 <- summary(lm(as.numeric(y) ~ x))
  df <- data.frame("tvalue" = temp3$coef[2,3],"pvalue" = temp3$coef[2,4])
  colnames(df) <- paste0(prefix,"_",colnames(df))
  return(df)
}
doParallel::registerDoParallel(cores = 4)

## RNAi 
centroids1 <- centroids[colnames(RNAi),]
all(colnames(RNAi) %in% rownames(centroids1))
RNAi.output <- plyr::adply(.data = RNAi,.margins = 1,.fun = function(row){
  cbind(
    get_tval_pval_lm(y = row, x = centroids1$BL1, prefix = "BL1"),
    get_tval_pval_lm(y = row, x = centroids1$BL2, prefix = "BL2"),
    get_tval_pval_lm(y = row, x = centroids1$M, prefix = "M"),
    get_tval_pval_lm(y = row, x = centroids1$LAR, prefix = "LAR")
  )
},.expand = F,.id = NULL, .progress = "time",.parallel = FALSE)
rownames(RNAi.output) <- rownames(RNAi)
write.csv(RNAi.output, file = file.path(dir.output,"depmap_RNAi_TNBC_result.csv"))

## sgRNA
centroids2 <- centroids[colnames(sgRNA),]
all(colnames(sgRNA) %in% rownames(centroids2))
sgRNA.output <- plyr::adply(.data = sgRNA,.margins = 1,.fun = function(row){
  cbind(
    get_tval_pval_lm(y = row, x = centroids2$BL1, prefix = "BL1"),
    get_tval_pval_lm(y = row, x = centroids2$BL2, prefix = "BL2"),
    get_tval_pval_lm(y = row, x = centroids2$M, prefix = "M"),
    get_tval_pval_lm(y = row, x = centroids2$LAR, prefix = "LAR")
  )
},.expand = F,.id = NULL, .progress = "time",.parallel = FALSE)
rownames(sgRNA.output) <- rownames(sgRNA)
write.csv(sgRNA.output, file = file.path(dir.output,"depmap_sgRNA_TNBC_result.csv"))
