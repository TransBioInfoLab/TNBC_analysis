#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# The following script will create Fig1C. TCGA heatmap with RNA-seq, CNA, mut,
# Clinical information and others
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
library(ggpubr)
library(ggplot2)
library(dplyr)
library(readr)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Folders where we was saved the data
# See script: 01_data_download.R
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data,"TCGA")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o Data                                                                    |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# Load supplemental table 2
S2 <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
)  %>% dplyr::filter(subtype != "UNS")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
# o plot                                                                    |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
p <- ggplot(
  data = S2,  
  mapping = aes(
    factor(subtype,levels = c("BL1","BL2","M","LAR")), 
    fill = factor(TIME_CLASS,levels = c("FI","SR","MR","ID"))
  )
) + geom_bar(position = "fill") + 
  ylab("Frequency") + 
  xlab("TNBC subtype") + 
  ggthemes::theme_base() + 
  labs(fill="TIME Class") + 
  scale_fill_manual(values= c(
    "FI" = "#FCFFA4FF",
    "ID" = "#420A68FF",
    "MR" = "#AE305CFF",
    "SR" = "#F8850FFF"
  )
  )
ggsave(file.path(dir.plots,"Barplot_Time_class_per_subtype.pdf"),width = 5,height = 5)