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
library(ggpubr)
library(plyr)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
S2.a <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS")  %>% dplyr::select(subtype,ImmuneScore...95) %>%
  dplyr::mutate(Xcell_immune_score = as.numeric(ImmuneScore...95)) %>% 
  dplyr::mutate(dataset = "TCGA")

S2.b <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "B-METABRIC"
) %>% dplyr::filter(!is.na(TNBCtype)) %>% 
  dplyr::filter(TNBCtype != "UNS") %>% dplyr::select(TNBCtype, Xcell_immune) %>% 
  dplyr::mutate(Xcell_immune_score = as.numeric(Xcell_immune)) %>% dplyr::rename(subtype = TNBCtype) %>% 
  dplyr::mutate(dataset = "Metabric")



S2.c <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "C-CPTAC"
)  %>% dplyr::filter(TNBCtype != "UNS")  %>% dplyr::select(TNBCtype, ImmuneScore) %>%
  dplyr::mutate(Xcell_immune_score = as.numeric(ImmuneScore)) %>% dplyr::rename(subtype = TNBCtype)%>% 
  dplyr::mutate(dataset = "CPTAC")


S2.d <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "D-MET500"
)  %>% dplyr::filter(subtype != "UNS")  %>% dplyr::select(subtype,Xcell_immune) %>% 
  dplyr::mutate(Xcell_immune_score = as.numeric(Xcell_immune))  %>% 
  dplyr::mutate(dataset = "MET500")

data <- plyr::rbind.fill(S2.d,S2.c,S2.b,S2.a)

data$dataset <- factor(data$dataset, levels = rev(c("Metabric", "CPTAC", "TCGA", "MET500")))
data$subtype <- factor(data$subtype, levels = c("BL1", "BL2", "M", "LAR","TNBC"))

my_comparisons <- list( c("M", "LAR"),c("M", "BL1"), c("M", "BL2"))

p <- ggboxplot(
  data, 
  x = "subtype", 
  y = "Xcell_immune_score",
  fill = "subtype",
  add = "jitter",
  palette = c(
    "#FF3300",
    "#FFCC00",
    "#3366FF",
    "#66FF33",
    "grey"
  )
) + facet_wrap( ~ dataset, ncol = 2,scales = "free") +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    label = "p.signif",
    vjust = 1.5
  ) + labs(y = "Xcell Immune", x = "") + theme(strip.background = element_rect(fill=alpha("white", 1.0)))

ggsave(
  plot = p,
  filename = file.path(
    dir.plots,
    "Boxplots_immune_score_xcell_by_TNBC_subtype_cohorts.pdf"
  ),
  height = 8,
  width = 7
)
