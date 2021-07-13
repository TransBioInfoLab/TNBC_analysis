#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Mutation burden Analysis
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
# Date: 28 June 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
library(beeswarm)
library(dplyr)
library(cowplot)
library(ggplot2)
library(gridGraphics)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots/TCGA"
dir.tcga <- file.path(dir.data,"TCGA")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

metadata <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% 
  dplyr::filter(subtype != "UNS") %>% 
  dplyr::mutate_at(.funs =  as.numeric,.vars = vars(dplyr::starts_with("Mut_Sig"))) 


beeswarm(
  Mut_Sig_1_APOBEC ~ subtype, data = metadata, 
  pch = 16,
  col = c("red","orange","green","blue"),
  main = 'beeswarm + bxplot')
bxplot(Mut_Sig_1_APOBEC ~ subtype, data = metadata,  add = TRUE)
p1 <- recordPlot()
plot.new() ## clean up device


beeswarm(
  Mut_Sig_2_DEAMiN ~ subtype, data = metadata, 
  pch = 16,
  col = c("red","orange","green","blue"),
  main = 'beeswarm + bxplot')
bxplot(Mut_Sig_2_DEAMiN ~ subtype, data = metadata,  add = TRUE)
p2 <- recordPlot()
plot.new() ## clean up device

beeswarm(
  Mut_Sig_3_MMR ~ subtype, data = metadata, 
  pch = 16,
  col = c("red","orange","green","blue"),
  main = 'beeswarm + bxplot')
bxplot(Mut_Sig_3_MMR ~ subtype, data = metadata,  add = TRUE)
p3 <- recordPlot()
plot.new() ## clean up device


beeswarm(
  Mut_Sig_4_DSB ~ subtype, data = metadata, 
  pch = 16,
  col = c("red","orange","green","blue"),
  main = 'beeswarm + bxplot')
bxplot(Mut_Sig_4_DSB ~ subtype, data = metadata,  add = TRUE)
p4 <- recordPlot()


p <- cowplot::plot_grid(p1,p2,p3,p4,nrow = 2,ncol = 2,scale = 0.8)
ggsave(filename = file.path(dir.plots,"mutation_sig_bewwswarm_S5E.pdf"),width = 10,height = 10,plot = p)
