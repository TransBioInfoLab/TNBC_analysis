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
library(RColorBrewer)
library(forcats)
library(reshape2)
library(ggplot2)

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
) %>% dplyr::filter(subtype != "UNS")

plotDF <- metadata %>% dplyr::select(c(subtype,Mut_Sig_1_APOBEC,Mut_Sig_2_DEAMiN,Mut_Sig_3_MMR,Mut_Sig_4_DSB)) %>% 
  dplyr::filter(Mut_Sig_1_APOBEC != "NA") %>%
  dplyr::mutate_at(.funs =  as.numeric,.vars = vars(dplyr::starts_with("Mut_Sig"))) 

plotDF.ordered <- plotDF %>% 
  dplyr::group_by(subtype)  %>% 
  dplyr::arrange(Mut_Sig_4_DSB,.by_group = T)

plotDF.final <- plotDF.ordered %>% dplyr::ungroup() %>% dplyr::mutate(subtype = NULL,x = 1:nrow(.)) %>% melt(id.vars = "x")
mycolors <- brewer.pal(4, "RdYlBu") %>% rev
names(mycolors) <-  c("Mut_Sig_1_APOBEC", "Mut_Sig_2_DEAMiN", "Mut_Sig_3_MMR", "Mut_Sig_4_DSB")

plot <- ggplot(data = plotDF.final, aes(x = x, y = value, fill = variable)) +
  geom_col(aes(y = value)) +
  scale_fill_manual(values = mycolors) +
  ggtitle(paste(unique(plotDF.ordered$subtype),collapse = ", ")) +
  ggthemes::theme_clean() +
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_blank()) + xlab(NULL)  

ggsave(plot = plot,filename = file.path(dir.plots,"Distribution_of_mutational_signatures_for_individual_TCGA_tumors_stratified_by_TNBC_subtype.pdf"))
