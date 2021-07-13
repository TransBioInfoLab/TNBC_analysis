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
library(ggpubr)
library(survminer)
library(survival)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data,"TCGA")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Aux function
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
survival_plot <- function(
  formula,
  palette = NULL,
  legend.labs = NULL,
  data = NULL,
  facet.by = NULL
) {
  
  p <- ggsurvplot(
    fit = formula,
    data = data, 
    size = 1,                 # change line size
    palette = palette,# custom color palettes
    conf.int = FALSE,          # Add confidence interval
    pval = TRUE, 
    pval.size = 4,# Add p-value
    risk.table = FALSE,        # Add risk table
    risk.table.col = "strata",# Risk table color by groups
    legend.labs = legend.labs,    # Change legend labels
    surv.plot.height = 2,
    test.for.trend = FALSE,
    facet.by = NULL,
    ggtheme = theme(
      aspect.ratio = 1, 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    ),
    tables.theme =  theme(aspect.ratio = 0.05)
  )
  #ggpar(p,font.x = c(14, "plain", "black"),font.tickslab=c(14, "plain", "black"), font.y=c(14, "plain", "black"),y.text.angle = 90 )
  return(p)
}

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Get data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
metadata <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS") %>% 
  dplyr::mutate(PFI.time.year = pmin((PFI.time / 365),10)) %>% 
  dplyr::mutate(subtype = factor(subtype,levels = c("BL1","BL2","M","LAR"))) %>%
  dplyr::mutate(MutationalLoad_0_low_1_high = factor(MutationalLoad_0_low_1_high,levels = c("0","1","NA")))

s <- survfit(formula = Surv(PFI.time.year, PFI) ~ subtype, data = metadata, type = "kaplan-meier",conf.type="log")
plot.subtype <- survival_plot(
  formula = s,   
  data = metadata,
  palette = c("subtype=BL1" = "red","subtype=M" = "blue", "subtype=LAR"= "green", "subtype=BL2"= "orange")
)

s <- survfit(
  formula = Surv(PFI.time.year, PFI) ~ MutationalLoad_0_low_1_high, 
  data = metadata %>% dplyr::filter(MutationalLoad_0_low_1_high != "NA"), 
  type = "kaplan-meier",
  conf.type="log"
)
plot.mutation.load.all <- survival_plot(
  formula = s,  
  data = metadata %>% dplyr::filter(MutationalLoad_0_low_1_high != "NA"),
  legend.labs = c("Low","High"),
  palette = c("Low" = "blue","High" = "red","NA" = "grey")
) + ggtitle("TNBC") + labs("color" = "Mutational Load")

plot.mutation.load.subtype <- plyr::llply(metadata$subtype %>% unique,.fun = function(x){
  s <- survfit(
    formula = Surv(PFI.time.year, PFI) ~ MutationalLoad_0_low_1_high, 
    data = metadata %>% dplyr::filter(subtype == x) %>% dplyr::filter(MutationalLoad_0_low_1_high != "NA"), 
    type = "kaplan-meier",
    conf.type = "log"
  )
  plot <- survival_plot(
    formula = s,   
    legend.labs = c("Low","High"),
    palette = c("Low" = "blue","High" = "red","NA" = "grey"),
    data = metadata %>% dplyr::filter(subtype == x) %>% dplyr::filter(MutationalLoad_0_low_1_high != "NA")
  )
  plot <- plot + ggtitle(x) + labs("color" = "Mutational Load")
  return(plot)
})
names(plot.mutation.load.subtype) <- metadata$subtype %>% unique
all <- survminer::arrange_ggsurvplots(c(list("All" = plot.mutation.load.all),plot.mutation.load.subtype),ncol = 1,nrow = 5)
ggsave(
  filename = file.path(dir.plots,"TCGA_Mutation_burden_survival.pdf"),
  width = 10,height = 20,
  plot = all
)
