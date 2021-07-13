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
library(survminer)
library(survival)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots/METABRIC"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
metabric <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "B-METABRIC"
) %>% dplyr::filter(!is.na(TNBCtype)) %>% 
  dplyr::filter(TNBCtype != "UNS")

metabric$TNBCtype <- factor( metabric$TNBCtype,levels = c("BL1","BL2","LAR","M"))

s <- survfit(formula = Surv(OS_YEAR_10YR, OS_EVENT_10YR) ~ TNBCtype, data = metabric, type = "kaplan-meier",conf.type="log")
print(s)

colors.tnbc <- c(
  "BL1" = "red",
  "BL2" = "orange",
  "LAR" = "green",
  "M" = "blue"
)

plot.subtype <- ggsurvplot(
  fit = s,
  data = metabric,
  size = 1,
  # change line size
  palette = colors.tnbc,
  # custom color palettes
  conf.int = FALSE,
  # Add confidence interval
  pval = TRUE,
  # Add p-value
  pval.size = 3,
  # Add p-value
  risk.table = FALSE,
  # Add risk table
  risk.table.col = "strata",
  # Risk table color by groups
  legend.labs = c("BL1", "BL2", "LAR", "M"),
  # Change legend labels
  surv.plot.height = 2,
  xlab="Time in Years",
  pval.coord = c(8, 0.0),
  legend.title="TNBC subtype",
  test.for.trend = TRUE,
  ggtheme = theme(
    aspect.ratio = 1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = c(0.1, 0.1),
    axis.line = element_line(colour = "black")
  ),
  tables.theme =  theme(aspect.ratio = 0.05),
) 


metabric$ESTIMATE_immune <- factor( metabric$ESTIMATE_immune,levels = c("low","high","med"))
s <- survfit(formula = Surv(OS_YEAR_10YR, OS_EVENT_10YR) ~ ESTIMATE_immune, data = metabric, type = "kaplan-meier",conf.type="log")
print(s)
plot.ESTIMATE_immune <- ggsurvplot(
  fit = s,
  data = metabric,
  size = 1,
  # change line size
  palette = c("high" = "red","med" = "orange","low" ="blue"),
  # custom color palettes
  conf.int = FALSE,
  # Add confidence interval
  pval = TRUE,
  # Add p-value
  pval.size = 3,
  # Add p-value
  risk.table = FALSE,
  # Add risk table
  risk.table.col = "strata",
  # Risk table color by groups  
  legend.labs = c("low","high","med"),
  # Change legend labels
  surv.plot.height = 2,
  test.for.trend = TRUE,
  xlab="Time in Years",
  pval.coord = c(8, 0.0),
  legend.title="Immune ESTIMATE",
  ggtheme = theme(
    aspect.ratio = 1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ),
  tables.theme =  theme(aspect.ratio = 0.05),
) 

all <- survminer::arrange_ggsurvplots(
  list("ESTIMATE_immune" = plot.ESTIMATE_immune,"Subtype" = plot.subtype),
  ncol = 1,
  nrow = 2,
  print = FALSE
)
ggsave(
  filename = file.path(dir.plots,"Metrabric_Survival_subtype_and_ESTIMATE_immune_S4_C_D.pdf"),
  width = 5,height = 8,
  plot = all
)
