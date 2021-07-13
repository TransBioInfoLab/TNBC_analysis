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
library(survival)
library(survminer)
library(TCGAbiolinks)
library(survMisc)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data"
dir.plots <- "plots/TCGA/survival/"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
S2.a <- readxl::read_xlsx(
  path = file.path(dir.data,"Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS" & Stage != "NA")  %>% 
  dplyr::select(patient, subtype,  OS.time, DSS.time, PFI.time, PFI, OS, DSS, Stage, TIME_CLASS, ImmuneScore...95, Age) %>%
  dplyr::mutate(Xcell_immune_score = as.numeric(ImmuneScore...95), DSS = as.numeric(DSS)) %>% 
  dplyr::mutate(dataset = "TCGA")

S2.a$Stage <- S2.a$Stage %>% {gsub("Stage_","",.)} %>% as.roman %>% as.numeric
S2.a$Age <- as.numeric(S2.a$Age)

threshHigh <- quantile(S2.a$Xcell_immune_score,2/3)
threshLow <- quantile(S2.a$Xcell_immune_score,1/3)

S2.a$ImmuneScore_Thresh <- "Med"
S2.a$ImmuneScore_Thresh[S2.a$Xcell_immune_score > threshHigh] <- "High"
S2.a$ImmuneScore_Thresh[S2.a$Xcell_immune_score < threshLow] <- "Low"
S2.a$ImmuneScore_Thresh <- factor(S2.a$ImmuneScore_Thresh, levels = c("High","Med","Low"))

S2.a$subtype <- factor(S2.a$subtype)
S2.a$TIME_CLASS <- factor(S2.a$TIME_CLASS)

aux <- S2.a
aux$ImmuneScore_Thresh <- aux$subtype <- aux$TIME_CLASS <- "TNBC all"
data <- S2.a %>% rbind(aux)
data$ImmuneScore_Thresh <- factor(data$ImmuneScore_Thresh, levels = c("TNBC all","High","Med","Low"))
data$subtype <- factor(data$subtype, levels = c("TNBC all",levels(S2.a$subtype)))
data$TIME_CLASS <- factor(data$TIME_CLASS, levels = c("TNBC all",levels(S2.a$TIME_CLASS)))


for(analysis in c("ImmuneScore_Thresh","TIME_CLASS","subtype")){
  
  df <- data.frame(
    survTime = c("PFI.time","OS.time","DSS.time"),
    survEvent = c("PFI","OS","DSS"),
    survStrata = c(analysis)
  )
  
  if (analysis == "ImmuneScore_Thresh"){
    palette = c("red","orange","green")
    legend.labs = c("High","Med","Low")  
  } else if (analysis == "TIME_CLASS"){
    palette = c(
      "FI" = "#FCFFA4FF",
      "ID" = "#420A68FF",
      "MR" = "#AE305CFF",
      "SR" = "#F8850FFF"
    )
    legend.labs = c("FI","ID","MR","SR")  
  } else if (analysis == "subtype"){
    palette = c(
      "BL1" = "red",
      "BL2" = "orange",
      "LAR" = "green",
      "M" = "blue"
    )
    
    legend.labs = c("BL1","BL2","LAR","M")  
  }
  
  plots <- plyr::alply(df,.margins = 1,.fun = function(row){
    
    survTitle <- paste0(row$survTime," ",row$survEvent," ~ ",row$survStrata)
    
    S2.a$Time <- S2.a[[row$survTime]]
    S2.a$Event <- S2.a[[row$survEvent]]
    S2.a$survStrata <- S2.a[[row$survStrata]]
    
    data$Time <- data[[row$survTime]]
    data$Event <- data[[row$survEvent]]
    data$survStrata <- data[[row$survStrata]]
    
    res.cox <- coxph(Surv(Time, Event) ~ survStrata + Age + Stage, data = data)
    p1 <- ggforest(
      res.cox, 
      data = data,
      main =  paste0(
        "Hazard ratio ",
        ifelse(row$survEvent == "PFI","progression-free interval", 
               ifelse(row$survEvent == "OS","Overall survival", "Disease-free survival")
        )
      )
    )
    
    s = survfit(
      formula = Surv(Time, Event) ~ survStrata,
      data = S2.a, 
      type = "kaplan-meier",
      conf.type = "log"
    )
    
    
    p2 <- ggsurvplot(
      s, 
      data = S2.a, 
      xlim = c(0,3652.5),
      break.time.by = 730.5,
      size = 1,                 # change line size
      palette = palette,
      legend.labs = legend.labs,   # Change legend labels
      conf.int = FALSE,          # Add confidence interval
      pval = TRUE, 
      pval.method = TRUE,
      xscale = "d_y",
      pval.size = 3,# Add p-value
      risk.table = FALSE,        # Add risk table
      risk.table.col = "strata",# Risk table color by groups
      surv.plot.height = 2,
      test.for.trend = FALSE,
      ggtheme = theme(
        aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
      ),
      tables.theme =  theme(aspect.ratio = 0.05)
    )
    
    pmerged <- ggarrange(
      p1,
      p2$plot,
      labels = NULL,
      ncol = 2, 
      nrow = 1,
      common.legend = FALSE
    )
  })
  pmerged <- ggarrange(
    plotlist = plots,
    ncol = 1, 
    nrow = 3,
    common.legend = FALSE
  )
  
  ggsave(
    pmerged , 
    filename = file.path(dir.plots,paste0("TNBC_HR_KM_",analysis,".pdf")),
    width = 16,
    height = 12
  )
}
