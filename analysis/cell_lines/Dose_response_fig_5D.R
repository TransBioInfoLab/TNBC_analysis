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
# Date: 29 June 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggsci)
library(plyr)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/drug_response/"
dir.plots <- "plots/drug_response/"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

files <- dir(dir.data,full.names = TRUE)

plots <- plyr::llply(files,.fun = function(f){
  data <- readr::read_csv(f)
  data_percent <- data * 100
  
  # remove E0771 and 231 and 148
  data_percent <- data_percent[,1:11]
  data_melted <- reshape2::melt(data_percent, id.var='DOSE')
  ###  ggplot the results  ###
  
  data_summarized <- ddply(
    data_melted, c("DOSE", "variable"), summarise,
    N    = 10,
    mean = mean(value,na.rm=TRUE),
    sd   = sd(value,na.rm=TRUE),
    se   = sd / sqrt(N)
  )
  
  library(RColorBrewer)
  myColors <- brewer.pal(length(levels(data_summarized$variable)),"RdYlBu")
  names(myColors) <- levels(data_summarized$variable)
  colScale <- scale_colour_manual(name = "grp",values = myColors)
  
  
  p1 <- ggplot(
    data = data_summarized,          # specify the data frame with data
    aes(x = DOSE, y = mean, col = variable)
  ) +  
    geom_smooth(method = "loess",se = FALSE)+  
    geom_point() +
    xlab(paste0(gsub("_drug.csv","",basename(f))," (uM)")) +   # label x-axis
    ylab("Viability (% control)") +    # label y-axis
    ylim(0, 151)+
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.05) +
    theme_bw() +  scale_colour_manual(name = "grp",values = myColors) + 
    scale_color_brewer(palette = "RdYlBu") + theme(aspect.ratio = 1)

  #theme(axis.text=element_text(size=10),axis.title=element_text(size=10))
  ##############
  

  data_percent$logDOSE <- log10(data_percent$DOSE)
  data_percent <- data_percent[,-1]
  data_melted_log <- reshape2::melt(data_percent, id.var = 'logDOSE')
  ###  ggplot the results  ###
  p2 <- ggplot(
    data=data_melted_log,          # specify the data frame with data
    aes(x = logDOSE, y = value, col = variable)) +   # specify x and y
    geom_smooth(method = "loess",se = FALSE) +  
    geom_point() +
    # make a scatter plot
    xlab(paste0(gsub("_drug.csv","",basename(f))," (uM)")) +   # label x-axis
    ylab("Viability (% control)") +    # label y-axis
    ylim(0, 151)+
    theme_bw() + 
    scale_color_brewer(palette = "RdYlBu")+ theme(aspect.ratio=1)
  list("mean" = p1,"log" = p2)
})
names(plots) <- gsub("_drug.csv","",basename(files))
merged <- ggpubr::ggarrange(
  plots$CPI$mean,
  plots$MAK$mean,
  plots$TAZ$mean,
  ncol = 3,common.legend = TRUE,
  legend = "right"
)

ggsave(plot = merged,filename = file.path(dir.plots,"Cell_line_drug_response_Figure_5D.pdf"),width = 10,height = 5)


merged.log <- ggpubr::ggarrange(
  plots$CPI$log,
  plots$MAK$log,
  plots$TAZ$log,
  ncol = 3,common.legend = TRUE,
  legend = "right"
)

ggsave(plot = merged.log,filename = file.path(dir.plots,"Cell_line_drug_response_log_Figure_5D.pdf"),width = 10,height = 5)
