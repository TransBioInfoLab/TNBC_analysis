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
library(ggplot2)
library(ggVennDiagram)
library(ggpubr)

#devtools::install_github("gaospecial/ggVennDiagram")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/cell_lines/"
dir.plots <- "plots/cell_lines"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
peaks <- readr::read_csv(file = file.path(dir.data,"Promoter_peaks.csv"))
tab <- readr::read_csv(file = file.path(dir.data,"H3K27_peak_union.csv"))

list <-  list(
  "CAL51" = na.omit(peaks$CAL51),
  "CAL120" = na.omit(peaks$CAL120),
  "BT549" = na.omit(peaks$BT549)
)
venn <- ggVennDiagram(
  x = list, 
  set_color = c("#edb477", "#699583", "#454d93"),
  label_alpha = 0,
  label_color = "white"
)

ggsave(plot = venn,filename = file.path(dir.plots,"figure_S11_B.pdf"))

tab$log10fdr <- -log10(tab$`FDR q-value`)
tab$Pathway <- tab$`Gene Set Name`
p <- ggbarplot(
  tab,
  x = "Gene Set Name",
  y = "log10fdr",
  orientation = "horiz",
  fill = "#75294c",
  sort.val = "asc",
  xlab = "Gene ontology pathway",
  ylab = expression(paste(-Log[10], " (FDR q-value)")),
) + theme(
  axis.line.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y = element_blank()
) + geom_text(aes(label=Pathway, y = 1), hjust = 'left', size = 4, color = "white") 
ggsave(plot = p,filename = file.path(dir.plots,"figure_S11_C.pdf"))

