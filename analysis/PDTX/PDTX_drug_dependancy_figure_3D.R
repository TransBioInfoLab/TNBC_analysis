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
# Date: 21 June 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(ggplot2)
library(reshape2)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.data <- "data/PDTX//"
dir.plots <- "plots/PDTX"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
TNBC_ballon <- read.csv(
  file.path(dir.data,"Caldas_samples_tvalue_selected.csv"),
  sep = ",",
  na.strings = "N/A",
  header = TRUE,
  comment.char = "#"
) %>% dplyr::mutate(X = NULL)
TNBC_ballon <- TNBC_ballon[grep("AZD8931",TNBC_ballon$DRUG_TARGET,invert = TRUE),]

# we will in the next stepts order the data by the drug sensitivity (negative t-value)

# First for each drug map it to the TNBC subtype it has the most negative effect
TNBC_ballon <- TNBC_ballon %>% 
  group_by(DRUG_TARGET) %>% 
  summarise("Lowest_drug_sentivity" = SUBTYPE[which.min(effect)]) %>% 
  merge(TNBC_ballon)

# Now we will get the subtype and drug ordered by effect
order <- TNBC_ballon %>% 
  group_by(Lowest_drug_sentivity) %>%  filter(SUBTYPE == Lowest_drug_sentivity) %>% 
  summarise("DRUG_TARGET" = DRUG_TARGET[order(effect, decreasing = F)])  %>% 
  dplyr::select(Lowest_drug_sentivity, DRUG_TARGET) 
order$order <- 1:nrow(order)

# Add new order to the table
TNBC_ballon <-  merge(TNBC_ballon,order)
TNBC_ballon <- TNBC_ballon[order(TNBC_ballon$order),]

# levels will be used to maintain the drugs in the correct order
TNBC_ballon$DRUG_TARGET <- factor(TNBC_ballon$DRUG_TARGET,levels = unique(TNBC_ballon$DRUG_TARGET) %>% as.character())

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plots
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
g <- ggplot(TNBC_ballon, aes(x = effect, y = DRUG_TARGET)) +
  geom_point(aes(
    alpha = 0.8,
    color = SUBTYPE,
    size = neg_log10p,
  )) +
  labs(color = "TNBC subtype") +
  labs(size = "-log10(p-val)") +
  geom_vline(xintercept = 0, linetype = "dashed",alpha = 0.4)  +
  scale_color_manual(values = c("red", "orange", "green", "blue","grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position = "right") + 
  xlim(-5, 5) +
  theme(axis.text = element_text(size = 4)) +
  ylab("Drug target") + 
  xlab("Drug sensitivity (T-value)") + 
  theme(legend.position = "bottom", legend.box = "horizontal") +
  guides(
    size = guide_legend(
      nrow = 1, 
      title.position = "top", 
      title.hjust = 0.5
    )
  ) + guides( 
    color = guide_legend(
      nrow = 1, 
      title.position = "top", 
      title.hjust = 0.5,
      override.aes = list(size = 3)
    )
  ) + theme(legend.key = element_rect(fill = "white", colour = "white")) 

p <- g + facet_grid(rows = vars(Lowest_drug_sentivity),scales = "free") + theme(strip.text.y = element_blank())

# we will fix the high of the each grid based on the number of genes
library(grid)
gt = ggplot_gtable(ggplot_build(p))
size <- gt$heights[7]
gt$heights[7] <- 0.4 * size
gt$heights[9] <- 0.15 * size
gt$heights[11] <- 0.3 * size
gt$heights[13] <- 0.2 * size

pdf(file.path(dir.plots,"caldas_balloon.pdf"),width = 8,height = 9)
grid.draw(gt)
dev.off()