#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Description
# The following script will create heatmap with GSVA results
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
# Date: 21 Jun 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VisualizeCNplot <-  function(
    scores_filt_call,
    titleplot = "test",
    tabAnno_conf_Amp,
    tabAnno_conf_Del,
    chrArms = TRUE,
    FDR.thresh,
    GeneToVisualize,
    anno_Amp_GeneWidePeak,
    anno_Del_GeneWidePeak,
    show.names = "all"
) {
    
    scores_filt_call$Chromosome <- paste0("chr", scores_filt_call$Chromosome)
    for (i in 1:9) {
        scores_filt_call[scores_filt_call$Chromosome == paste0("chr", i), "Chromosome"] <- paste0("chr0", i)
    }
    
    scores_filt_call <- cbind(scores_filt_call, ScoreSignif = rep(NA, nrow(scores_filt_call)))
    scores_filt_call <- cbind(scores_filt_call, ChrPosLine = rep(0, nrow(scores_filt_call)))
    scores_filt_call <- cbind(scores_filt_call, GeneAmp = rep("No", nrow(scores_filt_call)))
    scores_filt_call <- cbind(scores_filt_call, GeneDel = rep("No", nrow(scores_filt_call)))
    
    scores_filt_call$GeneAmp <- as.character(scores_filt_call$GeneAmp)
    scores_filt_call$GeneDel <- as.character(scores_filt_call$GeneDel)
    
    scores_filt_call_merged <- NULL
    
    # Significance threshold for q-values.  
    # Regions with q-values below this number are considered significant. (DEFAULT=0.25)
    # -log10(q-value) FDR in Gistic2 as default is 0.25
    # FDR.thresh <- 1 / 10^0.25
    # -log10(FDR.thresh)
    
    scores_filt_FDR_Amp <- scores_filt_call[scores_filt_call$`-log10(q-value)` > -log10(FDR.thresh), ]
    scores_filt_FDR_Del <- scores_filt_call[scores_filt_call$`-log10(q-value)` > -log10(FDR.thresh), ]
    
    scores_filt_FDR_Amp <- scores_filt_FDR_Amp[scores_filt_FDR_Amp$Type %in% "Amp", ]
    scores_filt_FDR_Del <- scores_filt_FDR_Del[scores_filt_FDR_Del$Type %in% "Del", ]
    
    transformationFormulaAxis2 <- rbind(scores_filt_FDR_Amp, scores_filt_FDR_Del)
    transformationFormulaAxis2_thresh <- mean(
        transformationFormulaAxis2$score / transformationFormulaAxis2$`-log10(q-value)`
    )
    
    scoreThreshAmp <- min(scores_filt_FDR_Amp$score)
    scoreThreshDel <- min(scores_filt_FDR_Del$score)
    
    if (chrArms == TRUE) {
        for (j in 1:23) {
            curChr <- paste0("chr0", j)
            if (j > 9) {
                curChr <- paste0("chr", j)
            }

            scores_filt_call_cur <- scores_filt_call[scores_filt_call$Chromosome %in% curChr, ]
            scores_filt_call_cur <- scores_filt_call_cur[order(scores_filt_call_cur$`Region Start [bp]`, decreasing = FALSE), ]
            
            if (nrow(scores_filt_call_cur) != 0) {
                scores_filt_call_cur[scores_filt_call_cur$`Region Start [bp]` ==  max(scores_filt_call_cur$`Region Start [bp]`), "ChrPosLine"] <- curChr
            }
            
            tabAnno_conf <- tabAnno_conf_Amp
           
            if (length(colnames(tabAnno_conf)) > 1) {
                curpos <- 1
                for (k in 2:(length(colnames(tabAnno_conf)) - 1)) {
                    curSignif_Chr_q <-
                        unlist(strsplit(as.character(colnames(
                            tabAnno_conf
                        )[k]), "q"))[1]
                    curSignif_Chr_p <-
                        unlist(strsplit(as.character(colnames(
                            tabAnno_conf
                        )[k]), "p"))[1]
                    
                    if (curSignif_Chr_p == j) {
                        if (length(intersect(scores_filt_call_cur$Type, "Amp")) == 1) {
                            scores_filt_call_cur[grep("Amp", scores_filt_call_cur$Type)[curpos], "ScoreSignif"] <- colnames(tabAnno_conf)[k]
                            
                            GeneSignif_Amp_Del <- as.matrix(tabAnno_conf[, colnames(tabAnno_conf)[k]])
                            if (length(intersect(GeneSignif_Amp_Del, GeneToVisualize)) != 0) {
                                scores_filt_call_cur[grep("Amp", scores_filt_call_cur$Type)[curpos], "GeneAmp"] <-
                                    as.character(paste(
                                        intersect(GeneSignif_Amp_Del, GeneToVisualize),
                                        collapse = ";"
                                    ))
                            }
                            
                            curpos <- curpos + 1
                         
                        }
                    }
                    
                   
                    if (curSignif_Chr_q == j) {
                        if (length(intersect(scores_filt_call_cur$Type, "Amp")) == 1) {
                            scores_filt_call_cur[grep("Amp", scores_filt_call_cur$Type)[curpos], "ScoreSignif"] <- colnames(tabAnno_conf)[k]
                            
                            GeneSignif_Amp_Del <- as.matrix(tabAnno_conf[, colnames(tabAnno_conf)[k]])
                          
                            if (length(intersect(GeneSignif_Amp_Del, GeneToVisualize)) != 0) {
                                scores_filt_call_cur[grep("Amp", scores_filt_call_cur$Type)[curpos], "GeneAmp"] <-
                                    as.character(paste(
                                        intersect(GeneSignif_Amp_Del, GeneToVisualize),
                                        collapse = ";"
                                    ))
                            }
                            
                            curpos <- curpos + 1
                           
                        }
                    }
                    
                }
            } #end amp
            
            # scores_filt_call_merged <- rbind(scores_filt_call_merged,scores_filt_call_cur )
            # working with Del
            
            
            tabAnno_conf <- tabAnno_conf_Del
            
            if (length(colnames(tabAnno_conf)) > 1) {
                curpos <- 1
                for (k in 2:(length(colnames(tabAnno_conf)) - 1)) {
                    curSignif_Chr_q <-
                        unlist(strsplit(as.character(colnames(
                            tabAnno_conf
                        )[k]), "q"))[1]
                    curSignif_Chr_p <-
                        unlist(strsplit(as.character(colnames(
                            tabAnno_conf
                        )[k]), "p"))[1]
                    
                    if (curSignif_Chr_p == j) {
                        if (length(intersect(scores_filt_call_cur$Type, "Del")) == 1) {
                            scores_filt_call_cur[grep("Del", scores_filt_call_cur$Type)[curpos], "ScoreSignif"] <- colnames(tabAnno_conf)[k]
                            
                            GeneSignif_Amp_Del <- as.matrix(tabAnno_conf[, colnames(tabAnno_conf)[k]])
                            if (length(intersect(GeneSignif_Amp_Del, GeneToVisualize)) !=
                                0) {
                                scores_filt_call_cur[grep("Del", scores_filt_call_cur$Type)[curpos], "GeneDel"] <-
                                    as.character(paste(
                                        intersect(GeneSignif_Amp_Del, GeneToVisualize),
                                        collapse = ";"
                                    ))
                            }
                            
                            curpos <- curpos + 1
                            
                        }
                    }
                    
                    if (curSignif_Chr_q == j) {
                        if (length(intersect(scores_filt_call_cur$Type, "Del")) == 1) {
                            scores_filt_call_cur[grep("Del", scores_filt_call_cur$Type)[curpos], "ScoreSignif"] <- colnames(tabAnno_conf)[k]
                            
                            GeneSignif_Amp_Del <- as.matrix(tabAnno_conf[, colnames(tabAnno_conf)[k]])
                            if (length(intersect(GeneSignif_Amp_Del, GeneToVisualize)) != 0) {
                                scores_filt_call_cur[grep("Del", scores_filt_call_cur$Type)[curpos], "GeneDel"] <-
                                    as.character(paste(
                                        intersect(GeneSignif_Amp_Del, GeneToVisualize),
                                        collapse = ";"
                                    ))
                            }
                            
                            curpos <- curpos + 1
                        
                        }
                    }
                    
                }
            } #end del
            
            scores_filt_call_merged <- rbind(scores_filt_call_merged, scores_filt_call_cur)
        }
    }
    
    if (chrArms == FALSE) {
        scores_filt_call_merged <- scores_filt_call
    }
    
    df <- subset(
        scores_filt_call_merged,
        select = c(
            "Type",
            "score",
            "Chromosome",
            "Region Start [bp]",
            "Region End [bp]",
            "-log10(q-value)",
            "ScoreSignif",
            "ChrPosLine",
            "GeneAmp",
            "GeneDel"
        )
    )
    
    colnames(df)[colnames(df) == "Region Start [bp]"] <- "start"
    colnames(df)[colnames(df) == "Region End [bp]"] <- "end"
    
    df[df$Type == "Del", "score"] <- -(df[df$Type == "Del", "score"])
    
    df <- cbind(Xpos = 1:nrow(df), df)
    
    # adding gene symbol to related chr start end position
    
    # 1. working with Amp Genes
    if (length(intersect(df$GeneAmp, GeneToVisualize)) == 0) {
        df_GR_ann_Amp <- 0
        length(df_GR_ann_Amp) <- 0
    }
    
    if (length(intersect(df$GeneAmp, GeneToVisualize)) != 0) {
        geneAmp_list <- NULL
        for (id in 1:nrow(df)) {
            curGene <- df$GeneAmp[id]
            
            if (curGene != "No") {
                geneAmp_list <- c(geneAmp_list, curGene)
            }
        }
        geneAmp_list_split <- NULL
        
        for (id in 1:length(geneAmp_list)) {
            curGene <- geneAmp_list[id]
            
            if (curGene != "No") {
                geneAmp_list_split <- c(geneAmp_list_split, unlist(strsplit(curGene, ";")))
            }
        }
        
        geneAmp_list_split <- sort(geneAmp_list_split)
        
        df_GR_ann_Amp <- anno_Amp_GeneWidePeak[anno_Amp_GeneWidePeak$GeneSymbol %in% geneAmp_list_split, ]
        
        # ending with Amp Genes
        df_GR_ann_Amp <- df_GR_ann_Amp[order(df_GR_ann_Amp$score, decreasing = TRUE), ]
        df_GR_ann_Amp <- df_GR_ann_Amp[!duplicated(df_GR_ann_Amp$GeneSymbol), ]
    }
    
    # 2. working with Del Genes
    if (length(intersect(df$GeneDel, GeneToVisualize)) != 0) {
        geneDel_list <- NULL
        for (id in 1:nrow(df)) {
            curGene <- df$GeneDel[id]
            
            if (curGene != "No") {
                geneDel_list <- c(geneDel_list, curGene)
            }
        }
        geneDel_list_split <- NULL
        
        for (id in 1:length(geneDel_list)) {
            curGene <- geneDel_list[id]
            
            if (curGene != "No") {
                geneDel_list_split <- c(geneDel_list_split, unlist(strsplit(curGene, ";")))
            }
        }
        
        geneDel_list_split <- unique(sort(geneDel_list_split))
       
        df_GR_ann_Del <- anno_Del_GeneWidePeak[anno_Del_GeneWidePeak$GeneSymbol %in% geneDel_list_split, ]
        
        # ending with Amp Genes
        df_GR_ann_Del <- df_GR_ann_Del[order(df_GR_ann_Del$score, decreasing = TRUE), ]
        df_GR_ann_Del <- df_GR_ann_Del[!duplicated(df_GR_ann_Del$GeneSymbol), ]
        # ending with Del Genes
    }
    # ending gene symbol to related chr start end position
    
    df_mod <- df
    df_mod$GeneAmp <- rep("No", nrow(df_mod))
    df_mod$GeneDel <- rep("No", nrow(df_mod))
  
    if (length(df_GR_ann_Del) != 0) {
        df_GR_ann_merged <- df_GR_ann_Del
    }
    if (length(df_GR_ann_Amp) != 0) {
        df_GR_ann_merged <- rbind(df_GR_ann_Amp, df_GR_ann_Del)
    }
    
    df_mod <- cbind(
        df_mod,
        ChrStartEnd = paste0(
            gsub("chr", "", gsub("chr0", "", df_mod$Chromosome)),
            "_", df_mod$start, "_", df_mod$end)
    )
    
    df_mod$ChrStartEnd <- as.character(df_mod$ChrStartEnd)
    
    df_GR_ann_merged <- cbind(
        df_GR_ann_merged,
        ChrStartEnd = paste0(
            df_GR_ann_merged$Chr,
            "_",
            df_GR_ann_merged$start,
            "_",
            df_GR_ann_merged$end
        )
    )
    df_GR_ann_merged$ChrStartEnd <- as.character(df_GR_ann_merged$ChrStartEnd)
    

    if (show.names == "all") {
        df_GR_ann_merged <- df_GR_ann_merged
    } else if (show.names == "significant") {
        df_GR_ann_merged <- df_GR_ann_merged[df_GR_ann_merged$`-log10(q-value)` > -log10(FDR.thresh), ]
    }

    for (curSelIew in c("Amp", "Del")) {
        tabGRcur_ampDel <- df_GR_ann_merged[df_GR_ann_merged$Type %in% curSelIew, ]
        
        for (iw in 1:nrow(tabGRcur_ampDel)) {
            curgene <- tabGRcur_ampDel[iw, ]
            tabGRcur <- tabGRcur_ampDel[tabGRcur_ampDel$GeneSymbol %in% curgene$GeneSymbol, ]
            df_mod[df_mod$ChrStartEnd == tabGRcur$ChrStartEnd, paste0("Gene", curSelIew)] <- curgene$GeneSymbol
        }
    }

    df <- df_mod
    
    df$Chromosome <- factor(df$Chromosome, levels = unique(sort(df$Chromosome)))
    require(ggrepel)
    
    p4 <- ggplot() + geom_bar(
        aes(y = score, x = Xpos, col = Type),
        data = df,
        stat = "identity"
    ) 
    p4 <- p4 + scale_y_continuous(limits = c(-0.75, 1.5))
    p4 <- p4 + geom_hline(
        yintercept = scoreThreshAmp,
        linetype = "dashed",
        color = "orange"
    )
    p4 <- p4 + geom_hline(
        yintercept = -scoreThreshDel,
        linetype = "dashed",
        color = "orange"
    )
    
    for (ichr in 1:23) {
        curPosChrsel <- df[df$Chromosome %in% paste0("chr0", ichr), ]
        if (ichr > 9) {
            curPosChrsel <- df[df$Chromosome %in% paste0("chr", ichr), ]
        }
        
        curPosChrsel <- curPosChrsel[curPosChrsel$ChrPosLine != 0, ]
        print(paste0("chr", ichr, " and ", curPosChrsel$Xpos))
        p4 <- p4 + geom_vline(
            xintercept = curPosChrsel$Xpos,
            linetype = "dashed",
            color = "black",
            alpha = 0.2)
        
    }
    # Add chr names
    chr.labels <- df %>% dplyr::filter(ChrPosLine != 0)
    p4 <- p4 + scale_x_continuous(
        breaks = c(chr.labels$Xpos),
        limits = c(0,max(df$Xpos)), 
        expand = c(0,0),
        labels = gsub("chr23","chrX",c(paste0(chr.labels$Chromosome)))
    )
    
    p4 <- p4 + labs(
        title = paste0("Amplification and Deletion with Gistic2 in ", titleplot),
        x = "Chromosomes",
        y = "G-score"
    )
    
    df_Amp_gene <- df[df$Type %in% "Amp", ]
    df_Amp_gene <- df_Amp_gene[df_Amp_gene$GeneAmp != "No", ]
    
    df_Del_gene <- df[df$Type %in% "Del", ]
    df_Del_gene <- df_Del_gene[df_Del_gene$GeneDel != "No", ]
    
    if (nrow(df_Amp_gene) != 0) {
        p4 <- p4 + geom_point(
            data = df_Amp_gene,
            shape = 1,
            color = "black",
            aes(
                #color = Type,
                x = Xpos, 
                y = score)
        )
        
        p4 <- p4 + scale_shape(solid = FALSE)
        
        p4 <- p4 + geom_label_repel(
            data = df_Amp_gene,
            aes(label = GeneAmp, x = Xpos, y = score),
            # size = 5,
            nudge_y      = 0.5,
            segment.size  = 0.2,
            nudge_x      = 0.5,
            box.padding   = 0.35,
            point.padding = 0.5,
            size = 3,
            label.size = 0.1,
            segment.color = 'grey50'
        )
    }
    
    
    if (nrow(df_Del_gene) != 0) {
        p4 <- p4 + geom_point(
            data = df_Del_gene,
            shape = 1,
            color = "black",
            aes(
                x = Xpos, 
                y = score
            )
        )
        
        p4 <- p4 + geom_label_repel(
            data = df_Del_gene,
            aes(label = GeneDel, x = Xpos, y = score),
            nudge_y      = -0.75,
            nudge_x      = -0.75,
            box.padding   = 0.35,
            point.padding = 0.5,
            segment.size  = 0.2,
            size = 3,
            label.size = 0.1,
            segment.color = 'grey50'
        )
        
    }
    
    return(p4)
}