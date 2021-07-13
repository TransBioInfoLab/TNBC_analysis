DMRbarplot <- function(
  methdiff,
  plot,
  qvalue.cutoff,
  meth.cutoff,
  exclude,
  subtitle,
  cMax,
  logFC = FALSE
) {
  
  if (logFC == FALSE) {
    diffPrint <-
      paste0(" & methylation diff. >=", meth.cutoff * 100, " %")
  } else {
    diffPrint <- paste0(" & methylation abs(logFC) >=", meth.cutoff)
  }
  
  x <- methdiff
  temp.hyper = x[x$qvalue < qvalue.cutoff & x$meth.diff >= meth.cutoff, ]
  temp.hypo = x[x$qvalue < qvalue.cutoff & x$meth.diff <= -meth.cutoff, ]
  
  dmc.hyper = 100 * nrow(temp.hyper) / nrow(x) # get percentages of hypo/ hyper
  dmc.hypo = 100 * nrow(temp.hypo) / nrow(x)
  
  all.hyper.hypo = data.frame(
    number.of.hypermethylated = nrow(temp.hyper),
    percentage.of.hypermethylated = dmc.hyper,
    number.of.hypomethylated = nrow(temp.hypo),
    percentage.of.hypomethylated = dmc.hypo
  )
  
  # plot barplot for percentage of DMCs per chr
  dmc.hyper.chr = merge(
    as.data.frame(table(temp.hyper$chr)),
    as.data.frame(table(x$chr)), by = "Var1"
  )
  dmc.hyper.chr = cbind(
    dmc.hyper.chr,
    perc = 100 * dmc.hyper.chr[, 2] / dmc.hyper.chr[, 3]
  )
  
  dmc.hypo.chr = merge(as.data.frame(table(temp.hypo$chr)),
                       as.data.frame(table(x$chr)), by = "Var1")
  dmc.hypo.chr = cbind(
    dmc.hypo.chr,
    perc = 100 * dmc.hypo.chr[, 2] / dmc.hypo.chr[, 3]
  )
  
  # merge hyper hypo per chromosome
  dmc.hyper.hypo = merge(
    dmc.hyper.chr[, c(1, 2, 4)],
    dmc.hypo.chr[, c(1, 2, 4)], by = "Var1"
  )
  
  names(dmc.hyper.hypo) = c(
    "chr",
    "number.of.hypermethylated",
    "percentage.of.hypermethylated",
    "number.of.hypomethylated",
    "percentage.of.hypomethylated"
  )
  if (plot) {
    if (!is.null(exclude)) {
      dmc.hyper.hypo = dmc.hyper.hypo[!dmc.hyper.hypo$chr %in% exclude, ]
    }
    
    
    matdmc <- dmc.hyper.hypo
    matdmc_hyper <- matdmc[, c(1, 3)]
    matdmc_hyper <- cbind(matdmc_hyper, DMR = rep("hyper", nrow(matdmc_hyper)))
    
    matdmc_hypo <- matdmc[, c(1, 5)]
    matdmc_hypo <- cbind(matdmc_hypo, DMR = rep("hypo", nrow(matdmc_hypo)))
    
    
    matdmc_hyper <- as.data.frame(matdmc_hyper)
    matdmc_hypo <- as.data.frame(matdmc_hypo)
    
    colnames(matdmc_hyper)[2] <- "percentage"
    colnames(matdmc_hypo)[2] <- "percentage"
    
    colnames(matdmc_hyper)[1] <- "chr"
    colnames(matdmc_hypo)[1] <- "chr"
    
    matdmc_DMR <- rbind(matdmc_hyper, matdmc_hypo)
    matdmc_DMR$chr <- as.factor(matdmc_DMR$chr)
    matdmc_DMR$DMR <- as.factor(matdmc_DMR$DMR)
    
    matdmc_DMR <-  matdmc_DMR[order(matdmc_DMR$DMR, decreasing = FALSE), ]
    
    p4 <- ggplot() + geom_bar(
      aes(y = percentage, x = chr, fill = DMR),
      data = matdmc_DMR,
      stat = "identity"
    )
    
    p4 <- p4 + scale_x_discrete(limits = paste0("chr", 22:1))
    p4 <- p4 + scale_fill_manual(values = c("hyper" = "#f7cb55", "hypo" = "#3e96bd"))
    p4 <- p4 + coord_flip(ylim = c(c(0, cMax)))
    p4 <- p4 + labs(y = paste("% (percentage) \n", subtitle)) +
      ggtitle(
        paste0(
          "% of hyper & hypo methylated regions per chromosome",
          "\n qvalue < ",
          qvalue.cutoff,
          diffPrint
        )
      )
    
    
    p4 <- p4 + theme(text = element_text(size = 10)) + ggthemes::theme_clean()
    
    return(p4)
    
  } else{
    list(diffMeth.per.chr = dmc.hyper.hypo, diffMeth.all = all.hyper.hypo)
  }
  
}
