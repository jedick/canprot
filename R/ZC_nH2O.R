# canprot/R/ZC_nH2O.R
# plot and summarize ZC and nH2O/residue of proteins
# 20170706 jmd

ZC_nH2O <- function(pdat, plot.it=TRUE) {
  nH2O <- pdat$pcomp$residue.basis[, "H2O"]
  ZC <- pdat$pcomp$ZC
  if(plot.it) {
    # set symbol shape and color
    # hollow red square for up, filled blue circle for down in cancer
    col <- ifelse(pdat$up2, cpcol$red, cpcol$blue)
    pch <- ifelse(pdat$up2, 0, 19)
    cex <- ifelse(pdat$up2, 1, 0.8)
    # shuffle the order of points to mitigate overplotting effects
    i <- sample(1:length(ZC))
    plot(ZC[i], nH2O[i], xlab=cplab$ZC, ylab=cplab$nH2O, col=col[i], pch=pch[i], cex=cex[i])
    title(pdat$description)
    mtext(pdat$dataset, side=4, cex=0.85, las=0, adj=0, line=-0.1)
  }
  # calculate and print sample size, difference of means, CLES, p-value
  ZC1 <- ZC[!pdat$up2]
  ZC2 <- ZC[pdat$up2]
  ZC.mean1 <- mean(ZC1)
  ZC.mean2 <- mean(ZC2)
  ZC.diff <- ZC.mean2 - ZC.mean1
  nH2O_1 <- nH2O[!pdat$up2]
  nH2O_2 <- nH2O[pdat$up2]
  nH2O.mean1 <- mean(nH2O_1)
  nH2O.mean2 <- mean(nH2O_2)
  nH2O.diff <- nH2O.mean2 - nH2O.mean1
  ZC.p.value <- stats::wilcox.test(ZC1, ZC2)$p.value
  ZC.CLES <- 100*CLES(ZC1, ZC2)
  nH2O.p.value <- stats::wilcox.test(nH2O_1, nH2O_2)$p.value
  nH2O.CLES <- 100*CLES(nH2O_1, nH2O_2)
  message(paste0(pdat$dataset, " (", pdat$description ,"): n1 ", length(ZC1), ", n2 ", length(ZC2)))
  message(paste0("ZC     MD ", format(round(ZC.diff, 3), nsmall=3, width=6),
               ", CLES ", round(ZC.CLES), "%",
               ", p-value ", format(round(ZC.p.value, 3), nsmall=3)))
  message(paste0("nH2O   MD ", format(round(nH2O.mean2 - nH2O.mean1, 3), nsmall=3, width=6),
               ", CLES ", round(nH2O.CLES), "%",
               ", p-value ", format(round(nH2O.p.value, 3), nsmall=3), "\n"))
  out <- data.frame(dataset=pdat$dataset, description=pdat$description,
    n1=length(ZC1), n2=length(ZC2),
    ZC.mean1, ZC.mean2, ZC.diff, ZC.CLES, ZC.p.value,
    nH2O.mean1, nH2O.mean2, nH2O.diff, nH2O.CLES, nH2O.p.value, stringsAsFactors=FALSE)
  return(invisible(out))
}
